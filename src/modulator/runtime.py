from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional


def find_project_root(start: str | Path) -> Path:
    start_path = Path(start).resolve()
    for candidate in [start_path, *start_path.parents]:
        if (candidate / "workflow" / "scripts" / "assemble_transcripts.py").exists() and (candidate / "config" / "config.yaml").exists():
            return candidate
    raise FileNotFoundError(
        f"Could not find a modulator project root from {start_path}; expected workflow/scripts/assemble_transcripts.py and config/config.yaml."
    )


def ensure_mpl_config_dir(project_root: Path) -> None:
    cache_root = project_root / ".cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    mpl_dir = cache_root / "matplotlib"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))


def resolve_path(project_root: Path, value: str | os.PathLike[str] | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return (project_root / path).resolve()


def is_set(value) -> bool:
    if value is None:
        return False
    if isinstance(value, str) and value.strip().lower() in {"", "none", "null", "na"}:
        return False
    return True


def as_bool(value, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value != 0
    if isinstance(value, str):
        s = value.strip().lower()
        if s in {"1", "true", "yes", "y", "on"}:
            return True
        if s in {"0", "false", "no", "n", "off"}:
            return False
        if s in {"", "none", "null", "na"}:
            return default
    return bool(value)


def normalize_pivot_mode(value, default: str = "auto") -> str:
    """Resolve a pivot-output setting to one of 'auto' | 'on' | 'off'.

    Accepts the tri-state strings directly, and maps legacy booleans / truthy-falsey strings
    (True/False, "true"/"false", "yes"/"no", 1/0) onto 'on'/'off' so old configs keep working.
    """
    if value is None:
        return default
    if isinstance(value, str):
        s = value.strip().lower()
        if s in {"auto", "smart", "adaptive"}:
            return "auto"
        if s in {"", "none", "null", "na"}:
            return default
    return "on" if as_bool(value, default=(default == "on")) else "off"


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


@dataclass
class TaskResult:
    label: str
    completed: bool


def format_command(cmd: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in cmd)


def _read_rss_kib(pid: int) -> int:
    """Instantaneous proportional memory (KiB) of one pid via /proc; 0 if gone or non-Linux.

    Prefers Pss (from /proc/<pid>/smaps_rollup): a shared page is split across the processes
    sharing it, so SUMMING Pss over a process tree gives the true unique physical footprint. This
    matters because modulator's multi-process stages (assemble's worker pool, the genotype modkit
    pool, parallel modkit_zn samples) all fork and share the ~3 GB reference / libs via copy-on-write
    -- summing VmRSS there double-counts those shared pages and massively over-reports (e.g. 18 GiB
    where the real peak is ~3). Falls back to VmRSS when smaps_rollup is unavailable (older kernels,
    permissions), which is exact for single-process stages and only over-counts multi-process ones.
    """
    try:
        with open(f"/proc/{pid}/smaps_rollup", "r") as fh:
            for line in fh:
                if line.startswith("Pss:"):
                    return int(line.split()[1])          # /proc reports KiB
    except (FileNotFoundError, ProcessLookupError, ValueError, PermissionError):
        pass
    try:
        with open(f"/proc/{pid}/status", "r") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError, ValueError, PermissionError):
        return 0
    return 0


def _descendant_pids(root_pid: int) -> list[int]:
    """root_pid plus all transitive descendants, from /proc PPID links (Linux)."""
    try:
        entries = os.listdir("/proc")
    except FileNotFoundError:
        return [root_pid]
    children: dict[int, list[int]] = {}
    for name in entries:
        if not name.isdigit():
            continue
        pid = int(name)
        try:
            with open(f"/proc/{pid}/stat", "rb") as fh:
                data = fh.read()
            # field 2 (comm) may contain spaces/')'; ppid is the field after state.
            rparen = data.rfind(b")")
            after = data[rparen + 2:].split()
            ppid = int(after[1])
        except (FileNotFoundError, ProcessLookupError, IndexError, ValueError):
            continue
        children.setdefault(ppid, []).append(pid)
    out, stack, seen = [], [root_pid], set()
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        out.append(pid)
        stack.extend(children.get(pid, ()))
    return out


def _tree_rss_kib(root_pid: int) -> int:
    return sum(_read_rss_kib(p) for p in _descendant_pids(root_pid))


class PeakRSSSampler:
    """Background thread tracking the peak of the summed proportional memory (Pss) across
    root_pid's whole process subtree. Robust to run_parallel's concurrent subprocesses and to
    ProcessPoolExecutor grandchildren (all are descendants of root_pid). Pss splits shared/COW
    pages among sharers, so the sum over forked workers is the true physical footprint (not the
    inflated VmRSS sum). Pure stdlib; no-ops on non-Linux (no /proc)."""

    def __init__(self, root_pid: int, interval: float = 0.2):
        self.root_pid = root_pid
        self.interval = interval
        self.peak_kib = 0
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._supported = os.path.isdir("/proc")

    def _run(self) -> None:
        while not self._stop.is_set():
            self.peak_kib = max(self.peak_kib, _tree_rss_kib(self.root_pid))
            self._stop.wait(self.interval)
        self.peak_kib = max(self.peak_kib, _tree_rss_kib(self.root_pid))

    def __enter__(self) -> "PeakRSSSampler":
        if self._supported:
            self._thread = threading.Thread(target=self._run, name="rss-sampler", daemon=True)
            self._thread.start()
        return self

    def __exit__(self, *exc) -> bool:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2 * self.interval + 1.0)
        return False

    @property
    def peak_gib(self) -> float:
        return self.peak_kib / (1024.0 * 1024.0)


def run_command(cmd: list[str], *, cwd: Path, label: str, verbose: bool = True,
                on_peak: Optional[Callable[[str, float], None]] = None,
                sample_interval: float = 0.2) -> None:
    if verbose:
        print(f"[modulator] {label}: {format_command(cmd)}", flush=True)
    if on_peak is None:
        subprocess.run(cmd, cwd=str(cwd), check=True)
        return
    # Measured path: Popen so we can sample the child's whole subtree peak RSS.
    proc = subprocess.Popen(cmd, cwd=str(cwd))
    with PeakRSSSampler(proc.pid, interval=sample_interval) as sampler:
        ret = proc.wait()
    on_peak(label, sampler.peak_gib)
    if ret != 0:
        raise subprocess.CalledProcessError(ret, cmd)


def require_tools(tools: Iterable[str]) -> None:
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        raise RuntimeError(f"Missing required external tools in PATH: {', '.join(sorted(missing))}")


def run_parallel(tasks: list[tuple[str, Callable[[], None]]], jobs: int) -> list[TaskResult]:
    if not tasks:
        return []
    if jobs <= 1 or len(tasks) == 1:
        results = []
        for label, fn in tasks:
            fn()
            results.append(TaskResult(label=label, completed=True))
        return results

    results: list[TaskResult] = []
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        future_map = {executor.submit(fn): label for label, fn in tasks}
        for future in as_completed(future_map):
            label = future_map[future]
            future.result()
            results.append(TaskResult(label=label, completed=True))
    return results
