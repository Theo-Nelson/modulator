from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


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


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


@dataclass
class TaskResult:
    label: str
    completed: bool


def format_command(cmd: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in cmd)


def run_command(cmd: list[str], *, cwd: Path, label: str, verbose: bool = True) -> None:
    if verbose:
        print(f"[modulator] {label}: {format_command(cmd)}", flush=True)
    subprocess.run(cmd, cwd=str(cwd), check=True)


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
