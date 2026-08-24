from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import yaml


AUTOPARSE_KEYS = (
    "assembler",
    "modkit",
    "aggregation",
    "toggles",
    "mods",
    "ref_bases",
    "test_diffs",
    "classify_diffs",
    "multigene_filter",
    "report",
    "genotype",
)


def _autoparse(value: Any) -> Any:
    if isinstance(value, str):
        s = value.strip()
        if (s.startswith("{") and s.endswith("}")) or (s.startswith("[") and s.endswith("]")):
            try:
                return json.loads(s)
            except Exception:
                try:
                    return yaml.safe_load(s)
                except Exception:
                    return value
    return value


def _normalize_dict(config: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(config)
    for key in AUTOPARSE_KEYS:
        if key in out:
            out[key] = _autoparse(out[key])
    return out


def load_config(path: str | Path) -> dict[str, Any]:
    cfg_path = Path(path)
    with cfg_path.open() as handle:
        loaded = yaml.safe_load(handle) or {}
    if not isinstance(loaded, dict):
        raise ValueError(f"Config at {cfg_path} must be a mapping.")
    return _normalize_dict(loaded)


def parse_override_value(raw: str) -> Any:
    try:
        return yaml.safe_load(raw)
    except Exception:
        return raw


def set_nested(config: dict[str, Any], dotted_key: str, value: Any) -> None:
    parts = dotted_key.split(".")
    cur = config
    for part in parts[:-1]:
        nxt = cur.get(part)
        if not isinstance(nxt, dict):
            nxt = {}
            cur[part] = nxt
        cur = nxt
    cur[parts[-1]] = value


def _warn_unknown_override_keys(base: dict[str, Any], keys: list[str]) -> None:
    """Warn (don't fail) when a `--set` key does not correspond to anything in the loaded config.

    A typo like ``genotype.enabel=true`` or ``genotpye.enable=true`` otherwise silently does nothing.
    We walk ``base`` segment by segment: a missing segment while the current node is still a dict is an
    unknown key. As soon as the current node is NOT a dict (a scalar/None leaf), we STOP -- everything
    below is user-defined structure (e.g. the ``nfail_score_k`` per-mod map, whose keys 'a'/'17802' are
    values, not schema), so we never false-positive on it.
    """
    import sys
    # Top-level keys the pipeline reads via .get() but that are NOT literals in the shipped config
    # (the user is expected to supply them via --set), so they must not be flagged as typos.
    extra_known_top = {"reference_fa", "reference_gtf"}
    for dotted in keys:
        node = base
        parts = dotted.split(".")
        if parts[0] in extra_known_top:
            continue
        for i, part in enumerate(parts):
            if not isinstance(node, dict):
                break                      # reached a scalar/None leaf -> user structure below; stop
            if part not in node:
                shown = ".".join(parts[: i + 1])
                print(f"[modulator] warning: --set {dotted!r} targets unknown config key {shown!r} "
                      f"(not in the loaded config); it will be created but likely does nothing. "
                      f"Check for a typo.", file=sys.stderr)
                break
            node = node[part]


def apply_overrides(config: dict[str, Any], overrides: list[str]) -> dict[str, Any]:
    out = copy.deepcopy(config)
    keys = []
    for item in overrides:
        if "=" not in item:
            raise ValueError(f"Override must look like key=value, got: {item}")
        key, raw_value = item.split("=", 1)
        keys.append(key)
        set_nested(out, key, parse_override_value(raw_value))
    # Validate against the ORIGINAL config (before overrides mutated it in) so a typo'd key that
    # set_nested just created is still reported as unknown.
    _warn_unknown_override_keys(config, keys)
    return _normalize_dict(out)
