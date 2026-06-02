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


def apply_overrides(config: dict[str, Any], overrides: list[str]) -> dict[str, Any]:
    out = copy.deepcopy(config)
    for item in overrides:
        if "=" not in item:
            raise ValueError(f"Override must look like key=value, got: {item}")
        key, raw_value = item.split("=", 1)
        set_nested(out, key, parse_override_value(raw_value))
    return _normalize_dict(out)
