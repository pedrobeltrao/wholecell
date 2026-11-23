"""Seeding: one master seed per build, named streams derived from it (never Python hash()).

    rng = seeds.rng(cfg, "pack.counts")            # numpy Generator for a named stream
    s = seeds.seed_int(cfg, "pack.tile", tile_id)   # 32-bit integer seed for numba kernels

Stream ids are hashed with sha256, so they are stable across processes, machines and PYTHONHASHSEED values.

Legacy seeds (exact reproduction of a published build) are declared in the cell config under
`seeds.legacy: {<stream>: <int or expression in 'k'>}`; when present they take precedence for that stream and are
recorded in the build stamp. Example (Jurkat v7): seeds.legacy = {pack.counts: 3, pack.tile: "k*7919+11"}.
"""
from __future__ import annotations

import hashlib

import numpy as np


def _stream_key(name: str) -> int:
    return int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "little")


def master_seed(cfg) -> int:
    return int(cfg.get("seeds", {}).get("master", 0))


def legacy(cfg, name: str, k: int | None = None):
    """legacy seed for stream `name` (and key k), or None"""
    leg = (cfg.get("seeds", {}) or {}).get("legacy", {}) or {}
    if name not in leg:
        return None
    v = leg[name]
    if isinstance(v, str):
        return int(eval(v, {"__builtins__": {}}, {"k": int(k or 0)}))   # tiny arithmetic expressions in k only
    return int(v)


def seed_int(cfg, name: str, *keys: int) -> int:
    lg = legacy(cfg, name, keys[0] if keys else None)
    if lg is not None:
        return lg
    ss = np.random.SeedSequence([master_seed(cfg), _stream_key(name), *[int(k) for k in keys]])
    return int(ss.generate_state(1, dtype=np.uint32)[0])


def rng(cfg, name: str, *keys: int) -> np.random.Generator:
    return np.random.default_rng(seed_int(cfg, name, *keys))
