"""Serial pre-placed builders: an ordered list of registered `preplaced` builders that share one context.

    python -m wholecell preplaced --cell C --build B --step serial

Some pre-placed layers depend on each other (the yeast chromosome arms start at the kinetochores on the microtubule
plus ends; the FG mesh threads around the transport receptors; mRNA threads avoid the polysomes ...). Such builders
run in one process, in the configured order, on one shared context:

    preplaced:
      serial:
        context: <registered preplaced_context>    # sets up cell geometry, library pack, pools (cell-specific)
        context_params: {...}
        out: preplaced/serial                      # build-relative output folder
        builders:                                  # order = build order (and the order of every pool deduction)
          - {name: npc, builder: core.npc_pores, params: {...}}
          - ...

A serial builder is `fn(ctx, params)` registered as kind `preplaced`. It reads and extends the context:
  ctx.cell, ctx.lp            geometry evaluator and library pack (types, sphere trees, bounding radii)
  ctx.exp3                    expected copies per (library type, ingredient, region): builders deduct what they place
  ctx.pre                     the pre-pass arrays (written as preplaced.npz; content contract = yeast s3 prepass.npz)
  ctx.obs                     named obstacle sphere sets (centres, radii) handed from builder to builder
  ctx.state                   other values handed on (e.g. the microtubule dict, chromosome start points)
  ctx.K                       the geometry kernels of the cell's evaluator (kern.py); ctx.seed(s), ctx.lib_path(id)
  ctx.log, ctx.consume(), ctx.emit(), ctx.pre_instances()
  (ctx.mod(<engine module>): the yeast context's vendored modules; used by yeast plugin builders only)

Core context (`core.cellgeom`, context.py): `exports` writes module files for the pack index from the instance groups
and bead arrays the builders produced (Jurkat).

Random numbers: every builder opens its own generator(s) from its stream number (params `seed`, ...); none draws
from a generator another builder used. Builders run in one process in the configured order, so the state of any
global RNG a numerical kernel touches (numba's np.random, re-seeded by the kernels at entry) is the same as in the
reference order as well.

Outputs (<bdir>/<out>/): preplaced.npz (all arrays), layers/<name>.npz (arrays each builder added or replaced),
expected_after_preplaced.parquet (ctx.exp3 after the last builder), lib_ids.npy (ctx.lp ids: the type order of
preplaced.npz), serial_manifest.json (order, keys, counts, seconds; written last).
"""
from __future__ import annotations

import json
import os
import time

import numpy as np


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


class Context:
    """shared state of one serial pre-pass run (filled by the cell's preplaced_context, extended by the builders)"""

    def __init__(self, cell=None, lp=None, exp3=None, log=log, **kw):
        self.cell, self.lp, self.exp3, self.log = cell, lp, exp3, log
        self.pre: dict = {}
        self.obs: dict = {}
        self.state: dict = {}
        self._mod = kw.pop("mod", None)
        self._pre_instances = kw.pop("pre_instances", None)
        self._emit = kw.pop("emit", None)
        for k, v in kw.items():
            setattr(self, k, v)

    def mod(self, name):
        """a numerical module of the cell's engine (yeast: the vendored engine pipeline module of that name)"""
        return self._mod(name)

    def emit(self, types, pos, q, lo=None, hi=None):
        """sphere-tree spheres (centres, radii) of instances (types, positions, quaternions) inside [lo, hi]"""
        return self._emit(types, pos, q, lo, hi)

    def pre_instances(self, pre=None):
        """types, positions, quaternions of every pre-placed instance group in pre (default ctx.pre)"""
        return self._pre_instances(self.pre if pre is None else pre)

    def consume(self, ing_ids, n):
        """remove n instances from the expected totals of these ingredients (all regions, proportionally); returns n taken"""
        # origin: yeast scripts/s3_pack.py:consume @ e3b7448 (unchanged)
        exp3 = self.exp3
        r = exp3.ing_id.isin(ing_ids); tot = float(exp3.loc[r, "expected"].sum()); take = min(max(n, 0), tot)
        if tot > 0:
            exp3.loc[r, "expected"] *= 1 - take / tot
        self.log(f"  {'/'.join(ing_ids)}: {take:.0f} of {tot:.0f} used ({n:.0f} needed)")
        return take


def _n(v):
    return int(len(v)) if np.ndim(v) else 1


def run_builders(ctx, builders, out_dir=None, write_layers=True):
    """run the serial builders in order on ctx; returns the manifest rows"""
    from wholecell import registry
    rows = []
    for b in builders:
        fn = registry.get("preplaced", b["builder"])
        before = dict(ctx.pre)
        t0 = time.time()
        fn(ctx, dict(b.get("params") or {}))
        sec = time.time() - t0
        keys = [k for k, v in ctx.pre.items() if k not in before or before[k] is not v]
        rows.append(dict(name=b["name"], builder=b["builder"], keys=keys, counts={k: _n(ctx.pre[k]) for k in keys},
                         seconds=round(sec, 1)))
        print(f"[preplaced] {b['name']} ({b['builder']}): {len(keys)} arrays, {sec:.0f} s", flush=True)
        if out_dir and write_layers and keys:
            os.makedirs(os.path.join(out_dir, "layers"), exist_ok=True)
            np.savez(os.path.join(out_dir, "layers", f"{b['name']}.npz"), **{k: ctx.pre[k] for k in keys})
    return rows


def run(cfg, bdir, out=None):
    """the configured serial pre-pass of this cell into <bdir>/<out> (default preplaced.serial.out or preplaced/serial)"""
    from wholecell import registry
    from . import builders, context  # noqa: F401  (register the core builders and the core context core.cellgeom)
    sc = (cfg.get("preplaced") or {}).get("serial")
    if not sc:
        raise SystemExit("preplaced --step serial: no preplaced.serial section in the cell config")
    out_dir = os.path.join(bdir, out or sc.get("out", "preplaced/serial"))
    done = os.path.join(out_dir, "serial_manifest.json")
    if os.path.exists(done):
        print(f"preplaced serial: {done} exists, step skipped", flush=True)
        return
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    ctx = registry.get("preplaced_context", sc["context"])(cfg, bdir, dict(sc.get("context_params") or {}))
    rows = run_builders(ctx, sc["builders"], out_dir)
    np.savez(os.path.join(out_dir, "preplaced.npz"), **ctx.pre)
    ctx.exp3.to_parquet(os.path.join(out_dir, "expected_after_preplaced.parquet"))
    # library order of the context's library pack (the type index of preplaced.npz; the pack's type order, A9)
    np.save(os.path.join(out_dir, "lib_ids.npy"), np.array([str(i) for i in ctx.lp["ids"]]))
    exports = []
    if sc.get("exports"):                  # module files for the pack index (core context: wholecell/preplaced/context.py)
        from .context import export_modules
        exports = export_modules(ctx, sc, out_dir)
    tmp = done + f".tmp{os.getpid()}"
    json.dump(dict(context=sc["context"], builders=rows, arrays=len(ctx.pre), seconds=round(time.time() - t0, 1),
                   **({"exports": exports, "frozen_obstacles": getattr(ctx, "frozen_layers", [])} if exports else {})),
              open(tmp, "w"), indent=1)
    os.replace(tmp, done)
    print(f"[preplaced] serial: {len(rows)} builders, {len(ctx.pre)} arrays, {time.time() - t0:.0f} s", flush=True)
    return ctx
