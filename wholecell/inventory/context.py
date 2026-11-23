"""Step context passed to adapters, builders and rules: config access, frozen inputs, stage paths, write-once outputs."""
from __future__ import annotations

import os
import shutil

from .. import build


class StepContext:
    """ctx.cfg (resolved CellConfig), ctx.inv (inventory section), ctx.dir (inventory stage folder),
    ctx.params(name) -> component parameters (inventory.params.<name>), ctx.frozen(name) -> declared frozen input path,
    ctx.path(rel) -> path inside the inventory folder, ctx.out(rel) -> write-once output path."""

    def __init__(self, cfg, bdir: str, step: str):
        self.cfg = cfg
        self.bdir = bdir
        self.step = step
        self.inv = cfg.get("inventory", {}) or {}
        self.dir = os.path.join(bdir, "inventory")
        self.report: list[str] = []

    # -- inputs
    def frozen(self, name: str) -> str:
        return self.cfg.frozen(name)

    def resolve(self, v):
        """a file reference in params: declared frozen input name | 'build:<path relative to the build dir>' | path"""
        if isinstance(v, str) and v in (self.cfg.get("frozen_inputs") or {}):
            return self.frozen(v)
        if isinstance(v, str) and v.startswith("build:"):
            return os.path.join(self.bdir, v[len("build:"):])
        return v

    def params(self, name: str) -> dict:
        return dict((self.inv.get("params", {}) or {}).get(name, {}) or {})

    def path(self, rel: str) -> str:
        return os.path.join(self.dir, rel)

    def need(self, rel: str) -> str:
        p = self.path(rel)
        if not os.path.exists(p):
            raise FileNotFoundError(f"inventory step '{self.step}' needs {p} (run the producing step first)")
        return p

    # -- outputs
    def out(self, rel: str) -> str:
        return build.write_once(self.path(rel))

    def write_tsv(self, df, rel: str, index: bool = False, **kw) -> str:
        p = self.out(rel)
        build.atomic_save(p, lambda t: df.to_csv(t, sep="\t", index=index, **kw))
        return p

    def copy(self, src: str, rel: str) -> str:
        p = self.out(rel)
        build.atomic_save(p, lambda t: shutil.copyfile(src, t))
        return p

    def say(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.report.append(s)

    def write_report(self, rel: str, lines=None):
        p = self.out(rel)
        txt = "\n".join(self.report if lines is None else lines) + "\n"
        build.atomic_save(p, lambda t: open(t, "w").write(txt))
        return p
