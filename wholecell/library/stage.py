"""Library stage: python -m wholecell library --cell C --build B [--step import] [nproc=N]

Two modes (config `library.mode`):

frozen (Jurkat v7, testcell): import a FROZEN structure library (one npz per species: residue beads + sphere tree; TM
species oriented):
    library:
      mode: frozen
      frozen: <frozen input name: {path: <folder with <species>.npz (+ auxiliary files)>,
                                   manifest: <tsv(.gz) with columns file, md5, bytes; committed in the repo>, md5: <of it>}>
      (or frozen_dir / frozen_manifest directly)
      copy: false        # false: the build references the frozen folder (verified file by file); true: copy into the build
Writes <bdir>/library/library_manifest.tsv (species, path, md5) and library_meta.json. Every file listed in the frozen
manifest must exist with the same md5 and the folder must hold no extra file, else the stage fails.
Decision (2026-11-15): Jurkat references the frozen v7 folder (176 MB, 16,246 files) instead of copying it per build;
the md5 check of every file makes the reference as safe as a copy, and the folder is part of the read-only archive.

build (yeast policy, unification plan §1 "Structure library", A3): a scripted, canonical library written to a NEW
folder per build (<bdir>/library/lib, never edited in place: built in lib.partial/, renamed when complete, files made
read-only), from declared frozen start entries plus scripted builders, and derived packing shapes (pack_shapes/,
§2b P8). Steps, one plain function each (STEPS), run in this order:
    start          copy the frozen start entries (md5-checked against their manifest)
    <steps>        the builders listed in library.steps, in order (a later write of the same id wins, as in the yeast
                   s0 build): complexes | af_models | pdb_oligomers | pdb_rna (yeast);
                   jurkat_assemblies (option tm_graft) | jurkat_monomers | reorient_tm | tm_route | adaptive_trees
                   (Jurkat, B3: wholecell/library/jurkat.py; inputs under library.jurkat)
    seal           lib.partial -> lib, library_manifest.tsv (species, path, md5, origin), library_meta.json
    pack_shapes    pack_shapes/pack_trees_<per>.npz + pack_meta.npz
    library:
      mode: build
      start: <frozen input name (folder + manifest)> | {frozen: <name>, files: <TSV file, md5, reason: a subset>}
      species: <TSV ing_id, kind (complex|protein|...), uniprot, pdb (from the assembly stage; yeast: plugin adapter)>
      pdb_dir: <RCSB assembly cache (read only; offline)>
      afdb_tar: <frozen input name or path of the AlphaFold DB proteome tar>
      topology: {path: <UniProt topology TSV>, columns: [acc, topo, tm, intramem]}
      steps: [{step: complexes, table: <TSV ing_id, orient, membrane, tm_uniprot, reason>, per: 40}, ...]
      pack_shapes: {per: 15, ids: <optional .npy / .txt list of type ids; default every entry>,
                    max_beads: <optional: larger species keep the library tree>}
      nproc: 8
`{build}` in a path is the build directory.

Frozen mode, `--step pack_shapes` (after the import): packing shapes of the imported species (pack.shapes, §2b P8),
`library.pack_shapes: {per: 15}`.

`python -m wholecell.library.stage make-manifest <dir> <out.tsv.gz>` writes a frozen manifest for a folder.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import sys


def _md5(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _open(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def make_manifest(folder, out):
    rows = []
    for f in sorted(os.listdir(folder)):
        p = os.path.join(folder, f)
        if os.path.isfile(p):
            rows.append((f, _md5(p), os.path.getsize(p)))
    if out.endswith(".gz"):
        with gzip.GzipFile(out, "wb", mtime=0) as fh:          # mtime=0: byte-reproducible manifest
            fh.write(("file\tmd5\tbytes\n" + "".join(f"{a}\t{b}\t{c}\n" for a, b, c in rows)).encode())
    else:
        with open(out, "w") as fh:
            fh.write("file\tmd5\tbytes\n" + "".join(f"{a}\t{b}\t{c}\n" for a, b, c in rows))
    return len(rows)


def read_manifest(path):
    with _open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]


def frozen_pack_shapes(lc, out_dir, nproc=None):
    """pack_shapes/ (pack_trees_<per>.npz + pack_meta.npz) for the npz species of the imported library manifest"""
    from . import pack_shapes as ps
    man = os.path.join(out_dir, "library_manifest.tsv")
    if not os.path.exists(man):
        raise RuntimeError("library pack_shapes: no library_manifest.tsv yet (run the library stage first)")
    pc = lc.get("pack_shapes") or {}
    with open(man) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(hdr, ln.rstrip("\n").split("\t"))) for ln in f if ln.strip()]
    ent = {r["species"]: r["path"] for r in rows if r["path"].endswith(".npz")}
    dst = os.path.join(out_dir, "pack_shapes")
    if os.path.exists(os.path.join(dst, "pack_meta.npz")):
        print("pack_shapes exist, step skipped", flush=True)
        return
    ps.run(ent, dst, per=int(pc.get("per", 15)), nproc=int(nproc or lc.get("nproc", 8)), log=lambda *a: print(*a, flush=True))


def run(cfg, bdir, step=None, **kw):
    lc = cfg.get("library") or {}
    mode = lc.get("mode", "frozen")
    if mode == "build":
        from . import build
        return build.run(cfg, bdir, step=step, **kw)
    if mode != "frozen":
        raise ValueError(f"library.mode '{mode}': expected 'frozen' or 'build'")
    out_dir = os.path.join(bdir, "library")
    out = os.path.join(out_dir, "library_manifest.tsv")
    if step == "pack_shapes":                            # packing shapes of a frozen library (pack.shapes, A7 P8)
        return frozen_pack_shapes(lc, out_dir, kw.get("nproc"))
    if os.path.exists(out):
        print("library manifest exists, step skipped", flush=True)
        return
    if lc.get("frozen"):                                 # name of a frozen input {path: dir, manifest: file}
        fi = cfg["frozen_inputs"][lc["frozen"]]
        lc = dict(lc, frozen_dir=fi["path"], frozen_manifest=fi["manifest"])
    src = lc["frozen_dir"]
    man = read_manifest(lc["frozen_manifest"])
    have = sorted(f for f in os.listdir(src) if not f.startswith(".") and os.path.isfile(os.path.join(src, f)))
    listed = [r["file"] for r in man]
    extra = sorted(set(have) - set(listed)); missing = sorted(set(listed) - set(have))
    if extra or missing:
        raise RuntimeError(f"library {src}: {len(missing)} files missing (e.g. {missing[:3]}), {len(extra)} unlisted (e.g. {extra[:3]})")
    copy = bool(lc.get("copy", False))
    dst = os.path.join(out_dir, "lib") if copy else src
    if copy:
        os.makedirs(dst, exist_ok=True)
    bad = []
    rows = []
    for r in man:
        p = os.path.join(src, r["file"])
        m = _md5(p)
        if m != r["md5"]:
            bad.append(r["file"])
            continue
        if copy:
            shutil.copy2(p, os.path.join(dst, r["file"]))
        stem = r["file"].rsplit(".", 1)[0]
        rows.append((stem, os.path.join(dst, r["file"]), m))
    if bad:
        raise RuntimeError(f"library {src}: {len(bad)} files differ from the frozen manifest (e.g. {bad[:5]})")
    os.makedirs(out_dir, exist_ok=True)
    tmp = out + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write("species\tpath\tmd5\n" + "".join(f"{a}\t{b}\t{c}\n" for a, b, c in rows))
    n_npz = sum(1 for _, b, _ in rows if b.endswith(".npz"))
    json.dump(dict(mode=mode, frozen_dir=src, frozen_manifest=lc["frozen_manifest"],
                   frozen_manifest_md5=_md5(lc["frozen_manifest"]), copied=copy, files=len(rows), species_npz=n_npz),
              open(os.path.join(out_dir, "library_meta.json"), "w"), indent=1)
    os.replace(tmp, out)
    print(f"library: {n_npz:,} species npz ({len(rows):,} files) verified against {lc['frozen_manifest']}"
          f"{' and copied' if copy else ' (referenced, not copied)'}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "make-manifest":
        print(make_manifest(sys.argv[2], sys.argv[3]), "files")
    else:
        raise SystemExit(__doc__)
