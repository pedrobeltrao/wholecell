"""Compare two build folders by CONTENT, ignoring what legitimately differs between two builds of the same cell.

    python scripts/ci_compare_builds.py <build A> <build B> [--out result.json]

Used by scripts/ci_testcell.sh for (1) determinism: the same commit built with 8 and 4 workers must be identical, and
(2) refactor checks: a commit that should not change outputs is compared with the previous CI build.
Exit code 0 = identical content, 1 = differences (listed in the json), 2 = usage error.

What is ignored, and why (every rule is listed here, nothing else is normalised):
  - folders/files that describe the run, not the cell: logs/, _code/, .snakemake/, _smk/, snakemake.log, launch.log,
    BUILD_STAMP.json, and the checksum lists MANIFEST.tsv / MANIFEST.sha256 (they hash the files compared here);
  - the build id, build folder and code commit wherever they appear in text (replaced by placeholders);
  - wall-clock fields: JSON keys in TIME_KEYS, TSV columns in TIME_COLUMNS, the last column of pack/pack_log.tsv
    (whose rows are also sorted: they are written in worker completion order);
  - YCT2 tiles (docs/reference/tile_format_v2.md): header word 15 = build_tag (bytes 60-63 of each tile), and in
    tiles/index.json the build_tag and the per-tile offset / gzip size / crc32 that follow from it;
  - renamed stages (A0c, "prepass" -> "preplaced"): a build made before a rename is compared with a newer one through
    the alias map wholecell.build.LEGACY_NAMES. File paths are matched in the old names (preplaced/... = prepass/...,
    pack/preplaced_pos.npy = pack/prepass_pos.npy, ...), and in text content (JSON, TSV, string arrays) each new name
    of the map is replaced by its old name before comparing (so a config key or path `preplaced` = `prepass`). The
    result lists the paths matched this way ("aliased").
"""
import argparse
import json
import os
import sys
import zlib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from wholecell.build import LEGACY_NAMES, legacy_rel  # noqa: E402

ALIAS_TEXT = sorted(LEGACY_NAMES.items(), key=lambda kv: len(kv[0]), reverse=True)   # longest first

SKIP_DIRS = {"logs", "_code", ".snakemake", "_smk"}
SKIP_FILES = {"snakemake.log", "launch.log", "BUILD_STAMP.json", "MANIFEST.tsv", "MANIFEST.sha256"}
TIME_KEYS = {"seconds", "finished", "started", "wall_s", "elapsed_s", "stamp_sha256", "created", "build_tag"}
TIME_COLUMNS = {"s", "seconds", "wall_s", "elapsed_s"}
YCT2_TAG_BYTES = slice(60, 64)          # header word 15 (u32 little endian) = first 4 bytes of sha256(build_id)
INDEX_VOLATILE_COLS = (4, 5, 12)        # tiles rows: off, gz_bytes, crc32


def build_identity(bdir):
    """strings that name this particular build: folder, id, commit (full and short)"""
    out = [os.path.abspath(bdir), os.path.basename(os.path.abspath(bdir))]
    stamp = os.path.join(bdir, "BUILD_STAMP.json")
    if os.path.exists(stamp):
        s = json.load(open(stamp))
        out += [s.get("commit", ""), s.get("commit7", "")]
    return sorted({x for x in out if x}, key=len, reverse=True)      # longest first, so the folder goes before the id


def list_files(bdir):
    files = []
    for root, dirs, names in os.walk(bdir):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS)
        for n in sorted(names):
            if n not in SKIP_FILES:
                files.append(os.path.relpath(os.path.join(root, n), bdir))
    return files


def norm_text(text, ident):
    for i, s in enumerate(ident):
        text = text.replace(s, f"<BUILD{i}>")
    for new, old in ALIAS_TEXT:                                      # renamed stages: new name -> old name
        text = text.replace(new, old)
    return text


def drop_time_keys(obj):
    if isinstance(obj, dict):
        return {k: drop_time_keys(v) for k, v in obj.items() if k not in TIME_KEYS}
    if isinstance(obj, list):
        return [drop_time_keys(v) for v in obj]
    return obj


def norm_json(path, ident):
    obj = drop_time_keys(json.loads(norm_text(open(path).read(), ident)))
    if os.path.basename(path) == "index.json" and isinstance(obj, dict) and obj.get("format") == "YCT2":
        obj["tiles"] = [[0 if i in INDEX_VOLATILE_COLS else v for i, v in enumerate(r)] for r in obj.get("tiles", [])]
    return obj


def norm_tsv(path, ident):
    rows = [ln.split("\t") for ln in norm_text(open(path).read(), ident).splitlines()]
    if path.endswith(os.path.join("pack", "pack_log.tsv")):          # no header; last column = seconds;
        return sorted(r[:-1] for r in rows)                          # rows are in worker completion order
    if not rows:
        return rows
    keep = [i for i, c in enumerate(rows[0]) if c not in TIME_COLUMNS]
    return [[r[i] for i in keep if i < len(r)] for r in rows]


def yct2_members(path):
    """decompressed gzip members of a YCT2 shard with the build_tag zeroed"""
    data, out = open(path, "rb").read(), []
    while data:
        d = zlib.decompressobj(31)
        raw = bytearray(d.decompress(data))
        if raw[:4] == b"YCT2":
            raw[YCT2_TAG_BYTES] = b"\0\0\0\0"
        out.append(bytes(raw))
        data = d.unused_data
    return out


def same_npz(pa, pb, ida, idb):
    a, b = np.load(pa, allow_pickle=False), np.load(pb, allow_pickle=False)
    if sorted(a.files) != sorted(b.files):
        return False, f"keys differ: {sorted(set(a.files) ^ set(b.files))}"
    bad = []
    for k in a.files:
        x, y = a[k], b[k]
        if x.dtype.kind == "U" and y.dtype.kind == "U":
            x = np.vectorize(lambda s: norm_text(s, ida))(x) if x.size else x
            y = np.vectorize(lambda s: norm_text(s, idb))(y) if y.size else y
        if x.shape != y.shape or not np.array_equal(x, y):
            bad.append(k)
    return (not bad), (f"arrays differ: {bad}" if bad else "")


def compare_file(rel, a, b, ida, idb, rel_b=None):
    pa, pb = os.path.join(a, rel), os.path.join(b, rel_b or rel)
    if open(pa, "rb").read() == open(pb, "rb").read():
        return True, ""
    ext = os.path.splitext(rel)[1]
    try:
        if ext == ".json":
            return norm_json(pa, ida) == norm_json(pb, idb), "json content"
        if ext in (".tsv", ".txt", ".yaml", ".csv"):
            return norm_tsv(pa, ida) == norm_tsv(pb, idb), "text content"
        if ext == ".npz":
            return same_npz(pa, pb, ida, idb)
        if ext == ".parquet":
            import pandas as pd
            return pd.read_parquet(pa).equals(pd.read_parquet(pb)), "parquet content"
        if ext == ".bin" and open(pa, "rb").read(2) == b"\x1f\x8b":
            return yct2_members(pa) == yct2_members(pb), "YCT2 tiles"
        if ext == ".gz" and rel.endswith("index.json.gz"):
            return True, "gzip copy of index.json (compared there)"
    except Exception as exc:                                         # unreadable: report, never pass silently
        return False, f"could not compare ({type(exc).__name__}: {exc})"
    return False, "bytes differ"


def compare(a, b):
    ida, idb = build_identity(a), build_identity(b)
    ma = {legacy_rel(r): r for r in list_files(a)}                  # canonical (old) name -> path in that build
    mb = {legacy_rel(r): r for r in list_files(b)}
    fa, fb = set(ma), set(mb)
    res = {"a": os.path.abspath(a), "b": os.path.abspath(b), "only_in_a": sorted(ma[k] for k in fa - fb),
           "only_in_b": sorted(mb[k] for k in fb - fa), "differ": {}, "n_compared": 0, "aliased": {}}
    for k in sorted(fa & fb):
        ra, rb = ma[k], mb[k]
        if ra != rb:
            res["aliased"][ra] = rb
        ok, why = compare_file(ra, a, b, ida, idb, rb)
        res["n_compared"] += 1
        if not ok:
            res["differ"][ra if ra == rb else f"{ra} ~ {rb}"] = why
    res["identical"] = not (res["only_in_a"] or res["only_in_b"] or res["differ"])
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("a"); ap.add_argument("b"); ap.add_argument("--out")
    args = ap.parse_args(argv)
    for d in (args.a, args.b):
        if not os.path.isdir(d):
            print(f"not a build folder: {d}", file=sys.stderr)
            return 2
    res = compare(args.a, args.b)
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(f"{'IDENTICAL' if res['identical'] else 'DIFFERENT'}: {res['n_compared']} files compared, "
          f"{len(res['differ'])} differ, {len(res['only_in_a'])} only in A, {len(res['only_in_b'])} only in B"
          + (f"; {len(res['aliased'])} matched through renamed-stage aliases (wholecell.build.LEGACY_NAMES)"
             if res["aliased"] else ""))
    for rel, why in list(res["differ"].items())[:20]:
        print(f"  differ: {rel} ({why})")
    for rel in (res["only_in_a"] + res["only_in_b"])[:20]:
        print(f"  missing on one side: {rel}")
    return 0 if res["identical"] else 1


if __name__ == "__main__":
    sys.exit(main())
