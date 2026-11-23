"""Stage the Hugging Face inputs dataset (local only; nothing is uploaded).

    python scripts/deposit/stage_hf.py --out <dir> [--repo-id pedrobeltrao/wholecell-inputs]

Reads data_manifest.tsv (one row per input), the cell configs (frozen inputs, md5) and the yeast read set
(plugins/yeast/data/pipeline_data_read_set.txt + scripts/deposit/yeast_read_set_terms.tsv) with the roots of the site
file (WHOLECELL_SITE), and writes into <dir>:
  stage/<licence folder>/<root>/<relpath>   symlinks to the input files (no copy)
  MANIFEST.tsv     every dataset file: hf path, kind (file | tar shard), bytes, md5, licence, target root and path, group
  TARS.tsv         members of each tar shard (directories with > 1,000 files are shipped as tar shards of <= 2 GB)
  EXCLUDED.tsv     inputs that are not deposited, with the reason and the route to get them
  SIZES.tsv        bytes per licence folder and group
The dataset card, upload and fetch scripts are copied next to it. The upload script builds the tar shards and uploads;
the fetch script downloads into the site roots (docs/REBUILD.md).
"""
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import os
import shutil
import sys

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
from wholecell.config import load, load_site  # noqa: E402

TAR_MIN_FILES = 1000
SHARD_BYTES = 2_000_000_000
FOLDERS = {   # licence folder -> licence of its contents
    "ours": "CC BY-NC 4.0 (this work: Beltrao lab and G. Magnusson; derived from the cited sources)",
    "cc_by": "CC BY 4.0 (third party: UniProt, AlphaFold DB, SIFTS, STRING, ...; attribute the source)",
    "cc0": "CC0 1.0 (third party: wwPDB, Complex Portal, PRIDE, GEO, EMPIAR)",
    "open": "source terms (EMBL-EBI Expression Atlas, NCBI GEO: public, no restriction stated)",
    "cc_by_sa": "CC BY-SA 4.0 (Human Protein Atlas tables and tables derived from them)",
    "af3_pool": "AlphaFold 3 output terms of use + the usage request of https://github.com/jurgjn/pooled-ppi-yeast",
}
# manifest name -> (route, licence folder, group). route: hf | hf_tar | source | excluded | builds | generated
# group: required (read by jurkat.yaml / yeast.yaml) | reference (reference presets, v7 / X5 reproduction) | optional
OURS = "ours"
R = {}
for n in ["inv_chromatin_v8_json;pre_chromatin_v8_rep5_json", "inv_cytoskeleton_v6_json", "inv_cytoskeleton_v7_json;pre_cytoskeleton_v7_json",
          "inv_fit_genes;val_loc_fit_genes", "inv_holdout_genes;val_loc_holdout_genes", "pre_lamina_v7_json", "inv_legacy_complexes_v0",
          "inv_legacy_ingredients_v7", "inv_membrane_lib_ids", "inv_mt_centerlines", "inv_sec_lab2026_jurkat_nt",
          "inv_sec_lab2026_jurkat_tcract24h", "inv_surface_jurkat_glyco2026", "inv_tier2_transfer", "inv_tm_route_libs",
          "inv_v4pack_loc_compare", "inv_v5pack_loc_compare", "lib_afdb_fragment_index_tmroute", "lib_afdb_human_tar_index",
          "pre_chromatin_drop_lamina_v7", "pre_chromatin_drop_lamina_v7_npy", "pre_chromatin_v8_rep5;val_chromatin_rep5",
          "pre_cytoskeleton_v7", "pre_lamina_v7", "pre_mrna_nascent_v8", "pre_mrna_v1", "val_chromatin_genomic", "val_chromatin_s1_em",
          "val_chromatin_s1_lit", "val_chromatin_v8_metrics", "val_chromatin_v8_metrics_table", "val_ribo_em", "geometry_jurkat",
          "lib_afdb_plddt_summary", "inv_sec_lab2026_yeast_nacl", "inv_sec_lab2026_yeast_ypd", "inv_sec_pang2020_yeast",
          "inv_sec_swamy2022_yeast", "yeast_x5_lib_start", "yeast_x5_work_seed", "geometry_yeast"]:
    R[n] = ("hf", OURS, "required")
R["inv_struct_library;lib_structures_v7"] = ("hf_tar", OURS, "required")
for n in ["inv_afdb_fetch_status", "inv_afdb_heterodimer_table", "inv_afdb_homodimer_table", "inv_itzhak2016", "inv_sifts",
          "inv_uniprot_domains", "inv_uniprot_full;val_uniprot_human", "inv_uniprot_reviewed", "inv_uniprot_topology",
          "lib_tm_route_uniprot", "lib_uniprot_human_table", "lib_uniprot_glyco_topo", "inv_string_links", "val_opencell"]:
    R[n] = ("hf", "cc_by", "required")
R["lib_afdb_complex_models"] = ("hf_tar", "cc_by", "required")
R["yeast_afdb_dimers"] = ("hf_tar", "cc_by", "required")
R["inv_complex_portal"] = ("hf", "cc0", "required")
R["inv_pdb_assemblies;lib_pdb_assemblies"] = ("hf_tar", "cc0", "required")
R["inv_npc_7r5k"] = ("within", "cc0", "required")                 # a file of the PDB assembly set
R["yeast_pdb_cache"] = ("hf", "cc0", "required")
for n in ["inv_sec_fossati_features", "inv_sec_fossati_fractions"]:
    R[n] = ("hf", "cc0", "required")
for n in ["inv_geiger2012_pg", "inv_rna_emtab2706", "inv_rna_gse145453_a1", "inv_rna_gse145453_a2", "inv_rna_gse276689_rep1",
          "inv_rna_gse276689_rep2", "inv_rna_gse301789"]:
    R[n] = ("hf", "open", "required")
for n in ["inv_hpa_celllines;val_hpa_celllines", "inv_hpa_jurkat_manifest", "inv_hpa_subcellular;val_hpa_subcellular",
          "inv_rna_hpa", "inv_tier2_v5", "val_chromatin_s1_dapi"]:
    R[n] = ("hf", "cc_by_sa", "required")
for n in ["lib_afdb_human_tar", "yeast_afdb_tar"]:
    R[n] = ("source", "", "required")
R["yeast_pipeline_data"] = ("read_set", "", "required")
for n in ["inv_yeast_complex_portal", "inv_yeast_lanz2024_s3", "inv_yeast_paxdb", "inv_yeast_uniprot", "loc_yeast_cyclops_wt1",
          "loc_yeast_cyclops_wt2", "loc_yeast_cyclops_wt3", "loc_yeast_morgenstern2017", "loc_yeast_uniprot_features",
          "loc_yeast_vogtle2017", "hidden_yeast_string_links", "hidden_yeast_uniprot_topology"]:
    R[n] = ("within_read_set", "", "required")
for n in ["geom_jgeom_v7", "pre_npc_v3", "pre_npc_v3_json", "pre_polysomes_v7", "pre_polysomes_v7_json", "pre_mrna_drop_v7",
          "inv_sec_assembly_thresholds", "val_yeast_buildW_anchors", "val_yeast_buildW_protein_features", "val_yeast_x5_build_json",
          "val_yeast_x5_build_log"]:
    R[n] = ("hf", OURS, "reference")
R["yeast_af3_pool_summary"] = ("excluded", "", "reference")
R["yeast_af3_pool_pairs"] = ("within_work_seed", "af3_pool", "reference")   # staged with yeast_x5_work_seed, in af3_pool/
R["hidden_yeast_inputs_dir"] = ("partial", "", "required")
R["hidden_legacy_pipeline_dir"] = ("hf_tar", OURS, "reference")
for n in ["hidden_v7_reference_tables", "ref_jurkat_v7_reference_validate", "ref_jurkat_step2_build", "ref_yeast_x5repro_build",
          "donor_yeast_pack_bedcea9"]:
    R[n] = ("builds", "", "optional")
R["geometry_testcell"] = ("generated", "", "optional")
R["env_yeast_x5_python"] = ("generated", "", "required")

EXCLUDE_WHY = {
    "yeast_af3_pool_summary": ("in-house AF3 pool predictions (444 MB, includes predictions not in the public release); read only "
                               "by yeast_x5_port.yaml (whole-engine X5 reproduction)",
                               "not available; yeast_x5_core.yaml uses the deposited af3_pool/ pool_pairs.parquet instead"),
    "lib_afdb_human_tar": ("AlphaFold DB proteome tar (5.2 GB), served by EBI", "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000005640_9606_HUMAN_v6.tar"),
    "yeast_afdb_tar": ("AlphaFold DB proteome tar (1.0 GB), served by EBI", "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000002311_559292_YEAST_v6.tar"),
}


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def read_yeast_manifest(path):
    out = []
    for line in open(path):
        parts = [x for x in line.rstrip("\n").split("  ") if x != ""]
        if not parts:
            continue
        if len(parts[0]) == 32:
            md5, size, _, rel = parts[:4]
        else:
            size, _, rel = parts[:3]; md5 = None
        out.append((rel[2:] if rel.startswith("./") else rel, int(size), md5))
    return out


def read_repo_manifest(path):
    T = pd.read_csv(path, sep="\t", dtype=str)
    return [(r.file, int(r.bytes), r.md5) for r in T.itertuples()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--repo-id", default="pedrobeltrao/wholecell-inputs")
    ap.add_argument("--no-md5", action="store_true", help="do not compute missing md5s (faster)")
    a = ap.parse_args()
    site = load_site()
    roots = {k: site.get(k) for k in ("data_root", "yeast_root", "afdb_root", "geom_root", "results_root")}
    M = pd.read_csv(os.path.join(REPO, "data_manifest.tsv"), sep="\t", dtype=str, keep_default_na=False)
    fi = {}
    for c in ("jurkat", "yeast", "jurkat_v7_port", "yeast_x5_core"):
        fi.update(load(c)["frozen_inputs"])
    unknown = [n for n in M.name if n not in R]
    if unknown:
        raise SystemExit(f"no route for manifest rows {unknown}: add them to stage_hf.R")
    stage = os.path.join(a.out, "stage")
    if os.path.exists(stage):
        shutil.rmtree(stage)
    files, tars, excl = [], [], []
    cur = {"cell": ""}

    def add_file(folder, root, rel, src, size, md5, licence, group, name):
        hf = f"{folder}/{root}/{rel}"
        dst = os.path.join(stage, hf)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        os.symlink(os.path.realpath(src), dst)
        if md5 is None and not a.no_md5:
            md5 = md5_file(src)
        files.append(dict(hf_path=hf, kind="file", bytes=size, md5=md5 or "", licence=FOLDERS[folder], target_root=root,
                          target_path=rel, group=group, cell=cur["cell"], input=name))

    def add_dir_tar(folder, root, rel, srcdir, members, group, name):
        """members: [(relpath in dir, bytes, md5)] -> tar shards (built at upload time)"""
        shard, size, k = [], 0, 0
        def flush():
            nonlocal shard, size, k
            if shard:
                hf = f"{folder}/{root}/{rel}.part{k:03d}.tar"
                files.append(dict(hf_path=hf, kind="tar", bytes=size, md5="", licence=FOLDERS[folder], target_root=root,
                                  target_path=rel, group=group, cell=cur["cell"], input=name))
                tars.extend(dict(tar=hf, source=os.path.join(srcdir, m), member=m, bytes=b, md5=h or "") for m, b, h in shard)
                k += 1; shard, size = [], 0
        for m in sorted(members):
            if size + m[1] > SHARD_BYTES:
                flush()
            shard.append(m); size += m[1]
        flush()

    def members_of(name, path):
        man = next((fi[n]["manifest"] for n in name.split(";") if (fi.get(n) or {}).get("manifest")), None)
        if man and man.endswith(".tsv.gz"):
            return read_repo_manifest(man)
        if man:
            return read_yeast_manifest(man)
        import subprocess                  # one find call: much faster than os.walk + stat on a parallel file system
        txt = subprocess.run(["find", "-L", path, "-type", "f", "-printf", "%P\t%s\n"], capture_output=True, text=True,
                             check=True).stdout
        return [(m, int(b), None) for m, b in (ln.split("\t") for ln in txt.splitlines() if ln)]

    for r in M.itertuples():
        route, folder, group = R[r.name]
        cur["cell"] = r.cell
        root = r.root
        src = os.path.join(str(roots.get(root) or f"<{root}>"), r.relpath)
        if route == "hf" and r.kind == "file":
            md5 = r.checksum if r.checksum_algo == "md5" and len(r.checksum) == 32 else None
            add_file(folder, root, r.relpath, src, os.path.getsize(src), md5, None, group, r.name)
        elif route == "hf" and r.kind == "dir":
            mem = members_of(r.name, src)
            for m, b, h in mem:
                if r.name == "yeast_x5_work_seed" and m == "organisation/pool_pairs.parquet":
                    add_file("af3_pool", root, f"{r.relpath}/{m}", os.path.join(src, m), b, h, None, "reference", "yeast_af3_pool_pairs")
                    continue
                add_file(folder, root, f"{r.relpath}/{m}", os.path.join(src, m), b, h, None, group, r.name)
            mf = next((fi[n]["manifest"] for n in r.name.split(";") if (fi.get(n) or {}).get("manifest")), None)
            if mf and mf.startswith(str(roots["data_root"])) and not mf.endswith(".tsv.gz"):    # the manifest file itself
                add_file(OURS, "data_root", os.path.relpath(mf, roots["data_root"]), mf, os.path.getsize(mf), None, None, group, r.name + " (manifest)")
        elif route == "hf_tar":
            mem = members_of(r.name, src)
            if r.name == "hidden_legacy_pipeline_dir":       # code + data of the legacy tree, without lib/ (deposited separately)
                mem = [m for m in mem if not m[0].startswith(("lib/", "__pycache__/")) and "/__pycache__/" not in m[0]]
            add_dir_tar(folder, root, r.relpath, src, mem, group, r.name)
            mf = next((fi[n]["manifest"] for n in r.name.split(";") if (fi.get(n) or {}).get("manifest")), None)
            if mf and mf.startswith(str(roots["data_root"])) and not mf.endswith(".tsv.gz"):
                add_file(OURS, "data_root", os.path.relpath(mf, roots["data_root"]), mf, os.path.getsize(mf), None, None, group, r.name + " (manifest)")
        elif route == "read_set":
            pats = [ln.strip() for ln in open(os.path.join(REPO, "plugins/yeast/data/pipeline_data_read_set.txt"))
                    if ln.strip() and not ln.startswith("#")]
            T = pd.read_csv(os.path.join(REPO, "scripts/deposit/yeast_read_set_terms.tsv"), sep="\t", comment="#", dtype=str,
                            keep_default_na=False)
            man = (fi.get("yeast_pipeline_data") or {}).get("manifest")
            for m, b, h in read_yeast_manifest(man):
                if not any(fnmatch.fnmatch(m, p) for p in pats):
                    continue
                t = next((t for t in T.itertuples() if fnmatch.fnmatch(m, t.pattern)), None)
                if t is None:
                    raise SystemExit(f"no terms row for read-set file {m}")
                if t.route == "hf":
                    add_file(t.licence_folder, root, f"{r.relpath}/{m}", os.path.join(src, m), b, h, None, group, r.name)
                else:
                    excl.append(dict(input=r.name, cell=r.cell, target_root=root, target_path=f"{r.relpath}/{m}", bytes=b, md5=h or "",
                                     why=t.licence, route=f"fetch from source: {t.source}"))
            mf = man
            add_file(OURS, "data_root", os.path.relpath(mf, roots["data_root"]), mf, os.path.getsize(mf), None, None, group, r.name + " (manifest)")
        elif route == "partial":                            # yeast inputs/: only the files the build reads
            for f_, folder_ in (("yeast_organelle_meshes.zip", "cc0"),):
                p = os.path.join(src, f_)
                add_file(folder_, root, f"{r.relpath}/{f_}", p, os.path.getsize(p), None, None, "optional", r.name)
            excl.append(dict(input=r.name, cell=r.cell, target_root=root, target_path=f"{r.relpath}/IntActYeast.txt", bytes="", md5="",
                             why="not read by the build", route="not needed"))
        elif route in ("source", "excluded"):
            why, how = EXCLUDE_WHY[r.name]
            excl.append(dict(input=r.name, cell=r.cell, target_root=root, target_path=r.relpath, bytes=r.size_bytes, md5=r.checksum,
                             why=why, route=how))
        elif route == "builds":
            excl.append(dict(input=r.name, cell=r.cell, target_root=root, target_path=r.relpath, bytes=r.size_bytes, md5=r.checksum,
                             why="comparison reference (optional; the build runs without it)",
                             route="builds dataset (https://huggingface.co/datasets/TODO/wholecell-builds)"))
        elif route == "generated":
            excl.append(dict(input=r.name, cell=r.cell, target_root=root, target_path=r.relpath, bytes="", md5="",
                             why="made locally", route="testcell geometry: toolkits/analytic/generate.py; env: envs/yeast_x5.pip.txt"))
        # within / within_read_set / within_work_seed: covered by the PDB tar, the read set or the work seed
    F = pd.DataFrame(files)
    dup = F[F.hf_path.duplicated()]
    F = F.drop_duplicates("hf_path")
    os.makedirs(a.out, exist_ok=True)
    F.to_csv(os.path.join(a.out, "MANIFEST.tsv"), sep="\t", index=False)
    pd.DataFrame(tars).to_csv(os.path.join(a.out, "TARS.tsv"), sep="\t", index=False)
    pd.DataFrame(excl).to_csv(os.path.join(a.out, "EXCLUDED.tsv"), sep="\t", index=False)
    F["folder"] = F.hf_path.str.split("/").str[0]
    S = F.groupby(["folder", "group"]).agg(files=("hf_path", "size"), bytes=("bytes", "sum")).reset_index()
    S.to_csv(os.path.join(a.out, "SIZES.tsv"), sep="\t", index=False)
    for f, src in (("upload_hf.py", "scripts/deposit/upload_hf.py"), ("fetch_inputs.py", "scripts/fetch_inputs.py"),
                   ("README.md", "scripts/deposit/dataset_card.md")):
        txt = open(os.path.join(REPO, src)).read().replace("pedrobeltrao/wholecell-inputs", a.repo_id)
        if f == "README.md":
            txt += "\n## Sizes (generated by stage_hf.py)\n\n| folder | group | files | GB |\n|---|---|---|---|\n" + "".join(
                f"| {x.folder} | {x.group} | {x.files:,} | {x.bytes / 1e9:.2f} |\n" for x in S.itertuples())
        open(os.path.join(a.out, f), "w").write(txt)
    print(S.to_string(index=False))
    print(f"{len(F)} dataset files ({F.bytes.sum() / 1e9:.2f} GB; {(F.kind == 'tar').sum()} tar shards with {len(tars)} members); "
          f"{len(excl)} excluded rows; {len(dup)} duplicate paths merged")


if __name__ == "__main__":
    main()
