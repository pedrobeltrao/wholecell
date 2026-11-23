# origin: yeast yeast_cell_dev_bundle/analysis/rcsb_coverage.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Best S. cerevisiae RCSB entry per unmodeled Complex Portal complex (subunit-type mass coverage)."""
import sys, re, json, requests, collections
import pandas as pd

S = sys.argv[1]
import os  # PORT PATCH 5 (release 2026-11-23): author-machine paths -> env (see PORT_PATCHES.md)
B = os.environ.get("YC_BUNDLE", "yeast_cell_dev_bundle")
sys.path.insert(0, f"{B}/pipeline")
import inventory
inventory.D = f"{B}/pipeline/data"
df, _ = inventory.build()
acc_of = dict(zip(df.gene, df.uniprot)); mass = dict(zip(df.gene, df.mass_da)); gene_of = {v: k for k, v in acc_of.items()}
cp = pd.read_csv(f"{B}/pipeline/data/complexportal_559292.tsv", sep="\t")
subs = {}
for r in cp.itertuples(index=False):
    sub = {gene_of[a]: max(1, int(s)) for a, s in re.findall(r"([A-Z0-9]{6,10})(?:-\w+)?\((\d+)\)", r[4]) if a in gene_of}
    subs[r[0]] = sub
C = pd.read_csv(f"{S}/complex_coverage.csv")
um = C[C.modeled == 0]
sel = pd.concat([um.head(70), um[um.tm].head(30)]).drop_duplicates("cpx")


def search(accs):
    nodes = [{"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": a}} for a in accs]
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes}, "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": 500}}}
    r = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=q, timeout=60)
    return [x["identifier"] for x in r.json().get("result_set", [])] if r.status_code == 200 else []


def entries(ids):
    q = '{ entries(entry_ids: %s) { rcsb_id rcsb_entry_info { resolution_combined experimental_method deposited_polymer_monomer_count } struct { title } polymer_entities { rcsb_polymer_entity_container_identifiers { uniprot_ids } } } }' % json.dumps(ids)
    return requests.post("https://data.rcsb.org/graphql", json={"query": q}, timeout=120).json()["data"]["entries"]


rows = []
for c in sel.itertuples():
    sub = subs[c.cpx]
    hits = collections.Counter()
    for g_ in sub:
        hits.update(search([acc_of[g_]]))
    ids = [i for i, _ in hits.most_common(40)]
    best = None
    for i in range(0, len(ids), 20):
        for e in entries(ids[i:i + 20]) or []:
            acc = set(a for pe in e["polymer_entities"] for a in (pe["rcsb_polymer_entity_container_identifiers"]["uniprot_ids"] or []))
            cov = sum(mass[g] * sub[g] for g in sub if acc_of[g] in acc) / sum(mass[g] * sub[g] for g in sub)
            res = (e["rcsb_entry_info"]["resolution_combined"] or [99])[0]
            key = (round(cov, 2), -res)
            if best is None or key > best[0]:
                best = (key, e["rcsb_id"], cov, res, e["rcsb_entry_info"]["experimental_method"], e["struct"]["title"][:60])
    rows.append(dict(cpx=c.cpx, name=c.name, tm=c.tm, nsub=c.nsub, copies=round(c.copies), mass_kda=round(c.mass_kda),
                     n_entries=len(ids), best=best[1] if best else "", cov=round(best[2], 2) if best else 0,
                     res=best[3] if best else None, method=best[4] if best else "", title=best[5] if best else ""))
    print(rows[-1], flush=True)
R = pd.DataFrame(rows)
R.to_csv(f"{S}/rcsb_coverage.csv", index=False)
print(R[["name", "tm", "copies", "mass_kda", "best", "cov", "res", "method"]].to_string(index=False))
