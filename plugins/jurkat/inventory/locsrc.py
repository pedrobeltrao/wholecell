# origin: jurkat/validate/loc_compare.py (legacy workspace, 2026-11-15): hpa_table, HPA_MAP, UP_RULES, up_comp,
# parse_uniprot, dist, REL_W, W_UNIPROT and the hpa/uniprot part of sources()
# changes: only what inventory rule v5 (C, P) reads: per gene {"hpa": {q, w}, "uniprot": {q, w}} with LOC_VERSION = v5
# semantics ('actin filaments' adds 0.3 x weight PM support). OpenCell and GO parts of sources() are not computed (v5
# does not read them). Paths are arguments. Numerics identical to legacy.
import re

import pandas as pd

HPA_MAP = {
    "nucleoplasm": "nucleus", "nuclear bodies": "nucleus", "nuclear speckles": "nucleus", "nucleoli": "nucleus",
    "nucleoli fibrillar center": "nucleus", "nucleoli rim": "nucleus", "mitotic chromosome": "nucleus", "kinetochore": "nucleus",
    "nuclear membrane": "NE", "cytosol": "cytosol", "cytoplasmic bodies": "cytosol", "rods & rings": "cytosol",
    "aggresome": "cytosol", "actin filaments": "cytosol", "microtubules": "cytosol", "microtubule ends": "cytosol",
    "intermediate filaments": "cytosol", "centrosome": "cytosol", "centriolar satellite": "cytosol", "basal body": "cytosol",
    "primary cilium": "cytosol", "primary cilium tip": "cytosol", "primary cilium transition zone": "cytosol",
    "mitotic spindle": "cytosol", "midbody": "cytosol", "midbody ring": "cytosol", "cytokinetic bridge": "cytosol",
    "cleavage furrow": "cytosol", "peroxisomes": "cytosol", "lipid droplets": "cytosol", "mitochondria": "mito",
    "endoplasmic reticulum": "ER", "golgi apparatus": "golgi", "plasma membrane": "PM", "cell junctions": "PM",
    "vesicles": "endolysosomal", "endosomes": "endolysosomal", "lysosomes": "endolysosomal"}
UP_RULES = [
    (r"nucle(us|ar) (inner |outer )?membrane|nucle(us|ar) envelope|nuclear pore|nuclear lamina", "NE"),
    (r"mitochondri", "mito"),
    (r"endoplasmic reticulum|microsome|sarcoplasmic reticulum", "ER"),
    (r"golgi", "golgi"),
    (r"peroxisom", "cytosol"), (r"lipid droplet", "cytosol"),
    (r"stress granule|p-body", "cytosol"),
    (r"endosom|lysosom|vesic|vacuol|phagosom|autophag|melanosom|granule|multivesicular|endomembrane|secretory", "endolysosomal"),
    (r"cell membrane|plasma membrane|cell surface|membrane raft|cell junction|focal adhesion|synap|uropod|immunological|"
     r"sarcolemma|apical|basolateral|projection.*membrane|ruffle membrane|lamellipodium membrane|microvillus membrane", "PM"),
    (r"nucleol", "nucleus"),
    (r"nucle|chromosom|kinetochore|centromere|telomere|pml body|cajal|speckle", "nucleus"),
    (r"cell cortex|cortical actin", "cytosol"),
    (r"cytoplasm|cytosol|cytoskeleton|centrosome|centriol|midbody|perikaryon|myofibril|spindle|cilium|flagell|sarcomere|"
     r"cell projection|axon|dendrit|lamellipod|filopod|ruffle|microvill|growth cone|actin|microtubule|intermediate filament", "cytosol"),
    (r"secreted|extracellular|cell wall", "secretory"),
    (r"^membrane$|^membrane;|^membrane ", "anymem"),
]
UP_RULES_C = [(re.compile(a), b) for a, b in UP_RULES]
REL_W = {"Enhanced": 1.0, "Supported": 0.8, "Approved": 0.5, "Uncertain": 0.0}
W_UNIPROT = 0.6


def up_comp(head):
    h = head.strip().lower()
    for rx, c in UP_RULES_C:
        if rx.search(h): return c
    return None


def parse_uniprot(text):
    if not isinstance(text, str): return []
    t = re.sub(r"\{[^}]*\}", "", text).replace("SUBCELLULAR LOCATION:", " ")
    t = re.sub(r"Note=.*?(?=SUBCELLULAR|$)", "", t, flags=re.S)
    out = []
    for st in re.split(r"\.\s+|\.$", t):
        st = re.sub(r"^\s*\[[^\]]*\]:\s*", "", st.strip())
        if not st: continue
        head = [x.strip() for x in st.split(";")][0]
        out.append((head, up_comp(head)))
    return out


def dist(terms):
    q = {}
    for c, w in terms:
        if c: q[c] = max(q.get(c, 0.0), w)
    return q


def sources(G, hpa_path, uniprot_path):
    """G: DataFrame indexed by gene with a `uniprot` column -> {gene: {"hpa": {q, w}, "uniprot": {q, w}}}"""
    hpa = pd.read_csv(hpa_path, sep="\t").drop_duplicates("Gene name").set_index("Gene name")
    up = pd.read_csv(uniprot_path, sep="\t", usecols=["Entry", "Subcellular location [CC]"]).set_index("Entry")
    S = {}
    for g in G.index:
        d = {}
        if g in hpa.index:
            h = hpa.loc[g]; rel = h["Reliability"]
            terms = []
            for col, w in (("Main location", 1.0), ("Additional location", 0.3)):
                for t in str(h[col]).split(";") if isinstance(h[col], str) else []:
                    c = HPA_MAP.get(t.lower())
                    if c: terms.append((c, w))
                    if t.lower() == "actin filaments": terms.append(("PM", 0.3 * w))
            d["hpa"] = dict(q=dist(terms), w=REL_W.get(rel, 0.0))
        u = G.at[g, "uniprot"]
        if isinstance(u, str) and u in up.index:
            st = parse_uniprot(up.loc[u]["Subcellular location [CC]"])
            q = dist([(c, 1.0) for _, c in st])
            d["uniprot"] = dict(q=q, w=W_UNIPROT if q else 0.0)
        S[g] = d
    return S
