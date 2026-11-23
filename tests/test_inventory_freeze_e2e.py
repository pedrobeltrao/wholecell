# new in wholecell (2026-11-20)
"""Integration: the Jurkat inventory stage of the CURRENT (freeze) config runs end to end (inventory only, ~6 min).
Regression for the freeze build failure at inventory_copies (surface 'future source' enabled). Reuses geometry /
preplaced of the accepted Step-2 build (same geometry package). Skipped when that build is absent."""
import json
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from wholecell.config import load_site  # noqa: E402
# the accepted Jurkat Step-2 build (upstream donor); env WHOLECELL_E2E_UPSTREAM overrides
UP = os.environ.get("WHOLECELL_E2E_UPSTREAM") or os.path.join(str(load_site().get("results_root") or ""), "jurkat",
                                                             "20261118-40d7967-s20261115-step2")


def test_jurkat_inventory_e2e():
    if not os.path.isdir(UP):
        print("skip: upstream build absent"); return
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "b")
        r = subprocess.run([sys.executable, os.path.join(REPO, "scripts", "inventory_e2e.py"), "--cell", "jurkat", "--out", out,
                            "--upstream-from", UP], capture_output=True, text=True, timeout=3600)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
        J = json.load(open(os.path.join(out, "inventory_e2e.json")))
        assert J["ok"] and any(s["step"] == "copies" and s["rc"] == 0 for s in J["steps"])
        import pandas as pd
        T = pd.read_csv(os.path.join(out, "inventory", "surface.tsv"), sep="\t")
        assert len(T) and T.pm_fraction.notna().sum() == 0 and T.uniprot.nunique() <= 596     # presence evidence only


if __name__ == "__main__":
    test_jurkat_inventory_e2e()
    print("ok test_jurkat_inventory_e2e")
