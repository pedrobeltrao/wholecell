"""Test of the X5-port fix for validation/vacuole_surface.py (2026-11-17), seconds on 1 CPU:

    conda_envs/yeast_x5/bin/python tests/test_yeast_vacuole_log.py   (or pytest; the env must import plugins.yeast)

vacuole_surface.py (vendored, never edited) parses "vacuole(s): N compartment(s), V um3, surface S um2" from
work/s2_geometry.log. In the port the geometry is made by the toolkit step, which logs those lines to
work/geometry_package.log; plugins/yeast/runner.append_generator_vacuole_lines appends them to s2_geometry.log after a
marker line. Checks: the vendored regex finds the model numbers in the result; other lines are not copied; a second
call appends nothing; the fallback <build>/logs/engine_geometry_package.log is used when the work copy is missing.
"""
import os
import re
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from plugins.yeast import runner  # noqa: E402

GEN_LOG = """$ python toolkits/analytic/yeast.py --engine-code ... (cwd ...)
18:39:10 cell outline: EMPIAR-11416 cell_06
18:39:16 vacuole: one real segmented vacuole (EMPIAR-11416 cell_06/vacuole_01), extent [2361. 2119. 2533.] nm
18:39:17   vacuole pass 1: scale 1.000, volume 5.836 um3 (target 5.848)
18:39:24 vacuole: 1 compartment, 5.84 um3, surface 16.2 um2, sphericity 0.97 (soft X-ray tomography of G1 cells: 0.98)
18:39:30 {"mc_label_volumes_um3": {"vacuole": 5.82}}
"""
S2_LOG = "18:40:02 MC volumes (um3): {'vacuole': 5.822}\n18:40:30 volumes (um3): {'vacuole': 5.834}\n"


def vendored_regex():
    """the pattern of plugins/yeast/engine/validation/vacuole_surface.py model() (read from the vendored source)"""
    src = open(os.path.join(REPO, "plugins", "yeast", "engine", "validation", "vacuole_surface.py")).read()
    return re.search(r're\.findall\(r"(.+?)", log\)', src).group(1)


def _fake(d, where="work"):
    work = os.path.join(d, "engine", "work")
    os.makedirs(work)
    os.makedirs(os.path.join(d, "logs"))
    gen = os.path.join(work, "geometry_package.log") if where == "work" else os.path.join(d, "logs", "engine_geometry_package.log")
    open(gen, "w").write(GEN_LOG)
    open(os.path.join(work, "s2_geometry.log"), "w").write(S2_LOG)
    return work


def test_append_and_parse():
    with tempfile.TemporaryDirectory() as d:
        work = _fake(d)
        log = open(os.path.join(work, "s2_geometry.log")).read()
        assert re.findall(vendored_regex(), log) == []                      # the port defect
        assert runner.append_generator_vacuole_lines(work, d) == 3
        log = open(os.path.join(work, "s2_geometry.log")).read()
        assert log.startswith(S2_LOG)                                       # the stage log itself is kept
        assert runner.VAC_MARK in log and "mc_label_volumes" not in log and "cell outline" not in log
        m = re.findall(vendored_regex(), log)
        assert m[-1] == ("1", "5.84", "16.2")
        assert runner.append_generator_vacuole_lines(work, d) == 0          # idempotent
        assert open(os.path.join(work, "s2_geometry.log")).read() == log


def test_fallback_build_log():
    with tempfile.TemporaryDirectory() as d:
        work = _fake(d, where="logs")
        assert runner.append_generator_vacuole_lines(work, d) == 3
        assert re.findall(vendored_regex(), open(os.path.join(work, "s2_geometry.log")).read())[-1] == ("1", "5.84", "16.2")


if __name__ == "__main__":
    test_append_and_parse()
    test_fallback_build_log()
    print("ok")
