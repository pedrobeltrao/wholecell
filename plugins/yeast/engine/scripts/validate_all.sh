#!/bin/bash
# origin: yeast yeast_cell_dev_bundle/scripts/validate_all.sh @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
# Re-run every comparison of the model with data after a rebuild (run_all.sh), regenerate the figures and the PDF of
# model_validation_and_tuning.md. Logs: $WORK/validate_all.log. Needs the GPUs (contact pass) and, for the PDF,
# the headless-browser libraries in .chromelibs (HANDOVER).
set -e
set -o pipefail
cd "$(dirname "$0")/.."
B=$(pwd)
PY="$B/.venv/bin/python"
WORK=$($PY -c "import sys; sys.path.insert(0, 'pipeline'); import config; print(config.WORK)")
export YC_NPROC=${YC_NPROC:-$(nproc)}
export PYTHONHASHSEED=0
log() { echo "$(date +%H:%M:%S) $*"; }
{
log "contact pass (decode tiles, GPU sphere-tree gaps)"
rm -f "$WORK/contacts/instances.npz"
# GPU steps take a lock shared by all builds on this machine: two validations at once ran out of GPU memory (2026-11-15)
GPU_LOCK=${YC_GPU_LOCK:-$B/.gpu.lock}     # build_isolated.sh sets one shared path for all runs; this fallback must be a folder that exists
(cd pipeline && flock "$GPU_LOCK" $PY contacts.py)
log "pair statistics for kinases"
(cd pipeline && $PY encounters.py --kinases)
cd validation
log "localization benchmark (dev, test)"
for s in dev test; do $PY localization.py --split $s --out "$WORK/validation/latest" | grep -v Warn; done
log "proximity tests"
[ -f "$WORK/validation/proximity_sets.tsv" ] || $PY proximity_sets.py
$PY proximity_test.py | grep -A12 "gap 10"
$PY ppiseq_test.py | tail -3
$PY singh2020_test.py
log "localization of signalling proteins"
$PY kinase_localization.py | tail -1
log "ribosomes vs cryo-ET, literature numbers, copy numbers"
$PY ribosome_cryoet.py | grep -E "model_uM|nn_median"
(cd ../pipeline && $PY genome_map.py > /dev/null)
$PY chromatin_stats.py
$PY hic_test.py
$PY nucleosome_test.py > /dev/null
$PY mrna_stats.py > /dev/null
$PY mrna_occupancy.py > /dev/null || echo "mrna_occupancy failed"
$PY mrna_localization.py > /dev/null
$PY trna_nucleolus.py | grep -E 'within_overlap'
$PY dnase_test.py 2>/dev/null | grep -E 'core_over_linker|protected'
$PY silent_chromatin_test.py > /dev/null 2>&1
$PY image_nc.py 2>&1 | grep -E '"spearman"|auroc|"status"|WARNING' || echo "WARNING: image_nc failed"   # wholecell port P7: runs without the raw archive (cache), reports its absence
$PY er_ribosome_density.py | grep -E '"per_um2"'
$PY literature_checks.py
$PY copy_numbers.py > /dev/null
log "diffusion sweeps of this build, 8 cytosol boxes pooled (simgpu/sweeps_pooled.py)"
(cd ../simgpu && flock "$GPU_LOCK" $PY sweeps_pooled.py --gpus ${YC_NGPU:-4}) || echo "sweeps_pooled failed: diffusion uses the copied sweeps"   # wholecell port P2: GPU count from the environment
$PY diffusion.py > /dev/null
$PY anchors.py > /dev/null
log "placement checks (who is missing, tile seams, membrane shells and layers next to membranes)"
$PY placement_checks.py
$PY vacuole_surface.py > /dev/null || echo "vacuole_surface failed"
$PY membrane_proteomics.py > /dev/null || echo "membrane_proteomics failed"
$PY metabolism.py > /dev/null || echo "metabolism failed"
log "negative control: no condensate in the log-phase cytosol (mean-field mixture, validation/cond_test.py)"
mkdir -p "$WORK/condensates"
$PY cond_test.py > "$WORK/condensates/cond_test.log" 2>&1 || echo "cond_test failed (see $WORK/condensates/cond_test.log)"
(grep -E "^log phase|^stress" "$WORK/condensates/cond_test.log" | cut -c1-160) || true
log "simulated microscopy"
$PY microscopy_gallery.py | grep -v Warn
$PY microscopy_idr0078.py | grep -v Warn | tail -1
log "figures"
$PY make_figures.py | grep -v Warn
cd ..
log "PDF"
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$B/.chromelibs FONTCONFIG_FILE=$B/.chromelibs/fonts/fonts.conf \
  $PY scripts/md_to_pdf.py model_validation_and_tuning.md docs/model_validation_and_tuning.pdf || echo "PDF failed (non-fatal in the wholecell port, P3)"
log "result files against the build stamp"
$PY scripts/build_stamp.py check
log "done"
} 2>&1 | tee "$WORK/validate_all.log"
