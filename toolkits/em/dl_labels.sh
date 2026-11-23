# origin: scratch_checks/dl_labels.sh (legacy workspace, 2026-11-15)
# ported, not yet restructured, not run: absolute paths replaced by toolkits/em/_paths (inputs inp(), outputs out() in the work folder); run from {data_root}: WHOLECELL_DATA_ROOT must be exported
#!/bin/bash
# Mirror jrc_jurkat-1 labels + groundtruth crops; 8 layers in parallel. Resumable.
# *_seg / contact / other layers: all scales (small, sparse). *_pred layers: s2-s4 only
# (a single *_pred layer is ~53 GB at s0; full-res predictions are fetched later as crops around the target cell).
cd "${WHOLECELL_DATA_ROOT:?export WHOLECELL_DATA_ROOT (site data_root)}"
AWS=conda_envs/vcell/bin/aws
S=s3://janelia-cosem-datasets/jrc_jurkat-1
D=data/em/jrc_jurkat-1
$AWS s3 ls --no-sign-request $S/jrc_jurkat-1.n5/labels/ | awk '/PRE/{print $2}' | sed 's#/##' | \
  xargs -P 8 -I{} sh -c "case {} in *_pred) EX='--exclude s0/* --exclude s1/*';; *) EX='';; esac; \
    $AWS s3 sync --no-sign-request --only-show-errors \$EX $S/jrc_jurkat-1.n5/labels/{} $D/jrc_jurkat-1.n5/labels/{} && echo done {}"
$AWS s3 sync --no-sign-request --only-show-errors $S/jrc_jurkat-1.zarr/recon-1/labels/groundtruth $D/jrc_jurkat-1.zarr/recon-1/labels/groundtruth && echo GT_DONE
echo ALL_DONE
