#!/bin/bash
# origin: yeast yeast_cell_dev_bundle/web/package.sh @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
# Assemble yeast_cell_viewer/ (built app + data + source + README) and zip it.
set -euo pipefail
W=/workspace/yc3d/web
S=$W/pkg/yeast_cell_viewer
rm -rf $W/pkg && mkdir -p $S/source/src $S/tests
cp -r $W/app/dist/index.html $W/app/dist/assets $W/app/dist/molstar $S/
cp -rL $W/data $S/data
cp $W/app/src/*.ts $S/source/src/
cp $W/app/index.html $W/app/package.json $W/app/package-lock.json $W/app/vite.config.ts $W/app/vite.single.config.ts $W/app/tsconfig.json $S/source/
cp $W/test_viewer.py $W/tests/dec_js.mjs $W/tests/cmp_raw.py $W/tests/decoder_roundtrip.json $S/tests/
cp /workspace/yc/webenc.py $S/tests/webenc.py
cp $W/README.md $S/README.md
cat > $S/serve.sh <<'EOF'
#!/bin/sh
# Serve the viewer at http://localhost:8000
cd "$(dirname "$0")" && python3 -m http.server "${1:-8000}"
EOF
cat > $S/serve.bat <<'EOF'
@echo off
cd /d "%~dp0"
python -m http.server 8000
EOF
chmod +x $S/serve.sh
cd $W/pkg
# gz tiles are already compressed: store them (-n suffixes)
zip -q -r -n .gz:.png:.ico yeast_cell_viewer.zip yeast_cell_viewer
ls -la yeast_cell_viewer.zip
unzip -l yeast_cell_viewer.zip | tail -1
