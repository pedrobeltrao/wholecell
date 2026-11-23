// Extract decodeTile + instPose from worker.ts, transpile with esbuild, decode tiles, dump poses.
import fs from 'fs'; import zlib from 'zlib';
import { transformSync } from '/workspace/yc3d/web/app/node_modules/esbuild/lib/main.js';
const src = fs.readFileSync('/workspace/yc3d/web/app/src/worker.ts', 'utf8');
const a = src.indexOf('const TILE'), b = src.indexOf('// ------------------------------------------------------------------ spatial hash');
let code = src.slice(a, b).replace(/^const cache[\s\S]*?^async function init[\s\S]*?^}\n/m, '');
code = transformSync(code, { loader: 'ts' }).code + '\nexport { decodeTile, instPose };';
fs.writeFileSync('/tmp/dec_extract.mjs', code);
const { decodeTile, instPose } = await import('/tmp/dec_extract.mjs');
const out = {};
for (const key of process.argv.slice(2)) {
  const gz = fs.readFileSync(`/workspace/yc3d/web/data/tiles/${key}.bin.gz`);
  const raw = zlib.gunzipSync(gz); const buf = raw.buffer.slice(raw.byteOffset, raw.byteOffset + raw.byteLength);
  const D = decodeTile(buf); const P = new Float64Array(12); const rows = [];
  for (let i = 0; i < D.nI; i++) { const r = instPose(D.inst, i, D.lo, P); rows.push([r.type, ...P]); }
  out[key] = { lo: D.lo, rows, nfree: D.free.length / 4, nV: D.verts.length / 3, nT: D.tris.length / 3 };
}
fs.writeFileSync('/tmp/dec_js.json', JSON.stringify(out));
console.log('decoded', Object.keys(out).length, 'tiles');
