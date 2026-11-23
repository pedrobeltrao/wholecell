# origin: yeast yeast_cell_dev_bundle/web/build_preview.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Build a single-file preview (yeast_column_preview.html): a contiguous block of tiles
(cell wall -> PM -> cortical ER -> cytoplasm -> nuclear envelope with an NPC -> nucleoplasm),
with the real viewer code and a structure library restricted to the types present.
All data are embedded as base64 in window.__YC_EMBED__ (same type indices as the full app)."""
import json, gzip, base64, os, sys, re, numpy as np
sys.path.insert(0, '/workspace/yc')
from webenc import decode_tile

DATA = '/workspace/yc3d/web/data'
SINGLE = '/workspace/yc3d/web/app/dist_single/index.html'
OUT = sys.argv[1] if len(sys.argv) > 1 else '/workspace/yc3d/web/yeast_column_preview.html'
BLOCK = [(ix, iy, iz) for ix in range(-11, -7) for iy in (0, 1) for iz in (0, 1)]
NPC = [-1856.0, 250.0, 122.0]

I = json.load(open(f'{DATA}/tiles/index.json'))
cols = I['tile_columns']; ci = {c: i for i, c in enumerate(cols)}
rows = [r for r in I['tiles'] if (r[ci['ix']], r[ci['iy']], r[ci['iz']]) in set(BLOCK)]
print('tiles in block', len(rows))
I2 = dict(I); I2['tiles'] = rows
emb = {}
present = set()
tiles_of_type = {}
for k, r in enumerate(rows):
    key = f"{r[ci['ix']]}_{r[ci['iy']]}_{r[ci['iz']]}"
    blob_gz = open(f'{DATA}/tiles/{key}.bin.gz', 'rb').read()
    emb[f'tiles/{key}.bin.gz'] = blob_gz
    d = decode_tile(gzip.decompress(blob_gz))
    ts = np.unique(d['inst'][0])
    if d['ghost'] is not None: ts = np.union1d(ts, np.unique(d['ghost'][0]))
    for t in d['inst'][0].tolist(): pass
    for t in np.unique(d['inst'][0]).tolist(): tiles_of_type.setdefault(t, []).append(k)
    present |= set(ts.tolist())
print('types present', len(present))

# library subset: keep all type records (same indices), but only beads of present types
tj = json.load(open(f'{DATA}/library/types.json'))
beads = np.fromfile(f'{DATA}/library/beads.bin', np.uint8).reshape(-1, 8)
parts, bo = [], 0
for i, t in enumerate(tj['types']):
    if i in present:
        parts.append(beads[t['bo']:t['bo'] + t['nb']]); t['bo'] = bo; bo += t['nb']
    else:
        t['bo'] = bo; t['nb'] = 0; t['nnt'] = 0
tj['n_beads'] = bo; tj['preview_subset'] = True
sub = np.concatenate(parts) if parts else np.zeros((0, 8), np.uint8)
print('beads kept', bo, 'of', len(beads))
emb['library/types.json'] = gzip.compress(json.dumps(tj, separators=(',', ':')).encode(), 9)
emb['library/beads.bin'] = gzip.compress(sub.tobytes(), 9)
emb['library/spheretree.bin'] = gzip.compress(open(f'{DATA}/library/spheretree.bin', 'rb').read(), 9)
emb['tiles/index.json'] = gzip.compress(json.dumps(I2, separators=(',', ':')).encode(), 9)
# types -> tiles (indices into the subset tile list)
NT = len(tj['types'])
off = np.zeros(NT + 1, np.uint32); lst = []
for t in range(NT):
    v = tiles_of_type.get(t, []); lst += v; off[t + 1] = off[t] + len(v)
emb['types_tiles.bin'] = gzip.compress(off.tobytes() + np.array(lst, np.uint16).tobytes(), 9)
emb['geometry.json'] = gzip.compress(open(f'{DATA}/geometry.json', 'rb').read(), 9)

# preview bookmarks centred on the block
n_z = [0, 0, 1]
bm = [
    dict(name='Overview', caption='A 800 x 400 x 400 nm block of the cell, cut through a nuclear pore complex. Left to right: cell wall, plasma membrane, cortical ER, a thin layer of cytoplasm, the nuclear envelope with an NPC, and nucleoplasm with chromatin.',
         target=[-1800, 200, 122], normal=n_z, up=[0, 1, 0], width=900, slab=60),
    dict(name='Cell wall', caption='Mannoproteins and procedural mannan brush (outer 35 nm) over the beta-glucan layer; plasma membrane below.',
         target=[-2080, 200, 122], normal=n_z, up=[1, 0, 0], width=260, slab=40),
    dict(name='Plasma membrane + cortical ER', caption='Plasma membrane with Pma1 and other transmembrane proteins; cortical ER sheet ~25 nm beneath it.',
         target=[-2010, 200, 122], normal=n_z, up=[1, 0, 0], width=240, slab=30),
    dict(name='Nuclear pore complex', caption='Nuclear pore complex (PDB 7N9F) spanning the double nuclear envelope.',
         target=NPC, normal=n_z, up=[1, 0, 0], width=260, slab=50),
    dict(name='Chromatin', caption='Nucleosome chains with linker DNA, RNA polymerases and nuclear proteins.',
         target=[-1560, 200, 122], normal=n_z, up=[1, 0, 0], width=220, slab=40),
    dict(name='Cytoplasm (thin layer)', caption='Ribosomes, tRNAs, glycolytic enzymes and mRNA between the cortical ER and the nucleus.',
         target=[-1935, 120, 122], normal=n_z, up=[1, 0, 0], width=160, slab=30),
]
emb['bookmarks.json'] = json.dumps(bm).encode()

html = open(SINGLE).read()
js = 'window.__YC_EMBED__={' + ','.join(f'"{k}":"{base64.b64encode(v).decode()}"' for k, v in emb.items()) + '};window.__YC_PREVIEW__=true;'
tag = '<script>' + js + '</script>'
html = html.replace('<head>', '<head>' + tag, 1)
open(OUT, 'w').write(html)
print('wrote', OUT, f'{os.path.getsize(OUT)/1e6:.1f} MB', {k: len(v) for k, v in emb.items() if not k.startswith('tiles/')})
