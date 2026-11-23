# Whole yeast cell: interactive 3D viewer

This is a GPU-accelerated browser viewer for a structural model of one complete haploid *Saccharomyces cerevisiae* cell (BY4741, YPD, log phase). The cell contains 25.2 million macromolecules of 5,483 species. Each molecule is drawn from its actual PDB or AlphaFold structure at one bead per residue or nucleotide.

## Run it

```bash
unzip yeast_cell_viewer.zip
cd yeast_cell_viewer
python3 -m http.server 8000
# open http://localhost:8000 in Chrome or Edge
```

You need a local web server because browsers block `fetch()` from `file://` pages. Any static server works, for example `npx serve` or nginx. The folder can also be copied as-is to any static host, such as GitHub Pages, S3 or Netlify.

### Browser requirements (WebGPU)
| Browser | Status |
|---|---|
| Chrome / Edge 113+ (Windows, macOS, ChromeOS, Android) | Works out of the box |
| Chrome on Linux | Start with `--enable-unsafe-webgpu --enable-features=Vulkan`, or enable `chrome://flags/#enable-unsafe-webgpu` |
| Safari 26+ (macOS / iOS) | Supported |
| Firefox 141+ (Windows); other platforms behind `dom.webgpu.enabled` | Experimental |

The viewer was developed against a desktop discrete or integrated GPU. A GPU with 2 GB or more of memory is recommended. The app allocates about 0.6–1.2 GB of GPU buffers, depending on the adapter limits.

The full data set is about 530 MB. It is streamed tile by tile, so only the region you look at is downloaded. On localhost this is instant.

## Controls
| Action | Input |
|---|---|
| Orbit | Left-drag |
| Pan | Right-drag, or Shift + left-drag |
| Zoom toward cursor | Mouse wheel |
| Move the cut plane | Shift + wheel, or the *Position* slider |
| Fly to a point | Double-click |
| Inspect a molecule | Click it; the panel shows gene, UniProt, PDB/AF source, pLDDT, copies per cell |
| Full-atom structure | *Open atomic structure* in the inspect panel (Mol*, fetched from RCSB / AlphaFold DB) |
| Fly | W A S D, Q / E (down / up), Z / X (roll) |
| Guided views | Keys 1–9, 0, or the buttons; T starts or stops the tour |
| Cut plane | C toggles the cut; V aligns the plane to the view; `,` / `.` step the plane; `[` / `]` change the slab depth |
| Organelle surfaces | O |
| Find a molecule | Type a gene, protein name or UniProt ID in the search box; Enter flies to the nearest copy and highlights all copies |
| Deselect | Esc |
| Screenshot (2x) | P |

**Display options**
- Color modes: Goodsell compartment palette, structure source (PDB vs AlphaFold), AlphaFold pLDDT, abundance (copies per cell), random per protein type, and compartment.
- Molecules / Organelles mode (organelle mode shows smooth membrane surfaces).
- Toggles for cell wall, lipids, polymers and individual compartments.
- Quality presets. *Auto* adapts the level of detail to keep the frame rate up.

## What you are looking at
- **Proteins and complexes.** 5,483 structure types: assembled PDB complexes (ribosome, nucleosome, NPC, RNA polymerases, FAS, 26S proteasome, CCT, V-ATPase, …), PDB homo-oligomers and AlphaFold monomers. AlphaFold residues with pLDDT < 50 are trimmed; the full model loads in Mol*. Copy numbers come from the PaxDb integrated proteome scaled to 4.2 × 10⁷ proteins per cell. Compartments come from UniProt localization, with GO as fallback.
- **Packing.** Sphere-tree random sequential addition in 200 nm tiles, packed in checkerboard phases with halos. Membrane proteins are oriented along the membrane normal.
  - Residue-level clash rate: 4.6 × 10⁻⁶ of sampled beads.
  - Crowding: cytosol 88.7 mg/mL, nucleoplasm 63.2, mitochondrial matrix 109.9.
- **Lipids, β-glucan and mannan.** These are generated procedurally near the camera. Generation is deterministic, so the same place always looks the same.
  - Lipids: individual molecules at 0.65 nm² each, in the measured class composition (Ejsing et al. 2009). They are removed wherever protein residues occupy the bilayer.
  - β-glucan and mannan: densities derived from the measured wall composition.
- **Chromatin, mRNA, linker DNA.** Chromatin is 73,152 nucleosomes on persistent walks. mRNA and linker DNA are schematic bead chains.

## Known limitations
- Organelle positions are illustrative. Volumes and counts come from the literature, but the arrangement is hand-designed, not taken from a tomogram.
- Mitochondria are under-filled (matrix about 55% of the expected copies placed) because of packing saturation.
- Copy numbers carry about 2-fold uncertainty.
- 63,238 transmembrane-protein copies and 6,888 copies with no usable structure were not placed.
- Not modelled: the cytoskeleton (actin cables and patches, microtubules), metabolites and water, the FG-repeat mesh inside the NPC, rDNA, and the 2-micron plasmid.
- The frame rate depends heavily on the GPU. The viewer was verified for correctness only on a software renderer (SwiftShader). Interactive speed on real hardware was not measured.

## Folder layout
```
index.html, assets/        viewer (TypeScript/WebGPU, compiled)
molstar/                   Mol* 4.18 viewer (MIT) for full-atom structures
data/library/              types.json (per-type metadata), beads.bin (residue beads), spheretree.bin (LOD)
data/tiles/                6,416 gzip tiles (instances, chromatin/mRNA polylines, membrane meshes)
data/geometry.json         analytic compartment geometry (SDF parameters, NPC positions, mesh grid)
data/organelle_meshes.bin.gz, bookmarks.json, types_tiles.bin
source/                    TypeScript source. Rebuild: cd source && npm install && mkdir -p public &&
                           ln -s ../../data public/data && cp -r ../molstar public/ && npx vite build
tests/                     headless test harness, JS-vs-Python tile decoder round-trip test and result
```
