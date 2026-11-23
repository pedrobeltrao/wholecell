# Plan: Interactive 3D whole-yeast-cell viewer (WebGPU, real structures)

## Summary
Build a browser app that lets you view one complete haploid *S. cerevisiae* cell (BY4741, YPD, same inventory as the existing model), fly from the whole cell (~4.3 µm) down to individual residues, and cut the cell along any plane in 3D. Every macromolecule is drawn from its actual PDB or AlphaFold coordinates, and clicking a molecule opens its full-atom file.

This needs two new pieces of work:
- **A real 3D packing of the entire cell.** About 32 M macromolecule instances. The current model has only a column and a 2D section.
- **A GPU-driven renderer.** It stores each of the ~5,700 structures once and draws instances of them, streaming 200 nm tiles and switching level of detail automatically.

### Decisions (you delegated these to me; recorded as assumptions)
| Topic | Decision |
|---|---|
| Fidelity | Every molecule uses its real coordinates at **one bead per residue/nucleotide** wherever it is large enough on screen to resolve. Far-away molecules switch to sphere-tree, then single-sphere stand-ins. This happens only below pixel scale, so it is not visible. **Full atoms on click:** the actual mmCIF is fetched from RCSB or AlphaFold DB and shown in an embedded Mol* panel. Full atoms for the whole cell (~1.2 × 10¹¹ atoms) is infeasible in any browser. |
| Lipids | All ~460 M lipids are drawn individually but **generated on the GPU** (deterministic seeds; measured class composition; 0.65 nm²/lipid). They are generated only for tiles close to the camera and removed wherever protein residues occupy the membrane. Far away, membranes are drawn as two-tone bilayer bands. Cell-wall β-glucan fibrils and mannan brushes are generated the same way. |
| 3D organelle layout | **Data-driven but illustrative.** Organelle volumes and counts come from the existing inventory and literature geometry; their positions are arranged by hand. No tomogram is used. |
| Hosting | (1) The full app as a static-site archive (`yeast_cell_viewer.zip`) that you unzip and serve with `python -m http.server`, or deploy to any static host. (2) A small single-file preview (`yeast_column_preview.html`, column model only) that may open directly in the results panel, if that panel allows WebGPU in its frame. |
| Graphics API | **WebGPU** (Chrome/Edge, Safari 26+, recent Firefox). If WebGPU is unavailable, the app shows a clear message. No WebGL fallback. |
| Disordered regions | AlphaFold residues with pLDDT < 50 stay trimmed in the packed cell, because extended "spaghetti" tails are unphysical and cause clashes. The inspect panel reports the trimmed fraction, and the full untrimmed model loads on click. |

## 1. Whole-cell 3D model (Python; extends existing code)

### 1a. Inventory v2
Reuse the existing inventory: 4.2 × 10⁷ proteins, compartment weights, 10 curated complexes, 18 homo-oligomers, tRNA 3.3 × 10⁶, mRNA 3 × 10⁴, nucleosomes 73,152, NPC 117. Add one extension:
- Rank Complex Portal complexes by total mass per cell.
- For the top 40 not already modeled, query the RCSB search API for an *S. cerevisiae* entry whose assembly covers at least 80% of the subunits by mass. Take the best resolution.
- Add each one found as an assembled ingredient and subtract its subunits from the monomer pool. Candidates include Sec61/translocon, OST, eIF2/eIF3, the exosome, Mediator, SAGA/TFIID, COPI, the CCR4–NOT core and the 20S proteasome.
- Save the result as `inventory/ingredients_v2.csv`, including a table of complexes that were tried but not found.

### 1b. 3D geometry
Each compartment is a signed-distance function (SDF), evaluated in numba:
- **Cell and wall.** Cell sphere r = 2,156 nm. Wall 110 nm: mannan 0–35 nm, glucan 35–110 nm. Plasma membrane below it.
- **Cortical ER.** Curved sheets covering 48% of the PM, 23 nm beneath it and 24 nm thick, with rounded rims.
- **Nucleus.** Sphere r = 889 nm. Nuclear envelope with a 34 nm perinuclear space. 117 NPCs (PDB 7N9F) placed on the envelope at least 150 nm apart, each sitting in a toroidal pore membrane. Nucleolus is a crescent of 20% of the nuclear volume, on the side facing away from the vacuole.
- **Vacuole.** 7% of cell volume.
- **Mitochondria.** 2% of cell volume. A branched tubular network of radius 150 nm (about 12 µm of total tubule length) with lamellar cristae.
- **Lipid droplets.** About 10 droplets, r = 215 nm, 1% of volume.
- **Golgi/endosomes.** About 30 cisternae and vesicles, 1% of volume.
- **Peroxisomes.** About 4, 0.2% of volume.

Membrane meshes (marching cubes at 4 nm) are exported for far-distance drawing and for generating lipids.

### 1c. Packing
- **Tiling.** The cell is divided into 200 nm cubes (~5,600 inside the cell).
- **Expected counts per tile.** Compartment copies × (tile ∩ compartment volume or membrane area) / total, stochastically rounded so cell totals are conserved.
- **Packing method.** The existing sphere-tree random-sequential-addition packer, generalized in three ways:
  - region tests use the SDFs;
  - membrane proteins are oriented along the SDF gradient, with the cytoplasmic side set by compartment;
  - tiles are packed in **8 checkerboard phases** with a halo from neighbouring tiles, so tile boundaries have no clashes.
- **Order within a tile:** NPCs → membrane proteins → chromatin → large complexes → the rest, largest first → tRNA.
- **Chromatin.** Self-avoiding persistent nucleosome walks with 12.5 nm steps in the nucleoplasm, excluding the nucleolus. Linker DNA is generated between consecutive nucleosomes.
- **mRNA.** 30,000 threads.
- **Execution.** Run in parallel on a 16 CPU / 64 GB worker. Run a 50-tile pilot first to measure throughput, then extrapolate the full runtime.
- **Checkpoints.** Each finished phase is saved to shared storage, and the save is verified before moving on.

### 1d. Model validation (`model/wholecell3d_validation.json`)
- Placed vs expected copies per compartment and per ingredient (Poisson z-scores).
- Residue-level clash rate, from 1 M sampled beads, including pairs that straddle tile boundaries.
- Containment: no residue in the wrong compartment.
- Orientation of membrane proteins.
- Crowding (mg/mL) per compartment.
- Fraction of mass that has PDB vs AlphaFold structures.
- Unplaced counts. Mitochondria are expected to saturate, as they did in 2D; this will be reported, not hidden.

## 2. Data format (static files, streamed)
- **`library/`** (~25–35 MB). `types.json` holds per-type metadata: gene, UniProt, name, copies, compartment, structure source (PDB ID or AF), mean pLDDT, trimmed fraction and palette class. Binary files hold:
  - residue beads in float16, plus kind and pLDDT as uint8;
  - beads grouped into **64-bead clusters** for GPU-driven drawing;
  - the sphere-tree LOD, the single-sphere LOD and a bounding sphere.
- **`tiles/`** holds one gzip binary per 200 nm tile. Each instance takes 16 bytes: type uint16; position as 3 × uint16 quantized within the tile (0.003 nm steps); quaternion packed smallest-three (int16 × 3 plus 2 bits). Each tile also stores seeds for lipids, glucan and mannan, and its chromatin and mRNA polylines. An `index.json` gives tile bounds, counts, compartments present and file sizes.
- **`meshes/`** holds simplified organelle and membrane meshes (~10 MB) for the overview, minimap and far-membrane drawing.
- **Expected total:** about 520 MB raw, 350–450 MB gzipped. Archive size will be similar.
- **Encoding tests.** A Python encoder and a JS decoder are tested against each other: position error < 0.01 nm and rotation error < 0.1°.

## 3. Web app (TypeScript + Vite, WebGPU; no heavy framework)

### Renderer
- **Streaming.** A web worker fetches and decodes tiles that intersect the view frustum and the visible slab behind the cut plane, prioritized by screen size. Tiles go into fixed GPU buffer pages with least-recently-used eviction under a GPU memory budget (default 1.2 GB, adjustable).
- **Per-frame compute pass.**
  - Culling: frustum culling; clip-plane culling (molecules entirely in front of the cut are dropped; molecules crossing the cut are sliced per pixel); optional Hi-Z occlusion culling.
  - LOD: chosen per instance by projected size. L0 = all residue clusters; L1 = sphere tree; L2 = one sphere.
  - Output: indirect-draw lists, so there are no per-type draw calls.
- **Drawing.** Ray-cast sphere impostors with correct depth. Molecules crossing the clip plane get flat **cut caps**, giving the Goodsell cross-section look. Lipids, glucan and mannan are generated by compute shaders only for L0 tiles, and lipids are masked by a GPU spatial hash of nearby protein residues.
- **Style passes** on a deferred G-buffer (color, ID, depth, normal):
  - Goodsell mode: flat palette color, per-type lightness jitter, outlines at ID and depth discontinuities, depth fog.
  - Optional shaded mode with ambient occlusion and lighting.
  - The palette matches the existing figures.
- **Picking** reads back the ID buffer to identify the instance and its type.

### Features
- **Navigation.** Orbit around the cell, a logarithmic zoom wheel from 5 µm to 2 nm, fly mode (WASD), and double-click to focus. Scale bar, overview minimap showing the slice plane, and FPS and draw-count stats.
- **Slicing.** A cut plane with a drag gizmo for orientation and offset. X/Y/Z and "through nucleus + pore" presets. Half-space cutaway or thin-slab mode with adjustable thickness. An animated sweep through the cell.
- **Search.** Autocomplete on gene, UniProt ID or protein name. Matching copies are highlighted and everything else is dimmed. Shows copies in the cell and in view, with "fly to nearest copy".
- **Inspect panel.**
  - Contents: gene, name, UniProt link, copies per cell, compartment, mass, and structure source (PDB ID or AF with pLDDT and trimmed fraction).
  - Atomic model: a **"Load atomic model"** button opens Mol* with the real file from RCSB or AlphaFold DB. If offline, it falls back to the residue model.
- **Color modes:** Goodsell compartment palette, PDB vs AlphaFold, per-residue pLDDT, log abundance, random per type.
- **Visibility toggles** by compartment and by class: ribosomes, tRNA, mRNA, chromatin, lipids, wall, NPCs, membranes.
- **Guided tour** of about 10 bookmarks, each a camera position plus a cut and a caption: cell wall, Pma1-rich PM, cortical ER, ribosome-crowded cytoplasm, nuclear pore, chromatin, nucleolus, mitochondrial cristae, vacuole, lipid droplet.
- **Screenshot** export at high resolution.
- **About panel** listing sources, assumptions and caveats.
- **Quality presets** plus automatic LOD-distance scaling to hold the target frame rate.

### Preview
The same renderer is bundled into one HTML file with the column model (27.6 k instances plus the library subset) embedded. Expected size is about 10–20 MB.

## 4. Testing, acceptance, compute

### Tests
- Python unit tests for the SDFs, tile expected counts, encoder round-trip and cross-tile clashes.
- Node test for the JS decoder.
- Headless Chromium (Playwright, WebGPU via SwiftShader) loads the app and checks:
  - no console errors;
  - screenshots at every tour bookmark;
  - clicking a known pixel returns the expected molecule type;
  - search for TDH3 highlights the expected count.
- CORS check that AlphaFold DB and RCSB files can be fetched from the browser.

### Acceptance
- At least 99% of expected instances placed outside mitochondria; shortfall reported per compartment.
- Clash rate < 10⁻⁵.
- All tour bookmarks render correctly in headless tests.
- **Design target, not measurable here:** ≥ 30 fps at 1080p on a mid-range GPU (Apple M1/M2, RTX 3060). The sandbox has no hardware GPU, so real-GPU performance cannot be verified here. The in-app stats overlay and quality presets let you check and tune it.

### Compute estimate
- **Basis.** The column packing measured ~24 k molecules in about 5 s on one core. Scaling to 32 M molecules is ~1.8 core-hours, or ~5 core-hours allowing for SDF overhead. On one 16 CPU / 64 GB worker that is about 30–60 min. Peak memory is ~1 GB per process.
- **Tile encoding and meshes:** under 20 min.
- **Disk:** about 1.5 GB of intermediates, ~0.5 GB of deliverables.
- **Development and headless testing** of the renderer take most of the effort.
- **No HPC is needed.**

### Deliverables
- `yeast_cell_viewer.zip`: the app, data and a README with the run command.
- `yeast_column_preview.html`
- `inventory/ingredients_v2.csv`
- `model/wholecell3d_validation.json`
- Headless screenshots
- A short `report_yeast_cell_viewer.md`: methods, assumptions, controls and limitations.

### Known limitations to state
- The organelle arrangement is illustrative.
- The far-distance LOD and the GPU-generated lipids/polysaccharides are generated, not stored. They are deterministic and reproducible.
- No cytoskeleton, metabolites, nuclear basket or FG mesh.
- Mitochondrial saturation.
- Abundances are uncertain about 2-fold.
