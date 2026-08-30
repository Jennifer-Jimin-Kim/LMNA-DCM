# LMNA-DCM Atlas — portal

Static site. Same template/engine as the Human PVR Cell Atlas.

## Run locally
    cd lmna-dcm-atlas
    python3 -m http.server 8000
    # http://localhost:8000

Opening index.html directly with file:// will NOT work — the page fetches
data/*.json, which browsers block on file://.

## Deploy (Vercel)
    npx vercel --prod
Framework preset: "Other". No build step; output directory is this folder.

## Contents
    index.html          page shell (PVR CSS + markup)
    js/atlas.js         snRNA explorer engine (UMAP, composition, states, markers)
    js/spatial.js       Visium viewer (spots over H&E) + QC table
    data/atlas.json     53,253 nuclei — UMAP, lineage, cell state, group, 35-gene panel
    data/spatial.json   4 Visium sections — spot coords + 20-gene panel + metrics
    img/*.png           tissue_lowres_image.png per section

## Regenerating data
    python3 export.py           # h5ad  -> data/atlas.json
    python3 export_spatial.py   # Space Ranger outs -> data/spatial.json + img/

## Data mapping used
- lineage      = obs['Integrated_celltype_Human_mouse_atlas_250709'] (12, full coverage)
- cell state   = obs['Celltype_integrated'] (35, full coverage)
- group (A/B)  = obs['group'] -> Control (24,300) / LMNA-DCM (28,953)
- sample rows  = donor_id > sample > sample_id > orig.ident > dataset (114 rows)
- expression   = X (log-normalised), scaled per gene to its 99.5th percentile

## Placeholders still to fill
- hero background image (currently a gradient)
- institution / funder logos (dashed boxes in Team & funding)
- all download hrefs, publication and code links
- named people in the Team grid
