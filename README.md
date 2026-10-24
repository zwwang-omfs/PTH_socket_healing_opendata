# PTH_socket_healing_opendata

Analysis code and numerical data for the manuscript **"Teriparatide Restores Socket Healing under Bone Marrow Suppression"**.

Raw and representative image datasets are deposited in Figshare (link in the manuscript's Data Availability statement).

## Repository structure

```
code/
  image_analysis/          Spatial distribution of EMCN+ vessels and Osx+ cells (Python)
    socket_margin/         Shared module: socket ROI, masks, distance features, KDE, plotting
    kde_spatial_distribution.py
    distance_heatmaps.py
    vessel_skeleton.py
  statistics/              Statistical models (R and Python)
    lmm_section_staining.R
    install_r_packages.R
    linear_model_hc3_microct_cbc.py
    two_way_anova_in_vitro.py
    build_figure_significance.py
data/
  images/                  Binary masks used by code/image_analysis
    emcn/d3, emcn/d7       <id>.png (traced socket outline), <id>_threshold.png (EMCN+ vessels, black)
    osx/d3, osx/d7         <id>_osx.png (Osx+ cells, white)
  statistics/              Input tables for code/statistics
    lmm_section_staining/  One row per root region: group, mouse_id, root, value
    microct_cbc/           One column per group, one row per mouse
    in_vitro/              CY_conc (4-HC, uM), P_trea (iPTH), outcome
numerical_data/            Per-figure value tables (first submission)
results/
  image_analysis/          Created by the image-analysis scripts
  statistics/              Model outputs used in the manuscript
    figure_significance/   Tables used for the significance labels in the figures
```

### Sample naming

Image sample IDs have the form `<day><group><mouse>_<n>`:

- `<day>`: `3` = day 3, `7` = day 7 after tooth extraction.
- `<group>`: `c` = control (Ctrl), `cy` = cyclophosphamide (Cy), `cyp` = cyclophosphamide + intermittent PTH (Cy + iPTH), `p` = intermittent PTH (iPTH).
- `<mouse>`: mouse number within the group. On day 3, the six mice of each group were split between immunofluorescence and light-sheet imaging, so the immunofluorescence mice carry a two-digit number (`11`, `12`, `13`). On day 7 the two cohorts were handled separately from the start, and a single digit is used (`1`, `2`, `3`).
- `_<n>`: serial number of the image from that mouse. It does not correspond to a fixed root or root region.

Examples: `3cyp12_2` = day 3, Cy + iPTH, mouse 12, image 2; `7c1_3` = day 7, Ctrl, mouse 1, image 3.

Exception: the day-7 iPTH socket outlines and EMCN masks (`data/images/emcn/d7`) use the two-digit form (`7p11_1` …), whereas the corresponding Osx masks use the one-digit form (`7p1_1` …). The code maps `7p<m>_<n>` to `7p1<m>_<n>` automatically.

Image acquisition, thresholding and quantification are described in the Methods and Supplementary Methods of the article.

### TRAP and HE images

The bright-field TRAP (days 3 and 7) and HE (day 7) images are deposited in Figshare, sorted by day (`d3`, `d7`) and group (`c`, `cy`, `cyp`, `p`). File names have the form `<day><group>_<mouse>_<magnification>_<roots>[_<n>]` (TRAP) or `<day><group><mouse>_<magnification>_<roots>` (HE), for example `7cyp_12_10x_d_1.jpg`:

- `<mouse>`: `11`, `12`, `13` on both days.
- `<magnification>`: `4x` or `10x`.
- `<roots>`: the extraction-socket root region(s) in the field: `d` = distal, `m` = mesial, `p` = palatal. A 4× field may cover two neighbouring roots (`m_d` or `md`, `m_p`).
- `_<n>`: serial number when several images were taken of the same field.

The sections were prepared with the Kawamoto film method. Some sections are slightly wrinkled and do not lie in a single focal plane, so several images focused at different depths were taken of the same root region; these carry the serial numbers `_1`, `_2`, …. Not every section contains all three roots, so some mice have no image of a given root.

### In vitro images

The in vitro images (Fig. 5) are not included in this repository. The raw and representative bright-field images are deposited in Figshare together with the other image datasets; this repository holds only the quantified values (`data/statistics/in_vitro`, `numerical_data/Figure5*`). Groups are named `cy<c>` (4-HC at `<c>` µM) and `cyp<c>` (4-HC at `<c>` µM + iPTH).

| Experiment | Cell line | 4-HC (µM) | Time points imaged | Representative images | Quantified value | Input table |
|---|---|---|---|---|---|---|
| Scratch (wound healing) assay, Fig. 5B/C | HUEhT-2 | 0, 2, 4, 8 | 0, 2, 4, 6 h | 0 h and 4 h | Healing rate (mm²/h) | `hueht2_scratch.csv` |
| Tube formation assay, Fig. 5D–F | HUEhT-2 | 0, 4, 8 | 2 h | 2 h | Number of meshes; total branching length | `hueht2_tube_Nb_meshes.csv`, `hueht2_tube_tot_branching_length_new.csv` |
| Scratch (wound healing) assay, Fig. 5H/I | MC3T3-E1 | 0, 4, 8, 16 | 0, 4, 6, 12, 24 h | 0 h and 12 h | Wound closure | `mc3t3_scratch.csv` |

Images were taken at 5× magnification. Wound areas (scratch assays) and tube networks were measured in ImageJ; these steps are not scripted.

### CCK-8 and RT-qPCR data

CCK-8 and RT-qPCR results are provided as tables only, in `data/statistics/in_vitro`. Each file has three columns: `CY_conc` (4-HC, µM), `P_trea` (iPTH: `No`/`Yes`) and the outcome; each row is one independent observation.

| Assay | Cell line | 4-HC (µM) | n per condition | Outcome column | Files |
|---|---|---|---|---|---|
| CCK-8 cell viability, Fig. 5A | HUEhT-2 | 0, 2, 4, 8, 16 | 5 | `Viability` | `cck8_hueht2.csv` |
| CCK-8 cell viability, Fig. 5G | MC3T3-E1 | 0, 2, 4, 8, 16 | 5 | `Viability` | `cck8_mc3t3.csv` |
| RT-qPCR, day 3 of osteogenic differentiation, Fig. 5J | MC3T3-E1 | 0, 4, 8, 16 | 3 | `ΔΔCT` | `mc3t3_pcr_d3_alp.csv`, `mc3t3_pcr_d3_col1.csv`, `mc3t3_pcr_d3_runx2.csv` (*Alp*, *Col1*, *Runx2*) |

- `Viability` is the absorbance normalized to the mean of the 0 µM / no-iPTH condition (mean = 1).
- `ΔΔCT` is the ΔCt of each sample minus the mean ΔCt of the 0 µM / no-iPTH condition (mean = 0); relative expression is 2^−ΔΔCt.

### Micro-CT data

The raw micro-CT images are not included in this repository or in the Figshare deposit. The maxillary bone mineral density analysis cannot be performed in ImageJ and requires dedicated micro-CT software (TRI/3D-BON); the raw image data in that software's format are too large for the available storage. This repository contains the morphometric parameters exported for each mouse, in `data/statistics/microct_cbc` (one column per group: Control, Cyclophosphamide, Cyclophosphamide + Teriparatide, Teriparatide; one row per mouse, n = 6 per group):

| Site | Parameters | Files |
|---|---|---|
| Maxilla (extraction socket), Fig. 1J | BMD, BMC/TV, BV/TV | `BMD.csv`, `BMC-TV.csv`, `BV-TV.csv` |
| Leg (trabecular bone) | BV/TV, Tb.N, Tb.Th, Tb.Sp, Conn.D, average thickness | `BVTV.csv`, `Tb.N.csv`, `Tb.Th.csv`, `Tb.Sp.csv`, `Conn.D.csv`, `thickness_average.csv` |

### Complete blood count data

The complete blood counts of the mice are in the same folder, `data/statistics/microct_cbc`, in the same layout (one column per group, one row per mouse, n = 3 per group). Blood collection, timing and units are described in the Methods of the article.

| Parameter | File |
|---|---|
| White blood cells (WBC), Fig. 1B | `WBC.csv` |
| Red blood cells (RBC), Fig. 1B | `RBC.csv` |
| Platelets (Plt), Fig. 1B | `Plt.csv` |

## Figure / table to code

| Manuscript item | Code | Input | Output |
|---|---|---|---|
| Fig. 1B, 1D, S1 (CBC, micro-CT) | `statistics/linear_model_hc3_microct_cbc.py` | `data/statistics/microct_cbc` | `results/statistics/microct_cbc_hc3` |
| Fig. 1F/H/J, 3B/C/E/F, 4C–F; Table S2 (section staining) | `statistics/lmm_section_staining.R` | `data/statistics/lmm_section_staining` | `results/statistics/lmm_section_staining` |
| Fig. 5 (in vitro) | `statistics/two_way_anova_in_vitro.py`, `statistics/build_figure_significance.py` | `data/statistics/in_vitro` | `results/statistics/in_vitro_two_way_anova`, `results/statistics/figure_significance` |
| Fig. 4A/G, S4A, S6A (heatmaps) | `image_analysis/distance_heatmaps.py` | `data/images` | `results/image_analysis/heatmaps` |
| Fig. 4B/H, S4B, S6B (KDE curves) | `image_analysis/kde_spatial_distribution.py` | `data/images` | `results/image_analysis/kde` |
| Fig. S5 (initial vessel centrelines) | `image_analysis/vessel_skeleton.py` | `data/images/emcn` | `results/image_analysis/skeleton` |

Steps performed interactively in Imaris or ImageJ are not scripted: image thresholding and socket-outline tracing, CD31/EMCN co-localization areas, micro-CT morphometry (TRI/3D-BON, ImageJ/BoneJ), the semi-manual correction of vessel centrelines (Fig. S5), and the 3D light-sheet renderings (Fig. 2).

## Running the analyses

Python 3.14 and R 4.6 were used. From the repository root:

```bash
python -m pip install -r requirements.txt

# Statistics
Rscript code/statistics/install_r_packages.R
Rscript code/statistics/lmm_section_staining.R
python code/statistics/linear_model_hc3_microct_cbc.py
python code/statistics/two_way_anova_in_vitro.py
python code/statistics/build_figure_significance.py

# Image analysis
python code/image_analysis/kde_spatial_distribution.py --marker emcn
python code/image_analysis/kde_spatial_distribution.py --marker osx
python code/image_analysis/distance_heatmaps.py --marker emcn --day d7 --samples 7p11_2
python code/image_analysis/distance_heatmaps.py --marker osx --day d3 --samples 3p11_1
python code/image_analysis/vessel_skeleton.py data/images/emcn/d3 --out-dir results/image_analysis/skeleton/d3
```

The Python scripts print their options with `--help`. The statistics scripts write to `results/statistics/` by default and therefore overwrite the reference outputs; use `--output-dir` (Python) or the `LMM_OUTPUT_DIR` environment variable (R) to write elsewhere.

## Reproducibility notes

- Planned contrasts are adjusted for multiplicity with the multivariate-t method. In the Python analyses this probability is computed by quasi-Monte Carlo integration with a fixed seed; on a different computer the adjusted p-values agree with the files in `results/statistics/` to within about 1e-3 (unadjusted p-values agree to machine precision), and no significance label changes. The R analysis reproduces `results/statistics/lmm_section_staining/` exactly.
- Image analysis: the per-pixel distances, KDE curves and heatmaps depend only on the masks in `data/images`. Library versions can shift the KDE curves slightly (in our checks by less than 0.02 percentage points of the area-zone proportions).
- One pixel = 1.243 µm in the analyzed maximum-intensity projections.

## Citation

If you use these data or code, please cite the associated article (citation details will be added upon publication).

## Contact

For questions about the data or analysis, please open an issue in this repository or contact the corresponding author.

## License

License terms will be added upon acceptance of the manuscript.
