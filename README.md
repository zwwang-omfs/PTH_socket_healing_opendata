# PTH_socket_healing_opendata

Raw numerical data and analysis resources for the manuscript **"Teriparatide Restores Socket Healing under Bone Marrow Suppression"**.

This repository makes available the quantitative data underlying the main figures of the manuscript, together with a place for the analysis notebooks used to generate the figures from the CSV files.

## Repository structure

- `numerical_data/`: CSV files containing the numerical values that underlie the main and selected supplementary figures.
- `code/`: Folder where Python scripts and/or Jupyter notebooks for data analysis will be stored.

## Data description

All datasets are stored as comma-separated values (`.csv`). Columns are labeled to indicate experimental groups, time points, and measured variables (see each file for exact headers and units).

The files in `numerical_data/` are organized by figure and panel as follows:

- **Figure 1B (systemic blood parameters)**
	- `Figure1B_Plt.csv` – platelet-related measurements underlying Fig. 1B.
	- `Figure1B_RBC.csv` – red blood cell–related measurements underlying Fig. 1B.
	- `Figure1B_WBC.csv` – white blood cell–related measurements underlying Fig. 1B.

- **Figure 1D (histology / woven bone)**
	- `Figure1D_HE_woven_bone.csv` – quantitative woven bone data from H&E staining underlying Fig. 1D.

- **Figure 1F & 1H (TRAP-related measurements)**
	- `Figure1F_d3_trap.csv` – day 3 TRAP-related quantitative data (Fig. 1F).
	- `Figure1H_d7_trap.csv` – day 7 TRAP-related quantitative data (Fig. 1H).

- **Figure 1J (micro-CT parameters)**
	- `Figure1J_BMC-TV.csv` – bone mineral content per tissue volume (BMC/TV) measurements.
	- `Figure1J_BMD.csv` – bone mineral density (BMD) measurements.
	- `Figure1J_BV-TV.csv` – bone volume per tissue volume (BV/TV) measurements.

- **Figure 3B–3F (co-signal analyses)**
	- `Figure3B_day3_cosignal-emcn.csv` – day 3 co-signal measurements involving EMCN (Fig. 3B).
	- `Figure3C_day3_cosignal-roi.csv` – day 3 co-signal measurements within regions of interest (Fig. 3C).
	- `Figure3E_day7_cosignal-emcn.csv` – day 7 co-signal measurements involving EMCN (Fig. 3E).
	- `Figure3F_day7_cosignal-roi.csv` – day 7 co-signal measurements within regions of interest (Fig. 3F).

- **Figure 4C–4F (vascular / filament analyses)**
	- `Figure4C_d3_emcnarea_roiarea.csv` – day 3 EMCN-positive area normalized to ROI area (Fig. 4C).
	- `Figure4C_d7_emcnarea_roiarea.csv` – day 7 EMCN-positive area normalized to ROI area (Fig. 4C).
	- `FIgure4D_d3_length_area.csv` – day 3 vessel/structure length per area (Fig. 4D).
	- `Figure4D_d7_length_area.csv` – day 7 vessel/structure length per area (Fig. 4D).
	- `Figure4E_d3_filament_area.csv` – day 3 filament area–related measurements (Fig. 4E).
	- `Figure4E_d7_filament_area.csv` – day 7 filament area–related measurements (Fig. 4E).
	- `Figure4F_d3_filament_length.csv` – day 3 filament length–related measurements (Fig. 4F).
	- `Figure4F_d7_filament_length.csv` – day 7 filament length–related measurements (Fig. 4F).

- **Figure 5A–5I (cellular assays and tube formation)**
	- `Figure5A_cck8_hueht2.csv` – CCK-8 assay data for HUEhT-2 cells (Fig. 5A).
	- `Figure5C_hueht2_scratch.csv` – scratch (wound-healing) assay data for HUEhT-2 cells (Fig. 5C).
	- `Figure5E_hueht2_tube_Nb_meshes.csv` – tube formation assay: number of meshes (HUEhT-2, Fig. 5E).
	- `Figure5F_hueht2_tube_tot_branching_length.csv` – tube formation assay: total branching length (HUEhT-2, Fig. 5F).

- **Figure 5G–5J (MC3T3-E1 assays and gene expression)**
	- `Figure5G_cck8_mc3t3.csv` – CCK-8 assay data for MC3T3-E1 cells (Fig. 5G).
	- `Figure5I_mc3t3_scratch.csv` – scratch assay data for MC3T3-E1 cells (Fig. 5I).
	- `Figure5J_mc3t3_pcr_d3_alp.csv` – day 3 qPCR data for Alp expression (Fig. 5J).
	- `Figure5J_mc3t3_pcr_d3_col1.csv` – day 3 qPCR data for Col1 expression (Fig. 5J).
	- `Figure5J_mc3t3_pcr_d3_runx2.csv` – day 3 qPCR data for Runx2 expression (Fig. 5J).

All data are aggregated quantitative measurements; no personally identifiable information is included.

The raw image data underlying these quantitative measurements will be made available via a DOI-linked repository (DOI to be added once assigned).



## Citation

If you use this dataset or the accompanying code in your own work, please cite the associated manuscript:

> Authors. *Teriparatide Restores Socket Healing under Bone Marrow Suppression*. Journal, year (details will be updated upon publication).

The exact citation (authors, journal, year, DOI) will be added here once the article is formally published.

## Contact

For questions about the data, analysis, or potential reuse, please open an issue in this GitHub repository or contact the corresponding author of the manuscript.

## License / data use

License details for this dataset will be finalized and documented here upon manuscript acceptance. Until then, the data are provided for research and educational purposes. For other uses (e.g., commercial applications or redistribution as part of another dataset), please contact the authors.
