# Co-cluster Smart heat meter data and analyse clusters
This repository includes the code, used for this publication: 

**Schaffer, M., Vera-Valdés, J. E., & Marszal-Pomianowska, A. (2024). Exploring smart heat meter data: A co-clustering driven approach to analyse the energy use of single-family houses. Applied Energy, 371, 123586. https://doi.org/10.1016/j.apenergy.2024.123586.**

If you use this code, please cite the above, mentioned publication. 

## Code Overview

This GitHub repository contains the code used in the analysis and generation of figures for our scientific publication. As the used data is not publicly available, thus the code for the initial data processing is not included. See also [Data](#data). The code is organized into several files, each serving specific purposes:

0. **[`00_get_weather_data.R`](00_get_weather_data.R)**: Fetch needed weather data from the Danish Meteorological Institute (DMI). Requires an free API key from DMI.
1. **[`01_visualise_shm_data.R`](01_visualise_shm_data.R)**: Visualize the SHM data of all 4798 selected buildings. This script produces Figure 2.
2. **[`02_analyse_bcs.R`](02_analyse_bcs.R)**: Analyse the Pearson correlation coefficient of the building characteristics (Figure A_17) and calculate the Generalised Variance Inflation Factor (GVIF) (Figure 3)
3. **[`03_co_cluster_shm_data.R`](03_co_cluster_shm_data.R)**: Performs the co-clustering and the analysis of the co-clustering of the SHM data. First the optimal number of basis functions is determined. Then the data is clustered before the temporal and energy use clusters are analysed.This script produces Figure 4 to Figure 8.
4. **[`04_bcs_vs_cluster.R`](04_bcs_vs_cluster.R)**: Visually analyse the distribution of the Building Characteristics in respect to the clusters.This script produces Figure 9 and the supplementary material figures. 
5. **[`05_mlrgl.R`](05_mlrgl.R)**: Multinomial logistic regression with group lasso penalty of the building characteristics against the clusters and analyses the result. This script produces Figure 10, Figure 11, Figure B.18 and the data for Table 3.
6. **[`06_vsurf.R`](06_vsurf.R)**: Runs VSURF and analyse result before optimising a random forest with the selected variables. This script produces Figure 12, Figure 13 and data for Table 4.
7. **[`07_vsurf_merged_cluster.R`](07_vsurf_merged_cluster.R)**: Identical to [`06_vsurf.R`](06_vsurf.R) with the exception that it uses the by domain knowledge merged clusters. This script produces Figure 14, Figure 15 and data for Table 5
8. **[`08_classification_visualisation.R`](08_classification_visualisation.R)**: Decision trees to understand what BCs and how differentiate the found clusters.This script produces Figure 16 and Figures for the supplementary materials

## Data

The data used in our publication is not included in this repository.
The data used is described in this publication and can be accessed via the specified conditions:

M. Schaffer, M. Veit, A. Marszal-Pomianowska, M. Frandsen, M.Z. Pomianowski, E. Dichmann, C.G. Sørensen, J. Kragh, Dataset of smart heat and water meter data with accompanying building characteristics, Data Br. 52 (2024) 109964. https://doi.org/10.1016/j.dib.2023.109964.

The run the above scripts three external data files, which are not included are required. Below the expected format of these three files is described.

**01_shm_data.fst**

   Smart meter data as a .fst file with hourly readings for each building. The file should contain the following columns:

   | heat_meter_id | time_rounded  | energy_area | 
   | ----------- | ------------- | ---------------- |
   | Heat meter ID uniquely identifying each building | Rounded hourly reading times | Total energy use processed with SPMS and normalised by the building area in kWh/h/m²|


**02_characteristics_bbr.csv**

   All used building characteristics originating from the BBR data. The heat_meter_id is used as identifier. For a detailed overview refer to [`Table 2`](https://www.sciencedirect.com/science/article/pii/S0306261924009693?via%3Dihub#tbl2) of the publication. 

**02_characteristics_epc.csv**

   All used building characteristics originating from the EPC data. The heat_meter_id is used as identifier. For a detailed overview refer to [`Table 2`](https://www.sciencedirect.com/science/article/pii/S0306261924009693?via%3Dihub#tbl2) of the publication. 
   
## Bibtext
```
@article{SCHAFFER2024123586,
abstract = {The ongoing digitalisation of the district heating sector, particularly the installation of smart heat meters (SHMs), is generating data with unprecedented extent and temporal resolution. This data offers potential insights into heat energy use at a large scale, supporting policymakers and district heating utility companies in transforming the building sector. Clustering is crucial for representing this wealth of data in human-understandable groups, necessitating consideration of seasonality. Advancing current research in clustering SHM data, this work applies an established co-clustering approach, FunLBM, considering seasonal variation without fixed season definitions. Furthermore, to enhance the understanding of differentiating factors between clusters, the possibility to understand cluster memberships based on 26 building characteristics was analysed using classification and variable selection methods. Applying FunLBM on a large-scale hourly dataset from single-family houses revealed six well-separated energy use clusters each distributed over six-temporal clusters, which are correlated with the exterior temperature, yet not following fixed seasons. Variable selection and classification showed that building characteristics describing the building with a high level of detail are insufficient to explain cluster membership (Matthew's correlation coefficient (MCC) ≈0.3). By merging the energy use clusters based on profile and magnitude similarities, classification performance significantly improved (MCC ≈0.5). In both cases, simple and readily available building characteristics yield similar insights to detailed ones, emphasising their cost-effectiveness and practicality.},
author = {Schaffer, Markus and Vera-Vald{\'{e}}s, J Eduardo and Marszal-Pomianowska, Anna},
doi = {10.1016/j.apenergy.2024.123586},
issn = {03062619},
journal = {Applied Energy},
keywords = {Classification,Co-clustering,District heating,Smart heat meter,Variable selection},
month = {oct},
pages = {123586},
title = {{Exploring smart heat meter data: A co-clustering driven approach to analyse the energy use of single-family houses}},
url = {https://www.sciencedirect.com/science/article/pii/S0306261924009693 https://linkinghub.elsevier.com/retrieve/pii/S0306261924009693},
volume = {371},
year = {2024}
}
```
   
