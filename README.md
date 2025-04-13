<!-- 
README for "Spatial Clustering of Citizen Science Data Improves Downstream Species Distribution Models"

# Nahian Ahmed
# December 10, 2024 
-->

# Spatial Clustering of Citizen Science Data Improves Downstream Species Distribution Models


<p align="center">
  <img src="site_clustering_process_flow.png"  width="800">
</p>


## Data 

Data available at <https://doi.org/10.5281/zenodo.14362178>.

`occupancy_feature_raster` contains raster file of occupancy/site features, and `checklist_data` contains eBird checklists of 31 species. 

Please place both directories in root directory before executing code.


## Requirements

Code was last run using `R` version 4.4.2, `gdal` version 3.10.0, `geos` version 3.13.0, and `proj` version 9.5.0.

Required `R` packages are listed in `dependencies.txt`.

## Instructions


Run the following command:

`Rscript run_species_experiments.R`

Replace `species_names` in `run_species_experiments.R` to run experiments for specific species. Abbreviations of species names are in `checklist_data` (e.g, `NOFL` = Northern Flicker, `COHA` = Cooper's Hawk, etc.).


@article{Ahmed_Roth_Hallman_Robinson_Hutchinson_2025, 
title={Spatial Clustering of Citizen Science Data Improves Downstream Species Distribution Models}, volume={39}, url={https://ojs.aaai.org/index.php/AAAI/article/view/34993}, DOI={10.1609/aaai.v39i27.34993}, 
number={27}, journal={Proceedings of the AAAI Conference on Artificial Intelligence}, author={Ahmed, Nahian and Roth, Mark and Hallman, Tyler A. and Robinson, W. Douglas and Hutchinson, Rebecca A.}, year={2025}, month={Apr.}, pages={27775-27783} }


## Citation

```
@inproceedings{Ahmed_Roth_Hallman_Robinson_Hutchinson_2025,
  title={Spatial Clustering of Citizen Science Data Improves Downstream Species Distribution Models},
  author={Ahmed, Nahian and Roth, Mark and Hallman, Tyler A and Robinson, W Douglas and Hutchinson, Rebecca A},
  journal={Proceedings of the AAAI Conference on Artificial Intelligence},
  volume={39}, 
  number={27},
  pages={27775-27783},
  url={https://ojs.aaai.org/index.php/AAAI/article/view/34993},
  DOI={10.1609/aaai.v39i27.34993},
  month={Apr.},
  year={2025},
}
````
