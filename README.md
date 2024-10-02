# Preleukemic Mouse Cell Atlas Pipeline
In [this study](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00252-5), we present a comprehensive analysis of the effects of mutations on hematopoiesis using single-cell RNA sequencing of 269,048 mouse bone marrow HSPCs from eight different preleukemic mutant mouse models (*Jak2* V617F, *Calr* 52-bp del, *Flt3* ITD, *Npm1c*, *Idh1* R132H, *Dnmt3a* R882H, *Ezh2* KO and *Utx* KO). 

This repository contains analysis notebooks and scripts that can be used to reproduce the main figures as well as to analyze new datasets.

An interactive web portal for this dataset is available [here](https://gottgens-lab.stemcells.cam.ac.uk/preleukemia_atlas/).<br>
*If you cannot access the web portal, you can download the data [here](https://drive.google.com/drive/folders/1vDyf3N17QraIbJ_RgLdrhQg9uWVnxUs1).

## Notebooks
The analysis notebooks demonstrate how you can reproduce our results and apply our methods to your own data:
  - <ins>**[Reference-based data integration](https://github.com/TomoyaIsobe/PMCA/tree/main/1_data_integration.ipynb)**</ins>
  - <ins>**[Differential abundance analysis](https://github.com/TomoyaIsobe/PMCA/tree/main/2_differential_abundance.ipynb)**</ins>
  - <ins>**[Differential fate probability analysis](https://github.com/TomoyaIsobe/PMCA/tree/main/3_fate_probability.ipynb)**</ins>
  - <ins>**[Differential metabolic flux analysis](https://github.com/TomoyaIsobe/PMCA/tree/main/4_metabolic_flux.ipynb)**</ins>
  - <ins>**[Cell type-wise differential expression analysis](https://github.com/TomoyaIsobe/PMCA/tree/main/5_pseudobulk_differential_expression.ipynb)**</ins>
  - <ins>**[Pseudo-temporal differential expression analysis](https://github.com/TomoyaIsobe/PMCA/tree/main/6_pseudotemporal_differential_expression.ipynb)**</ins>
