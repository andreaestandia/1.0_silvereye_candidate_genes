# Candidate gene length polymorphisms are linked to dispersive behaviour: a mechanism behind the “paradox of the great speciator”?

This repository contains the code necessary to reproduce all analyses in Estandía et al. (2023) Candidate gene polymorphisms are linked to dispersive and migratory behaviour: searching for a mechanism behind the “paradox of the great speciators” *bioRxiv*
DOI: https://doi.org/10.1101/2023.01.19.524190

To download this repository just open a terminal and paste:

```git clone https://github.com/andreaestandia/1.0_silvereye_candidate_genes.git```

The structure is the following:

```src``` contains `ae-candidate-gene-analysis_source.R` all the functions used in the R notebooks and my_scripts includes all the python functions I wrote to do data wrangling with BEAGLE, BA3, and VCF files 

`docs` contains the newest version of the MS

`notebooks` contains all the R notebooks and bash code to generate the final BEAGLE file that I used to estimate population structure with PCAngsd and NGSadmix.

You should download: i) the data and integrate it within the 1.0_silvereye_candidate_genes folder as a folder called ```data```, and ii) ```reports``` folder.

You can find both here: https://doi.org/10.5061/dryad.63xsj3v6f

If you have any questions, feel free to send me an email andrea.estandia@biology.ox.ac.uk
