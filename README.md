# Candidate gene length polymorphisms are linked to dispersive behaviour: a mechanism behind the “paradox of the great speciator”?

This repository contains the code necessary to reproduce all analyses in Estandía et al. (2022) Candidate gene length polymorphisms are linked to dispersive behaviour: a mechanism behind the “paradox of the great speciator”? *bioRxiv*

To download this repository just open a terminal and paste:

```git clone https://github.com/andreaestandia/1.0_silvereye_candidate_genes.git```

The structure is the following:

```docs``` contains the:

* Preprint
* Supplementary material

```src``` contains `ae-candidate-gene-analysis_source.R` all the functions used in the R notebooks and `subset_beagle.py` to subset any BEAGLE file by providing a list of individuals

`notebooks` contains all the R notebooks and bash code to generate the final BEAGLE file that I used to estimate population structure with PCAngsd and NGSadmix.

You should download: i) the data and integrate it within the 1.0_silvereye_candidate_genes folder as a folder called ```data```, and ii) ```reports``` folder.

You can find both here: [dryad link]

If you have any questions, feel free to send me an email andrea.estandia@biology.ox.ac.uk