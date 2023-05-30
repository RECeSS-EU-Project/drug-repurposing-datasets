#!/bin/bash

## Documentation: https://github.com/DGRC-PT/HPOSim_Helper_Scripts (Version May 2023)
if [ ! -d "HPOSim_Helper_Scripts/" ]; then git clone https://github.com/DGRC-PT/HPOSim_Helper_Scripts.git; fi

helper_folder="HPOSim_Helper_Scripts/hpodata.rda update/"

## hpodata.rda creation
if [ ! -f $helper_folder"genes_to_phenotype.txt" ]; then wget -O $helper_folder"genes_to_phenotype.txt" "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"; fi
if [ ! -f $helper_folder"phenotype_to_genes.txt" ]; then wget -O $helper_folder"phenotype_to_genes.txt" "http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"; fi
## Take into account the update of HPO from April 2023 (disease id no longer in the 8th column, but the 3rd)
sed -i 's/8/3/g' $helper_folder"convert_to_prepare_rda.py"
cd $helper_folder && python3 convert_to_prepare_rda.py genes_to_phenotype.txt phen_R disease_hpo disease_gene terms term_disease phenotype_to_genes.txt; cd ../../; 
cd $helper_folder && Rscript "make_RDA.R"; cd ../../;

## HPO.db creation (R 4.1.2)
helper_folder="HPOSim_Helper_Scripts/HPO.db update/"
if [ ! -f $helper_folder"hp.obo" ]; then wget -O $helper_folder"hp.obo" "https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo"; fi
if [ ! -d $helper_folder"HPO.db/" ]; then cd $helper_folder && wget "https://downloads.sourceforge.net/project/hposim/HPO.db_1.9.tar.gz" && gzip -d "HPO.db_1.9.tar.gz" && tar -xvf "HPO.db_1.9.tar"; cd ../..; fi
cd $helper_folder && python3 obo2sqlite.py hp.obo && mv HPO.sqlite HPO.db/inst/extdata/ && tar -czvf HPO.db_2.0.tar.gz HPO.db; cd ../../;

## Use new HPO.db+HPOSim package in R
R -e 'install.packages("'$helper_folder'HPO.db_2.0.tar.gz", repos = NULL, type = "source")'
R -e 'install.packages("igraph");library(devtools);install_github("cran/HPOSim")'