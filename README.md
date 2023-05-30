# Drug Repurposing Datasets for Collaborative Filtering Methods

This repository is a part of the EU-funded [RECeSS project](https://recess-eu-project.github.io) (#101102016), and hosts the notebooks and code files which helped to build the following novel drug repurposing datasets:

- TRANSCRIPT dataset [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7982969.svg)](https://doi.org/10.5281/zenodo.7982969)

- PREDICT dataset [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7982964.svg)](https://doi.org/10.5281/zenodo.7982964)

The latest versions (v2.0.0 and v2.0.1 of May 29th, 2023) of those datasets correspond to the current notebooks in folder *notebooks/*. The notebooks corresponding to the v1.0.0 and v1.0.1 (December 28th, 2022) are hosted in folder *notebooks/v1.0.0/*. To get more details about the datasets themselves, please refer to the corresponding Zenodo pages.

## Running the notebooks

### Global variables

Update to your liking the following file *paths_global.py*. The corresponding files and how to obtain them are detailed in notebooks/FEATURELESS_dataset.ipynb.

### Environment

In order to run the notebooks:

```
conda create -n drug_repurposing python=3.8.5 -y
conda activate drug_repurposing
python3 -m pip install -r requirements.txt
conda deactivate
conda activate drug_repurposing
cd notebooks/ && jupyter notebook
```

## Licence

These datasets and code are under an [OSI-approved](https://opensource.org/licenses/) [MIT license](https://raw.githubusercontent.com/RECeSS-EU-Project/drug-repurposing-datasets/main/LICENSE). Note that the final matrices for PREDICT are not available, due to the usage rules of the DrugBank dataset. However, one can run the corresponding notebook PREDICT_dataset.ipynb in order to get the final matrices. A partial dataset is downloadable on Zenodo.

## Community guidelines with respect to contributions, issue reporting, and support

[Pull requests](https://github.com/RECeSS-EU-Project/drug-repurposing-datasets/pulls) and [issue flagging](https://github.com/RECeSS-EU-Project/drug-repurposing-datasets/issues) are welcome, and can be made through the GitHub interface. Support can be provided by reaching out to recess-project-at-proton.me.