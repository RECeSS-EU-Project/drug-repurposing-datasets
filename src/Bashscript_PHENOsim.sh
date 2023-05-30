#!/bin/bash

## Installation: http://pymeshsim.systemsgenetics.cn/install.html (needs Python 3.6.5)
conda create -n pysim python=3.6.5
conda activate pysim
apt-get install -y python3-bcolz git

CURRENT_DIR=$(pwd)

## Download file (requires (academic) sign-in)
if [ -f "~/public_mm_linux_main_2020.tar.bz2" ]; then echo "Download file at https://data.lhncbc.nlm.nih.gov/umls-restricted/ii/tools/MetaMap/download/public_mm_linux_main_2020.tar.bz2 and place it in folder ~/ (/home/<user>/)" && exit; fi

expected_sha="de166603003a6fa6cb2fc1eba34f84443e57229bb757397b906200b21eb14428"
found_sha=$(sha256sum ~/public_mm_linux_main_2020.tar.bz2 | cut -d" " -f1)
if [ "$found_sha" != "$expected_sha" ]; then echo "Wrong file" && exit; fi

if [ -d "~/public_mm/" ]; then cd ~/public_mm/ && tar -xvjf ~/public_mm_linux_main_2020.tar.bz2 && bash ./bin/install.sh && export PATH=$PATH":~/public_mm/bin"; fi

cd $CURRENT_DIR
python -m pip install Cython==0.26 --install-option='--no-cython-compile' -y 
python -m pip install pandas==0.22 pymysql==1.0.2 -y 
conda install -c conda-forge bcolz -y 
git clone https://github.com/luozhhub/pyMeSHSim.git
cd pyMeSHSim && git checkout 7771d29428c0b61f664ce6cec90211b97e9ef9da && python ./setup.py install

sed "s/df = pd.DataFrame.from_items/from collections import OrderedDict\n        df = pd.DataFrame.from_dict(OrderedDict/" ~/miniconda3/envs/pysim/lib/python3.6/site-packages/bcolz/ctable.py | sed "s/.* for key in keys.*/            (key, self[key][:]) for key in keys),/" > ctable.py
mv ctable.py ~/miniconda3/envs/pysim/lib/python3.6/site-packages/bcolz/

## https://pymeshsim.readthedocs.io/en/latest/tutorial.html#sim
## http://pymeshsim.systemsgenetics.cn/tutorial.html
DISEASES=$(cat diseases.txt | sed 's/,/","/g')

~/public_mm/bin/skrmedpostctl start
~/public_mm/bin/wsdserverctl start 
python3 -c 'diseases_=["'$DISEASES'"];import numpy as np;from pyMeSHSim.Sim.similarity import termComp;simCom = termComp();calcTerm=lambda xi:exec("try:print(simCom.getMeSHConcept(cui=xi))\nexcept:print(\"this is an invalid UMLS ID\")");[calcTerm(x) for ix,x in enumerate(diseases_)]' | sed '/this is an invalid UMLS ID/d' | cut -d' ' -f2,4 | cut -d',' -f1,2 > disease_phenotype_MeshTERMS.out
~/public_mm/bin/wsdserverctl stop
~/public_mm/bin/skrmedpostctl stop

rm -f diseases.txt
conda deactivate