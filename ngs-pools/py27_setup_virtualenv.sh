#!/usr/bin/env bash
virtualenv -p python env2.7
source env2.7/bin/activate
pip install pandas
pip install ipython
pip install biopython
git clone https://github.com/Twistbioscience/twistdb.git
pip install -e ./twistdb