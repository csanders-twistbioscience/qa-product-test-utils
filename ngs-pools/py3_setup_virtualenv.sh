virtualenv -p python3 env
source env/bin/activate
pip3 install pandas
pip3 install ipython
pip3 install biopython
git clone https://github.com/Twistbioscience/twistdb.git
pip3 install -e ./twistdb
