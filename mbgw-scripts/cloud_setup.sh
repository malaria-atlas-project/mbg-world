#!/bin/bash

#apt-get install screen

cd /usr/lib/python2.5/site-packages/
rm -r -f pymc*
cd
svn checkout http://pymc.googlecode.com/svn/trunk/ pymc
cd pymc
python setupegg.py install
cd 

apt-get install python-boto

rm -r -f mbg-world
git clone git://github.com/malaria-atlas-project/mbg-world.git
cd mbg-world
ln -s ../datafiles datafiles
python setup.py develop
cd

rm -r -f mbg-world
git clone git://github.com/malaria-atlas-project/generic_mbg.git
cd generic_mbg 
ln -s ../datafiles datafiles
python setup.py develop
cd

rm -r -f st-cov-fun
git clone git://github.com/malaria-atlas-project/st-cov-fun.git
cd st-cov-fun
f2py -c fst_cov_fun.f -m fst_cov_fun
python setup.py install
cd

rm -r -f map_utils
git clone git://github.com/malaria-atlas-project/map_utils.git
cd map_utils
python setup.py install
cd

rm -r -f pr-incidence
git clone git://github.com/malaria-atlas-project/pr-incidence.git
cd pr-incidence
python setup.py install
cd