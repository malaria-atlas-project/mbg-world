# Check out latest version of pymc only:
svn export http://pymc.googlecode.com/svn/trunk/ pymc

# Check out latest revision of MBGW only:
git clone git://github.com/malaria-atlas-project/mbg-world.git --depth 1

# File transfers:
datafiles/auxiliary_data/landSea-e.hdf5
datafiles/auxiliary_data/age_dist_model.hdf5
datafiles/auxiliary_data/parameter_model.hdf5
datafiles/auxiliary_data/periurb.hdf5
datafiles/auxiliary_data/urb.hdf5
#Also any traces.


# install pymc
/share/apps/python/bin/python setup.py install --prefix ~
export PYTHONPATH=~/lib/python2.6/site-packages:$PYTHONPATH