# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from setuptools import setup
from numpy.distutils.misc_util import Configuration
config = Configuration('mbgw',parent_package=None,top_path=None)

# config.add_extension(name='st_cov_fun.fst_cov_fun',sources=['mbgw/st_cov_fun/fst_cov_fun.f'])
config.add_extension(name='cf_helper',sources=['mbgw/cf_helper.f'])
config.add_subpackage(subpackage_name='st_cov_fun',standalone=True)

import os
# agepr_files = []
# for f in os.listdir('mbgw/agepr'):
#     if f[-5:]=='agepr':
#         agepr_files.append('mbgw/agepr/%s'%f)
#     
# aux_dfiles = ['./datafiles']    
# for aux in os.listdir('./datafiles/auxiliary_data'):
#     if aux[-4:]=='hdf5':
#         aux_dfiles.append('datafiles/auxiliary_data/%s'%aux)

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(  version="1.0",
            description="The Malaria Atlas Project's model-based geostatistical software.",
            author="Anand Patil and Peter Gething", 
            author_email="anand.prabhakar.patil@gmail.com",
            url="www.map.zoo.ox.ac.uk",
            license="Creative Commons License",
            requires=['NumPy','PyMC','PyTables','SciPy','RPy','Matplotlib'],
            long_description="""
            blablabla
            """,
            # packages=["mbgw/google_earth","mbgw/joint_simulation","mbgw/povray","mbgw/st_cov_fun","mbgw/agepr"],
            # data_files=[('mbgw/auxiliary_data',aux_dfiles), ('mbgw/agepr',agepr_files)],
            # packages=["mbgw/google_earth","mbgw/joint_simulation","mbgw/povray","mbgw/st_cov_fun","mbgw/agepr","mbgw/auxiliary_data"],
            # data_files=[('mbgw/auxiliary_data',aux_dfiles), ('mbgw/agepr',agepr_files)],
            **(config.todict()))

