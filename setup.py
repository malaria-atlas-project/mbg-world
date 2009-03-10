# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from setuptools import setup
from numpy.distutils.misc_util import Configuration
config = Configuration('mbgw',parent_package=None,top_path=None)

config.add_extension(name='cf_helper',sources=['mbgw/cf_helper.f'])
config.add_subpackage(subpackage_name='st_cov_fun',standalone=True)

config.packages.extend(["mbgw","mbgw/google_earth","mbgw/joint_simulation","mbgw/povray","mbgw/agepr"])
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
            **(config.todict()))

