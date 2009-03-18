import tables
import os, sys
import mbgw

__root__ = mbgw.__path__[0] + '/../datafiles/auxiliary_data'
data_dict = {}
for fname in os.listdir(__root__):
    if fname[-4:]=='hdf5':
        try:
            data_dict[fname[:-5]] = tables.openFile('%s/%s'%(__root__,fname)).root
        except:
            cls, inst, tb = sys.exc_info()
            print 'Warning: unable to import auxiliary data. Error message: \n'+inst.message
locals().update(data_dict)
age_dist_model = age_dist_model.chain1.PyMCsamples
parameter_model = parameter_model.chain1.PyMCsamples
