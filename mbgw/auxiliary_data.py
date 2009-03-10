import tables
import os, sys
import mbgw
try:
    __root__ = mbgw.__path__[0] + '/../datafiles/auxiliary_data'
    data_dict = {}
    for fname in os.listdir(__root__):
        if fname[-4:]=='hdf5':
            data_dict[fname[:-5]] = tables.openFile('%s/%s'%(__root__,fname)).root
    locals().update(data_dict)
    age_dist_model = age_dist_model.chain1.PyMCsamples
    parameter_model = parameter_model.chain1.PyMCsamples
except:
    cls, inst, tb = sys.exc_info()
    print 'Warning: unable to import auxiliary data. Error message: \n'+inst.message
