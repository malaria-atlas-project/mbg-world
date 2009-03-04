import tables
import os, sys
import mbgw
try:
    __root__ = mbgw.__path__[0] + '/../datafiles/auxiliary_data'
    data_dict = {}
    for fname in os.listdir(__root__):
        if fname[-4:]=='hdf5':
            # print fname
            data_dict[fname[:-5]] = tables.openFile('%s/%s'%(__root__,fname)).root
    locals().update(data_dict)
    age_dist_model = age_dist_model.chain1.PyMCsamples
    parameter_model = parameter_model.chain1.PyMCsamples
except:
    cls, inst, tb = sys.exc_info()
    print 'Warning: unable to import auxiliary data. Error message: \n'+inst.message
        # exec('%s=tables.openFile(%s/%s)'%(fname[:-5],__root__,fname))
# urb = tables.openFile(__root__+'/urb.hdf5').root
# periurb = tables.openFile(__root__+'/periurb.hdf5').root
# parameter_model = tables.openFile(__root__+'/parameter_model.hdf5').root.chain1.PyMCsamples
# age_dist_model = tables.openFile(__root__+'/age_dist_model.hdf5').root.chain1.PyMCsamples
# ndvi = tables.openFile(__root__+'/mbg-ndvi.hdf5').root
