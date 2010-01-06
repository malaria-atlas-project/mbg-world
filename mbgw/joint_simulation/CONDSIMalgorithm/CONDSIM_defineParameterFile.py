import sys
import shutil
PARAMFILE = str(sys.argv[1])
shutil.copy(PARAMFILE,'CONDSIM_params.py')
