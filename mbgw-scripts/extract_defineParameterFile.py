import sys
import shutil
PARAMFILE = str(sys.argv[1])
shutil.copy(PARAMFILE,'extract_params.py')
