from pylab import *
from pymc import *
from tables import *
from numpy import *
import os,sys
from mbgw import master_grid
from mbgw.master_grid import missing_val
import gc
from IPython.kernel import client
from animate_realization import targ, chunk_to_str

"""
python animate_realization.py fname missing_val t_chunk
"""

fname = sys.argv[1]
missing_val = float32(sys.argv[2])
t_chunk = int(sys.argv[3])


t_start = 0
t_end = -1


hf = openFile(fname)
r = hf.root.realizations
if t_end<0:
    t_end = r.shape[-1] + t_end -1
if t_start < 0:
    t_start = r.shape[-1] + t_start + 1
nr = r.shape[0]
hf.close()

chunk_str = chunk_to_str(t_chunk)

os.system('rm -r anim-scratch')
os.mkdir('anim-scratch')

mec = client.MultiEngineClient()
mec.reset()
mec.execute('from animate_realization import *')
mec.push({'fname':fname,
            'missing_val':missing_val,
            't_start':t_start,
            't_end':t_end,
            't_chunk':t_chunk,
            'chunk_str':chunk_str})
            
mec.scatter('T',range(t_start, t_end, t_chunk))
print mec.execute('print fname, missing_val, t_start, t_end, t_chunk, chunk_str')
print mec.execute('print T')
print mec.execute('[targ(t,fname,missing_val,t_start,t_end,t_chunk,chunk_str) for t in T]')