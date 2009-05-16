from pylab import *
from pymc import *
from tables import *
from numpy import *
import os,sys
from mbgw import master_grid
from mbgw.master_grid import missing_val
import gc
from IPython.kernel import client
from animation_support import targ, chunk_to_str

"""
python realization.py fname missing_val t_chunk
Note: you need to create an ipython cluster before using this. Such as:

ipcluster -n 4

for a 4-engine cluster.
"""

fname = sys.argv[1]
continent = sys.argv[2]
missing_val = float32(sys.argv[3])
t_chunk = int(sys.argv[4])


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
mec.execute('from animation_support import *')
mec.push({'fname':fname,
            'continent':continent,
            'missing_val':missing_val,
            't_start':t_start,
            't_end':t_end,
            't_chunk':t_chunk,
            'chunk_str':chunk_str})
            
mec.scatter('T',range(t_start, t_end, t_chunk))

print mec.execute('[targ(t,continent,fname,missing_val,t_start,t_end,t_chunk,chunk_str) for t in T]')