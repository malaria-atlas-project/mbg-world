from pylab import *
from pymc import *
from tables import *
from numpy import *
import os,sys
from mbgw import master_grid
from mbgw.master_grid import missing_val
import matplotlib
import gc
matplotlib.interactive('False')

def chunk_to_str(c):
    yrs = c/12
    if yrs==0:
        yr_str = ''
    elif yrs==1:
        yr_str = '1 year'
    else:
        yr_str = '%i years'%yrs
        
    mos = int(remainder(c,12))
    if mos==0:
        mo_str = ''
    elif mos==1:
        mo_str = '1 month'
    else:
        mo_str = '%i months'%mos

    if yrs>0 and mos>0:
        join_str = ', '
    else:
        join_str = ''
    
    return yr_str + join_str + mo_str
    
moname = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']    


fname = sys.argv[1]
# fname = 'test_sim.hdf5'
# t_chunk = 12
t_chunk = 1
t_start = 0
t_end = -1

chunk_str = chunk_to_str(t_chunk)

hf = openFile(fname)
r = hf.root.realizations
V = hf.root.PyMCsamples.cols.V[:]
def time_to_str(i):
    t = hf.root.t_axis[i] + 2009
    mo = int(remainder(t,1)*12)
    yr = int(t)
    return moname[mo] + ' %i'%yr

if t_end<0:
    t_end = r.shape[-1] + t_end -1
if t_start < 0:
    t_start = r.shape[-1] + t_start + 1

os.system('rm -r anim-scratch')
os.mkdir('anim-scratch')
sl = np.empty((r.shape[1], r.shape[2], t_chunk))
sh = sl.shape
for i in xrange(r.shape[0]):
    os.mkdir('anim-scratch/%i'%i)
    for j in xrange(t_start, t_end, t_chunk):
        # hf = openFile(fname)
        # r = hf.root.realizations
        # V = hf.root.PyMCsamples.cols.V[:]
        print i,j/float(t_end)
        for k in xrange(t_chunk):
            sl[:,:,k] = r[i,:,:,j + k]
            
        out = sl + np.random.normal(size=sh)*np.sqrt(V[i])
        out = invlogit(out.ravel()).reshape(sh)
        out = np.mean(out, axis=2)
        out = ma.masked_array(out, out==missing_val)
        imshow(out.T)
        colorbar()
        axis('off')
        title('%s starting %s'%(chunk_str, time_to_str(j)))
        savefig('anim-scratch/%i/%i.png'%(i,j))
        close('all')
        
        gc.collect()
        hf.flush()