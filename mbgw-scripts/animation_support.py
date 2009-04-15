from pylab import *
from pymc import *
from tables import *
from numpy import *
import os,sys
from mbgw import master_grid
from mbgw.master_grid import missing_val
import matplotlib
import processing
import gc
matplotlib.interactive(False)

"""
python animate_realization.py fname missing_val t_chunk
"""

def chunk_to_str(c):
    yrs = c/12
    if yrs==0:
        yr_str = ''
    elif yrs==1:
        yr_str = '1 year'
    else:
        yr_str = '%i years'%yrs
        
    mos = int(rem(c,12))
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


def targ(inner,fname,missing_val,t_start,t_end,t_chunk,chunk_str):
    outer=0
    print fname, missing_val, t_start, t_end, t_chunk, chunk_str
    print outer,inner/float(t_end)
    hf = openFile(fname)
    r = hf.root.realizations
    V = hf.root.PyMCsamples.cols.V[:]
    
    sl = np.empty((r.shape[1], r.shape[2], t_chunk))
    sh = sl.shape
    
    t = hf.root.t_axis[outer] + 2009
    mo = int(rem(t,1)*12)
    yr = int(t)
    time_str = moname[mo] + ' %i'%yr
    
    
    for k in xrange(t_chunk):
        sl[:,:,k] = r[outer,:,:,inner + k]
    sl[sl==missing_val]=NaN
    out = 0
    
    out = sl + np.random.normal(size=sh)*np.sqrt(V[outer])
    out = invlogit(out.ravel()).reshape(sh)
    out = np.mean(out, axis=2)
    out = ma.masked_array(out, isnan(out))
    
    imshow(out.T)
    colorbar()
    axis('off')
    title('%s starting %s'%(chunk_str, time_str))
    savefig('anim-scratch/%i.png'%(inner))
    close('all')
    
    del sl, out
    gc.collect()
    hf.close()
    
    return sh
    
    
