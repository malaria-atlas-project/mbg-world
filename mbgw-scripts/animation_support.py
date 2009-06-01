from pylab import *
from pymc import *
from tables import *
from numpy import *
import os,sys
from mbgw import master_grid
from mbgw.master_grid import missing_val, xllc, yllc, cellsize, nrows, ncols
from map_utils import grid_convert
import matplotlib
from mpl_toolkits import basemap
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


def targ(inner,continent,fname,missing_val,t_start,t_end,t_chunk,chunk_str):
    grid_lims = getattr(master_grid, continent+'_lims')
    xllc_here = (xllc + cellsize * (grid_lims['leftCol']-1))
    xurc_here = (xllc + cellsize * (grid_lims['rightCol']-1))
    yllc_here = (yllc + cellsize * (nrows-grid_lims['bottomRow']))
    yurc_here = (yllc + cellsize * (nrows-grid_lims['topRow']))

    
    outer=0


    hf = openFile(fname)
    r = hf.root.realizations
    V = hf.root.PyMCsamples.cols.V[:]
    
    sl = np.empty((r.shape[1], r.shape[2], t_chunk))
    sh = sl.shape
    
    t = hf.root.t_axis[inner] + 2009
    mo = int(round(remainder(t*12,12)))
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
    
    b = basemap.Basemap(llcrnrlon=xllc_here, llcrnrlat=yllc_here,
                        urcrnrlon=xurc_here, urcrnrlat=yurc_here)
    
    b.imshow(grid_convert(out,'y-x+','x+y+'), cmap=matplotlib.cm.hot)
    # b.drawcoastlines(color='.8',linewidth='2')
    b.drawcountries(color='.8',linewidth='1')
    colorbar()
    # axis('off')
    # title('%s starting %s'%(chunk_str, time_str))
    title(str(yr))
    savefig('anim-scratch/%i.png'%(inner))
    close('all')
    
    del sl, out
    gc.collect()
    hf.close()
    
    
