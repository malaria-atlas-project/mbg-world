import tables as tb
from mbgw import master_grid
from mbgw.master_grid import deg_to_rad
import numpy as np 
import pymc as pm

origin = np.array([[0,0,0]])

def realization2acf(fname, index, nugget=False):
    """
    Converts an hdf5 file containing a realization to a theoretical autocovariance function.
    """
    hf=tb.openFile(fname)
    C = hf.root.group0.C[index]
    if nugget:
        V = hf.root.PyMCsamples.cols.V[index]
    else:
        V = 0
    hf.close()
    def acf(delta, C=C, V=V):
        d=np.array(delta, copy=True)
        is_origin = np.array([(row==0).all() for row in d])
        d[:,:2] *= deg_to_rad
        return np.asarray(C(origin,d)).squeeze() + V*is_origin
    acf.covdoc ="""
    Covariance function wrapped: %s.
    
    Nugget: %i.

    Covariance function parameters: \n\t- %s
        """%(C.eval_fun.__name__,V, '\n\t- '.join([item[0]+': '+ str(item[1]) for item in C.params.iteritems()]))
    
    acf.__doc__ ="""
    An autocovariance function. This function takes one argument: 
    an array of dimension (n,3). The columns should correspond to the
    lon (degrees), lat (degrees), and time (years) of the difference.
    """+acf.covdoc
    return acf

def acf2semivariance(acf):
    """
    Converts an autocovariance function to a semivariance function.
    """
    def semivar(delta, acf=acf):
        return acf(origin) - acf(delta)
    semivar.covdoc = acf.covdoc
    semivar.__doc__ ="""
    A semivariance function. This function takes one argument: 
    an array of dimension (n,3). The columns should correspond to the
    lon (degrees), lat (degrees), and time (years) of the difference.
    """ + acf.covdoc
    return semivar

def realization2semivariance(fname, index, nugget=False):
    return acf2semivariance(realization2acf(fname,index, nugget))
    
def semivar2directional(semivar):
    """
    Three directional semivariance functions: one of longitude (in degrees),
    one of latitude (in degrees), one of time (in years),
    """
    out = ()
    dir_string = ['longitude: degrees', 'latitude: degrees', 'time: years']
    for i in xrange(3):
        def dirsem(x,semivar=semivar,i=i):
            zerovec = 0.*x
            delta = np.vstack((zerovec,)*(i)+(x,)+(zerovec,)*(2-i)).T
            return semivar(delta)
        dirsem.__doc__ =     """
    A one-directional semivariance function. This function takes one argument: 
    an array of dimension (n,3) indicating lag in %s.
    """%dir_string[i] + semivar.covdoc
        out = out + (dirsem,)
    return out
        
if __name__ == '__main__':
    hf = tb.openFile('/Users/anand/renearch/mbg-world/mbgw-scripts/test_sim_100000000.000000.hdf5')
    
    acf = realization2acf('/Users/anand/renearch/mbg-world/mbgw-scripts/test_sim_100000000.000000.hdf5',0)
    semivar = realization2semivariance('/Users/anand/renearch/mbg-world/mbgw-scripts/test_sim_100000000.000000.hdf5',0)
    dirsems = semivar2directional(semivar)
    x=np.array([[0,0,0],
        [1e-6,1e-6,1e-6],
        [1,1,1],
        [2,2,2],
        [1e6,1e6,1e6]])
        
    print 'Autocovariance: ', acf(x)
    print 'Semivariance: ',semivar(x)
    print 'Directional semivariances: ', np.array([dirsems[i](x[:,i]) for i in xrange(3)])