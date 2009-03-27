from numpy import *
from numpy.testing import *
import mbgw
from mbgw.EP import frontend_interface as fi
from pymc import *
from copy import copy, deepcopy
import tables as tb
import nose

spatial = True
with_urban=True
with_stukel=False
chunk=2
covariates=False

disttol=10/6378.
ttol=1/12.

N_param_vals = 500
N_per_param = 500
N_nearest = 40
N_age_samps = 1000

trace_thin = 10
trace_burn = 1000

# Store an array of age-correction factors for future use
try:
    hf.close()
except:
    pass

root = mbgw.__path__[0]
tracefile = tb.openFile(root+'/../datafiles/good-traces/QRYPFPR101108Ken_KenyaThin_Run_11.2.2009_urb_periurb.hdf5')

rad_to_km = 6378.1/pi
km_to_rad = 1./rad_to_km
rad_to_deg = 180./pi
deg_to_rad = 1./rad_to_deg

n_digits = 4
iter = 10000
tol = .0001

# class test_EP_frontend(TestCase):
class test_EP_frontend(object):
        
    def test_robust_observe_and_eval(self):
        """Makes sure that robust_observe_and_eval works with:
            - The standard Covariance object
            - FullRankCovariance, which will fail 'natively' and be observed 'manually' within the function
            - FullRankCovariance with small nugget, which will succeed
            - FullRankCovariance with small negative nugget.
            
        All these results should agree."""
        fcov = gp.cov_funs.gaussian.euclidean
        M1 = gp.Mean(lambda x: 0.*x)
        C1 = gp.Covariance(fcov, amp=1, scale=1)
        
        M2 = gp.Mean(lambda x: 0.*x)
        C2 = gp.FullRankCovariance(fcov, amp=1, scale=1)
        
        M3 = gp.Mean(lambda x: 0.*x)
        C3 = gp.FullRankCovariance(fcov, amp=1, scale=1)
        
        M4 = gp.Mean(lambda x: 0.*x)
        C4 = gp.FullRankCovariance(fcov, amp=1, scale=1)
        
        lm = ones(3)
        sm = arange(3)
        pm = arange(-3,0)
        
        mo1, co1 = fi.robust_observe_and_eval(lm, M1, C1, sm, 0., pm)
        mo2, co2 = fi.robust_observe_and_eval(lm, M2, C2, sm, .0000001, pm)
        mo3, co3 = fi.robust_observe_and_eval(lm, M3, C3, sm, 0., pm)
        mo4, co4 = fi.robust_observe_and_eval(lm, M4, C4, sm, -.0000001, pm)
        
        for i in xrange(3):
            for j in xrange(i+1,4):
                exec('assert_almost_equal(mo%i,mo%i)'%(i+1,j+1))
                exec('assert_almost_equal(co%i,co%i)'%(i+1,j+1))
    
    def test_simulate_on_pred(self):
        """Makes sure the simulations, sans age correction, are distributed as they should be."""
        nsamp = 100000
        Vp = [.5]
        M_mesh = [1]
        nmonths = [1]
        fac_array = ones((1000,2))
        lo_age = [0]
        up_age = [1]

        samps = fi.simulate_on_pred(nsamp, Vp, M_mesh, nmonths, fac_array, lo_age, up_age).squeeze()

        obs_dist, bin_edges = histogram(samps, 50, normed=True)
        x = (bin_edges[:-1] + bin_edges[1:])*.5
        theor_dist = 1./sqrt(2.*pi*Vp[0])/x/(1-x) * exp(-1./2/Vp[0]*(log(x/(1-x)) - M_mesh[0])**2)
        
        delta = (obs_dist - theor_dist)/obs_dist.max()
        n_differences = sum(abs(delta)>.005)
        assert(n_differences<25)
        assert(abs(delta).max()<.05)
        
    def test_resample(self):
        """Makes sure that the resampling procedure produces samples from the posterior."""
        mu=3.
        V=1.
        
        obs_mu = 2.
        obs_V = .05
        
        samps = random.normal(size=100000)*sqrt(V) + mu
        mps = array([normal_like(x,obs_mu,1./obs_V) for x in samps])
        ui, nc = fi.resample_with_mp(mps)
        
        resamps = []
        for i, n in zip(ui,nc):
            resamps.extend([samps[i]]*n)
        resamps = array(resamps)    
        
        new_samps = random.normal(size=100000) * sqrt(V*obs_V/(V+obs_V)) + (mu*obs_V + obs_mu*V)/(V+obs_V)
        
        rehist, reedge = histogram(resamps,50,normed=True)
        newhist, newedge = histogram(new_samps,reedge,normed=True)
        x=(reedge[1:]+reedge[:-1])*.5
    
        delta = (rehist-newhist)/rehist.max()
        n_differences = sum(abs(delta)>.005)
        
        assert(n_differences<25)
        assert(abs(delta).max()<.05)
        
    def test_inputs(self):
        """Makes sure the input dictionaries are wrangled correctly."""
        
        ip = [{'lon': 1,
                'lat': 3,
                'month': 3,
                'year': 1985,
                'lo_age': 5,
                'up_age': 15, 
                'n': 1000},
                {'lon': -1,
                'lat': -3,
                'month': 0,
                'year': 0,
                'lo_age': 0,
                'up_age': 100, 
                'n': 1e6}]
                
        op = [{'lon': 2,
                'lat': 4,
                'month': 9,
                'year': 2009,
                'lo_age': 2,
                'up_age': 10,
                'nmonths': 3},
                {'lon': 5,
                'lat': 7,
                'month': 5,
                'year': 2000,
                'lo_age': 1,
                'up_age': 8,
                'nmonths': 2}]
        
        din = fi.dicts2rec(ip)
        dout =  fi.dicts2rec(op)        

        pred_mesh, samp_mesh = fi.ra_to_mesh(din, dout)
        
        correct_samp_mesh = array([[1.*deg_to_rad, 3.*deg_to_rad, 1985+3./12-2009],
                                    [-1*deg_to_rad, -3*deg_to_rad, -2009.]])
                                    
        correct_pred_mesh = array([[2*deg_to_rad, 4*deg_to_rad, 2009+9./12-2009],
                                    [2*deg_to_rad, 4*deg_to_rad, 2009+10./12-2009],
                                    [2*deg_to_rad, 4*deg_to_rad, 2009+11./12-2009],
                                    [5*deg_to_rad, 7*deg_to_rad, 2000+5./12-2009],
                                    [5*deg_to_rad, 7*deg_to_rad, 2000+6./12-2009]])
                                    
        assert_equal(correct_samp_mesh, samp_mesh)
        assert_almost_equal(correct_pred_mesh, pred_mesh,10)
        
    def test_dryrun(self):
        """Does a dry run, simulating input data from the user."""    
        
        Np = 4    
        Ns = 5
        
        N_exam = array([2]*Ns)
        
        alpl = random.randint(1,12,size=Np)
        alph = alpl+random.randint(30, size=Np)
    
        alsl = random.randint(1,12,size=Ns)
        alsh = alsl+random.randint(30, size=Ns)
        
        latp = np.atleast_1d(runiform(-5., 5., size=Np))
        lonp = np.atleast_1d(runiform(33., 40., size=Np))
        yrp = np.atleast_1d(array([2007]*Np))
        mop = np.atleast_1d(random.randint(0,12,size=Np))
        
        lats = np.atleast_1d(runiform(-5., 5., size=Ns))
        lons = np.atleast_1d(runiform(33., 40., size=Ns))
        yrs = np.atleast_1d(array([2007]*Ns))
        mos = np.atleast_1d(random.randint(0,12,size=Ns))
    
        ip = [{ 'lon': lons[i], 
                'lat': lats[i], 
                'month': mos[i], 
                'year': yrs[i], 
                'lo_age': alsl[i], 
                'up_age': alsh[i], 
                'n': N_exam[i]} for i in range(Ns)]
        op = [{ 'lon': lonp[i], 
                'lat': latp[i], 
                'year': yrp[i], 
                'month': mop[i], 
                'lo_age': alpl[i], 
                'up_age': alph[i], 
                'nmonths': 4} for i in range(Np)]
    
        fi.update_posterior(ip, op, tracefile, trace_thin, trace_burn, N_param_vals, N_per_param, N_nearest, utilities=[np.std, np.mean])
    
if __name__ == '__main__':
    # nose.runmodule()
    tester = test_EP_frontend()
    # tester.test_simulate_on_pred()
    # tester.test_robust_observe_and_eval()
    # tester.test_resample()
    tester.test_dryrun()
    # tester.test_inputs()