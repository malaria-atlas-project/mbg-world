from numpy.testing import *
import mbgw
from mbgw import EP
import nose,  warnings
from numpy import *
from scipy import integrate

# FIXME: Known failure modes: 
#  - negative variance is too small, so there is no posterior density. You can protect against this. 

n_digits = 4
iter = 10000
tol = .0001

def normal_like(x, mu, t):
    like = - 0.5 * t * (x-mu)**2
    # like = like + 0.5*log(0.5*t/pi)
    return like

def matsqrt(K):
    e,v = linalg.eigh(K)
    out= v*sqrt(e.astype('complex'))
    return out

def geto_observe(M, C, V, val):
    try:
        sig = linalg.cholesky(C+diag(V))
    except:
        sig = matsqrt(C+diag(V))
    C_sig = linalg.solve(sig, C)
    C_post = C - dot(C_sig.T, C_sig)
    M_post = M + dot(linalg.solve(sig.T, C_sig).T, val-M)
    return M_post.astype('float'), C_post.astype('float')

def standard_EP_t(M_pri, C_pri, nugs, obs_mus, obs_Vs, mu_guess=None, V_guess=None):    
    N = len(M_pri)
    lps = [lambda x, m=obs_mus[i], v=obs_Vs[i]: normal_like(x,m,1./v) for i in xrange(N)]

    # Do EP algorithm
    E = EP.EP(M_pri, C_pri, lps, nugs, mu_guess, V_guess)
    E.fit(iter, tol=tol)        
    
    # Independently observe with the inferred 'mean' and 'variance'.
    M_post, C_post = geto_observe(M_pri, C_pri, obs_Vs + nugs, obs_mus)
    
    # Make sure the observing arithmetic is going right.
    assert_almost_equal(M_post, E.M, n_digits)
    assert_almost_equal(C_post/E.C, E.C*0+1., n_digits-2)
    
    # Make sure it's correctly finding mu and V
    assert_almost_equal(E.mu, obs_mus, n_digits)
    assert_almost_equal(E.V*0+1., (obs_Vs+nugs)/(E.V+nugs), n_digits-2)
        
    return E, M_post, C_post
    
# class test_mbgw(TestCase):
class test_EP(TestCase):
    
    N = 3
    M_pri = random.normal(size=N)
    sig_pri = random.normal(size=(N, N))
    C_pri = dot(sig_pri.T, sig_pri)
    
    def test_expectations(self):
        f = lambda x: normal_like(x, -pi, 1./2)
        lo, hi = EP.estimate_envelopes(f, -pi, sqrt(10), 13., ftol=.000000001)
        
        funs=[lambda x:x, lambda x:x**2]
        
        x = linspace(lo, hi, 1000)
        dx=x[2]-x[1]
        post_vec = exp([f(xi) for xi in x])

        p = integrate.simps(post_vec, dx=dx)
        
        # Return E_pri [like_fn(x)] and the posterior expectations of funs(x).        
        moments = []
        for fun in funs:
            moments.append(integrate.simps(post_vec * fun(x), dx=dx) / p)
        
        # print moments[1]-moments[0]**2-2
        assert_almost_equal(moments[0], -pi, 3)
        assert_almost_equal(moments[1]-moments[0]**2, 2, 2)
            
    def test_low_V(self):
        
        # Moderate-sized positive variance and nugget.
        obs_Vs = random.normal(size=self.N)**2 * .1
        nugs = random.normal(size=self.N)**2 * .3
        
        # 'mean' and log-probability functions.
        obs_mus = random.normal(size=self.N)
        
        standard_EP_t(self.M_pri, self.C_pri, nugs, obs_mus, obs_Vs, obs_mus, obs_Vs)
        
    def test_neg_V(self):
        # Moderate-sized negative variance and nugget.
        obs_Vs = -diag(self.C_pri)*3
        nugs = diag(self.C_pri)*.2
        
        # 'mean' and log-probability functions.
        obs_mus = random.normal(size=self.N)
        
        standard_EP_t(self.M_pri, self.C_pri, nugs, obs_mus, obs_Vs)
        
    def test_small_neg_V(self):
        # Negative variance that is smaller than prior variance + nugget.
        obs_Vs = -diag(self.C_pri)
        nugs = diag(self.C_pri)*.2
        
        # 'mean' and log-probability functions.
        obs_mus = random.normal(size=self.N)

        try:
            standard_EP_t(self.M_pri, self.C_pri, nugs, obs_mus, obs_Vs)
            raise AssertionError, 'Improper posterior failed to trigger error in estimate_envelopes'
        except RuntimeError:
            pass
        
    def test_hi_V(self):
        # High variance and nugget.
        obs_Vs = ones(self.N)*100000
        nugs = ones(self.N)*100000
        
        # 'mean' and log-probability functions.
        obs_mus = random.normal(size=self.N)
        
        E, M_post, C_post = standard_EP_t(self.M_pri, self.C_pri, nugs, obs_mus, obs_Vs)
        
        assert_almost_equal(M_post, self.M_pri, n_digits)
        assert_almost_equal(C_post, self.C_pri, n_digits-2)
        
    def test_tiny_V(self):
        # Moderate-sized positive variance and nugget.
        obs_Vs = random.normal(size=self.N)**2 * .00001
        nugs = random.normal(size=self.N)**2 * .00003
        
        # 'mean' and log-probability functions.
        obs_mus = random.normal(size=self.N)
        
        standard_EP_t(self.M_pri, self.C_pri, nugs, obs_mus, obs_Vs)
    

def call_and_check(meth):
    try:
        meth()
    except (RuntimeError, AssertionError):
        import sys
        cls, inst, tb = sys.exc_info()
        print 'Assertion or runtime error in %s: \n%s'%(meth.__name__, inst.message)
        
        
if __name__ == '__main__':
    
    # while True:
    #     print 'Running all tests'
    #     tester = test_EP()
    #     call_and_check(tester.test_expectations)
    #     call_and_check(tester.test_low_V)
    #     call_and_check(tester.test_hi_V)
    #     call_and_check(tester.test_small_neg_V)
    #     call_and_check(tester.test_neg_V)
    #     call_and_check(tester.test_tiny_V)
    nose.runmodule()

