# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

import pylab as pl
import pymc as pm
from numpy import *
from agepr import *

def binom_deviance(n,y,p):
    y_hat = n*p
    return 2.*sum((y*log(y/y_hat))[where(y>0)]) + 2.*sum(((n-y) * log(n-y) / (n-y_hat))[where(y<n)])

# c: .02-1.
# alph: 18-30
# b: .02-1.

# alpha = array([10., 10., 10., 3.])
sig_mean = array([.5, .4, .5, 1, .4, .5, 1])

# sigma = pm.Gamma('sigma', alpha, alpha/sig_mean)
sigma = pm.OneOverX('sigma', value=sig_mean)

p_mean_mu = array([1.3, 3.7, 2.3, 0, 3.7, 2.3, 0])
# p_mean = pm.MvNormal('p_mean', p_mean_mu, diag([10., 10., 10., 1.]))
p_mean = pm.Uninformative('p_mean',value=p_mean_mu)

R1 = pm.Uninformative('R1', zeros(6,dtype=float))
R2 = pm.Uninformative('R2', zeros(3,dtype=float))
R3 = pm.Uninformative('R3', zeros(3,dtype=float))

# For debugging
# R1.value = arange(1,7)
# R2.value = arange(7,10)
# R3.value = arange(10,13)

@pm.deterministic
def cholfac(R1=R1, R2=R2, R3=R3, sigma = sigma):
    """Cholesky factor of the covariance matrix."""

    cov = np.zeros((7,7),dtype=float)

    cov[0,0] = 1
    cov[1:,0] = R1
    cov[1:4,1:4] = pm.flib.expand_triangular(ones(3), R2)
    cov[4:7,4:7] = pm.flib.expand_triangular(ones(3), R3)    
    
    try:
        cholfac = asmatrix(diag(sigma)) * linalg.cholesky(cov)
    except linalg.LinAlgError:
        raise pm.ZeroProbability
    return cholfac
    
Fs = []
Ps= []
data_list = []
pred_data_list = []
deviance_list = []
pred_deviance_list = []
p_vec_list = []
P_prime_list = []
for name in datasets.iterkeys():
    
    safe_name = name.replace('.','_')
    
    P_prime_now = pm.Beta('P_prime_%s'%safe_name,3.,3.)
    p_vec_now = pm.MvNormalChol('p_vec_%s'%safe_name, p_mean, cholfac)
    p_vec_list.append(p_vec_now)
    P_prime_list.append(P_prime_now)
    
    b = pm.lam_dtrm('b', lambda p_vec = p_vec_now: 1./exp(p_vec[0]))
    
    if methods[name] == 'Microscopy':        
        
        # alpha, s and c depend on p_vec[1:4]
        c = pm.lam_dtrm('c', lambda p_vec = p_vec_now: 1./exp(p_vec[1]))
        alph = pm.lam_dtrm('alph', lambda p_vec = p_vec_now: exp(p_vec[2]))
        s = pm.lam_dtrm('s', lambda p_vec = p_vec_now: pm.invlogit(p_vec[3]))
    
    elif methods[name] == 'RDT':
        
        # alpha, s and c depend on p_vec[4:7]
        c = pm.lam_dtrm('c', lambda p_vec = p_vec_now: 1./exp(p_vec[4]))
        alph = pm.lam_dtrm('alph', lambda p_vec = p_vec_now: exp(p_vec[5]))
        s = pm.lam_dtrm('s', lambda p_vec = p_vec_now: pm.invlogit(p_vec[6]))        
        
    @pm.dtrm
    def this_F(c=c, alph=alph, a=age_bin_ctrs[name], s=s):
        """
        The function F, which gives detection probability.
        """
        out = empty(len(a))
        out[where(a<alph)] = 1.
        where_greater = where(a>=alph)
        out[where_greater] = (1.-s*(1.-exp(-c*(a-alph))))[where_greater]
        return out
    this_F.__name__ = 'F_%s'%safe_name
        
    
    @pm.dtrm
    def this_P(P_prime=P_prime_now, b=b, a=age_bin_ctrs[name], F=this_F):
        """
        The function P, which gives probability of a detected infection.
        """
        return P_prime * (1.-exp(-b*a)) * F
    this_P.__name__ = 'P_%s'%safe_name
    
    Fs.append(this_F)
    Ps.append(this_P)

    #Data
    this_data = pm.Binomial('data_%s'%safe_name, n=datasets[name].N, p=this_P, value=datasets[name].pos, isdata=True)
    data_list.append(this_data)
    
    # # Data from sample-space (for model checking)
    # this_pred_data = pm.Binomial('pred_data_%i'%i, n=N[:,i], p=this_P)
    # pred_data_list.append(this_pred_data)
    # 
    # # Deviance
    # this_deviance = pm.lam_dtrm('dev_%i'%i, lambda n=N[:,i], p=this_P, y=this_data: binom_deviance(n,y,p))
    # deviance_list.append(this_deviance)
    # 
    # # Deviance from sample-space (for model checking)
    # this_pred_deviance = pm.lam_dtrm('pred_dev_%i'%i, lambda n=N[:,i], p=this_P, y=this_pred_data: binom_deviance(n,y,p))
    # pred_deviance_list.append(this_pred_deviance)    
    
    # Initialize chain at MAP estimates to shorten burnin.
    M_start = pm.MAP([p_vec_now,P_prime_now,this_data])
    M_start.fit()

# Samples from predictive distribution of parameters.
p_pred = pm.MvNormalChol('p_predictive',p_mean,cholfac)
b_pred = pm.lam_dtrm('b', lambda p_vec = p_pred: 1./exp(p_vec[0]))
c_pred = pm.lam_dtrm('c', lambda p_vec = p_pred: 1./exp(p_vec[1]))
alph_pred = pm.lam_dtrm('alph', lambda p_vec = p_pred: exp(p_vec[2]))
s_pred = pm.lam_dtrm('s', lambda p_vec = p_pred: pm.invlogit(p_vec[3]))

# a_pred = (a[:-1] + a[1:])*.5

@pm.dtrm
def F_pred(c=c_pred, alph=alph_pred, a=a, s=s_pred):
    """
    A sample from the predictive distribution of F.
    """
    out = empty(a.shape[0])
    out[where(a<alph)] = 1.
    where_greater = where(a>=alph)
    out[where_greater] = 1.-s*(1.-exp(-c*(a[where_greater]-alph)))
    return out

@pm.dtrm
def P_pred(P_prime=1., b=b_pred, a=a, F=F_pred):
    """
    A sample from the predictive distribuiton of P.
    """
    return P_prime * (1.-exp(-b*a)) * F

# dev_pvals = empty(N_pops)
# hist_dev = empty(500)
# for j in xrange(N_pops):
#     for i in xrange(500):
#         q=pred_data_list[j].random()
#         hist_dev[i] = pred_deviance_list[j].value
# 
#     dev_pvals[j] = sum(hist_dev > deviance_list[j].value)/500.

all_variables = [sigma, p_mean, R1, R2, R3, cholfac, p_vec_list, P_prime_list, Fs, Ps, p_pred, b_pred, c_pred, alph_pred, s_pred, F_pred, P_pred]

M = pm.MCMC(all_variables, name='parameter_model', db='hdf5', comp_level=5)
M.use_step_method(pm.AdaptiveMetropolis, [p_mean, R1, R2, R3, sigma], 
    scales={p_mean: .01*ones(7), R1: .01*ones(6), R2: .01*ones(3), R3: .01*ones(3), sigma: .01*ones(7)})

for i in xrange(len(datasets)):
    M.use_step_method(pm.AdaptiveMetropolis, [p_vec_list[i],P_prime_list[i]], scales={p_vec_list[i]: .001*ones(7), P_prime_list[i]: [.001]}, delay=10000)