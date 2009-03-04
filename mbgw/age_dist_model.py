# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

# b has something to do with force of infection, also acquisition of immunity.
# Need to be faster than 10 yrs, slower than 1 yr. Between 1 and 5 years.
# 
# The bigger you get, the more you get bitten.
# 
# c can't asymptote any faster than about 10 years. No slower than 50. Probably 
# 
# Only thing different btween rdt and micro is alpha changes.
#
# All these records are probably microscopy but check with Carlos. Dave has some papers from Bob
# with age-stratified RDT records that he can send you.


import pymc as pm
from numpy import *
from agepr import *

N_pops = len(datasets)
S_guess = exp(-a/20.)
S_guess /= sum(S_guess)

# Model for central tendency of age distribution.
alph_pri = 10.

scg = []
for i in xrange(len(a)):
    asg_now = pm.lam_dtrm('asg_%i'%i, lambda alpha=alph_pri, S=S_guess, i=i : alpha*S[i])
    scg.append(pm.Gamma('scg_%i'%i, alpha=asg_now, beta=1., value=asg_now.value))
sc = pm.lam_dtrm('sc',lambda gams=scg: array(gams)/sum(gams))

# Model for age distributions of individual populations
# alph = pm.Gamma('alph',2.,2./1000.)
alph = pm.OneOverX('alph', value=50.)

asc = pm.lam_dtrm('asc',lambda alph=alph,S=sc: alph*S)

S = []
S_steps = []
data_list = []
for name in datasets.iterkeys():
    
    this_dataset = datasets[name]
    
    @pm.deterministic
    def this_asc_slice(asc=asc, slices = age_slices[name], dataset=this_dataset):
        out = []
        for i in xrange(len(dataset)):
            out.append(sum(asc[slices[i]]))
        out = array(out)
        return out
    
    S_now = pm.Dirichlet(('S_%s'%name).replace('.','_'), this_asc_slice)
    S.append(S_now)
    data_list.append(pm.Multinomial('data',n=sum(this_dataset.N),p=S_now,value=this_dataset.N,isdata=True))
    
S_pred = pm.Dirichlet('S_pred',asc)
    
M = pm.MCMC({'variables': [sc, scg, alph, asc, S, data_list, S_pred],
            '__name__': 'age_dist_model',
            'step_methods': S_steps},
            db='hdf5',comp_level = 5)

M.use_step_method(pm.Metropolis, alph, sig=.05)

for i in xrange(len(datasets)):    
    M.use_step_method(pm.DirichletMultinomial, S[i])
