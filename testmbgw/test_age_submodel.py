# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

from MAPdata import csv2recarray
from pylab import *
from tables import *
from numpy import *
# from Training import *
from agepr import *

# a_lo = range(15) + range(15,70,5)
# a_hi = range(15) + range(19,74,5)

__root__ = mbgw.__path__[0] + '/../datafiles/auxiliary_data/'

hfp = openFile(__root__+'parameter_model.hdf5')
proot = hfp.root.chain1.PyMCsamples
    
F_trace, P_trace = proot.col('F_pred'), proot.col('P_pred')

hfa = openFile(__root__+'age_dist_model.hdf5')
aroot = hfa.root.chain1.PyMCsamples

S_trace = aroot.col('S_pred')

T = csv2recarray('Combined_Testing.csv')
T = T[where(T.PF2>0)]
T = T[:len(T)/2]

low_ages = zip((T.L1, T.L2))
up_ages = zip((T.U1, T.U2))

figure(1)
clf()
# plot(T.PR1, T.PR2, 'ko', markersize=4)
"L1","U1","N1","PF1","L2","U2","N2","PF2"
obs_ratios = (array(T.PF1,dtype=float) * T.N2 / T.PF2 / T.N1)
hist(obs_ratios,50)
title('Observed PfPR pairs')

def plot_empirical_hist(a, label, style='k-'):
    x_plot = np.linspace(log(.5),log(10),500)
    y_plot = np.empty(500)
    for i in xrange(500):
        y_plot[i] = np.sum(log(a) < x_plot[i])/float(len(a))
    plot(x_plot, y_plot, style, label=label)

p_indices = random.randint(P_trace.shape[0], size=5000)
S_indices = random.randint(S_trace.shape[0], size=5000)

pred_pf_1 = []
pred_pf_2 = []
for i in xrange(5000):
    
    P_now = P_trace[p_indices[i]]
    S_now = S_trace[S_indices[i],0,:]
    
    j = random.randint(len(T))
    
    a_index_min = np.where(a<=T.L1[j])[0][-1]
    a_index_max = np.where(a>=T.U1[j])[0]
    if len(a_index_max)==0:
        a_index_max = len(a)-1
    else:
        a_index_max = a_index_max[0]
        
    PfPR_1_now = sum(P_now[a_index_min:a_index_max] * S_now[a_index_min:a_index_max])\
        / sum(S_now[a_index_min:a_index_max]) 
    
    PfPR_1_now = np.random.binomial(T.N1[j], PfPR_1_now) / float(T.N1[j])

    a_index_min = np.where(a<=T.L2[j])[0][-1]
    a_index_max = np.where(a>=T.U2[j])[0]
    if len(a_index_max)==0:
        a_index_max = len(a)-1
    else:
        a_index_max = a_index_max[0]        

    PfPR_2_now = sum(P_now[a_index_min:a_index_max] * S_now[a_index_min:a_index_max]) \
        / sum(S_now[a_index_min:a_index_max])

    PfPR_2_now = np.random.binomial(T.N2[j], PfPR_2_now) / float(T.N2[j])
    
    pred_pf_1.append(PfPR_1_now)
    pred_pf_2.append(PfPR_2_now)        
    
figure(2)
clf()
hist(array(pred_pf_1)/array(pred_pf_2),50)
pred_ratios = (array(pred_pf_1)/array(pred_pf_2))
title('Predicted dataset')
figure()
plot_empirical_hist(pred_ratios, 'Predicted', 'k-.')
plot_empirical_hist(obs_ratios, 'Observed', 'k-')
legend(loc=0)
xlabel('Log of ratio of PR values')
ylabel('Probability')
title('Empirical CDF of testing set and posterior predictive CDF')
savefig('cdf_comparison.pdf')
