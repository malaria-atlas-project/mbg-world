



NB: THIS FUNCTION IS NOW HOUSED IN supervisor.py

import tables as tb
import numpy as np
import time
import copy as cp
from rpy import *

r.source('/home/pwg/mbg-world/mbgw/joint_simulation/CONDSIMalgorithm/MVRNORM.R')
mvrnormPY=r['MVRNORM']

def predictPointsFromBlock(XYT_in,z_in, XYT_out,C,VERBOSE=FALSE):

    '''
    params to pass:
    XYT_in    : (2d numpy array)   three column array housing absolute x,y,t locations of thinned sample from unconditoned block
    XYT_out   : (2d numpy array)   as above but for data locations we are predicting to outside of block
    z_in      : (1d numpy array)   values of unconditoied block at XYZ_in locations
    C         : (method)           covariance function method obtained from mcmc output hdf5 file (hf.root.group0.C)

    returns:
    z_out     : (1d numpy array)    vector of simulated values
    ''' 

    # define and check lengths
    n_in = len(XYT_in[:,0])
    if (len(z_in)!=n_in): raise ValueError 'Length of XYT_in ('+str(len(XYT_in[:,0]))') does not match length of z_in ('+str(len(z_in))')'

    n_out = len(XYT_out[:,0])
    MaxToSim=float(n_in)

    # define seperate x,y,t vectors for ease
    x_in = XYT_in[:,0]
    y_in = XYT_in[:,1]
    t_in = XYT_in[:,2]

    x_out = XYT_out[:,0]
    y_out = XYT_out[:,1]
    t_out = XYT_out[:,2]

    # define zero-mean vectors
    mean_in = np.zeros(n_in)
    mean_out = np.zeros(n_out)
    
    # initialise vector for output values
    z_out = np.ones(n_out)*-9999

    ## populate full covariance matrices 
    tstart=time.time()
    t1=time.time()
    C_out_out = C((np.vstack((x_out,y_out,t_out)).T),(np.vstack((x_out,y_out,t_out)).T))
    tfor_C_out_out = time.time()-t1
    if VERBOSE: print '\ttfor_C_out_out ('+str(len(x_out))+'x'+str(len(x_out))+' to '+str(len(x_out))+'x'+str(len(x_out))+') was : '+str(tfor_C_out_out)

    t1=time.time()
    C_in_out = C(np.vstack((x_in,y_in,t_in)).T,np.vstack((x_out,y_out,t_out)).T)
    tfor_C_in_out = time.time()-t1
    if VERBOSE: print '\ttfor_C_in_out ('+str(len(x_in))+'x'+str(len(x_in))+' to '+str(len(x_out))+'x'+str(len(x_out))+') was : '+str(tfor_C_in_out)

    t1=time.time()
    C_in_in = C(np.vstack((x_in,y_in,t_in)).T,np.vstack((x_in,y_in,t_in)).T)
    tfor_C_in_in = time.time()-t1
    if VERBOSE: print '\ttfor_C_in_in ('+str(len(x_in))+'x'+str(len(x_in))+' to '+str(len(x_in))+'x'+str(len(x_in))+') was : '+str(tfor_C_in_in)

    # perform matrix operations
    t1=time.time()
    C_in_in.inv=np.linalg.inv(C_in_in)
    if VERBOSE: print '\ttime for inverting C_in_in (shape:'+str(np.shape(C_in_in))+') was : '+str(time.time()-t1)

    t1=time.time()
    PostMeanInterim=np.linalg.linalg.dot(C_in_out.T,C_in_in.inv)
    PostVar=C_out_out - np.linalg.linalg.dot((np.linalg.linalg.dot(C_in_out.T,C_in_in.inv)),C_in_out)
    PostMean=mean_out + np.linalg.linalg.dot(PostMeanInterim , (z_in-mean_in) )
    if VERBOSE: print '\ttime for remaining matrix operations was :'+str(time.time()-t1)

    t1=time.time()
    z_out=mvrnormPY(1,MU=PostMean,COV=PostVar)
    if VERBOSE: print '\ttime for joint simulation over '+str(len(sim))+' points was :'+str(time.time()-t1)

    if verbose: print 'Total time was '+str(tstart-time.time())

    # check z_out contains no unsimulated values (-9999)
    if(np.any(z_out==-9999)): raise Warning str(np.sum(z_out==-9999))+'unpredicted values found in z_out'

    # return 1d array of simulated values    
    return z_out

    ###### THIS SECTION IS REDUNDANT: IMPLEMENTS SEQUENTIAL SIMULATION OF SUBSETS OF PREDICTION LOCATIONS ########################### 
    # define sets of prediction locations to predict
    #
    ## if doing all predictions in one simulation, define a single set of predictions
    #if (len(x_out) <= MaxToSim):
    #    predSetList = []
    #    predSetList.append(np.arange(len(x_out)))
    #
    #
    ## if MaxToSim is less than total amount of prediction locations, then we need to subset the prediction into groups of MaxToSim
    #if (len(x_out) > MaxToSim):
    #
    #    # divide set of prediction locations into approximately equal sized groups, none exceeding MaxToSim
    #    nsetstemp=int(np.ceil(len(x_out)/MaxToSim))
    #    setsize=int(len(x_out)/nsetstemp)
    #    predlist = np.arange(len(x_out))
    #    np.random.shuffle(predlist)
    #    
    #    predSetList = []
    #    for i in arange(nsetstemp):
    #        set = predlist[(i*setsize):((i*setsize)+setsize):1]
    #        predSetList.append(set)
    #    if ((i*setsize)+setsize)< len(predlist):
    #        predSetList.append(predlist[((i*setsize)+setsize)::])
    #
    #    nsets = len(predSetList)
    #
    #    # loop through prediction sets and generate joint realisation for each
    #    for i in arange(nsets):
    #
    #        print 'On prediction set '+str(i)+' of '+str(nsets)
    #
    #        # define subset index 
    #        predIndex = predSetList[i]
    #    
    #        # define new subsetted covariance matrices incorporating this prediction set
    #        temp= C_out_out[predIndex,:]
    #        C_out_out_subset = temp[:,predIndex]
    #        C_in_out_subset = C_in_out[:,predIndex]
    #        mean_out_subset = mean_out[predIndex]
    #
    #        # perform matrix operations
    #        t1=time.time()
    #        C_in_in.inv=np.linalg.inv(C_in_in)
    #        print '\ttime for inverting C_in_in (shape:'+str(np.shape(C_in_in))+') was : '+str(time.time()-t1)
    #        t1=time.time()
    #        PostMeanInterim=np.linalg.linalg.dot(C_in_out_subset.T,C_in_in.inv)
    #        PostVar=C_out_out_subset - np.linalg.linalg.dot((np.linalg.linalg.dot(C_in_out_subset.T,C_in_in.inv)),C_in_out_subset)
    #        PostMean=mean_out_subset + np.linalg.linalg.dot(PostMeanInterim , (z_in-mean_in) )
    #        print '\ttime for remaining matrix operations was :'+str(time.time()-t1)
    #        t1=time.time()
    #        sim=mvrnormPY(1,MU=PostMean,COV=PostVar)
    #        print '\ttime for joint simulation over '+str(len(sim))+' points was :'+str(time.time()-t1)
    #        z_out[predIndex] = sim
    #
    #
    #        # update matrices and data vectors to incorporate these prediction as data in next round of predictions 
    #        
    #        ## predictions now added to data vector
    #        z_in = np.concatenate((z_in,z_out[predIndex]))
    #        mean_in = np.concatenate((mean_in,mean_out[predIndex]))
    #
    #        ## data-to-data covariance appended with previous data-to-prediction and pred-to-pred
    #        temp1 = hstack((C_in_in,C_in_out_subset))
    #        temp2 = hstack((C_in_out_subset.T,C_out_out_subset))
    #        C_in_in= vstack((temp1,temp2))
    #        
    #        ## data-to-prediction covariance appended with previous data-to-prediction 
    #        C_in_out = vstack((C_in_out,C_out_out[predIndex,:]))
    #####################################################################################################################################

