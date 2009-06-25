# Author: Anand Patil
# Date: 6 Feb 2009
# License: Creative Commons BY-NC-SA
####################################

import numpy as np
import pymc as pm
import mbgw
from mbgw.joint_simulation import *
import tables as tb
from fast_krige import *
import st_cov_fun
from parse_and_check import *
import time
import copy as cp
from mbgw.master_grid import *
from mbgw import auxiliary_data
import os
curpath = os.getcwd()
import mbgw
mbgw_root = __root__ = mbgw.__path__[0]
r_path = mbgw_root+'/joint_simulation/CONDSIMalgorithm'
os.chdir(r_path)
from examineRealization import *
from rpy import r
r.source("CONDSIMpreloop.R")
r.source("CONDSIMmonthloop.R")
r.source('MVRNORM.R')
mvrnormPY=r['MVRNORM']

os.chdir(curpath)
import scipy
from scipy import ndimage, mgrid
from map_utils import grid_convert

#####################TEMP
#from map_utils import getEmpiricalCovarianceFunction_STmarginals
#from map_utils import plotEmpiricalCovarianceFunction
#####################TEMP

from IPython.Debugger import Pdb

__all__ = ['create_realization', 'create_many_realizations','reduce_realizations','getThinnedBlockXYTZlists','array3d_2_XYTZlist','gridParams_2_XYTmarginallists','predictPointsFromBlock']

def get_covariate_submesh(name, grid_lims):
    """
    Matches the specified spatial mesh to the 'master' covariate mesh.
    """

    try:
        order = getattr(mbgw.auxiliary_data, name)._v_attrs.order
    except:
        order = 'y-x+'
    
    raw_shape = getattr(mbgw.auxiliary_data, name).data.shape
    raw = getattr(mbgw.auxiliary_data, name).data[grid_lims['topRow']:grid_lims['bottomRow']+1, grid_lims['leftCol']:grid_lims['rightCol']+1]
    out = grid_convert(raw, order, 'x+y+').copy()
    targ_shape = (grid_lims['rightCol']-grid_lims['leftCol']+1, grid_lims['bottomRow']-grid_lims['topRow']+1)

    if out.shape != targ_shape:
        raise ValueError, "Grid %s's shape is %s with order %s. Cannot be sliced to shape %s. \n\tGrid limits: %s\n\tnrows: %i\n\tncols: %i"\
        %(name, raw_shape, order, targ_shape, grid_lims, nrows, ncols)

    return out

def create_many_realizations(burn, n, trace, meta, grid_lims, start_year, nmonths, outfile_name, memmax, relp=1e-3, mask_name=None, n_in_trace=None, thinning=10,paramfileINDEX=0,NinThinnedBlock=0):
    """
    Creates N realizations from the predictive distribution over the specified space-time mesh.
    """
    
    
    
    # Establish grids
    xllc_here = (xllc + cellsize * (grid_lims['leftCol']-1))*deg_to_rad
    xurc_here = (xllc + cellsize * (grid_lims['rightCol']-1))*deg_to_rad
    yllc_here = (yllc + cellsize * (nrows-grid_lims['bottomRow']))*deg_to_rad
    yurc_here = (yllc + cellsize * (nrows-grid_lims['topRow']))*deg_to_rad    
    grids = [(xllc_here, xurc_here, grid_lims['rightCol']-grid_lims['leftCol']+1),
            (yllc_here, yurc_here, grid_lims['bottomRow']-grid_lims['topRow']+1),
            (start_year-2009., start_year-2009+(nmonths-1)/12., nmonths)]
    axes = [np.linspace(*grids[i]) for i in xrange(3)]
    grid_shape = (grids[0][2], grids[1][2], grids[2][2])
    
    if mask_name is not None:
        mask = get_covariate_submesh(mask_name, grid_lims)
    else:
        mask = np.ones(grid_shape[:2])
        
    if not mask.shape == grid_shape[:2]:
        raise ValueError, 'You screwed up the shapes.'

    # Check that all data are in bounds
    data_locs = meta.logp_mesh[:]    
    in_mesh = np.ones(data_locs.shape[0],dtype=bool)
    for i in xrange(len(data_locs)):
        l = data_locs[i]
        for j in xrange(3):
            if l[j] <= grids[j][0] or l[j] >= grids[j][1]:
                in_mesh[i]=False

    print '****',np.sum(in_mesh)
    # from pylab import plot,clf,show
    # 
    # print np.asarray(grids)
    # 
    # plot(data_locs[:,0], data_locs[:,1], 'k.')
    # dl = data_locs[np.where((data_locs[:,-1]<=grids[-1][1])*(data_locs[:,-1]>=grids[-1][0]))]
    # plot(dl[:,0],dl[:,1],'r.')
    # 
    # plot(grids[0][:2],grids[1][:2],'b.',markersize=16)
    # plot(grids[0][:2],grids[1][:2],'b.',markersize=16)
    # # show()
    # from IPython.Debugger import Pdb
    # Pdb(color_scheme='Linux').set_trace()   

    # Find the mesh indices closest to the data locations
    data_mesh_indices = np.empty(data_locs.shape, dtype=np.int)

    for i in xrange(len(data_locs)):
        for j in xrange(3):


            ####??????????
            data_mesh_indices[i,j] = np.argmin(np.abs(data_locs[i,j] - axes[j]))
            ####??????????


    if n_in_trace is None:
        n_in_trace=len(trace.group0.C) 
    spacing = (n_in_trace-burn)/n
    indices = np.arange(burn, n_in_trace, spacing)
    N = len(indices)    
    
    outfile = tb.openFile(outfile_name, 'w')
    outfile.createArray('/','lon_axis',axes[0],title='Longitude in radians')
    outfile.createArray('/','lat_axis',axes[1],title='Longitude in radians')
    outfile.createArray('/','t_axis',axes[2],title='Time in years since 2009')        
    outfile.createCArray('/','realizations',
                        tb.Float32Atom(), 
                        shape=(N,grid_shape[1],grid_shape[0],grid_shape[2]), 
                        filters=tb.Filters(complevel=1), 
                        chunkshape = (1,grid_shape[0],grid_shape[1],1))
    
    #Store information to help access the parameter samples that generated the realizations.
    outfile.createArray('/','indices',indices,
        title='Indices in the trace file that correspond to the realizations here.')
    new_table = outfile.createTable('/','PyMCsamples',
                        description=trace.PyMCsamples.description,
                        expectedrows=len(indices),
                        title='Trace of numeric-valued variables in model, thinned to be relevant to realizations')
    new_table.append(trace.PyMCsamples[slice(burn,n_in_trace,spacing)])
    outfile.createGroup('/','group0',title='Trace of object-valued variables in model, thinned to be relevant to realizations')
    for node in trace.group0._f_iterNodes():
        new_node = outfile.createVLArray('/group0',node.name,tb.ObjectAtom())
        [new_node.append(node[index]) for index in indices]
    outfile.root._v_attrs.orig_filename = trace._v_file.filename
    

    # Total number of pixels in month.
    npix = grid_shape[0]*grid_shape[1]/thinning**2
    # Maximum number of pixels in tile.
    npixmax = memmax/4./data_locs.shape[0]
    # Minimum number of tiles needed.
    ntiles = npix/npixmax
    # Blocks.
    n_blocks_x = n_blocks_y = np.ceil(np.sqrt(ntiles))

    print 'I can afford %i by %i'%(n_blocks_x, n_blocks_y)

    
    # Scatter this part to many processes
    for i in xrange(len(indices)):
        print 'Realization %i of %i'%(i,N)
        
        # Pull mean information out of trace
        this_M = trace.group0.M[indices[i]]
        mean_ondata = this_M(data_locs)
        covariate_mesh = np.zeros(grid_shape[:2])
        for key in meta.covariate_names[0]:
            try:
                this_coef = trace.PyMCsamples.col(key+'_coef')[indices[i]]
            except KeyError:
                print 'Warning, no column named %s'%key+'_coef'
                continue
            mean_ondata += getattr(meta, key)[:][meta.ui[:]] * this_coef
            this_pred_covariate = get_covariate_submesh(key+'5km-e_y-x+', grid_lims) * this_coef
            covariate_mesh += this_pred_covariate        

        # Pull covariance information out of trace
        this_C = trace.group0.C[indices[i]]
        this_C = pm.gp.NearlyFullRankCovariance(this_C.eval_fun, relative_precision=relp, **this_C.params)

        #data_vals = trace.PyMCsamples[i]['f'][in_mesh]
        #create_realization(outfile.root.realizations, i, this_C, mean_ondata, this_M, covariate_mesh, data_vals, data_locs, grids, axes, data_mesh_indices, n_blocks_x, n_blocks_y, relp, mask, thinning, indices)

        data_vals = trace.PyMCsamples[i]['f'][:]
        create_realization(outfile.root, i, this_C,trace.group0.C[indices[i]], mean_ondata, this_M, covariate_mesh, data_vals, data_locs, grids, axes, data_mesh_indices, np.where(in_mesh)[0], np.where(True-in_mesh)[0], n_blocks_x, n_blocks_y, relp, mask, thinning,indices,paramfileINDEX,NinThinnedBlock)
        outfile.flush()
    outfile.close()

def normalize_for_mapcoords(arr, max):
    arr /= arr.max()
    arr *= max

# def create_realization(out_arr,real_index, C, mean_ondata, M, covariate_mesh, tdata, data_locs, grids, axes, data_mesh_indices, n_blocks_x, n_blocks_y, relp, mask, thinning, indices):
def create_realization(outfile_root,real_index, C,C_straighfromtrace, mean_ondata, M, covariate_mesh, tdata, data_locs, grids, axes, data_mesh_indices, where_in, where_out, n_blocks_x, n_blocks_y, relp, mask, thinning,indices,paramfileINDEX,NinThinnedBlock):

    # define only realizations chunk of hdf5 realization file
    out_arr = outfile_root.realizations


    """
    Creates a single realization from the predictive distribution over specified space-time mesh.
    """
    grid_shape = tuple([grid[2] for grid in grids])


    #r.X11(width=8,height=4)
    #r.par(mfrow=(1,2))
    #r.plot(data_mesh_indices[:,0],data_mesh_indices[:,1],xlab="",ylab="",main="",cex=0.5)
    #r.plot(data_locs[:,0],data_locs[:,1],xlab="",ylab="",main="",cex=0.5)    
 
    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()



    thin_grids = tuple([grid[:2]+(grid[2]/thinning,) for grid in grids])    
    thin_grid_shape = tuple([thin_grid[2] for thin_grid in thin_grids])
    thin_axes = tuple([np.linspace(*thin_grid) for thin_grid in thin_grids])

    mapgrid = np.array(mgrid[0:grid_shape[0],0:grid_shape[1]], dtype=float)
    for i in xrange(2): normalize_for_mapcoords(mapgrid[i], thin_grid_shape[i]-1)

    thin_mapgrid = np.array(mgrid[0:thin_grid_shape[0], 0:thin_grid_shape[1]], dtype=float)
    for i in xrange(2): normalize_for_mapcoords(thin_mapgrid[i], grid_shape[i]-1)

    def thin_to_full(thin_row):
        return ndimage.map_coordinates(thin_row, mapgrid)
    def full_to_thin(row):
        return ndimage.map_coordinates(row, thin_mapgrid)
    thin_mask = np.array(np.round(full_to_thin(mask)), dtype='bool')

    # Container for x
    thin_x = np.empty(thin_grid_shape[:2] + (3,))
    mlon,mlat = np.meshgrid(*thin_axes[:2])
    thin_x[:,:,0] = mlon.T
    thin_x[:,:,1] = mlat.T
    thin_x[:,:,2] = 0

    x = np.empty(grid_shape[:2] + (3,))
    mlon,mlat = np.meshgrid(*axes[:2])
    x[:,:,0] = mlon.T
    x[:,:,1] = mlat.T
    x[:,:,2] = 0
    
    del mlon, mlat

    Cp = C.params
    
    # Prepare input dictionaries
    covParamObj = {'Scale': Cp['scale'][0],
                    'amp': Cp['amp'][0], 
                    'inc': Cp['inc'][0], 
                    'ecc': Cp['ecc'][0], 
                    't.lim.corr': Cp['tlc'][0], 
                    'scale.t': Cp['st'][0], 
                    'sin.frac': Cp['sf'][0]}
    gridParamObj = {'YLLCORNER': grids[1][0]*rad_to_deg, 
                    'CELLSIZE': (grids[1][1]-grids[1][0])/(grids[1][2]-1.)*rad_to_deg, 
                    'NROWS':grid_shape[1],
                    'NCOLS':grid_shape[0]}
    monthParamObj = {'Nmonths':grid_shape[2],'StartMonth':grids[2][0]}
    
    
    ################################~TEMP
    # Call R preprocessing function and check to make sure no screwy re-casting has taken place.
    #t1 = time.time()
    #os.chdir(r_path)
    #preLoopObj = r.CONDSIMpreloop(covParamObj,gridParamObj,monthParamObj,indices.min(), indices.max(),paramfileINDEX)
    #tree_reader = reader(file('listSummary_preLoopObj_original_%i_%i.txt'%(indices.min(), indices.max())),delimiter=' ')

    #preLoopClassTree, junk = parse_tree(tree_reader)
    #preLoopObj = compare_tree(preLoopObj, preLoopClassTree)
    
    #OutMATlist = preLoopObj['OutMATlist']
    #tree_reader = reader(file('listSummary_OutMATlist_original_%i_%i.txt'%(indices.min(), indices.max())),delimiter=' ')
    #OutMATClassTree, junk = parse_tree(tree_reader)
    #OutMATlist = compare_tree(OutMATlist, OutMATClassTree)
    #os.chdir(curpath)
    #preLoop_time = time.time()-t1
    #print "preLoop_time :"+str(preLoop_time)
        
    ## Create and store unconditional realizations
    #print '\tGenerating unconditional realizations.'
    #t1 = time.time()
    #for i in xrange(grid_shape[2]):
    #    print 'On month :'+str(i)
    #    #print 'OutMATlist:'
    #    #print OutMATlist
    #    os.chdir(r_path)
    #    monthObject = r.CONDSIMmonthloop(i+1,preLoopObj,OutMATlist, indices.min(), indices.max(),paramfileINDEX)
    #    #monthObject = r.CONDSIMmonthloop(i+1,preLoopObj,OutMATlist,paramfileINDEX)
    #    os.chdir(curpath)
    #    OutMATlist= monthObject['OutMATlist']
    #    MonthGrid = monthObject['MonthGrid']
    #    out_arr[real_index,:,:,i] = MonthGrid[:grid_shape[1],:grid_shape[0]]
    #t2 = time.time()
    #print '\t\tDone in %f'%(t2-t1)
    #print "monthloop_time :"+str(t2-t1)+" for "+str(grid_shape[2])+" months" 
    
    # delete unneeded R products
    #del OutMATlist, preLoopObj, MonthGrid, monthObject
    ################################~TEMP

    ################################~TEMP DIRECTLY JOIN SIMULATE UNCODITIONED BLOCK FOR TESTING   
    getUnconditionedBlock(out_arr,real_index,grids,C_straighfromtrace,NinThinnedBlock=None,relp=None,FULLRANK=False)
    print 'variance of unconditioned block = '+str(round(np.var(out_arr[np.isnan(out_arr)==False]),10))
    examineRealization(outfile_root,0,6,15,None,None,conditioned=False,flipVertical="FALSE",SPACE=True,TIME=True)
    ################################~TEMP

    
    # Figure out pdata
    pdata = np.empty(tdata.shape)
    for i in xrange(len(where_in)):
        pdata[where_in[i]] = out_arr[real_index, grid_shape[1]-1-data_mesh_indices[i,1], data_mesh_indices[i,0], data_mesh_indices[i,2]]




    # jointly simulate at data points conditional on block    

    ## first get XYZT list of locations of a regular thinned sample from block

    #array3d = out_arr[real_index,:,:,:]
    ThinnedBlockXYTZlists = getThinnedBlockXYTZlists (out_arr,real_index,grids,NinThinnedBlock)
    xyt_in = ThinnedBlockXYTZlists['xyt_in']
    z_in = ThinnedBlockXYTZlists['z_in']

    ## get locations of data outside block (referenced by grid location)
    #XYT_out_gridlocs = data_mesh_indices[where_out]
    # 
    ### convert these grid lcoations into actual position in radians and months
    #coordsDict = gridParams_2_XYTmarginallists(grids)
    #xcoords = coordsDict['xcoords']
    #ycoords = coordsDict['ycoords']
    #tcoords = coordsDict['tcoords']
    #
    #x_out = xcoords[XYT_out_gridlocs[:,0]]
    #y_out = ycoords[XYT_out_gridlocs[:,1]]
    #t_out = tcoords[XYT_out_gridlocs[:,2]]
    #   
    #xyt_out = np.vstack((x_out,y_out,t_out)).T
    
    # get locations of data outside block 
    xyt_out = data_locs[where_out]

    # now we have locations and values of thinned sample from the block, and locations we want to predict at outside the block,go ahead and 
    # get simulated values of the latter, conditonal on the former

    print '\tsimulating over '+str(len(xyt_out[:,0]))+' locations outside block using thinned block sample of '+str(len(z_in))+' points'
    t1=time.time()
    z_out = predictPointsFromBlock(xyt_in,z_in, xyt_out,C_straighfromtrace,relp)
    print '\ttime for simulation: '+str(time.time()-t1)

    #########################################CHECK COVARIANCE STRUCTURE
    #r.X11(width=12,height=12)
    #r.par(mfrow=(3,2))

    # test marginal space and time covariance structures of points outside block
    #cfdict_out = getEmpiricalCovarianceFunction_STmarginals(xyt_out,z_out,mu=0,margTol_S=0.05,margTol_T=0.9/12,nbins=20, cutoff = 0.8)
    #plotEmpiricalCovarianceFunction(cfdict_out['space'],CovModelObj=C_straighfromtrace,spaceORtime="space", cutoff = 0.8, title="Points outside (S) "+str(paramfileINDEX))
    #plotEmpiricalCovarianceFunction(cfdict_out['time'],CovModelObj=C_straighfromtrace,spaceORtime="time", cutoff = 0.8, title="Points outside (T)"+str(paramfileINDEX))

    # test marginal space and time covariance structures of points inside block
    #cfdict_in = getEmpiricalCovarianceFunction_STmarginals(xyt_in,z_in,mu=0,margTol_S=0.05,margTol_T=0.9/12,nbins=20, cutoff = 0.8)
    #plotEmpiricalCovarianceFunction(cfdict_in['space'],CovModelObj=C_straighfromtrace,spaceORtime="space", cutoff = 0.8, title="Points inside (S)"+str(paramfileINDEX))
    #plotEmpiricalCovarianceFunction(cfdict_in['time'],CovModelObj=C_straighfromtrace,spaceORtime="time", cutoff = 0.8, title="Points inside (T)"+str(paramfileINDEX))

    # test marginal space and time covariance structures of points inside and outside block
    #cfdict_inout = getEmpiricalCovarianceFunction_STmarginals(np.vstack((xyt_in,xyt_out)),np.hstack((z_in,z_out)),mu=0,margTol_S=0.05,margTol_T=0.9/12,nbins=20, cutoff = 0.8)
    #plotEmpiricalCovarianceFunction(cfdict_inout['space'],CovModelObj=C_straighfromtrace,spaceORtime="space", cutoff = 0.8, title="Points inside (S)"+str(paramfileINDEX))
    #plotEmpiricalCovarianceFunction(cfdict_inout['time'],CovModelObj=C_straighfromtrace,spaceORtime="time", cutoff = 0.8, title="Points inside (T)"+str(paramfileINDEX))

    #########################################CHECK COVARIANCE STRUCTURE

    

       
    # assign these values to pdata    
    pdata[where_out] = z_out

    ###############################~~TEMP     
    #return()
    #####################################
    
    
    # Bring in data.
    print '\tKriging to bring in data.'    
    print '\tPreprocessing.'
    t1 = time.time()    
    dev_posdef, xbi, ybi, dl_posdef = preprocess(C, data_locs, thin_grids, thin_x, n_blocks_x, n_blocks_y, tdata, pdata, relp, mean_ondata)
    t2 = time.time()
    print '\t\tDone in %f'%(t2-t1)

    ###############################~~TEMP
    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()
    #####################################
    
    thin_row = np.empty(thin_grid_shape[:2], dtype=np.float32)
    print '\tKriging.'
    t1 = time.time()
    for i in xrange(grid_shape[2]-1,-1,-1):
        thin_row.fill(0.)
        
        thin_x[:,:,2] = axes[2][i]
        x[:,:,2] = axes[2][i]
        
        krige_month(C, i, dl_posdef, thin_grid_shape, n_blocks_x, n_blocks_y, xbi, ybi, thin_x, dev_posdef, thin_row, thin_mask)
        row = ndimage.map_coordinates(thin_row, mapgrid)
        
        row += covariate_mesh
        row += M(x)
        row += grid_convert(out_arr[real_index,:,:,i], 'y-x+', 'x+y+')
 
        # NaN the oceans to save storage
        row[np.where(1-mask)] = missing_val
        
        out_arr[real_index,:,:,i] = grid_convert(row, 'x+y+','y-x+')
    
    ####################################TEMP
    print 'variance of conditioned block = '+str(round(np.var(out_arr[np.isnan(out_arr)==False]),10))
    examineRealization(outfile_root,0,6,15,None,None,conditioned=True,flipVertical="FALSE",SPACE=True,TIME=True)
    ########################################
        
            
    print '\t\tDone in %f'%(time.time()-t1)        
        

def reduce_realizations(filename, reduce_fns, slices, a_lo, a_hi, n_per):
    """
    Generates n_per * len(filename.root.realizations) realizations, 
    on the space-time slice defined by slice (a tuple of three slices) 
    and reduces them according to the function reduce. Reduce_fns should 
    be a list of Python functions of the form
    
    reduce(this_PR_chunk, product_sofar=None)
    
    and incorporate this_realization into product_sofar in the desired
    way. It should be robust to the product_sofar=None case, of course.
    a_lo and a_hi are the limits of the age range.
    """
    slices = tuple(slices)
    hf = tb.openFile(filename)
    hr = hf.root
    n_realizations = len(hr.realizations)
    products = dict(zip(reduce_fns, [None]*len(reduce_fns)))
    
    N_facs = int(1e5)
    
    # Get nugget variance and age-correction factors
    V = hr.PyMCsamples.col('V')[:]
    facs = mbgw.correction_factors.age_corr_factors_from_limits(a_lo, a_hi, N_facs)
    
    for i in xrange(n_realizations):
        # Pull out parasite rate chunk
        tot_slice = (slice(i,i+1,1),) + slices
        f_chunk = hr.realizations[tot_slice].squeeze()
        for j in xrange(n_per):
            chunk = f_chunk + np.random.normal(loc=0, scale=np.sqrt(V[i]), size=f_chunk.shape)
            chunk = pm.invlogit(chunk)
            chunk *= facs[np.random.randint(N_facs, size=np.prod(chunk.shape))]
            chunk = chunk.reshape(f_chunk.shape)
            
            for f in reduce_fns:
                product_sofar = products[f]
                products[f] = f(chunk, product_sofar)
    
    return products
    
def getThinnedBlockXYTZlists(relblock4d,real_index,grids,NinThinnedBlock=None):

    # extract grid parameters for ease
    ncols = grids[1][2]
    nrows = grids[0][2]
    nmonths = grids[2][2]

    # if thinning, do dirty calculation to approximately evenly spread sample of size NinThinnedBlock accross ST unconditioned block
    if NinThinnedBlock is not None:
        Nmonthstosample = int(np.ceil((nmonths/12.)*6))
        Tthinrate=nmonths/Nmonthstosample
        Npermonth = NinThinnedBlock/Nmonthstosample
        XYthinRate = int(np.ceil(np.sqrt(nrows*ncols)/np.sqrt(Npermonth)))
        data_footprint=(slice(real_index,real_index+1,1),slice(0,nrows,XYthinRate),slice(0,ncols,XYthinRate),slice(0,nmonths,Tthinrate))

    # if not thinning just define a null slice that will include everything
    if NinThinnedBlock is None:    
        data_footprint=(slice(real_index,real_index+1,1),slice(None,None,None),slice(None,None,None),slice(None,None,None))

    # extract unconditoned values from block at these locations and convert to 3d matrix (remove realisation dimension)
    z_cube = np.squeeze(relblock4d[data_footprint])
    z_cube = np.atleast_3d(z_cube)
    
    # now need to define correposnding long,lat, and time values for these locations
    
    ## first make cubes of lon,lat,time corresponding to z_cube
    ### get grid's marginal coordinates
    coordsDict = gridParams_2_XYTmarginallists(grids)
    xcoords = coordsDict['xcoords']
    ycoords = coordsDict['ycoords']
    tcoords = coordsDict['tcoords']

    ### thinned vectors
    xcoords = xcoords[slice(0,ncols,XYthinRate)]
    ycoords = ycoords[slice(0,nrows,XYthinRate)]
    tcoords = tcoords[slice(0,nmonths,Tthinrate)]

    ### get XYT and Z lists from this extracted block
    XYTZdict = array3d_2_XYTZlist(xcoords,ycoords,tcoords,z_cube)
    
    return(XYTZdict)

def array3d_2_XYTZlist(xcoords,ycoords,tcoords,z_cube=None, as4dcoordblock=False):

    '''
    array as a 3d numpy array
    xcoords,ycoords,tcoords are 1d numpy arrays containing coordinate positions of marginal axes
    (these can be obtained by a call to gridParams_2_XYTmarginallists)
    '''
    
    ### vectors converted to cubes
    xarray = np.vstack(((xcoords,))*len(ycoords))
    yarray = np.vstack(((ycoords,))*len(xcoords)).T
    x_cube = np.dstack(((xarray,))*len(tcoords))
    y_cube = np.dstack(((yarray,))*len(tcoords))
    t_cube = np.ones(np.product((len(ycoords),len(xcoords),len(tcoords)))).reshape((len(ycoords),len(xcoords),len(tcoords)))
    t_cube = t_cube*tcoords

    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()

    # if what we ant is a 4-d block containing coordinates (a3d block where every node is a 3-element array containg x,y,z):
    if as4dcoordblock is True:
        cubeshape = [t_cube.shape[0],t_cube.shape[1],t_cube.shape[2]]
        x_cube.resize(cubeshape+[1])
        y_cube.resize(cubeshape+[1])
        t_cube.resize(cubeshape+[1])
        xyt_cube = np.concatenate((x_cube,y_cube,t_cube),axis=3)
        return({"xyt_cube":xyt_cube})     
    
    ### check shapes of x,y,t,and z cubes are identical (if we passed a z cube)
    if z_cube is not None:
        if (np.shape(x_cube)!=np.shape(z_cube)): raise ValueError, 'shape of x_cube ('+str(np.shape(x_cube))+') does not match shape of z_cube ('+str(np.shape(z_cube))+')'
        if (np.shape(y_cube)!=np.shape(z_cube)): raise ValueError, 'shape of y_cube ('+str(np.shape(y_cube))+') does not match shape of z_cube ('+str(np.shape(z_cube))+')'
        if (np.shape(t_cube)!=np.shape(z_cube)): raise ValueError, 'shape of t_cube ('+str(np.shape(t_cube))+') does not match shape of z_cube ('+str(np.shape(z_cube))+')'

    # collapse x,y,t,and z cubes to 1d arrays  
    if z_cube is not None: z_in = np.ravel(z_cube)
    x_in = np.ravel(x_cube)
    y_in = np.ravel(y_cube)
    t_in = np.ravel(t_cube)
    xyt_in  = np.vstack((y_in,x_in,t_in)).T
    
    if z_cube is not None: return({"xyt_in":xyt_in,"z_in":z_in})
    if z_cube is None: return({"xyt_in":xyt_in})

def gridParams_2_XYTmarginallists(grids):

    '''
    grids is list of three tuples: [(min lon,max lon,n lon),(min lat,max lat,n lat),(min t,max t,n t)]
    returns coordinate values in same units for marginal x,y,t axes (dictionary of three 1-d arrays)
    '''
    
    # extract grid parameters for ease
    ncols = grids[0][2]
    nrows = grids[1][2]
    nmonths = grids[2][2]
    
    cellsize=((grids[0][1]-grids[0][0])/(grids[0][2]-1))
    tsize=((grids[2][1]-grids[2][0])/(grids[2][2]-1))
    xmin_centroid = grids[0][0]
    ymax_centroid = grids[1][1]
    tmin_centroid = grids[2][0]

    # full vectors of grid locations
    xcoords = xmin_centroid + (np.arange(ncols) * cellsize)    
    ycoords = ymax_centroid - (np.arange(nrows) * cellsize) 
    tcoords = tmin_centroid + (np.arange(nmonths)*tsize)
    
    return({"xcoords":xcoords,"ycoords":ycoords,"tcoords":tcoords})

def predictPointsFromBlock(XYT_in,z_in, XYT_out,C,relp,VERBOSE=False):

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
    if (len(z_in)!=n_in): raise ValueError, 'Length of XYT_in ('+str(len(XYT_in[:,0]))+') does not match length of z_in ('+str(len(z_in))+')'

    n_out = len(XYT_out[:,0])
    MaxToSim=float(n_in)

    M = pm.gp.Mean(lambda x:np.zeros(x.shape[:-1]))
    C = pm.gp.NearlyFullRankCovariance(C.eval_fun, relative_precision=relp, **C.params)
    pm.gp.observe(M,C,obs_mesh=XYT_in,obs_vals=z_in)
    f = pm.gp.Realization(M,C)
 
    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()
 
    # return 1d array of simulated values    
    return f(XYT_out)
    
def getUnconditionedBlock(relblock4d,real_index,grids,C,NinThinnedBlock=None,relp=None,FULLRANK=False):

    # get grid's marginal coordinates
    coordsDict = gridParams_2_XYTmarginallists(grids)
    xcoords = coordsDict['xcoords']
    ycoords = coordsDict['ycoords']
    tcoords = coordsDict['tcoords']

    # convert coordinates to 3d blocks  
    XYTdict =array3d_2_XYTZlist(xcoords,ycoords,tcoords, as4dcoordblock=True)
    xyt_cube=XYTdict['xyt_cube']

    # define latent mean function and optionally covert covariance function to sub full rank
    M = pm.gp.Mean(lambda x:np.zeros(x.shape[:-1]))
    if FULLRANK is True:
        C = pm.gp.NearlyFullRankCovariance(C.eval_fun, relative_precision=relp, **C.params)

    # define function for realisation
    f = pm.gp.Realization(M,C)

    # realise at specified locations (return same shape as sd coordinate block - should be 3d)
    simVector = f(xyt_cube)
    
    # insert these values into this realisation of the main output array
    relblock4d[real_index,:,:,:]=simVector[:,:,:]
 
    #from IPython.Debugger import Pdb
    #Pdb(color_scheme='Linux').set_trace()
 
    # return 3d array of simulated values    
    #return simVector

    
    
    
