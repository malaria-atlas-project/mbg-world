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
from mbgw.master_grid import *
from mbgw import auxiliary_data
import os
curpath = os.getcwd()
import mbgw
mbgw_root = __root__ = mbgw.__path__[0]
r_path = mbgw_root+'/joint_simulation/CONDSIMalgorithm'
os.chdir(r_path)
from rpy import r
r.source("CONDSIMpreloop.R")
r.source("CONDSIMmonthloop.R")
os.chdir(curpath)

__all__ = ['create_realization', 'create_many_realizations','reduce_realizations']

def get_covariate_submesh(name, grid_lims):
    """
    Matches the specified spatial mesh to the 'master' covariate mesh.
    """
    return getattr(mbgw.auxiliary_data, name).data[nrows-grid_lims['bottomRow']:nrows-grid_lims['topRow']+1,
                                                    grid_lims['leftCol']-1:grid_lims['rightCol']][::-1,:].T

def create_many_realizations(burn, n, trace, meta, grid_lims, start_year, nmonths, outfile_name, N_nearest, memmax, relp=1e-3, mask_name=None, n_in_trace=None):
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
    bad = []
    in_mesh = np.ones(data_locs.shape[0],dtype=bool)
    for l in data_locs:
        for j in xrange(3):
            if l[j] <= grids[j][0] or l[j] >= grids[j][1]:
                in_mesh[i]=False
                bad.append(l)
    if len(bad) > 0:
        bad = np.array(bad)
        bad[:,0] *= rad_to_deg
        bad[:,1] *= rad_to_deg
        bad[:,2] += 2009
        print 'Warning: The following data locations [lon,lat,t] are out of bounds: \n'+str(bad)

    # Find the mesh indices closest to the data locations
    data_mesh_indices = np.empty(data_locs.shape, dtype=np.int)

    for i in xrange(len(data_locs)):
        for j in xrange(3):
            data_mesh_indices[i,j] = np.argmin(np.abs(data_locs[i,j] - axes[j]))

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
                        shape=(N,) + grid_shape, 
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
    
    data_locs = data_locs[in_mesh]
    data_mesh_indices = data_mesh_indices[in_mesh]
    
    # Total number of pixels in month.
    npix = grid_shape[0]*grid_shape[1]
    # Maximum number of pixels in tile.
    npixmax = memmax/4./data_locs.shape[0]
    # Minimum number of tiles needed.
    ntiles = npix/npixmax
    # Blocks.
    n_blocks_x = n_blocks_y = np.ceil(np.sqrt(ntiles))
    
    # from IPython.Debugger import Pdb
    # Pdb(color_scheme='Linux').set_trace()
    
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
            mean_ondata += getattr(meta, key)[:][in_mesh] * this_coef
            this_pred_covariate = get_covariate_submesh(key, grid_lims) * this_coef
            covariate_mesh += this_pred_covariate

        # Pull covariance information out of trace
        this_C = trace.group0.C[indices[i]]
        this_C = pm.gp.NearlyFullRankCovariance(this_C.eval_fun, relative_precision=relp, **this_C.params)

        data_vals = trace.PyMCsamples[i]['f'][in_mesh]
        create_realization(outfile.root.realizations, i, this_C, mean_ondata, this_M, covariate_mesh, data_vals, data_locs, grids, axes, data_mesh_indices, n_blocks_x, n_blocks_y, relp, mask, N_nearest)
        outfile.flush()
    outfile.close()

def create_realization(out_arr,real_index, C, mean_ondata, M, covariate_mesh, tdata, data_locs, grids, axes, data_mesh_indices, n_blocks_x, n_blocks_y, relp, mask, N_nearest):
    """
    Creates a single realization from the predictive distribution over specified space-time mesh.
    """
    grid_shape = tuple([grid[2] for grid in grids])
    
    # Container for x
    x = np.empty(grid_shape[:2] + (3,))
    mlon,mlat = np.meshgrid(*axes[:2])
    x[:,:,0] = mlon.T
    x[:,:,1] = mlat.T
    x[:,:,2] = 0    
    del mlon, mlat

    Cp = C.params
    
    # # Prepare input dictionaries
    # covParamObj = {'Scale': Cp['scale'][0]*rad_to_km,
    #                 'amp': Cp['amp'][0], 
    #                 'inc': Cp['inc'][0], 
    #                 'ecc': Cp['ecc'][0], 
    #                 't.lim.corr': Cp['tlc'][0], 
    #                 'scale.t': Cp['st'][0], 
    #                 'sin.frac': Cp['sf'][0]}
    # gridParamObj = {'YLLCORNER': grids[1][0]*rad_to_deg, 
    #                 'CELLSIZE': (grids[1][1]-grids[1][0])/(grids[1][2]-1.)*rad_to_deg, 
    #                 'NROWS':grid_shape[1],
    #                 'NCOLS':grid_shape[0]}
    # monthParamObj = {'Nmonths':grid_shape[2],'StartMonth':grids[2][0]}
    # 
    # # Call R preprocessing function and check to make sure no screwy re-casting has taken place.
    # os.chdir(r_path)
    # preLoopObj = r.CONDSIMpreloop(covParamObj,gridParamObj,monthParamObj)
    # tree_reader = reader(file('listSummary_preLoopObj_original.txt'),delimiter=' ')
    # preLoopClassTree, junk = parse_tree(tree_reader)
    # preLoopObj = compare_tree(preLoopObj, preLoopClassTree)
    # 
    # OutMATlist = preLoopObj['OutMATlist']
    # tree_reader = reader(file('listSummary_OutMATlist_original.txt'),delimiter=' ')
    # OutMATClassTree, junk = parse_tree(tree_reader)
    # OutMATlist = compare_tree(OutMATlist, OutMATClassTree)
    # os.chdir(curpath)
    # 
    # # Create and store unconditional realizations
    # print '\tGenerating unconditional realizations.'
    # t1 = time.time()
    # for i in xrange(grid_shape[2]):
    #     os.chdir(r_path)
    #     monthObject = r.CONDSIMmonthloop(i+1,preLoopObj,OutMATlist)
    #     os.chdir(curpath)
    #     OutMATlist= monthObject['OutMATlist']
    #     MonthGrid = monthObject['MonthGrid']
    #     out_arr[real_index,:,:,i] = MonthGrid[::-1,:].T[:grid_shape[0], :grid_shape[1]]
    # t2 = time.time()
    # print '\t\tDone in %f'%(t2-t1)
    # 
    # # delete unneeded R products
    # del OutMATlist, preLoopObj, MonthGrid, monthObject
    
    # Figure out pdata
    pdata = np.empty(tdata.shape)
    for i in xrange(len(pdata)):
        pdata[i] = out_arr[(real_index,) + tuple(data_mesh_indices[i,:])]
    
    # Bring in data.
    print '\tKriging to bring in data.'    
    print '\tPreprocessing.'
    t1 = time.time()  
    dev, xbi, ybi, rel_data_ind = preprocess(C, data_locs, grids, x, n_blocks_x, n_blocks_y, tdata, pdata, relp, mean_ondata, N_nearest)   
    t2 = time.time()
    print '\t\tDone in %f'%(t2-t1)
    
    row = np.empty(grid_shape[:2], dtype=np.float32)
    print '\tKriging.'
    t1 = time.time()  
    for i in xrange(grid_shape[2]-1,-1,-1):    
        # print '\t Month %i of %i'%(i,grid_shape[2])
        row.fill(0.)
        
        x[:,:,2] = axes[2][i]

        C_eval = C(data_locs, data_locs)      
        krige_month(C, C_eval, i, data_locs, grid_shape, n_blocks_x, n_blocks_y, rel_data_ind, xbi, ybi, x, dev, row, mask, relp)
                
        row += covariate_mesh[:,::-1]
        row += M(x)
        row += out_arr[real_index,:,:,i]

        import pylab as pl
        import matplotlib
        matplotlib.interactive(True)
        pl.close('all')
        pl.figure(figsize=(8,14))
        pl.subplot(1,2,1)
        pl.imshow(row.T, interpolation='nearest', extent=[grids[0][0],grids[0][1],grids[1][0],grids[1][1]],cmap=matplotlib.cm.hot)
        pl.colorbar()
        pl.plot(data_locs[:,0],data_locs[:,1],'b.',markersize=2)        
        pl.axis('off')
        pl.subplot(1,2,2)
        row = pm.invlogit((row + np.random.normal(size=row.shape)*np.sqrt(2)).ravel()).reshape(row.shape)
        row[np.where(1-mask[:,::-1])] = 0
        pl.imshow(row.T, interpolation='nearest', extent=[grids[0][0],grids[0][1],grids[1][0],grids[1][1]],cmap=matplotlib.cm.hot,vmin=-.5,vmax=1.5)
        pl.colorbar()        
        pl.plot(data_locs[:,0],data_locs[:,1],'b.',markersize=2)                
        pl.axis('off')     
        pl.savefig('row%i.pdf'%i)   
        from IPython.Debugger import Pdb
        Pdb(color_scheme='Linux').set_trace()   

        # NaN  the oceans to save storage
        row[np.where(1-mask[:,::-1])] = missing_val
        
        out_arr[real_index,:,:,i] = row          
    t2 = time.time()
    print '\t\tDone in %f'%(t2-t1)        
        

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