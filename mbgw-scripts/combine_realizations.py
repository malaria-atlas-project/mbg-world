from tables import *
from numpy import *
import os, sys

path = 'oxgrid-outputs'
trace_id = 'QRYPFPR010708_Africa_Run_9.10.2008'

fnames = os.listdir(path)
fnames = filter(lambda n: n.find(trace_id) >= 0, fnames)

# Index ranges
iters_lo = map(lambda n: int(n.split('.')[-2].split('_')[-2]), fnames)
iters_hi = map(lambda n: int(n.split('.')[-2].split('_')[-1]), fnames)
iters = zip(iters_lo, iters_hi)
lowest = min(iters_lo)
highest = max(iters_hi)

# Create output file
outfile = openFile(path+'/'+fnames[0].replace('%i_%i'%(iters_lo[0],iters_hi[0]), '%i_%i'%(lowest, highest)),'w')

# Open a single realization file to use as a template
example = openFile(path+'/'+fnames[0])
sh = example.root.realizations.shape[1:]

# Record lat, lon, t meshes from template
for dir in ['lon','lat','t']:
    outfile.createArray('/','%s_axis'%dir, getattr(example.root, '%s_axis'%dir)[:])
outfile.createGroup('/','group0')

# Create containers for means and covariances from template.
for child in example.root.group0._f_listNodes():
    n = outfile.createVLArray('/group0', child.name, ObjectAtom())
        
# Create table to hold PyMC samples
outfile.createTable('/','PyMCsamples',example.root.PyMCsamples.description)

# Create compressed array to hold realizations
outfile.createCArray('/', 'realizations', Float32Atom(), (highest-lowest,) + sh, chunkshape=example.root.realizations.chunkshape, filters=Filters(complevel=1, complib='zlib'))

for il, ih in iters:
    in_fname = path+'/'+'_'.join(['realizations', trace_id, 'iterations', str(il), str(ih)])+'.hdf5'
    print 'Drawing from %s'%in_fname
    infile = openFile(in_fname)
    for j in xrange(ih-il):
        # Append mean and covariance objects
        for oa in outfile.root.group0._f_listNodes():
            oa.append(getattr(infile.root.group0, oa.name)[j])
        # Append scalar samples
        outfile.root.PyMCsamples.append(infile.root.PyMCsamples[j:j+1])
        # Append realizations
        for t in xrange(sh[-1]):
            outfile.root.realizations[il+j,:,:,t] = infile.root.realizations[j,:,:,t]