import os

lines = {'py': 0, 'f': 0,'r': 0}
topdirs = ['mbgw','mbgw-scripts','testmbgw']

for topdir in topdirs:
    for dirpath, dirnames, filenames in os.walk('../%s'%topdir):
        for ext in lines.iterkeys():
            for name in filter(lambda n: n.split('.')[-1].lower() == ext, filenames):
                f = file(dirpath+'/'+name)
                for line in f:
                    lines[ext] += 1
                
for ext, n in lines.iteritems():
    print '%i lines of code in files ending in %s'%(n,ext)
print '--------------------------------'
print '%i total lines of code'%sum(lines.itervalues())