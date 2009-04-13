import boto
from subprocess import PIPE, STDOUT, Popen
import time

__all__ = ['send_work', 'spawn_engines', 'stop_all_engines', 'map_jobs']

init_ssh_str = 'ssh -o "StrictHostKeyChecking=no" -i /amazon-keys/MAPbase.pem'
init_scp_str = 'scp -o "StrictHostKeyChecking=no" -i /amazon-keys/MAPbase.pem'

def send_work(e, cmd):
    """
    e: a boto Instance object (an Amazon EC2 instance)
    cmd: a string
    
    After the call, e.job will be a Popen instance. e.job.cmd will be cmd.
    """
    command_str = init_ssh_str + ' root@%s %s'%(e.dns_name, cmd)
    e.job = Popen(command_str, shell=True, stdout=PIPE, stderr=STDOUT)
    e.job.cmd = cmd

def stop_all_engines():
    """
    Shuts down ALL Amazon EC2 instances we own!
    """
    conn = boto.connect_ec2()
    for r in conn.get_all_instances():
        r.stop_all()

def spawn_engines(N_engines):
    """
    Starts up N_engines engines using the image with scipy, etc. on it.
    """
    conn = boto.connect_ec2()
    r = conn.run_instances('ami-88e80fe1', min_count=N_engines, max_count=N_engines, security_groups=['MAP'], instance_type='m1.lcmde')
    print 'Starting %s'%r.instances

def map_jobs(cmds, init_cmds=None, upload_files=None, interval=10, shutdown=True):    
    """
    cmds: list of strings that can be executed from the shell.
    
    Optional arguments:
        upload_files: list of paths.
          Will be uploaded to each instance before any init_cmds or cmds
          are sent to it.
        init_cmds: list of strings that can be executed from the shell.
          Will be executed once, in order, on each instance before any cmds
          are sent to it.
        interval: integer. Delay in seconds between polling instances.
        shutdown: boolean. Whether to shut down instances that are not needed.
          If True, all instances will be terminating by the time function returns.
          
    Returns a list of (cmd, output) tuples. Output is everything the process wrote
    to stdout and stderr... but there won't be anything in stderr, because only
    successful exits get written out. 
    """
    returns = []
    
    conn = boto.connect_ec2()
    r = conn.get_all_instances()[-1]
    print 'Extant engines are %s'%r.instances
        
    spawning_engines = [e for e in r.instances]
    running_engines = []
    done = False
    retcode = None
    
    while True:
        print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        print "Gotta check my fly-paper... see what's stickin."
    
        # Watch engines to see when they come alive.
        for e in spawning_engines:
            if e.update()==u'running':
                print '\n%s is up'%e
                e.job = None
                spawning_engines.remove(e)
                running_engines.append(e)
                if upload_files is not None:
                    print '\n\tUploading files to %s'%e
                    for upload_file in upload_files:
                        print '\t'+upload_file
                        p = Popen(init_scp_str +  ' %s root@%s:%s'%(upload_file, e.dns_name, upload_file), shell=True, stdout=PIPE,stderr=STDOUT)
                        retval = p.wait()
                        if retval:
                            raise ValueError, 'Upload failed! Output:\n' + p.stdout.read()
                if init_cmds is not None:
                    print '\n\tExecuting initial commands on %s'%e
                    for init_cmd in init_cmds:
                        print '\t$ %s'%init_cmd
                        p = Popen(init_ssh_str + ' root@%s '%e.dns_name+ init_cmd, shell=True, stdout=PIPE, stderr=STDOUT)
                        while p.poll() is None:
                            print '\t\tWaiting for %i...'%p.pid
                            time.sleep(10)
                        retval = p.poll()                    
                        if retval:
                            raise ValueError, 'Initial command failed! Output:\n' + p.stdout.read()
                        print '\tSuccessful.'    

        N_running = 0
        for e in running_engines:
            # See if previous work is done
            if e.job is not None:
                retcode = e.job.poll()

                if retcode is not None:
            
                    print '\n\t%s has completed\n\t$ %s\n\twith code %i'%(e,e.job.cmd,retcode)
                    if retcode>0:
                        print '\n\tAdding\n\t$ %s\n\tback onto queue because %s fucked it up, code %i. Message:\n'%(e.job.cmd, e, retcode)
                        for line in e.job.stdout:
                            print line
                        cmds.append(e.job.cmd)
                    else:
                        returns.append((e.job.cmd, e.job.stdout.read()))
                        
                else:
                    N_running += 1
        
            # In event of fault, move e back to spawning list.
            if e.update()!=u'running':
                running_engines.remove(e)
                spawning_engines.append(e)
        
            # Send work    
            elif (e.job is None or retcode is not None) and not done:
                if len(cmds) == 0:
                    print 'No jobs remaining in queue'
                    done = True
                else:
                    cmd=cmds.pop(0)
                    print '\n\tSending\n\t$ %s\n\tto %s'%(cmd,e)
                    send_work(e, cmd)
                    N_running += 1

            # Kill the engine
            if done and shutdown:
                print 'Stopping %s'%e
                e.stop()

        if done and N_running == 0:
            print 'All jobs complete!'
            break
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'        
        time.sleep(interval)
    
    if shutdown:
        r.stop_all()
    return returns
    
if __name__ == '__main__':
    # How many processes to spawn?
    # N_engines = 2
    N_jobs = 2
    iter_per_job = 2
    cmds = ['screen ipython amazon_joint_sim.py %i %i %i'%(i,iter_per_job,N_jobs) for i in xrange(N_jobs)]
    returns = map_jobs(cmds, 
                shutdown=True, 
                upload_files=['amazon_joint_sim.py','cloud_setup.sh'], 
                init_cmds=['bash /root/cloud_setup.sh'], 
                interval=20)    
    
    # cmds = ['python boto_upload.py %i'%i for i in xrange(N_jobs)]
    # returns = map_jobs(cmds, shutdown=False, upload_files=['boto_upload.py','cloud_setup.sh'], init_cmds=['bash /root/cloud_setup.sh', 'apt-get install python-boto'], interval=2)
