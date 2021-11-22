import logging
from subprocess import call
from os.path import isfile, basename, dirname


def SubmitJob(job_string, log_file = 'job_log.txt', memory = '2G', cwd = None, num_cores = 1, depends = None, wd = None, timeLim = '0-02:00', extra_flags = '',queue = 'short'):
    """Submits a job to the LSF queueing system. The string should be
    exactly what would be typed at the Unix shell in the current
    directory.
    
    Returns the job ID given by the LSF system.
    """
    
    from subprocess import Popen, PIPE
    import re
    from string import replace

    #KW $ and " have a special meaning in python, hence the escape character \
    # Escape shell command and double quotes 
    job_string = replace(job_string, '$','\$')
    job_string = replace(job_string, '"','\"')
    
    #ADDED from General.py
    if depends is not None:
        if type(depends) is str or type(depends) is int:
            depends = [str(depends)]
        
        #sbatch --dependency=afterok:1:2,afterany:3:4,?afternotok:5 submit.sh
        extra_flags += ' --dependency=afterok:%s --kill-on-invalid-dep=yes' %(':'.join(depends))
    
    if num_cores>1:
        extra_flags+=' -c %s'%(num_cores)
    
    if wd is not None:
        extra_flags+=' -D %s'%(wd)
    
    if memory>2:
        extra_flags+=' --mem=%s'%(memory)
    
    extra_flags+=' '
    
    
    
    #Would this work? or else just job? I would add the name argument, but I don't want to change the actual script if I don't absolutely have to do so.    
    name = job_string.split(" ")[1].split("/")[-1].rsplit(".",1)[0]
    out_file = name+"-%j.out"
    err_file = name+"-%j.err"
    
    
    #KW change
    #memory should be in KBytes.
    #no projects in Cerebellum (qconf -sprjl)
    #just chose a 8-node queue in Cerebellum - may want to make this a user input option
    #Don't need parallel environment, num_cores is only ever 1, so not used here
    #Maybe include job name: -N flag.
    queue_call = 'sbatch -p %s -t %s --job-name=%s -o %s -e %s --exclude=compute-f-17-[09-25]%s--wrap="%s"' \
                 %(queue, timeLim, name, out_file, err_file, extra_flags, job_string)
    
    print queue_call
    
    
    job_process = Popen(queue_call, shell = True, stdout = PIPE, stderr = PIPE, cwd = None)
    out_stream, err_stream = job_process.communicate()
    job_process.wait()
    
    logging.info(out_stream)
    print str(out_stream).strip()
    
    if len(err_stream) > 0:
        raise Exception('Error submitting job: %s \nException: %s' %(job_string, err_stream))
    
    #Find Job ID and return it
    #The immediate terminal output is like "Submitted batch job 44045936"
    job_id = re.search(r'Submitted batch job ([0-9]+)', out_stream)
    
    #KW originally, job_id.group()[1:-1].
    if not job_id is None:
        return str(job_id.group(1))
    else:
        raise Exception('Job ID could not be retrieved from slurm system. This is likely an error.')
    

def SubmitJobsAndWait(list_of_jobs, num_cores = 1, print_jobs = True, wd = None, memory = '2G', timeLim = '0-02:00', queue= 'short', **kwargs):
    """Submits a list of string specifiying jobs. These should be written as they'd be typed
    into the shell.
    kwargs go straight to the 'SubmitJob" function."""    
    
    from time import sleep

    if len(list_of_jobs) == 0:
        return 0

    submitted_jobs = set()
    #KW I don't think this function needs to change since it's calling the previous function.
    for job_string in list_of_jobs:
        print 'Submitting to O2: \n %s' %job_string
        submitted_jobs.add(SubmitJob(job_string, num_cores = num_cores, memory = memory, wd = wd, timeLim=timeLim, queue=queue, **kwargs))
    
    
    sleep(10)
    while True:
        
        running_jobs = WhichJobsAreRunning()
        
        #This intersects the set of running and submitted jobs
        #If the set is empty, this means that all the jobs have finished
        if len(running_jobs & submitted_jobs) == 0:
            break
        else:
            sleep(10)
    

#KW Changed ERISone    
def WhichJobsAreRunning():
    
    from subprocess import Popen, PIPE
    import re
    
    
    job_process = Popen('O2squeue', shell = True, stdout = PIPE, stderr = PIPE)
    
    out_stream, err_stream = job_process.communicate()
    
    
    job_process.wait()
    
    job_ids = []
    
    
    for l in out_stream.split('\n'):
        parse = re.search(r"^([0-9]+)\s+",l)
        if parse is not None:
            job_ids.append(parse.group(1))
    
    return set(job_ids)

#ADDED from General.py
def WaitForJobs(job_id_list):
    """Similar to SubmitJobsAndWait, but gets a list of job numbers for jobs that are assumed to have
    already been submitted."""
    
    from time import sleep
    
    if job_id_list is None:
        return 0

    submitted_jobs = set(job_id_list)
    
    sleep(10)
    while True:

        running_jobs = WhichJobsAreRunning()
        
        #This intersects the set of running and submitted jobs
        #If the set is empty, this means that all the jobs have finished
        if len(running_jobs & submitted_jobs) == 0:
            break
        else:
            sleep(10)

