import logging
from subprocess import call
from os.path import isfile, basename, dirname


def SubmitJob(job_string, log_file = 'job_log.txt', memory = 8, cwd = None, num_cores = 1, depends = None, wd = None, extra_flags = ''):
    """Submits a job to the SunGrid queueing system. The string should be
    exactly what would be typed at the Unix shell in the current
    directory.
    
    Returns the job ID given by the SunGrid system.
    """
    #KW This needs to change to a SunGrid queueing system using qsub for Cerebellum
    
    #KW string is a python module installed on Cerebellum
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

        #KW Comma Separated List! http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html
        job_depend_str=''.join(['%s,' %(_) for _ in depends])
        
        #get rid of the last comma
        job_depend_string = job_depend_str[:-1]

        extra_flags += '-hold_jid %s' %job_depend_string
    if wd is not None:
        working_dir='-wd %s'%wd
    else:
        working_dir='-cwd'
    
    #KW change
    #memory should be in MegaBytes - converted from KB: 8000KB=8MB Turns out I didn't need to, as I could just use K, but too late.
    #no projects in Cerebellum (qconf -sprjl)
    #just chose a 8-node queue in Cerebellum - may want to make this a user input option
    #Don't need parallel environment, num_cores is only ever 1, so not used here
    #Maybe include job name: -N flag.
    #1/17/18 KW changed the queues; 1-4 are down
    queue_call = 'qsub -q all.q@compute-0-6.local,all.q@compute-0-7.local,all.q@compute-0-8.local -l mem_total=%sM -o %s -b y %s -S /bin/sh %s "%s"' \
                 %(memory, log_file, working_dir, extra_flags, job_string)
    
    #print queue_call
    
    job_process = Popen(queue_call, shell = True, stdout = PIPE, stderr = PIPE, cwd = None)
    out_stream, err_stream = job_process.communicate()
    job_process.wait()
    
    logging.info(out_stream)
    print str(out_stream).strip()
    
    if len(err_stream) > 0:
        raise Exception('Error submitting job: %s \nException: %s' %(job_string, err_stream))
    
    #Find Job ID and return it
    #On Cerebellum, the job id is a number. The output files are structured JobName.oJobID & JobName.eJobID.
    #The immediate terminal output is "Your job JobID ("JobName") has been submitted"
    job_id = re.search(r'Your job ([0-9]*)', out_stream)
    
    #KW originally, job_id.group()[1:-1].
    if not job_id is None:
        return str(job_id.group(1))
    else:
        raise Exception('Job ID could not be retrieved from Sun Grid system. This is likely an error.')
    

def SubmitJobsAndWait(list_of_jobs, num_cores = 1, print_jobs = True, wd = None, memory = 12, **kwargs):
    """Submits a list of string specifiying jobs. These should be written as they'd be typed
    into the shell.
    kwargs go straight to the 'SubmitJob" function."""    

    from time import sleep

    if len(list_of_jobs) == 0:
        return 0

    submitted_jobs = set()
    #KW I don't think this function needs to change since it's calling the previous function.
    for job_string in list_of_jobs:
        print 'Submitting to Cerebellum: \n %s' %job_string
        submitted_jobs.add(SubmitJob(job_string, num_cores = num_cores, memory = memory, wd = wd, **kwargs))
    
    
    while True:
        
        running_jobs = WhichJobsAreRunning()
        
        #This intersects the set of running and submitted jobs
        #If the set is empty, this means that all the jobs have finished
        if len(running_jobs & submitted_jobs) == 0:
            break
        else:
            sleep(10)
    
def WhichJobsAreRunning():
    
    from subprocess import Popen, PIPE
    import re
    #KW Changed bjobs to qstat that only shows running jobs
    job_process = Popen(['qstat','-s','r'], shell = True, stdout = PIPE, stderr = PIPE, cwd = None)
    
    out_stream, err_stream = job_process.communicate()
    
    
    job_process.wait()
    
    #KW change
    #Assuming the current user is the one submitting the job.
    uid_job = Popen(['id','-u','-n'],stdout=PIPE, stderr=PIPE)
    uid, uid_err = uid_job.communicate()
    #job_ids = (line.split(' ')[0] for line in out_stream.split('\n') if uid[:-1] in line)
    job_ids = []
    
    #print "uid: %s" %uid[:-1]
    
    for l in out_stream.split('\n'):
        #print "%s\n" %l
        if uid[:-1] in l:
            #this will become an issue if 500K more jobs are submitted - then make it a 7 instead of 6
            val = re.search(r'^ ([0-9]{6}) ', l)
            if val is not None:
                #print "job id: %s\n" %val.group(0)[1:]
                job_ids.append(val.group(1))
    
    return set(job_ids)

#ADDED from General.py
def WaitForJobs(job_id_list):
    """Similar to SubmitJobsAndWait, but gets a list of job numbers for jobs that are assumed to have
    already been submitted."""

    from time import sleep


    if job_id_list is None:
        return 0

    submitted_jobs = set(job_id_list)

    while True:

        running_jobs = WhichJobsAreRunning()
        
        #This intersects the set of running and submitted jobs
        #If the set is empty, this means that all the jobs have finished
        if len(running_jobs & submitted_jobs) == 0:
            break
        else:
            sleep(10)

