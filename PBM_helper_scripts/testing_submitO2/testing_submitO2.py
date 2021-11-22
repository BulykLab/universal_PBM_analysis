from submitO2 import WhichJobsAreRunning, SubmitJobsAndWait, SubmitJob, WaitForJobs
import time

job = "python /n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_helper_scripts/waitFile.py"

print WhichJobsAreRunning()

print time.strftime('%c')
job_id = SubmitJob(job,timeLim='0-00:05',memory='10M')
print WhichJobsAreRunning()
print time.strftime('%c')
WaitForJobs([job_id])
print time.strftime('%c')
SubmitJobsAndWait([job],timeLim='0-00:05',memory='10M')
print WhichJobsAreRunning()
print time.strftime('%c')
