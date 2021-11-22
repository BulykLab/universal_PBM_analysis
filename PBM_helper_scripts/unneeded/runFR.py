import argparse
from argparse import ArgumentParser
from os.path import isdir
from os import listdir
from subprocess import call
import submitCerebellum
from submitCerebellum import SubmitJobsAndWait

base_dir = "/net/data/Out/kweinand/"

arg_obj = ArgumentParser()
arg_obj.add_argument('dir', metavar = 'dir', help='dir holding combinatorial files - IT IS EXTREMELY IMPORTANT TO GIVE THE FULL PATH. START WITH /net/... or ~/... Note: if protein name contains a period, it will be changed to a hypen')
arg_obj.add_argument('-init', metavar = 'init', help = 'init file - FYI - currently overwritting this variable; please contact Katy Weinand or current techie if variables within init file need to change')
#arg_obj.add_argument('env', metavar = 'env', help = 'R env bash file')
#arg_obj.add_argument('envDisp', metavar = 'envDisp', help = 'R env bash file with display info')
arg_obj.add_argument('notes', metavar = 'notes', help = 'notes about experiment like v number')
arg_obj.add_argument('person', metavar = 'person', help = 'person who ran the experiment')
args = arg_obj.parse_args()

if args.dir[-1] !='/':
    args.dir+='/'

args.init=base_dir+'dream5.init'
#newInit=args.dir+'copied.init'

#call('cp %s %s'%(oldInit,newInit),shell=True)
print 'Using %s as init file'%args.init

env = base_dir+'.envrc'
envDisp = base_dir+'.envrcWdisp'
print 'Using %s as env file'%env

jobs = []

if isdir(args.dir):
    for f in listdir(args.dir):
        if "_combinatorial.txt" in f and "original" not in f:
            PN = f.split("_combinatorial.txt")
            PN=PN[0].replace('.','-')
            newComboFile=args.dir+f[:-4]+'_FR.txt'
            
            call('sed "s/^/%s\\t/" %s > %s'%(PN, args.dir+f, newComboFile),shell=True)
            
            singleCmd = 'python26 %sanalyzeSingleFR.py %s'%(base_dir,' '.join([str(args.dir), str(args.init), str(env), str(envDisp), str(PN), str(newComboFile)]))
            jobs.append(singleCmd)
            
    print 'FeatureREDUCE jobs submitted'
    SubmitJobsAndWait(jobs, wd=args.dir)
    
    #print "Cleaning up logos that didn't finish correctly"
    #errorCmd = 'python %slogoErrors-LSF.py %s'%(base_dir,args.dir)
    #call(errorCmd,shell=True) #SubmitJobsAndWait already within script
    
    #note - this outputDir is specified in the init file, but I did hardcode the init file, so...
    outputDir = args.dir+'FR_results/'
    
    print 'Generating reverse complement logos.'
    RCcmd = 'python26 %sxmlRC.py %s %s'%(base_dir,outputDir, envDisp)
    call(RCcmd,shell=True) #SubmitJobsAndWait already within script
    
    print 'Generating precurser report.'
    brainCmd = 'python26 %smakeBrainReport3Types.py %s %s %s FR'%(base_dir,outputDir,args.notes,args.person)
    call(brainCmd,shell=True) #doesn't need to be run on Cerebellum
    
    print 'All files from this run have been generated except the final report. \nPlease visit http://thebrain.bwh.harvard.edu/KW_logos/buildReport.php in order to finalize which logos to put in the report. \nThe directory to use is %s and the type of logo is FR. \nAn email will be sent to you containing the report.'%outputDir
    
    #Added 7/5/16
    print '\nGenerating all logo report. Mostly for Badis runs, but may come in handy otherwise.'
    allCmd = 'python26 %sBadis09/makeReport.py %s %s %s'%(base_dir,outputDir,args.notes, args.person)
    call(allCmd,shell=True) #doesn't need to be run on Cerebellum
    
    print 'Done.'
    
else:
    raise Exception('Given directory is not a directory')

