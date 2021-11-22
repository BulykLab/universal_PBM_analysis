import argparse
from argparse import ArgumentParser
from subprocess import call
import submitCerebellum
from submitCerebellum import SubmitJob,WaitForJobs

arg_obj = ArgumentParser()
arg_obj.add_argument('dir', metavar = 'dir', help='dir holding combinatorial files - this one is only doing the initial run')
arg_obj.add_argument('init', metavar = 'init', help = 'init file')
arg_obj.add_argument('env', metavar = 'env', help = 'R env bash file')
arg_obj.add_argument('envDisp', metavar = 'envDisp', help = 'R env bash file with display info')
arg_obj.add_argument('PN', metavar = 'PN', help = 'Protein Name')
arg_obj.add_argument('input', metavar = 'input', help = 'input probe file')
args = arg_obj.parse_args()


runningJobs=[]


primaryScript = args.dir+args.PN+'_primary.sh'
with open(primaryScript,'w') as prim:
    prim.write('#!/bin/sh\n\n')
    prim.write('source %s\n\n'%args.env)
    prim.write('java -Xmx4000M FeatureReduce %s -l primary -i %s -c 0 2 1 -ids %s -displayMotifs No\n\n'%(args.init,args.input,args.PN))

primaryID = SubmitJob("bash %s"%primaryScript,wd=args.dir)
runningJobs.append(primaryID)
#primaryID='123'

expectedPrimFile=args.dir+'FR_results/primary.'+args.PN+'.final.psam.10nt.xml'

secondaryScript = args.dir+args.PN+'_secondary.sh'
with open(secondaryScript,'w') as secon:
    secon.write('#!/bin/sh\n\n')
    secon.write('source %s\n\n'%args.env)
    secon.write('java -Xmx4000M FeatureReduce %s -l secondary -i %s -c 0 2 1 -ids %s -residuals xml %s -displayMotifs No\n\n'%(args.init,args.input,args.PN,expectedPrimFile))
  
secondaryID = SubmitJob("bash %s"%secondaryScript,depends=primaryID,wd=args.dir)
runningJobs.append(secondaryID)
#secondaryID='234'

expectedSecFile = args.dir+'FR_results/secondary.'+args.PN+'.final.psam.10nt.xml'

primaryLogoScript = args.dir+args.PN+'_primary_logo.sh'
with open(primaryLogoScript,'w') as prim:
    prim.write('#!/bin/sh\n\n')
    prim.write('source %s\n\n'%args.envDisp)
    #prim.write('DISPLAY=:1.0\nexport DISPLAY\n\n')
    #prim.write('Xvfb :1 -fp /usr/share/X11/fonts/100dpi -screen 0 1280x1024x24 2>/dev/null &\npidToKill=$!\n\n')
    prim.write('java FeatureReduce -psam xml %s -saveLogo %s -batch\n\n'%(expectedPrimFile,expectedPrimFile[:-3]+'png'))
    #prim.write("echo 'Logo is made.'\n\nkill $pidToKill\n\necho 'Xvfb is no longer running'\n\n")
    prim.write("echo 'Logo is made.'\n\npkill 'Xvfb'\n\necho 'Xvfb is no longer running'\n\n")
    

primaryLogoID = SubmitJob("bash %s"%primaryLogoScript,depends=primaryID,wd=args.dir)
runningJobs.append(primaryLogoID)


secondaryLogoScript = args.dir+args.PN+'_secondary_logo.sh'
with open(secondaryLogoScript,'w') as secon:
    secon.write('#!/bin/sh\n\n')
    secon.write('source %s\n\n'%args.envDisp)
    #secon.write('DISPLAY=:1.0\nexport DISPLAY\n\n')
    #secon.write('Xvfb :1 -fp /usr/share/X11/fonts/100dpi -screen 0 1280x1024x24 2>/dev/null &\npidToKill=$!\n\n')
    secon.write('java FeatureReduce -psam xml %s -saveLogo %s -batch\n\n'%(expectedSecFile,expectedSecFile[:-3]+'png'))
    #secon.write("echo 'Logo is made.'\n\nkill $pidToKill\n\necho 'Xvfb is no longer running'\n\n")
    secon.write("echo 'Logo is made.'\n\npkill 'Xvfb'\n\necho 'Xvfb is no longer running'\n\n")

secondaryLogoID = SubmitJob("bash %s"%secondaryLogoScript,depends=secondaryID,wd=args.dir)
runningJobs.append(secondaryLogoID)

print '--- CHECKPOINT: waiting for all submitted jobs to finish... ---'
WaitForJobs(runningJobs)

print 'Jobs done.'
