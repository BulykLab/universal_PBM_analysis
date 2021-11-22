import argparse
from argparse import ArgumentParser
from os.path import isdir
from os import listdir
from subprocess import call
import submitCerebellum
from submitCerebellum import SubmitJobsAndWait


arg_obj = ArgumentParser()
arg_obj.add_argument('dir', metavar = 'dir', help='dir holding script files')
args = arg_obj.parse_args()


if args.dir[-1] !='/':
    args.dir+='/'


if isdir(args.dir):
    for f in listdir(args.dir):
        if "logo.sh.e" in f:
            cont = False
            with open(args.dir+f, 'r') as FILE:
                if 'Fatal server error' in FILE.read():
                    cont=True
            if cont==True:
                tabs = f.split('.e')
                #print 'bash %s'%tabs[0]
                SubmitJobsAndWait(["bash %s"%tabs[0]],wd=args.dir)
