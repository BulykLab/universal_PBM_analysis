import re
import argparse
from argparse import ArgumentParser
import numpy as np
from os.path import isdir
from os import listdir
import submitCerebellum
from submitCerebellum import SubmitJobsAndWait

arg_object = ArgumentParser()
arg_object.add_argument('dir', metavar = 'dir', help='dir with psam xml files')
arg_object.add_argument('envDisp', metavar = 'envDisp', help = 'R env bash file with display info')
args = arg_object.parse_args()


if args.dir[-1] !='/':
    args.dir+='/'

scriptFile = args.dir+'RClogos.sh'
SCRIPT = open(scriptFile,'w')
SCRIPT.write('#!/bin/sh\n\n')
SCRIPT.write('source %s\n\n'%args.envDisp)


if isdir(args.dir):
    for f in listdir(args.dir):
        if "final.psam.10nt.xml" in f:
            
            Apos = []
            Gpos = []
            Cpos = []
            Tpos = []
            
            IN = open(args.dir+f,'r')
            
            for i,line in enumerate(IN):
                colIdx = re.search(r'<col indx="([0-9]+)">',line)
                A = re.search(r'<weight sym="adenine" prob="([0-9,.,E,-]+)"',line)
                G = re.search(r'<weight sym="guanine" prob="([0-9,.,E,-]+)"\/>',line)
                C = re.search(r'<weight sym="cytosine" prob="([0-9,.,E,-]+)"\/>',line)
                T = re.search(r'<weight sym="thymine" prob="([0-9,.,E,-]+)"\/>',line)
                
                #if colIdx is not None:
                    #print 'new position: %s'%colIdx.group(1)
                if A is not None:
                    #print 'A'
                    Apos.append(A.group(1))
                if G is not None:
                    #print 'G'
                    Gpos.append(G.group(1))
                if C is not None:
                    #print 'C'
                    Cpos.append(C.group(1))
                if T is not None:
                    #print 'T'
                    Tpos.append(T.group(1))
            
            IN.close()
            
            Tnp=np.flipud(np.array(Apos))
            Gnp=np.flipud(np.array(Cpos))
            Cnp=np.flipud(np.array(Gpos))
            Anp=np.flipud(np.array(Tpos))
            
            if len(Tnp) != len(Gnp) or len(Tnp) != len(Cnp) or len(Tnp) != len(Anp):
                print len(Tnp)
                print len(Gnp)
                print len(Cnp)
                print len(Anp)
                raise Exception('%s: not rectangular'%f)
            
            outputFile=args.dir+'RC.'+f
            
            OUT = open(outputFile,'w')
            OUT.write('<MarkovModel>\n  <alphabet name="DNA"/>\n')
            for i in xrange(len(Tnp)):
                OUT.write('  <col indx="%s">\n'%(i+1))
                OUT.write('    <weight sym="adenine" prob="%s"/>\n'%Anp[i])
                OUT.write('    <weight sym="guanine" prob="%s"/>\n'%Gnp[i])
                OUT.write('    <weight sym="cytosine" prob="%s"/>\n'%Cnp[i])
                OUT.write('    <weight sym="thymine" prob="%s"/>\n'%Tnp[i])
                OUT.write('  </col>\n')

            OUT.write('</MarkovModel>\n')
            OUT.close()
            
            SCRIPT.write('java FeatureReduce -psam xml %s -saveLogo %s -batch\n\n'%(outputFile,outputFile[:-3]+'png'))


SCRIPT.write("echo 'Logos are made.'\n\npkill 'Xvfb'\n\necho 'Xvfb is no longer running'\n\n")
SCRIPT.close()

#print "bash %s"%scriptFile
SubmitJobsAndWait(["bash %s"%scriptFile],wd=args.dir)
