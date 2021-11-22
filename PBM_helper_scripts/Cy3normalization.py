import logging
import re
import argparse
from argparse import ArgumentParser
from subprocess import call 
from os.path import isfile, basename, dirname, isdir 
from os import makedirs
from time import sleep
from submitO2 import WhichJobsAreRunning, SubmitJobsAndWait, SubmitJob


#KW change paths
PBM_suite_path = '/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_analysis_suite/'

analysis_tools_path = '/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_helper_scripts/'

#################
### FUNCTIONS ###
#################

class ArgParser(ArgumentParser):

    import logging

    def AddStandardArgs(parser_object):
        """Adds standard arguments to an argument parser with specific
        variables (e.g., -verbose)."""
        
        #-verbose is the correct flag for qsub
        parser_object.add_argument('-v', '--verbose',
                                   help='Verbose mode (-vv for even more verbose)',
                                   action='count')
        
        #I didn't change anything here
        parser_object.add_argument('--debuglog', metavar = 'log.txt',
                                   help = 'If given, redirects debugging output to a file',
                                   default = None)
        
        return parser_object

    def Parse(self):
        
        parsed_args = self.parse_args()
        
        if parsed_args.verbose == 2: 
            logging.basicConfig(level=logging.DEBUG)
        elif parsed_args.verbose == 1:
            logging.basicConfig(level=logging.INFO)
        else:
            logging.basicConfig(level=logging.WARNING)
        
        return parsed_args
    
    @classmethod
    def default(cls, *args, **kwargs):
        
        #Somehow this should be a "self" rather than a return (or both)
        
        return cls.AddStandardArgs(cls(*args, **kwargs)) 
    
    @classmethod
    def fromList(cls, args_list, descriptions = None):
        
        parser_object = cls()
          
        if descriptions is None:
            [parser_object.add_argument(arg, required= True)
             for arg in args_list]
        
        else:
            if len(args_list) != len(descriptions):
                raise Exception('Different number of arguments \
                                and descriptions')
            else:
            
                [parser_object.add_argument(arg, help=descr, required = True)
                for arg,descr in zip(args_list, descriptions)]
            
        
        parser_object = cls.AddStandardArgs(parser_object)        
        
        return parser_object.parse_args()
    
def AllFilesInDir(directory_path):
    
    import os
    from os.path import isdir
    
    if not isdir(directory_path):
        raise Exception('The directory "%s" does not exist.' %directory_path)
    
    for filename in os.listdir(directory_path):
        yield filename


arg_object = ArgParser.default()
arg_object.add_argument('directory', metavar = 'PBM_directory',
                        help='Path to directory with array data', type = str)
arg_object.add_argument('sequence_file', metavar = 'array_format_file',
                        help='File containing info about the array format \
                             (e.g., 8x60k_v14_amadid_30265_analysis.txt)', type = str)
args = arg_object.Parse()


if not args.directory[-1] == '/':
    args.directory += '/'

files_in_dir = AllFilesInDir(args.directory)


files_to_use = sorted([filename for filename in files_in_dir if re.search('[1-8]-[1-8].gpr', filename) ])

if len(files_to_use) == 0:
    raise Exception('The given directory does not contain any .gpr files!')

matched_files = {}
experiment_name = None



for filename in files_to_use:
    
    block_match = re.search('_[1-9]-[1-9].gpr', filename)
    if not block_match:
        raise Exception('The file %s does not contain a description of the block it belongs to.' %filename)
    else:
        block_number = block_match.group()[1:4]
    
    if not block_number in matched_files:
        matched_files[block_number]={'Alexa488':[]}
    
    findCy3 = re.search('_[cC]y3_lp[0-9]*pg[0-9]*_([1-9])-[1-9].gpr', filename)
    if findCy3 is not None:
        #print filename+"\t"+findCy3.group(1)
        matched_files[block_number]['Cy3']=filename
        
    else:
        #print "bad filename: "+filename
        raise Exception('File does not contain properly formatted intensity information: %s' %filename)
    
        
#Iterate over all blocks and ask for names
protein_names = []
blocks = sorted(matched_files.keys())

for block_number in blocks:
    protein_names.append("block%s"%block_number[0])
    #print "Using name: %s" %protein_names[-1]



######################
###  Normalization ###
######################

base_normalize_call = 'perl '+PBM_suite_path+'normalizeCy3only.pl -i %s -c %s -s %s -o %s'

jobs_queue = []

for block_number, protein_name in zip(blocks, protein_names):
    masliner_file = "/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_helper_scripts/Alexa488_Cy3R2_dummyfile.gpr"
    
    if 'Cy3' in matched_files[block_number]:    
        cy3_file = args.directory + matched_files[block_number]['Cy3']
    
        normalize_call = base_normalize_call %(masliner_file, cy3_file, args.sequence_file, args.directory+protein_name)
        #print normalize_call
    else:
        raise Exception("There are no Cy3 files in this directory. Are you sure this is the program you want to use?")

    expected_file = args.directory + protein_name + '_regression.txt'
    if not isfile(expected_file):
        jobs_queue.append(normalize_call)
    else:
        print 'Normalization files already exist for %s. Using existing files.' %protein_name
        

print 'Submitting normalization jobs...'
SubmitJobsAndWait(jobs_queue)


#give the files a chance to show up in the directory.
sleep(20)

if not isfile(expected_file):
    raise Exception('The normalization procedure did not terminate successfuly.\n At least file \
                    %s is missing' %expected_file)

