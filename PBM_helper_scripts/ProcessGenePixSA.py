"""
Takes the GenePix files from a single PBM experiment and automates the rest of the analysis.
Author: Luis Barrera
"""

#Cerebellum edition

#test data in /n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_anaylsis_suite/PBM_test_data/
#to run: python /n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_helper_scripts/ProcessGenePixSA.py /n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_analysis_suite/PBM_test_data/ /n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_analysis_suite/8x60k_v14_amadid_30265_analysis.txt

import logging
import re
import argparse
from argparse import ArgumentParser
from subprocess import call 
from os.path import isfile, basename, dirname, isdir 
from os import makedirs
from time import sleep
from submitO2 import WhichJobsAreRunning, SubmitJobsAndWait, SubmitJob

import Genomics


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
    

def FlatFile(input_filename,skip=0, split=True, delimiter='\t', comment_mark = '#'):

    in_handle=open(input_filename,'r')
    

    for i,line in enumerate(in_handle):
        
        if i>=skip and line[0] != comment_mark:            
            line=line.strip()
            if split==True:
                line=line.split(delimiter)        
            #Process line
            yield line
            
    in_handle.close()
    
#ADDED    
def FixCombinatorial(combinatorial_file, low_thresh):

    old_combinatorial_file = combinatorial_file.replace('combinatorial','original_combinatorial')

    if not isfile(old_combinatorial_file):
        call('cp %s %s' %(combinatorial_file, old_combinatorial_file) , shell = True)

        print 'Copying raw combinatorial file to: %s' %(old_combinatorial_file)


    f_out = open(combinatorial_file, 'w')

    min_signal = min([float(fields[0]) for fields in FlatFile(old_combinatorial_file)])

    # Make sure no array values are negative
    if min_signal < 0:
        floor = min_signal
    else:
        floor = 0.0

    for fields in FlatFile(old_combinatorial_file): 

        mi = float(fields[0]) - floor
        seq = fields[1]

        if mi >= low_thresh:

            f_out.write('%s\t%s\n' %(mi, seq))


    f_out.close()


arg_object = ArgParser.default()
arg_object.add_argument('directory', metavar = 'PBM_directory',
                        help='Path to directory with array data', type = str)
arg_object.add_argument('sequence_file', metavar = 'array_format_file',
                        help='File containing info about the array format \
                             (e.g., 8x60k_v14_amadid_30265_analysis.txt)', type = str)
arg_object.add_argument('--trim',
                        help = 'If given, the PWMs will be trimmed before making logos.',
                        action = 'store_true')
arg_object.add_argument('--beeml2',
                        help = 'If given, will generate a secondary BEEML motif. IN TESTING - DON\'T USE',
                        action = 'store_true')
arg_object.add_argument('--low_filter_thresh',
                        help = 'Any probes with an intensity below this value will be thrown out. By default, all probes included.',
                        type = float, default = 0.0)
arg_object.add_argument('-p', help = 'Which k-mer pattern to use: 8of10, 8of12 (not configured for 7mers)', default = '8of10', metavar = '8of10',
                        choices = ['8of10', '8of12'])
arg_object.add_argument('-noSW', help = 'If given, Seed and Wobble will not run and no report will be generated.', action = 'store_true')
arg_object.add_argument('-email', help = 'If given, report with be send to this email.', type=str)
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
    
    if 'Alexa488' in filename:
        if experiment_name is None:
            experiment_name = filename.split('Alexa488')[0]
            print 'Assuming all of the Alexa488 files have the prefix: %s' %experiment_name
        else:
            assumed_prefix = filename.split('Alexa488')[0]
            if experiment_name != assumed_prefix:
                raise Exception('The prefix of the file (before Alexa488) does not match: %s vs. %s' \
                                %(assumed_prefix, experiment_name))
        
    block_match = re.search('_[1-9]-[1-9].gpr', filename)
    if not block_match:
        raise Exception('The file %s does not contain a description of the block it belongs to.' %filename)
    else:
        block_number = block_match.group()[1:4]
    
    if not block_number in matched_files:
        matched_files[block_number]={'Alexa488':[]}
    
    if re.search('Alexa488_lp[0-9]*pg[0-9]*_', filename):
        matched_files[block_number]['Alexa488'].append(filename)
    elif re.search('_[cC]y3_', filename):
        matched_files[block_number]['Cy3']=filename
    else:
        raise Exception('File does not contain properly formatted intensity information: %s' %filename)
    
        
#Iterate over all blocks and ask for names
protein_names = []
blocks = sorted(matched_files.keys())

person = raw_input('Enter the name of the person that performed the experiment (e.g., Anastasia): ')
annotations = raw_input('Enter any relevant information about the array (e.g., 032912_v11_273_2, AMADID # 030236): ')

useNames = raw_input('Are all 8 protein names in the filename? If yes, input the number of underscores before the first protein name. If no, press enter: ')

if len(useNames)>=1:
    offset = int(useNames)-1

counter = 1
fnames = experiment_name.split('_')
for block_number in blocks:
    index = int(block_number[0])
        
    print "=========================\nProcessing block %s" %(block_number)
    
    if 'Cy3' in matched_files[block_number]:    
        print "-- Cy3 file -- \n%s" %matched_files[block_number]['Cy3']
    else:
        print '-- No Cy3 file found --'
    print "-- Alexa488 file(s) (should be sorted by intensity) --\n%s" %'\n'.join(matched_files[block_number]['Alexa488'])
    
    if len(useNames)>=1:
        name = raw_input('\nCurrent name of block %s: %s\nIf this is wrong, input correct name (or else hit enter): \n'%(index, fnames[index+offset]))
        
        """original version
        if(len(name)>=1):
            if name in protein_names:
                raise Exception('Cannot use duplicate protein names. Please include an additional character to distinguish them.')
            protein_names.append(name)
        else:
            if fnames[index+2] in protein_names:
                raise Exception('Cannot use duplicate protein names. Please include an additional character to distinguish them.')
            protein_names.append(fnames[index+2])
        """
        #condensed version
        if(len(name)<1):
            name=fnames[index+offset]
        if name in protein_names:
            raise Exception('Cannot use duplicate protein names. Please include an additional character to distinguish them.')
        protein_names.append(name)
        
    else:
        name = raw_input('\nEnter the protein name (no slashes):\n') 
        
        if name in protein_names:
            raise Exception('Cannot use duplicate protein names. Please include an additional character to distinguish them.')
    
        if len(name) >= 1:
            protein_names.append(name)                                                       
        else:
            protein_names.append('block_%s' %counter)
    
    counter += 1
    print "Using name: %s" %protein_names[-1]
    
#stop = raw_input('stop here? y or n. REMOVE AFTER TESTING!!! ')

#if stop=='y':
#    raise Exception("You wanted this. Don't come crying to me.")
    
###############################################
###    Exclude files not in linear range    ###
###############################################

# This will exclude files that are problematic for Masliner,
# which happens when all ADJBSI values are outside
# its linear range. If this happens, Masliner will crash.

lr_lower, lr_upper = 200, 40000 #Lower (-ll) and upper (-lh) bound for what is considered the "linear range" by Masliner
min_linear_values = 20 #Minimum of lower values to require for a file to be used

for block_number in blocks:
    blacklist = []

    for a488_file in matched_files[block_number]['Alexa488']:

        linear_values_found = 0
    
        for line in FlatFile(args.directory + a488_file, split = False, skip = 2):

                if len(line) > 1:
                        if line[0] != '"': #Skip header lines

                            fields = line.strip().split('\t')

                            a488_val = float(fields[signal_ix])
                                       
                            if a488_val >= lr_lower and a488_val <= lr_upper:
                                if fields[flag_ix] == '0':
                                    linear_values_found += 1
                                    if linear_values_found >= min_linear_values:
                                        break
                                                
                        else:
                            if re.match('"Block"\t"Column"', line):
                                fields = line.strip().split('\t')
                                if '"F488 Median - B488"' in fields:
                                    signal_ix = fields.index('"F488 Median - B488"')
                                else:
                                    #KW added 10/27/17 - was just the AssertionError
                                    regex = [x for x in fields if re.match('"?F[0-9]+ Median - B[0-9]+"?',x)]
                                    if len(regex)==1:
                                        signal_ix = fields.index(regex[0])
                                    else:
                                        raise AssertionError('Found no A488 signal field in the GPR header.')
                                print "Using the following column as the signal value: %s"%(fields[signal_ix])
                                flag_ix = fields.index('"Flags"')
                            
                            """ #KW changed 10/27/17
                            if re.match('"Block"\t"Column"', line):
                                print a488_file #KW added
                                print line #KW added
                                fields = line.strip().split('\t')
                                #signal_ix = fields.index('"F488 Median"')
                                signal_ix = fields.index('"F488 Median - B488"')
                                flag_ix = fields.index('"Flags"')
                            """

        
        if linear_values_found >= min_linear_values:
            continue
            #print 'File %s has enough values in linear range.' %a488_file
        else:
            print 'File %s does not have enough values in linear range. Therefore, it will not be used.' %a488_file
            blacklist.append(a488_file)


    for bad_file in blacklist:
        matched_files[block_number]['Alexa488'].remove(bad_file)

######################
###    Masliner    ###
######################

master_masliner_command = "perl "+PBM_suite_path+"masliner.txt -g1 %s -g2 %s -f1 488 -f2 488 -o %s -ll 200 -lh 40000"
job_queue, masliner_files = [], []    

for block_number in blocks:
    
    print 'Now running masliner on block %s .' %block_number
    
    #What happens if all intensity files were bad files (not in linear range), so num_scans would equal 0, and there wouldn't be any existing masliner files...
    num_scans = len(matched_files[block_number]['Alexa488']) 
    
    masliner_filename = '%s%sAlexa488_MaslinerOutput_%s.GPR' %(args.directory, experiment_name, block_number)
    masliner_files.append(masliner_filename)
    
    if not isfile(masliner_filename):
        
        if num_scans > 2:
            
            intermediate_files = ['%s%s_output%s.GPR' %(args.directory, block_number, i+1) for i in range(num_scans-2) ]
            remaining_files = [args.directory+x for x in matched_files[block_number]['Alexa488'][2:]]
            output_files = intermediate_files[1:] + [masliner_files[-1]]
            
            assert len(intermediate_files) == len(remaining_files)
            assert len(output_files) == len(remaining_files)
            
            #First call is always the same
            masliner_call = master_masliner_command \
                            %(args.directory+matched_files[block_number]['Alexa488'][0],
                              args.directory+matched_files[block_number]['Alexa488'][1],
                              intermediate_files[0])
            
            
            #Rest of calls
            for intermediate_file, data_file, out_file in zip(intermediate_files, remaining_files, output_files):
                masliner_call += '; '+ master_masliner_command %(intermediate_file, data_file, out_file)
                
                
        elif num_scans == 2:
            
            masliner_call = master_masliner_command \
                            %(args.directory + matched_files[block_number]['Alexa488'][0],
                              args.directory + matched_files[block_number]['Alexa488'][1],
                              masliner_files[-1])
    
        elif num_scans == 1:
            
            masliner_call = 'cp %s%s %s' %(args.directory, matched_files[block_number]['Alexa488'][0],
                                       masliner_files[-1])
        
        if num_scans>0:
            if not isfile(masliner_filename):    
                call(masliner_call, shell = True)
                #ADDED
                fixed_file = masliner_filename + '.fixed'
                if not isfile(masliner_filename):
                    raise Exception('Masliner did not create the expected file: %s' %masliner_filename)
            else:
                print 'Masliner file already exists. Using %s' %masliner_filename
            
            if num_scans > 2:
                print 'Deleting intermediate files: \n' + '\n'.join(intermediate_files)
                for file_to_delete in intermediate_files:
                    if isfile(file_to_delete):            
                        call ('rm %s' %file_to_delete, shell = True)
        else:
            print "There are not enough valid scans for %s; it will be ignored in further analysis."%(block_number)
        
    else:
        print 'File %s already exists. Skipping masliner.' %masliner_filename
    
#ADDED
####################################################
###  Adjust BSI values to remove negative points ###
####################################################

# Increases ADJBSI values so that the smallest value is 0

adj_masliner_files = []
for masliner_file in masliner_files:

    #First pass: find the minimum value
    
    adj_masliner_file = masliner_file.replace('MaslinerOutput','AdjMaslinerOutput')
    
    adj_masliner_files.append(adj_masliner_file)
    
    if not isfile(adj_masliner_file) and isfile(masliner_file):
        print 'Adjusting %s to prevent negative ADJBSI values.' %masliner_file

        min_adjbsi_value = 0.0
        for line in FlatFile(masliner_file, split = False, skip = 2):
            
            if len(line) > 1:
                if line[0] != '"': #Skip header lines

                    fields = line.strip().split('\t')
                    adjbsi_value = float(fields[signal_ix])
                       
                    if adjbsi_value < min_adjbsi_value:
                        if fields[flag_ix] == '0':
                            min_adjbsi_value = adjbsi_value

                else:
                    if re.match('"Block"\t"Column"', line):
                        fields = line.strip().split('\t')
                        if 'ADJBSI' in fields:
                            signal_ix = fields.index('ADJBSI')
                        elif '"F488 Median - B488"' in fields:
                            signal_ix = fields.index('"F488 Median - B488"')
                        else:
                            #KW added 10/27/17 - was just the AssertionError
                            regex = [x for x in fields if re.match('"?F[0-9]+ Median - B[0-9]+"?',x)]
                            if len(regex)==1:
                                signal_ix = fields.index(regex[0])
                            else:
                                raise AssertionError('Found no A488 signal field in the GPR header.')
                        print "Using the following column as the signal value: %s"%(fields[signal_ix])
                        flag_ix = fields.index('"Flags"')


        

        mas_adj_out_handle = open(adj_masliner_file, 'w')

        #Second pass: write new values
        for line in FlatFile(masliner_file, split = False, skip = 0):
            
            if len(line) > 50 and line[0] != '"': #Skip header lines

                    fields = line.strip().split('\t')
                    fields[signal_ix] = int(float(fields[signal_ix]) - min_adjbsi_value)
            
                    Genomics.WriteTSV(fields, mas_adj_out_handle)

            else:
                    mas_adj_out_handle.write(line + '\n')
                                
        mas_adj_out_handle.close()


######################
###  Normalization ###
######################

base_normalize_call = 'perl '+PBM_suite_path+'normalize_agilent_array.pl -i %s -c %s -s %s -o %s'
mcy3_normalize_call = 'perl '+PBM_suite_path+'normalize_agilent_array-cy3.pl -i %s -s %s -o %s'

jobs_queue, comb_files_to_fix = [], []

for block_number, masliner_file, protein_name in zip(blocks, adj_masliner_files, protein_names):
    
    if 'Cy3' in matched_files[block_number]:    
        cy3_file = args.directory + matched_files[block_number]['Cy3']
    
        normalize_call = base_normalize_call %(masliner_file, cy3_file, args.sequence_file, args.directory+protein_name)
    else:
        
        normalize_call = mcy3_normalize_call %(masliner_file, args.sequence_file, args.directory+protein_name)

    expected_file = args.directory + protein_name + '_combinatorial.txt'
    if not isfile(expected_file) and isfile(masliner_file):
        jobs_queue.append(normalize_call)
        #ADDED
        comb_files_to_fix.append(expected_file)
    else:
        print 'Normalization files already exist for %s. Using existing files.' %protein_name
    
    #call(normalize_call, shell = True)    

print 'Submitting normalization jobs...'
SubmitJobsAndWait(jobs_queue)

if len(jobs_queue)!=0:
    #give the files a chance to show up in the directory.
    sleep(20)

#KW changed getting rid of the Exception and putting it in the combinatorial fixing step if file doesn't exist.
#if not isfile(expected_file):
#    raise Exception('The normalization procedure did not terminate successfuly.\n At least file \
#                    %s is missing' %expected_file)

#ADDED
for comb_file in comb_files_to_fix:
    if isfile(comb_file):
        FixCombinatorial(comb_file, args.low_filter_thresh)
    else:
        print "The normalization procedure did not produce the combinatorial file: %s"%(comb_file)
    
    
if args.noSW is True:
    print "You have chosen not to run the Seed and Wobble portion of the pipeline. No png files or report will be generated."
else:        
    ######################
    ###  Seed & Wobble ###
    ######################
    
    combinatorial_files = [args.directory + x + '_combinatorial.txt' for x in protein_names]
    
    jobs_queue = []
    
    single_analysis_flags = ''
    
    if args.trim is True:
        single_analysis_flags += ' --trim '
    if args.beeml2 is True:
        single_analysis_flags += ' -secbeeml '
    if args.p is not None:
        single_analysis_flags += ' -p %s ' %args.p
    
    expected_file = '' #placeholder
    for protein_name,filename in zip(protein_names, combinatorial_files):
                
        #KW added isfile 10/27/17
        if isfile(filename):
            expected_file = args.directory+'primary_'+protein_name+'_8mers_top_enrichment.txt'
            
            if not isfile(expected_file):
                snw_call = 'python '+ analysis_tools_path + 'AnalyzeSinglePBM.py -vv %s %s' %(single_analysis_flags, filename)
                
                #KW moved expected_file definition to copy with isfile 10/27/17
                jobs_queue.append(snw_call)
        
        
    print 'Submitting S&W + BEEML-PBM jobs...'
    SubmitJobsAndWait(jobs_queue)
    
    if len(jobs_queue)!=0:
        #Trim pngs are not showing up as existing...
        sleep(20)
    
    if not isfile(expected_file):
        raise Exception('Seed & Wobble did not generate the expected files. Check logs for debugging information.')
    
    ######################
    ###  Write report  ###
    ######################
    
    print "Finished S&W jobs. Writing report."
    
    insert_image = '<img src="%s.png" width=1100><br>'
    
    def Count8mersAboveThreshold(filename, threshold = 0.45):
        #Takes a S&W filename and counts the number of contiguous 8-mers that exceed a certain threshold.
        if not isfile(filename):
            return 0
        else:
            return sum((1 for line in FlatFile(filename, split = True, skip = 1) if line[2]!="NA" and float(line[2]) > threshold))
    
    def MakeTable(f_out, protein_names, motif_type, height=''):
        
        f_out.write('<table>')
        
        names_to_use = protein_names[:]
    
        if len(names_to_use) % 2 != 0:
            names_to_use.append('__blank__')
    
        for i in range(0, len(names_to_use), 2):
            pngName1 = args.directory + names_to_use[i] + motif_type + '.png'
            pngName2 = args.directory + names_to_use[i+1] + motif_type + '.png'
            #pngName1 = names_to_use[i] + motif_type + '.png'
            #pngName2 = names_to_use[i+1] + motif_type + '.png'
            
            if not isfile(pngName1) and '_trim' in motif_type:
                logging.warning('Trimmed motif file does not exist: %s' %pngName1)
                pngName1 = "%sNA.png"%(analysis_tools_path)
            if not isfile(pngName2) and '_trim' in motif_type:
                logging.warning('Trimmed motif file does not exist: %s' %pngName2)
                pngName2 = "%sNA.png"%(analysis_tools_path)
    
            f_out.write('<tr>\n')
            
            f_out.write('<td>'+'<figure><figcaption>%s</figcaption>'%("8-mers > 0.45: "+str(Count8mersAboveThreshold(args.directory + '%s_%s_8mers_11111111.txt'%(motif_type[1:] if "trim" not in motif_type else motif_type[1:-5],names_to_use[i]), 0.45)) if 'BEEML' not in motif_type and names_to_use[i]!="__blank__" else "8-mers > 0.45: N/A")+'<img src="%s" ' %pngName1 +'width=450 %s>' %height+'</figure></td>\n')
            f_out.write('<td>'+'<figure><figcaption>%s</figcaption>'%("8-mers > 0.45: "+str(Count8mersAboveThreshold(args.directory + '%s_%s_8mers_11111111.txt'%(motif_type[1:] if "trim" not in motif_type else motif_type[1:-5],names_to_use[i+1]), 0.45)) if 'BEEML' not in motif_type and names_to_use[i+1]!="__blank__" else "8-mers > 0.45: N/A")+'<img src="%s" ' %pngName2 +'width=450 %s></figure>' %height+ '</td>\n')
            
            f_out.write('</tr>\n')
        f_out.write('</table>')
        f_out.write('<hr>\n')
        
    num_above_thresh_primary = dict((protein,Count8mersAboveThreshold(args.directory + 'primary_%s_8mers_11111111.txt'%(protein), 0.45))
                        for protein in protein_names)
    
    num_above_thresh_secondary = dict((protein,Count8mersAboveThreshold(args.directory + 'secondary_%s_8mers_11111111.txt'%(protein), 0.45))
                        for protein in protein_names)
            
            
    with open(args.directory + 'report.html', 'w') as f_out:
        
        f_out.write('<html>\n')
        f_out.write('<center><font size = 13>')
        f_out.write('Array: %s<br>'% annotations)
        f_out.write('Done by: %s<br>'% person)
        f_out.write('</center></font>')
        
        MakeTable(f_out, protein_names, '_primary','height=110')
        MakeTable(f_out, protein_names, '_secondary','height=110')
        f_out.write('<br><br><br><br>')
        MakeTable(f_out, protein_names, '_BEEML')
        
        for protein in protein_names:
            
            for motif_type in ['_primary','_secondary','_BEEML']:
                f_out.write('\n<hr>')        
                f_out.write(insert_image %(protein + motif_type))
                f_out.write('<br>')
                if motif_type == '_primary':
                    f_out.write('8-mers > 0.45: <b>%s</b>  <br>' %(num_above_thresh_primary[protein]))
                elif motif_type == '_secondary':
                    f_out.write('8-mers > 0.45: <b>%s</b>  <br>' %(num_above_thresh_secondary[protein]))
        
        if args.trim is True:
            f_out.write('<hr><br><h3>Trimmed versions.</h3><br><br><br><br><br><br>')
            MakeTable(f_out, protein_names, '_primary_trim','height=200')
            f_out.write('<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>')
            MakeTable(f_out, protein_names, '_secondary_trim', 'height=200')
            f_out.write('<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>')
            MakeTable(f_out, protein_names, '_BEEML_trim', 'height=200')
    
        f_out.write('</html>')
    
    print 'Generating .pdf report'
    
    call('wkhtmltopdf %sreport.html %sreport.pdf' %(args.directory, args.directory), shell = True)
    
    if args.email is not None:
        print 'Emailing .pdf report'
        
        call('python %sKWsendPDF.py -to %s -pw Bulyk060115 -attach %sreport.pdf'%(analysis_tools_path,args.email,args.directory), shell=True)

    """
    #not an option on O2, or at least, not an easy option.
    print 'Generating precurser reports.'
    brainCmd1 = 'python %smakeBrainReport3Types.py %s %s %s SW'%(analysis_tools_path,args.directory,annotations,person)
    call(brainCmd1,shell=True) #doesn't need to be run on Cerebellum
    brainCmd2 = 'python %smakeBrainReport3Types.py %s %s %s BL'%(analysis_tools_path,args.directory,annotations,person)
    call(brainCmd2,shell=True) #doesn't need to be run on Cerebellum
    
    print 'All files from this run have been generated including a report. \nHowever, you can visit http://thebrain.bwh.harvard.edu/KW_logos/buildReport.php to determine which logo orientations to put in a final report. \nThe directory to use is %s and the type of logo can be either SW or BEEML. \nAn email will be sent to you containing the report.'%args.directory
    """

print "Finished Program."
    
