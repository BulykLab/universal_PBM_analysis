## Automates the analysis of a single PBM dataset
## Generates primary, secondary and BEEML PWMs and logos
## Author: Luis Barrera

#called from ProcessGenePixSA.py Seed-n-Wobble section

import Genomics
import logging
from subprocess import call
from os.path import isfile, basename, dirname
import argparse
from argparse import ArgumentParser
from submitO2 import WhichJobsAreRunning, SubmitJobsAndWait, SubmitJob, WaitForJobs

## Paths

#KW Change paths
PBM_suite_path = '/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_analysis_suite/'

BEEML_path = '/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/BEEML_PBM_Example/'
pattern_suffix = PBM_suite_path + 'pattern_files/patterns_'

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
        
        #-v SGE_DEBUG_LEVEL=? is an environment variable in Cerebellum
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


## Parse arguments

arg_obj = ArgParser.default()
arg_obj.add_argument('probe_file', metavar = 'TF_combinatorial.txt', help='File with probes and intensities')
arg_obj.add_argument('-p', help = 'Which k-mer pattern to use', default = '8of10', metavar = '8of10',
                        choices = ['8of10', '8of12'])
arg_obj.add_argument('-t', help = 'Pattern type (e.g., 4x44k_all_8mer)', default = '4x44k_all_8mer',
                        choices = ['4x44k_all_8mer'])
arg_obj.add_argument('-k', help = 'K-mer length to use (typically 8)', default = 8, type = int)
arg_obj.add_argument('-o','--override', help = 'If given, will override existing files.',
                        default = False, action = 'store_true')
arg_obj.add_argument('--trim', help = 'If given, will trim PWMs before making logos.', action='store_true')
arg_obj.add_argument('-skipbeeml', help = 'Will skip running BEEML-PBM', action='store_true')
arg_obj.add_argument('-secbeeml', help = 'Will perform a "secondary motif" analysis with BEEML PWMs. IN TESTING - DON\'T USE', action='store_true')
arg_obj.add_argument('--skip_ikaros', help = 'Will skip running IKAROS', action='store_true')
args = arg_obj.parse_args()

#Configure logging

if not args.debuglog is None:
    logging.basicConfig(filename = args.debuglog)

if args.verbose == 1:
    logging.basicConfig(level = logging.INFO)
elif args.verbose == 2:
    logging.basicConfig(level = logging.DEBUG)

## Run different S&W commands

## Include an option for the current directory instead of the directory in which the data file is stored
running_jobs = []

if '_combinatorial.txt' not in args.probe_file:
    raise Exception('For proper labeling of output files, the input file must be named *_combinatorial.txt')

output_dir = dirname(args.probe_file)+'/'

file_prefix = basename(args.probe_file).split('_combinatorial')[0]

primary_SNW_call = 'perl %sseed_and_wobble.pl %s %s %s %s %s'  \
                    %(PBM_suite_path, args.probe_file, args.k, pattern_suffix+args.p+'.txt',
                      pattern_suffix+args.t+'.txt', output_dir+'primary_'+file_prefix)
logging.info('Generating primary Seed & Wobble files...')
logging.debug('Command called: %s' %primary_SNW_call)

primary_pwm_file = output_dir+'primary_'+file_prefix+'_8mers_pwm.txt'



if not isfile(primary_pwm_file) or args.override is True:
    #call(primary_SNW_call, shell = True, stdout = None, stderr = None)
    #KW Change KB->MB for all memory in script
    snw_primary_job_id = SubmitJob(primary_SNW_call, memory = '1G', timeLim="0-00:45")
    running_jobs.append(snw_primary_job_id)

else:
    logging.warning('Warning: the primary file already exists. Using the existing file.')
    snw_primary_job_id = None





reranked_file = output_dir + file_prefix + '_reranked.txt'

rerank_call = 'perl %srerank.pl %s %s %s' %(PBM_suite_path, args.probe_file, primary_pwm_file,
                                            reranked_file)

logging.info('Re-ranking probe intensities...')
logging.debug(rerank_call)



if not isfile(reranked_file) or args.override is True:
    #call(rerank_call, shell = True, stdout = None)

    rerank_job_id = SubmitJob(rerank_call, depends = snw_primary_job_id, memory = '500M',timeLim='0-00:10')
    running_jobs.append(rerank_job_id)

else:
    logging.warning('Warning: the re-ranked file already exists. Using the existing file.')
    rerank_job_id = None

secondary_SNW_call = 'perl %sseed_and_wobble.pl %s %s %s %s %s' %(PBM_suite_path, reranked_file, args.k,
                                                              pattern_suffix+args.p+'.txt',
                                                              pattern_suffix+args.t+'.txt',
                                                              output_dir+'secondary_'+file_prefix)
logging.info('Generating secondary S&W files...')
logging.debug(secondary_SNW_call)

secondary_pwm_file = output_dir + 'secondary_'+file_prefix+'_8mers_pwm.txt'

if not isfile(secondary_pwm_file) or args.override is True:
    #call(secondary_SNW_call, shell = True, stdout = None)


    secondary_job_id = SubmitJob(secondary_SNW_call, depends = rerank_job_id, memory = '1G', timeLim='0-00:45')
    running_jobs.append(secondary_job_id)

else:
    logging.warning('Warning: the secondary file already exists. Using the existing file.')
    secondary_job_id = None


##################
### Run IKAROS ###
##################
#KW commented out on Luis's Instruction

#ikaros_filename = '/home/unix/barrera/Projects/2013/new_PBM_analysis_tool/IKAROS.py' 

#ikaros_out_prefix = output_dir + 'ikaros_' + file_prefix


#ikaros_call = 'python %s %s %s %s --predicted_combinatorial --trim_bottom 1000' \
 #             %(ikaros_filename, output_dir+'primary_'+file_prefix + '_8mers_11111111.txt', 
 #               output_dir + file_prefix + '_combinatorial.txt' , ikaros_out_prefix)

#print 'Need to check that IKAROS output is already present.'

#logging.info('Submitting IKAROS job...')
#logging.debug(ikaros_call)

#if not isfile(ikaros_out_prefix + '_alpha_df.csv') or args.override is True:
 #   if args.skip_ikaros is False:
  #      ikaros_job_id = SubmitJob(ikaros_call, depends = snw_primary_job_id, memory = 4)
   #     running_jobs.append(ikaros_job_id)
    #else:
     #   print 'Skippping IKAROS jobs because of --skip_ikaros flag.'


########################
### BEEML processing ###
########################

if args.skipbeeml is False:

    logging.info('Selecting the top seed primary motif for BEEML optimization')
    
    # BEEML doesn't like having multiple motifs in one file (which S&W gives)
    # This eliminates the rest of the matrices from the file
    
    start_reading, num_rows = False, 0
    top_pwm_file = output_dir + 'primary_'+file_prefix+'_top_pwm.txt'
    
    WaitForJobs([snw_primary_job_id]) # Need S&W primary before file is generated
    
    with open(top_pwm_file, 'w') as f_out:
        first_line = True
        for line in open(primary_pwm_file, 'r'):
            if first_line:
                f_out.write(line)
                first_line = False
            
            #print line
            if len(line) > 1 and start_reading is True:
                f_out.write(line)
                num_rows += 1
            
            if 'Probability matrix' in line and num_rows == 0:
                start_reading = True
                    
            if num_rows == 4:
                start_reading = False
    
    # Now that we have a seed matrix, BEEML can be run
    
    
    BEEML_pwm_file = output_dir + file_prefix + '_BEEML_pwm.txt'
    BEEML_call = "Rscript %sGeneratePWM.R %s %s %s" %(BEEML_path, args.probe_file, top_pwm_file,
                                                                   output_dir + file_prefix)
    logging.info(' Now processing files with BEEML-PBM')
    logging.debug(BEEML_call)
    
    if not isfile(BEEML_pwm_file) or args.override is True:
        
        #call(BEEML_call, shell = True)

        beeml_job_id = SubmitJob(BEEML_call, depends = snw_primary_job_id, memory = '2G',timeLim='0-00:30')
        running_jobs.append(beeml_job_id)

    else:
        logging.warning('Warning: the BEEML-PBM file already exists. Using the existing file.')
        beeml_job_id = None

#KW Secondary BEEML-PWMs are not fully tested, so we are skipping them for now.
args.secbeeml = False
if args.secbeeml: #If True, will try to generate a BEEML secondary motif
    
    logging.info('Generate BEEML secondary')
    
    BEEML_pfm_file = output_dir + 'pfm_'+ basename(BEEML_pwm_file)
    
    BEEML_pwm.SavePFM(BEEML_pfm_file)
    
    beeml_reranked_file = output_dir + file_prefix+'_BEEML_reranked.txt'
    beeml_rerank_call = 'perl %srerank.pl %s %s %s' %(PBM_suite_path, args.probe_file, BEEML_pfm_file,
                                            beeml_reranked_file)
    if not isfile(beeml_reranked_file) or args.override is True:
        call(beeml_rerank_call, shell = True)
    else:
        logging.info('BEEML reranked file already exists. Using existing file.')
        
    secondary_BEEML_call = "R CMD BATCH '--args %s %s %s' %sGeneratePWM.R" \
                           %(beeml_reranked_file, BEEML_pfm_file,
                             output_dir+file_prefix+'_2nd', BEEML_path)
    
    BEEML_2nd_file = output_dir + file_prefix+'_2nd_BEEML_pwm.txt'
    if not isfile(BEEML_2nd_file) or args.override is True:
        #call(secondary_BEEML_call, shell = True)

        beeml_sec_job_id = SubmitJob(secondary_BEEML_call, depends = beeml_job_id, memory = '2G',timeLim='0-00:30')
        running_jobs.append(beeml_sec_job_id)

    else:
        logging.info('BEEML secondary PWM file already exists. Using existing file.')
        beeml_sec_job_id = None
    

print '--- CHECKPOINT: waiting for all submitted jobs to finish... ---'
WaitForJobs(running_jobs)

#Generate primary logo
logging.info('Generating primary logo... (warnings can be ignored if they occur)')
primary_pwm = Genomics.PWM.fromSNW(primary_pwm_file)
primary_pwm.MakeLogo(output_dir + file_prefix+'_primary.png', '%s - Primary (%s - %s)' \
                         %(file_prefix, primary_pwm.top_seed, primary_pwm.top_escore) )
primary_pwm.MakeLogo(output_dir + file_prefix+'_primary_rc.png', '%s - Primary (%s - %s)' \
                         %(file_prefix, primary_pwm.top_seed, primary_pwm.top_escore), reverse_complement = True)
if args.trim is True:
    cont = primary_pwm.Trim()
    if cont:
        primary_pwm.MakeLogo(output_dir + file_prefix+'_primary_trim.png', '%s P'%file_prefix)
        primary_pwm.MakeLogo(output_dir + file_prefix+'_primary_rc_trim.png', '%s P_rc'%file_prefix)
    

#Generate secondary logo
logging.info('Generating secondary logo... (warnings can be ignored if they occur)')
secondary_pwm = Genomics.PWM.fromSNW(secondary_pwm_file)
secondary_pwm.MakeLogo(output_dir + file_prefix+'_secondary.png', '%s - Secondary (%s - %s)' \
                         %(file_prefix, secondary_pwm.top_seed, secondary_pwm.top_escore) )
secondary_pwm.MakeLogo(output_dir + file_prefix+'_secondary_rc.png', '%s - Secondary (%s - %s)' \
                         %(file_prefix, secondary_pwm.top_seed, secondary_pwm.top_escore), reverse_complement = True)
if args.trim is True:
    cont = secondary_pwm.Trim()
    if cont:
        secondary_pwm.MakeLogo(output_dir + file_prefix+'_secondary_trim.png', '%s S'%file_prefix)
        secondary_pwm.MakeLogo(output_dir + file_prefix+'_secondary_rc_trim.png', '%s S_rc'%file_prefix)
    

if args.skipbeeml is False:
    #Generate BEEML logo
    logging.info('Generating BEEML logo... (warnings can be ignored if they occur)')
    BEEML_pwm = Genomics.PWM.fromBEEML(BEEML_pwm_file)
    BEEML_pwm.MakeLogo(output_dir + file_prefix+'_BEEML.png', '%s - BEEML' %file_prefix)
    BEEML_pwm.MakeLogo(output_dir + file_prefix+'_BEEML_rc.png', '%s - BEEML' %file_prefix, reverse_complement = True)
    if args.trim is True:
        cont = BEEML_pwm.Trim()
        if cont:
            BEEML_pwm.MakeLogo(output_dir + file_prefix+'_BEEML_trim.png', '%s B'%file_prefix)
            BEEML_pwm.MakeLogo(output_dir + file_prefix+'_BEEML_rc_trim.png', '%s B_rc'%file_prefix, reverse_complement = True)


if args.secbeeml:
    BEEML_2nd_pwm = Genomics.PWM.fromBEEML(BEEML_2nd_file)
    if args.trim is True:
        BEEML_2nd_pwm.Trim()
    
    BEEML_2nd_pwm.MakeLogo(output_dir + file_prefix+'_BEEML_2nd.png', '%s - BEEML (2nd)' %file_prefix)
    BEEML_2nd_pwm.MakeLogo(output_dir + file_prefix+'_BEEML_2nd_rc.png', '%s - BEEML (2nd)' %file_prefix, reverse_complement = True)
