"""
Takes the GenePix files from a single PBM experiment and automates the rest of the analysis.
Author: Luis Barrera
"""

from barrera import Parsers, General
import logging
import re
from subprocess import call
from os.path import isfile, basename, dirname, isdir
from os import makedirs

lr_lower, lr_upper = 200, 40000 #Lower and upper bound for what is considered the "linear range" by Masliner
min_linear_values = 20 #Minimum of lower values to require for a file to be used

PBM_suite_path = '/home/unix/barrera/Apps/PBM_analysis_suite/'
#PBM_suite_path = '/home/unix/barrera/Apps/AnastasiaPBMsuite/'

analysis_tools_path = '/home/unix/barrera/Apps/Analyze_PBM/'

arg_object = General.ArgParser.default()
arg_object.add_argument('directory', metavar = 'PBM_directory',
                        help='Path to directory with array data', type = str)
arg_object.add_argument('sequence_file', metavar = 'array_format_file',
                        help='File containing info about the array format \
                             (e.g., 8x60k_v14_amadid_30265_analysis.txt)', type = str)
arg_object.add_argument('--trim',
                        help = 'If given, the PWMs will be trimmed before making logos.',
                        action = 'store_true')
arg_object.add_argument('--beeml2',
                        help = 'If given, will generate a secondary BEEML motif.',
                        action = 'store_true')
arg_object.add_argument('--no_submit',
                        help = 'If given, will not submit any sub-jobs. Mostly useful for debugging.',
                        action = 'store_true')
arg_object.add_argument('--low_filter_thresh',
                        help = 'Any probes with an intensity below this value will be thrown out. By default, all probes included.', type = float, default = 0.0)
args = arg_object.Parse()

if not args.directory[-1] == '/':
    args.directory += '/'

files_in_dir = General.AllFilesInDir(args.directory)


files_to_use = sorted([filename for filename in files_in_dir if re.search('[1-8]-[1-8].gpr', filename) ])

if len(files_to_use) == 0:
    raise Exception('The given directory does not contain any .gpr files!')

matched_files = {}
experiment_name = None

def FixCombinatorial(combinatorial_file, low_thresh):

    old_combinatorial_file = combinatorial_file.replace('combinatorial','original_combinatorial')

    if not isfile(old_combinatorial_file):
        call('cp %s %s' %(combinatorial_file, old_combinatorial_file) , shell = True)

        print 'Copying raw combinatorial file to :%s' %(old_combinatorial_file)


    f_out = open(combinatorial_file, 'w')

    min_signal = min([float(fields[0]) for fields in Parsers.FlatFile(old_combinatorial_file)])

    # Make sure no array values are negative
    if min_signal < 0:
        floor = min_signal
    else:
        floor = 0.0

    for fields in Parsers.FlatFile(old_combinatorial_file): 

        mi = float(fields[0]) - floor
        seq = fields[1]

        if mi >= low_thresh:

            f_out.write('%s\t%s\n' %(mi, seq))


    f_out.close()



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

counter = 1
using_cy3 = True
for block_number in blocks:
        
    print "=========================\nProcessing block %s" %(block_number)
    
    if 'Cy3' in matched_files[block_number]:    
        print "-- Cy3 file -- \n%s" %matched_files[block_number]['Cy3']
    else:
        print '-- No Cy3 file found --'
        using_cy3 = False
    print "-- Alexa488 file(s) (should be sorted by intensity) --\n%s" %'\n'.join(matched_files[block_number]['Alexa488'])
    
    name = raw_input('\nDouble check the files were inferred correctly and enter the protein name (no slashes):\n') 
    
    if name in protein_names:
        raise Exception('Cannot use duplicate protein names. Please include an additional character to distinguish them.')

    if len(name) >= 1:
        protein_names.append(name)                                                       
    else:
        protein_names.append('block_%s' %counter)
    
    counter += 1
    print "Using name: %s" %protein_names[-1]

###############################################
###    Exclude files not in linear range    ###
###############################################

# This will exclude files that are problematic for Masliner,
# which happens when all ADJBSI values are outside
# its linear range. If this happens, Masliner will crash.


for block_number in blocks:
    blacklist = []

    for a488_file in matched_files[block_number]['Alexa488']:

        linear_values_found = 0
    
        for line in Parsers.FlatFile(args.directory + a488_file, split = False, skip = 2):

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
                                #signal_ix = fields.index('"F488 Median"')
                                signal_ix = fields.index('"F488 Median - B488"')
                                flag_ix = fields.index('"Flags"')

        
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

master_masliner_command = "perl "+PBM_suite_path+"masliner.txt -g1 %s -g2 %s -f1 488 -f2 488 -o %s -ll " + str(lr_lower) + " -lh " + str(lr_upper)
job_queue, masliner_files = [], []    

for block_number in blocks:
    
    print 'Now running masliner on block %s .' %block_number
    print matched_files[block_number]['Alexa488']
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
        
        
        #Check that Masliner output file has too many negative values
        
        
        
        
        if not isfile(masliner_filename):    
            call(masliner_call, shell = True)
            fixed_file = masliner_filename + '.fixed'
        
            if not isfile(masliner_filename):
                raise Exception('Masliner did not create the expected file: %s' %masliner_filename)
        
            #call('python %sFixADJBSI.py %s %s' %(analysis_tools_path, masliner_filename, fixed_file), shell = True)
        
            # if isfile(fixed_file):
                
            #     call('mv %s %s' %(fixed_file, masliner_filename), shell = True)
        
            
        else:
            print 'Masliner file already exists. Using %s' %masliner_filename
        
        if num_scans > 2:
            print 'Deleting intermediate files: \n' + '\n'.join(intermediate_files)
            for file_to_delete in intermediate_files:
                if isfile(file_to_delete):            
                    call ('rm %s' %file_to_delete, shell = True)
        
    else:
        print 'File %s already exists. Skipping masliner.' %masliner_filename

####################################################
###  Adjust BSI values to remove negative points ###
####################################################

# Increases ADJBSI values so that the smallest value is 0

adj_masliner_files = []
for masliner_file in masliner_files:

    #First pass: find the minimum value
    
    adj_masliner_file = masliner_file.replace('MaslinerOutput','AdjMaslinerOutput')
    
    adj_masliner_files.append(adj_masliner_file)
    
    if not isfile(adj_masliner_file):
        print 'Adjusting %s to prevent negative ADJBSI values.' %masliner_file

        min_adjbsi_value = 0.0
        for line in Parsers.FlatFile(masliner_file, split = False, skip = 2):
            
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
                            raise AssertionError('Found no A488 signal field in the GPR header.')
                        flag_ix = fields.index('"Flags"')


        

        mas_adj_out_handle = open(adj_masliner_file, 'w')

        #Second pass: write new values
        for line in Parsers.FlatFile(masliner_file, split = False, skip = 0):
            
            if len(line) > 50 and line[0] != '"': #Skip header lines

                    fields = line.strip().split('\t')
                    fields[signal_ix] = int(float(fields[signal_ix]) - min_adjbsi_value)
            
                    General.WriteTSV(fields, mas_adj_out_handle)

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

    if not isfile(expected_file):
        jobs_queue.append(normalize_call)
        comb_files_to_fix.append(expected_file)
    else:
        print 'Normalization files already exist for %s. Using existing files.' %protein_name
    
    #call(normalize_call, shell = True)    

if not args.no_submit:
    print 'Submitting normalization jobs...'
    General.SubmitJobsAndWait(jobs_queue, memory = 1000)

if not isfile(expected_file):
    raise Exception('The normalization procedure did not terminate successfuly.\n At least file \
                    %s is missing' %expected_file)
    
for comb_file in comb_files_to_fix:
    FixCombinatorial(comb_file, args.low_filter_thresh)


######################
###  Seed & Wobble ###
######################

combinatorial_files = [args.directory + x + '_combinatorial.txt' for x in protein_names]

jobs_queue = []

single_analysis_flags = ''

if args.trim is True:
    single_analysis_flags += ' --trim '
if args.beeml2 is True:
    single_analysis_flags += ' --secbeeml '

for protein_name,filename in zip(protein_names, combinatorial_files):
    
    snw_call = 'python '+ analysis_tools_path + 'AnalyzeSinglePBM.py -vv %s %s' %(single_analysis_flags, filename)
    
    expected_file = args.directory+'primary_'+protein_name+'_8mers_top_enrichment.txt'
    jobs_queue.append(snw_call)
    
    
print 'Submitting S&W + BEEML-PBM jobs...'

if not args.no_submit:

    General.SubmitJobsAndWait(jobs_queue, memory = 1000)

if not isfile(expected_file):
    raise Exception('Seed & Woble did not generate the expected files. Check logs for debugging information.')

###################################
###  Generate diagnostic plots  ###
###################################

kmer_df_file = args.directory + 'kmer_df_combined.txt'
probe_df_file = args.directory + 'probe_df_combined.txt'
ikaros_df_file = args.directory + 'ikaros_df_combined.txt'
positional_df_file = args.directory + 'positional_df_combined.txt'
combinatorial_df_file = args.directory + 'combinatorial_df_combined.txt'
cy3_r2_df_file = args.directory + 'cy3_r2_df_combined.txt'
# Make data frames for use in R

#These files are generated by the BEEML jobs called above


kmer_files_to_combine = [args.directory + _ + '_kmer_df.txt' for _ in protein_names]
probe_files_to_combine = [args.directory + _ + '_probe_df.txt' for _ in protein_names]
ikaros_files_to_combine = [args.directory + 'ikaros_' + _ + '_alpha_df.csv' for _ in protein_names]
positional_files_to_combine = [args.directory + 'ikaros_' + _ + '_positional_df.csv' for _ in protein_names] 
combinatorial_files_to_combine = [args.directory + _ + '_original_combinatorial.txt' for _ in protein_names] 
regression_files_to_combine = [args.directory + _ + '_regression.txt' for _ in protein_names] 

call("for file in %s; do awk -vOFS='\t' -vNAME=$(basename $file | sed 's/_original_combinatorial.txt//' ) '{print NAME,$1,$2}' $file; done > %s" \
    %(' '.join(combinatorial_files_to_combine), combinatorial_df_file), shell = True)


# Make Cy3 plots
call('for file in ' + ' '.join(regression_files_to_combine) + '; do printf "%s\t%s\n" $(basename $file | sed "s/_regression.txt//") \
    $(grep R.2 $file | cut -f 2 -d " "); done > ' + cy3_r2_df_file, shell = True)
#KW should be able to make this without ikaros
call('Rscript ' + analysis_tools_path + 'Make_Cy3_R2_Plot.R ' + cy3_r2_df_file + ' ' + args.directory + 'cy3_r2_summary.pdf', shell = True)


# Exception('Stop')

call('cat %s > %s' %(' '.join(kmer_files_to_combine), kmer_df_file), shell = True)
call('cat %s > %s' %(' '.join(probe_files_to_combine), probe_df_file), shell = True)

call("cat %s | grep -v 'Experiment' | sed 's/ikaros_//' > %s" %(' '.join(ikaros_files_to_combine), ikaros_df_file), shell = True)
call("cat %s | grep -v 'Experiment' | sed 's/ikaros_//' > %s" %(' '.join(positional_files_to_combine), positional_df_file), shell = True)


# Make R plots of probe intensities
#KW should be able to make this without ikaros
R_call = 'Rscript --verbose /home/unix/barrera/Apps/Analyze_PBM/Make_Combinatorial_Scatter.R %s %s %0.5f' \
         %(combinatorial_df_file, args.directory + 'probe_intensity_comparison', args.low_filter_thresh)

call(R_call, shell = True)



# Make R plots of BEEML output
#KW can't make this one - ikaros
R_call = 'Rscript --verbose /home/unix/barrera/Apps/Analyze_PBM/Make_BEEML_Scatter.R %s %s %s %ssummary %s' \
            %(probe_df_file, kmer_df_file, ikaros_df_file, args.directory, positional_df_file)

print 'Making BEEML scatter plots in R:'
print R_call

call(R_call, shell = True)

# if not isfile(args.directory + 'summary_r2_plot.pdf'):
#     call(R_call, shell = True)
# else:
#     print 'Found BEEML plots. Skipping generation.'


######################
###  Write report  ###
######################

insert_image = '<img src="%s.png", width = 1000><br>'

def Count8mersAboveThreshold(filename, threshold = 0.45):
    """Takes a S&W filename and counts the number of contiguous 8-mers that exceed a certain threshold."""

    if not isfile(filename):
        return 0
    else:
        return sum((1 for line in Parsers.FlatFile(filename, split = True, skip = 1) if float(line[2]) > threshold))

def MakeTable(f_out, protein_names, motif_type):
    
    f_out.write('<table>')
    
    names_to_use = protein_names[:]

    if len(names_to_use) % 2 != 0:
        names_to_use.append('__blank__')

    for i in range(0, len(names_to_use), 2):

        f_out.write('<tr>\n')
    
        f_out.write('<td>'+'<img src="%s.png"' %(names_to_use[i] + motif_type) +'width=550><br>'+'</td>\n')
        f_out.write('<td>'+'<img src="%s.png"' %(names_to_use[i+1] + motif_type) +'width=550><br>'+ '</td>\n')
        
        f_out.write('</tr>\n')
    f_out.write('</table>')
    f_out.write('<hr>\n')
    
num_above_thresh_primary = {protein:Count8mersAboveThreshold(args.directory + 'primary_%s_8mers_11111111.txt'%(protein), 0.45)
                    for protein in protein_names}

num_above_thresh_secondary = {protein:Count8mersAboveThreshold(args.directory + 'secondary_%s_8mers_11111111.txt'%(protein), 0.45)
                    for protein in protein_names}
        
        
with open(args.directory + 'report.html', 'w') as f_out:
    
    f_out.write('<html>\n')
    f_out.write('<center><font size = 13>')
    f_out.write('Array: %s<br>'% annotations)
    f_out.write('Done by: %s<br>'% person)
    f_out.write('</center></font>')
    
    MakeTable(f_out, protein_names, '_primary')
    MakeTable(f_out, protein_names, '_secondary')
    MakeTable(f_out, protein_names, '_BEEML')
    
    for protein in protein_names:
        
        for motif_type in ['_primary','_secondary','_BEEML']:
            f_out.write('<hr>')        
            f_out.write(insert_image %(protein + motif_type))
            f_out.write('<br>')
            if motif_type == '_primary':
                f_out.write('8-mers > 0.45: <b>%s</b>  <br>' %(num_above_thresh_primary[protein]))
            elif motif_type == '_secondary':
                f_out.write('8-mers > 0.45: <b>%s</b>  <br>' %(num_above_thresh_secondary[protein]))

    f_out.write('</html>')
    
print 'Generating .pdf report'

call('wkhtmltopdf-amd64 %sreport.html %sreport.pdf' %(args.directory, args.directory), shell = True)

#This is defined in the Linear Range section
if using_cy3:
    cy3_file = '-a %scy3_r2_summary.pdf' %args.directory
else:
    cy3_file = ''

if isfile(args.directory + 'report.pdf'):
    print 'E-mailing report...'
    call(' echo "%s" | mutt -s "PBM Report: %s" -a %sreport.pdf -a %ssummary_r2_plot.pdf -a %ssummary_position_plot.pdf %s lbarrera@gmail.com' \
         %(annotations, experiment_name, args.directory, args.directory, args.directory, cy3_file), shell = True)
    
    #This would be called for me if I was using mutt - but I'm not - it was being dumb.
    #The attachments would differ depending on the diagnostic plots computed.
    #call(' echo "%s" | mutt -s "PBM Report: %s" -a %sreport.pdf %s katy.weinand@gmail.com' \
     #    %(annotations, experiment_name, args.directory, cy3_file), shell = True)
    
    
