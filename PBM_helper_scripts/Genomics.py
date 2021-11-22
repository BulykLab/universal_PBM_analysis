### Imports from default Python libraries
import shlex
import subprocess
import re
import math
import numpy as np
import logging


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

def GenerateTempFilename(temp_dir = '/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_tmp_files/', base_string = None):
    """Returns a full name for a temporary file that is hashed such that it is ~guaranteed
    to be unique.
    The file itself exists in the Broad's temporary high-performance storage.

    If base_string is specified, it will use that instead of a randomly generated one.
    I.e., the same temporary file will be returned given the same base_string.

    """
    
    import time
    import hashlib
    from random import random
    from os.path import isdir
    from os import makedirs

    if not isdir(temp_dir):
        makedirs(temp_dir)
    
    if base_string is None:
    
        return temp_dir + hashlib.md5(time.ctime() + str(random())).hexdigest()

    else:

        return temp_dir + hashlib.md5(base_string).hexdigest()

def WriteTSV(list_to_print, file_handle, rounding_decimals = 4):
    
    #Pre-process list by joining nested lists with commas
    list_to_print = [','.join(x) if type(x) is list else x for x in list_to_print]
    
    #Round the floating point values and convert everything to strings
    list_to_print = [str(round(x, rounding_decimals)) if type(x) is float else str(x) for x in list_to_print ]
    file_handle.write('\t'.join(list_to_print)+'\n')

def RCmatrix(matrix):
        #Reverse columns
        rc_matrix = np.fliplr(matrix)
        
        #Reverse rows
        return np.flipud(rc_matrix)


class PWM():
    """Class for dealing with nucleotide position weight matrices."""
    
    def __init__(self, matrix, background = np.array([0.25, 0.25, 0.25, 0.25]), suppress_pwm = False, name = None, reverse = False, **kwargs):
        
        if reverse is True:
            matrix = np.flipud(np.fliplr(matrix))
        
        self.pfm = np.array(matrix)
        self.rc_pfm = RCmatrix(self.pfm)
        
        if self.pfm.shape[0] != 4:
            raise Exception('Position weight matrices must have exactly 4 rows')
        
        self.background = np.array(background)
        if len(self.background) != 4:
            raise Exception('Length of background frequencies must be exactly 4.')
        
        self.ic = np.array([2+sum([self.pfm[j,i] * np.log2(self.pfm[j,i]) if self.pfm[j,i] > 0 else 0.0
                          for j in range(self.pfm.shape[0]) ]) 
                          for i in range(self.pfm.shape[1]) ])
        
        
       
        if not name is None:
            self.name = name
        else:
            self.name = 'Unspecified protein'
       
        if suppress_pwm is False:
            self.pwm = np.vstack( (1.0/self.background * self.pfm[:,i] \
                                   for i in range(self.pfm.shape[1]))).T 
            self.pwm = np.log2(self.pwm)
            
            self.rc_pwm = RCmatrix(self.pwm)


            self.max_score = sum(self.pwm.max(0)) # Sum of  maximum possible log odds score
            
        for additional_arg in kwargs:
            setattr(self, additional_arg, kwargs[additional_arg])
            
        
    def __len__(self):

        return self.pwm.shape[1]

        
    def toFile(self, filename):
        """Saves the matrix (as plain text) to a specified file.""" 
        
        with open(filename, 'w') as f_out:
            header = '> Name: %s Consensus: %s\n' %(self.name, self.asIUPAC())
            f_out.write(header)
            for i in range(4):
                WriteTSV(self.pfm[i,:], f_out)

    @classmethod
    def fromBEEML(cls, filename, **kwargs):
        """Reads an energy matrix as output by BEEML-PBM."""
        
        energy_matrix = []
        
        for line in FlatFile(filename, skip = 1):
            energy_matrix.append([-float(x) for x in line[1:]])
        
        energy_matrix = np.array(energy_matrix)
        
        pfm = np.array([[np.exp(energy_matrix[i,j])/sum(np.exp(energy_matrix[:,j]))
                         for j in range(energy_matrix.shape[1])]
                         for i in range(4)])
        
        return cls(pfm, energy_matrix = energy_matrix, **kwargs)
            
    @classmethod
    def fromSNW(cls, filename, **kwargs):
        """Reads the PFM from a standard S&W output file.
        It passes the additional kwargs for the top seed and the top e-score."""
            
        matrix = []
        start_reading, num_rows, first_row = False, 0, True
        
        for line in FlatFile(filename, split = True):
            if first_row is True:
                
                top_seed = line[1]
                top_escore = round(float(line[2]),3)
                first_row = False
            
            #print line
            if len(line) > 1 and start_reading is True:
                matrix.append([max(float(x), 1e-15) for x in line[1:]])
                num_rows += 1
            
            if line[0] == 'Probability matrix' and num_rows == 0:
                start_reading = True
                    
            if num_rows == 4:
                start_reading = False
        
        return cls(matrix, top_seed = top_seed, top_escore = top_escore, **kwargs)
    
    def Trim(self, threshold = 0.5):
        """Trims the matrix, eliminating positions that do not meet an IC
        threshold. Trimming stops when two consecutive positions above the
        threshold are encountered.  """
        left_bound=None
        right_bound=None
        
        
        IC_mask = [1 if self.ic[i] > threshold else 0 for i,_ in enumerate(self.ic)]
        N = len(self.ic)
        
    
        if not 1 in IC_mask:
            print 'The IC threshold given would trim all PWM positions.'
            logging.warning('The IC threshold given would trim all PWM positions.')
            return False
    
        for i in range(N-1):
            if IC_mask[i] == 1 and IC_mask[i+1] == 1:
                left_bound = i
                break
            
        for i in range(N-1, 0, -1):
            if IC_mask[i] == 1 and IC_mask[i-1] == 1:
                right_bound = i+1
                break
        
        if left_bound is None or right_bound is None:
            print "There are no consecutive positions that meet the IC threshold given, so no trimmed logo could be made."
            logging.warning('There are no consecutive positions that meet the IC threshold given, so no trimmed logo could be made.')
            return False
        
        self.pfm = self.pfm[:, left_bound:right_bound]
        self.pwm = self.pwm[:, left_bound:right_bound]
        self.ic = self.ic[left_bound:right_bound]
        
        self.rc_pfm = RCmatrix(self.pfm)
        self.rc_pwm = RCmatrix(self.pwm)
        
        return True
    
    def MakeLogo(self, filename, title = "", reverse_complement = False, num_seqs = 1000,  **kwargs):
        """Creates a webLogo png file corresponding to the PWM.
        Optionally, a title can be displayed at the top of the logo."""
        #from Bio import Motif
        import StringIO
        
        temp_file = GenerateTempFilename()
        #create new tempfile
        subprocess.call('touch %s'%temp_file, shell=True)
        #give it read-write permissions for all users
        subprocess.call('chmod 660 %s'%temp_file, shell=True)

        nucleotides = ['A', 'C', 'G', 'T']
        
        #For some reason, the BioPython logo maker does not work with StringIO
        #Therefore, use a temporary file, as inelegant as it may be.

        if reverse_complement:
            output_matrix = self.rc_pfm
        else:
            output_matrix = self.pfm
        
        if '.png' in filename.lower():
            logo_format = 'PNG'
        elif '.pdf' in filename.lower():
            logo_format = 'PDF'
        else:
            raise Exception('File selected for making logo does not contain an appropriate extension (png/pdf)')

        # Create a cumulative PWM

        cum_pfm = np.cumsum(output_matrix, axis = 0)

        # Create a sequence list that matches the PWM

        pfm_buffer = open(temp_file, 'w')

        for i in range(num_seqs):
            seq_list = []
            fraction = float(i)/num_seqs
            for j in range(len(self)):
                for k, nucl in enumerate(nucleotides):
                    if fraction < cum_pfm[k, j]:
                        seq_list.append(nucl)
                        break

            pfm_buffer.write(''.join(seq_list) + '\n')

        pfm_buffer.close()

        label_positions = ','.join([str(_) for _ in range(1,len(self)+1)])
        
        #subprocess.call(['/usr/bin/weblogo', '-f', '%s'%temp_file, '-o', '%s'%filename, '-F', '%s'%logo_format, '-P', '""', '-t', '"%s"'%title, '--resolution', '600', '--errorbars', 'NO', '-c', 'classic', '-s', 'large', '--annotate', '%s'%label_positions])
        
        subprocess.call('weblogo -f %s -o %s -F %s -P "" -t "%s" --resolution 600 --errorbars NO -c classic -s large --annotate %s' \
        %(temp_file, filename, logo_format, title, label_positions), shell = True)
        
        subprocess.call("rm %s"%(temp_file),shell=True)
