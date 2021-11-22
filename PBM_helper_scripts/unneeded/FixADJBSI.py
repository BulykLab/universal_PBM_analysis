from sys import argv


input_file = argv[1]
output_file = argv[2]


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


#Count number of negative values


found_index = False
counter = 0 
for line in FlatFile(input_file, split = False):

    line = line.strip()
    split_line = line.split('\t')
    
    if len(split_line) < 5:
        pass
        
    else:
        if 'ADJBSI' in line and found_index is False:
            index = split_line.index('ADJBSI')
            found_index = True
            
        elif len(line) > 1:
            
            bsi_value = float(split_line[index])
            
            if bsi_value < 0:
                counter +=1
                
found_index = False               
                    
if counter > 200:
    print '--- There are over 200 negative values in the Masliner file. OMG. ---'

    with open(output_file, 'w') as f_out:
            
                
        for line in FlatFile(input_file, split = False):
        
            line = line.strip()
            split_line = line.split('\t')
            
            if len(split_line) < 5:
                f_out.write(line + '\n')
                
            else:
                if 'ADJBSI' in line and found_index is False:
                    index = split_line.index('ADJBSI')
                    found_index = True
                    f_out.write(line + '\n')
                
                elif len(line) > 1:
                    
                    bsi_value = float(split_line[index])
                    
                    if bsi_value < 0:
                        split_line[index] = "1"
                    
                    f_out.write('\t'.join(split_line) + '\n')
        

