import argparse
from argparse import ArgumentParser
from os.path import isdir, isfile
from os import listdir
from subprocess import call

analysis_tools_path = '/n/data2/bch/medicine/bulyk/shared_software/universal_PBM_analysis/PBM_helper_scripts/'

arg_obj = ArgumentParser()
arg_obj.add_argument('pngs', metavar = 'pngFile', help='file with full paths to png files to use in report (output from BrainReport3)')
arg_obj.add_argument('notes', metavar = 'notes', help='notes about experiment, like v number')
arg_obj.add_argument('person', metavar = 'person', help='person who ran the experiment')
arg_obj.add_argument('type', metavar = 'type', help='FR,SW,BL')
args = arg_obj.parse_args()

if args.type == 'FR':
    args.type=''

directory = args.pngs[:args.pngs.rfind('/')]+'/'
output = directory+'chosenLogoReport%s.html'%args.type
print output

IN = open(args.pngs,'r')
f_out = open(output,'w')
f_out.write('<html>\n')
f_out.write('<center><font size = "6">\n')
f_out.write('Array: %s<br>\n'% args.notes)
f_out.write('Done by: %s<br>'% args.person)
if args.type=='FR':
    f_out.write('<br>')
f_out.write('</center></font>\n')

proteins = [] #protein names
primary = {} #PNG files for primary logo analysis
secondary = {} #PNG files for secondary logo analysis
pFN = {} #8mer_11111111.txt files for primary 8mer>0.45
sFN = {} #8mer_11111111.txt files for secondary 8mer>0.45

for i,line in enumerate(IN):
    pngName = line[line.rfind('/')+1:].strip()
    oriDir = line[:line.rfind('/')+1]
    
    if args.type=='':
        if pngName[:2]=="RC":
            protein = pngName.split('.')[2]
            motifType = pngName.split('.')[1]
        else:
            protein = pngName.split('.')[1]
            motifType = pngName.split('.')[0]
    
    if args.type=='SW' or args.type=='BL':
        protein = pngName.split(".")[0].split('_')[0]
        motifType = pngName.split(".")[0].split('_')[1]
    
    
    if protein not in proteins:
        proteins.append(protein)
    if motifType == "primary" or motifType == "BEEML":
        #primary.append(line.strip())
        primary[protein]=line.strip() 
        pFN[protein] = oriDir+motifType+"_"+protein+"_8mers_11111111.txt"
    if motifType == "secondary":
        #secondary.append(line.strip())
        secondary[protein]=line.strip()
        sFN[protein] = oriDir+motifType+"_"+protein+"_8mers_11111111.txt"
    
    
IN.close()

def getOneShellOutput(cmd):
    from subprocess import call, PIPE, Popen
    job_process = Popen(cmd, shell=True, stdout = PIPE, stderr = PIPE)
    out_stream, err_stream = job_process.communicate()
    if len(err_stream)!=0:
        print err_stream
    return out_stream.strip()

def Count8mersAboveThreshold(filename, threshold = 0.45):
    #Takes a S&W filename and counts the number of contiguous 8-mers that exceed a certain threshold.
    if not isfile(filename):
        return 0
    else:
        return getOneShellOutput("sed '1d' %s | awk '$3>%s' | wc -l"%(filename,threshold))

def MakeTable(f_out, proteins, pngNames, eightMers, height=''):
    
    f_out.write('<table>')
    
    names_to_use = proteins[:]

    if len(names_to_use) % 2 != 0:
        names_to_use.append('__blank__')

    for i in range(0, len(names_to_use), 2):
        pngName1 = pngNames[names_to_use[i]]
        if names_to_use[i+1] in pngNames:
            pngName2 = pngNames[names_to_use[i+1]]
        else:
            pngName2="NA"
        
        
        if not isfile(pngName1):
            pngName1 = "%splaceholder.png"%(analysis_tools_path)
        if not isfile(pngName2):
            pngName2 = "%splaceholder.png"%(analysis_tools_path)
        
        f_out.write('<tr>\n')
        
        caption1 = pngName1[pngName1.rfind('/')+1:].strip()
        caption2 = pngName2[pngName2.rfind('/')+1:].strip()
        
        if args.type=='SW':
            caption1 += " - 8mers > 0.45: "+str(eightMers[names_to_use[i]])
            if names_to_use[i+1] in eightMers.keys():
                caption2 += " - 8mers > 0.45: "+str(eightMers[names_to_use[i+1]])
        
        f_out.write('<td>'+'<figure><figcaption>%s</figcaption>'%caption1+'<img src="%s" ' %pngName1 +'width=450 %s>' %height+'</figure></td>\n')
        f_out.write('<td>'+'<figure><figcaption>%s</figcaption>'%caption2+'<img src="%s" ' %pngName2 +'width=450 %s></figure>' %height+ '</td>\n')
        
        f_out.write('</tr>\n')
    f_out.write('</table>')
    f_out.write('<hr>\n')


if args.type=='SW':
    num_above_thresh_primary = dict((protein,Count8mersAboveThreshold(pFN[protein], 0.45))
                        for protein in proteins)
    
    num_above_thresh_secondary = dict((protein,Count8mersAboveThreshold(sFN[protein], 0.45))
                        for protein in proteins)
else:
    num_above_thresh_primary={}
    num_above_thresh_secondary={}
        
insert_image = '<figure><img src="%s" width=1100 height=275><figcaption>%s</figcaption></figure>'


with open(output, 'w') as f_out:
    
    f_out.write('<html>\n')
    f_out.write('<center><font size = 13>')
    f_out.write('Array: %s<br>'% args.notes)
    f_out.write('Done by: %s<br>'% args.person)
    f_out.write('</center></font>')
    
    if len(primary)>0:
        MakeTable(f_out, proteins, primary, num_above_thresh_primary, 'height=110')
    if len(secondary)>0:
        MakeTable(f_out, proteins, secondary, num_above_thresh_secondary, 'height=110')
    f_out.write('<br><br><br><br><br>') #for Cerebellum, only 3 <br>
    
    for protein in proteins:
        
        if len(primary)>0:
            f_out.write('\n<hr>')
            PNGname = primary[protein]
            caption = PNGname[PNGname.rfind('/')+1:].strip()
            if args.type=='SW':
                caption+=' - 8-mers > 0.45: <b>%s</b>'%(num_above_thresh_primary[protein])
            f_out.write("<br>\n"+insert_image %(PNGname,caption)+"\n<br>")
        if len(secondary)>0:
            f_out.write('\n<hr>')
            PNGname = secondary[protein]
            caption = PNGname[PNGname.rfind('/')+1:].strip()
            if args.type=='SW':
                caption+=' - 8-mers > 0.45: <b>%s</b>'%(num_above_thresh_secondary[protein])
            f_out.write("<br>\n"+insert_image %(PNGname,caption)+"\n<br><hr><br>") #for Cerebellum, take out the last <br>
    
    f_out.write('\n<hr>\n</html>')

print 'Generating .pdf report'

call('wkhtmltopdf %s %s' %(output, output[:-4]+'pdf'), shell = True)

print 'OK'
