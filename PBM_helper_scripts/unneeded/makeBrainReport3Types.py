import argparse
from argparse import ArgumentParser
from os.path import isdir
from os import listdir
from subprocess import call


arg_obj = ArgumentParser()
arg_obj.add_argument('directory', metavar = 'dir', help='FULL PATH directory with png files - IT IS EXTREMELY IMPORTANT TO GIVE THE FULL PATH. START WITH /net/... or ~/...')
arg_obj.add_argument('annotations', metavar = 'notes', help='array notes')
arg_obj.add_argument('person', metavar = 'person', help='person who did the array')
arg_obj.add_argument('-names', metavar = 'proteins', help='if want in specific order, comma separted list of protein names')
arg_obj.add_argument('type', metavar = 'type', help='FeatureREDUCE => FR, Seed and Wobble => SW, BEEML => BL')
args = arg_obj.parse_args()

if args.type != 'FR' and args.type != 'SW' and args.type != 'BL':
    raise Exception('The type can only be FR, SW, or BL. Please try again.')


def MakeTable(f_out, protein_names, motif_type, height=''):
    
    f_out.write('<table>')
    
    names_to_use = protein_names[:]

    if len(names_to_use) % 2 != 0:
        names_to_use.append('__blank__')

    for i in range(0, len(names_to_use)):
        if args.type == "SW" or args.type == 'BL':
            pngName1 = names_to_use[i] + '_' + motif_type + '.png' #type change
            pngName2 = names_to_use[i] + '_' + motif_type + '_rc.png' #type change
        if args.type == 'FR':
            pngName1 = motif_type + '.' + names_to_use[i] + '.final.psam.10nt.png'
            pngName2 = 'RC.' + motif_type + '.' + names_to_use[i] + '.final.psam.10nt.png'
            
        
        varName = 'jkl-'+motif_type[0]+'-'+names_to_use[i]
        
        f_out.write('<tr><td colspan="2"><b>Block:</b> %s <b>Logo Type:</b> %s</td></tr>\n'%(names_to_use[i],motif_type))
        f_out.write('<tr><td colspan="2">\n')
        f_out.write('<input type="radio" name="%s" value="%s" <?php if (isset($_POST["%s"]) && $_POST["%s"]=="forward") { ?> checked="checked" <?php } ?> > forward\n'%(varName,args.directory+pngName1,varName,varName))
        f_out.write('<input type="radio" name="%s" value="%s" <?php if (isset($_POST["%s"]) && $_POST["%s"]=="RC") { ?> checked="checked" <?php } ?> > RC\n'%(varName,args.directory+pngName2,varName,varName))
        
        f_out.write('</td></tr>\n')
        
        f_out.write('<tr>\n')
        
        f_out.write('<td><figure>\n<figcaption>%s</figcaption>\n<img src="%s" width=450 %s>\n</figure></td>\n'%(pngName1,args.directory.replace('/','_-')+'/'+pngName1,height))
        f_out.write('<td><figure>\n<figcaption>%s</figcaption>\n<img src="%s" width=450 %s>\n</figure></td>\n'%(pngName2,args.directory.replace('/','_-')+'/'+pngName2,height))
        
        f_out.write('</tr>\n')
    f_out.write('</table>')
    f_out.write('<hr>\n')
    

def getNames(direct,PN):
    if isdir(direct):
        for f in listdir(direct):
            if args.type=='FR':
                if "final.psam.10nt.png" in f and "RC.primary" in f:
                    tabs = f.split('.')
                    PN.append(tabs[2])
            if args.type=='SW':
                if "primary_rc.png" in f: #type change
                    tabs = f.split('_') #type change
                    PN.append(tabs[0]) #type change
            if args.type=='BL':
                if "BEEML_rc.png" in f: #type change
                    tabs = f.split('_') #type change
                    PN.append(tabs[0]) #type change

if args.names is not None:
    protein_names = [x for x in args.names.split(',')]
else:
    protein_names = []
    getNames(args.directory,protein_names)
print protein_names

if len(protein_names) == 0:
    raise Exception('There are no PNG files in the directory you gave for the type you gave. Please try again.')


if args.type=='FR':
    suffix = ''
else:
    suffix = args.type

insert_image = '<img src="%s%s.%s.final.psam.10nt.png" width=1100 height=500><br>'  
with open(args.directory + 'brainReport3%s.php'%suffix, 'w') as f_out: #type change
    
    f_out.write('<html>\n')
    f_out.write('<center><font size = 13>')
    f_out.write('Array: %s<br>'% args.annotations)
    f_out.write('Done by: %s<br>'% args.person)
    f_out.write('</center></font>\n')
    f_out.write('<form method="post" action="brainReport3Submit.php">\n')
    f_out.write('<input type="hidden" name="FN" value="<?php echo isset($_POST["FN"]) ? $_POST["FN"] : "%s" ?>">\n'%(args.directory.replace('/','_-')+'/'+'chosenLogos%s.txt'%suffix)) #type change
    f_out.write('<input type="hidden" name="htmlFN" value="<?php echo isset($_POST["htmlFN"]) ? $_POST["htmlFN"] : "%s" ?>">\n'%(args.directory.replace('/','_-')+'/'+'chosenLogoReport%s.html'%suffix)) #type change
    f_out.write('<input type="hidden" name="pdfFN" value="<?php echo isset($_POST["pdfFN"]) ? $_POST["pdfFN"] : "%s" ?>">\n'%(args.directory.replace('/','_-')+'/'+'chosenLogoReport%s.pdf'%suffix)) #type change
    f_out.write('<input type="hidden" name="directory" value="<?php echo isset($_POST["directory"]) ? $_POST["directory"] : "%s" ?>">\n'%(args.directory.replace('/','_-')+'/'))
    f_out.write('<input type="hidden" name="notes" value="<?php echo isset($_POST["notes"]) ? $_POST["notes"] : "%s" ?>">\n'%(args.annotations))
    f_out.write('<input type="hidden" name="person" value="<?php echo isset($_POST["person"]) ? $_POST["person"] : "%s" ?>">\n'%(args.person))
    f_out.write('<input type="hidden" name="type" value="<?php echo isset($_POST["type"]) ? $_POST["type"] : "%s" ?>">\n'%args.type) #type change
    
    
    
    
    if args.type == 'BL':
        MakeTable(f_out, protein_names, 'BEEML')#,height=600)
        f_out.write('<br><hr><br>')
        #the SW pipe says BEEML_2nd is in testing and shouldn't be used. So not including it
        #MakeTable(f_out, protein_names, 'BEEML_2nd')#,height=600)
        #f_out.write('<br><hr><br>')
    else:
        MakeTable(f_out, protein_names, 'primary')#,height=600)
        f_out.write('<br><hr><br>')
        MakeTable(f_out, protein_names, 'secondary')#,height=600)
        f_out.write('<br><hr><br>')
    
    
    f_out.write('email: <input type="text" name="email" value="<?php echo isset($_POST["email"]) ? $_POST["email"] : "" ?>"><br><br>\n') #new
    f_out.write('password to email server: <input type="text" name="pw" value="<?php echo isset($_POST["pw"]) ? $_POST["pw"] : "" ?>"><br><br>\n') #new
    f_out.write('<input type="submit" name="submit" value="Submit">\n')
    f_out.write('</form>\n')
    
    
    f_out.write('\n</html>')


