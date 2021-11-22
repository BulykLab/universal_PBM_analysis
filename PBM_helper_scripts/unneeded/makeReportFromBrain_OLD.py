import argparse
from argparse import ArgumentParser
from os.path import isdir
from os import listdir
from subprocess import call


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

proteins = {}

for i,line in enumerate(IN):
    pngName = line[line.rfind('/')+1:].strip()
    
    if args.type=='':
        if pngName[:2]=="RC":
            protein = pngName.split('.')[2]
        else:
            protein = pngName.split('.')[1]
    
    if args.type=='SW' or args.type=='BL':
        protein = pngName.split('_')[0]
        
    if protein in proteins.keys():
        tmp = proteins[protein]
        tmp.append(line.strip())
        proteins[protein]=tmp
    else:
        proteins[protein]=[line.strip()]
    
    if i%8==0:
        if i!=0:
            f_out.write('</table>\n')
            f_out.write('<hr>\n')
            if args.type == '':
                f_out.write('<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>\n')
            if args.type=='BL':
                f_out.write('<br><br><br><br><br><br><br><br><br><br><br><br><br>\n')
        f_out.write('<table>\n')
    if i%2==0:
        f_out.write('<tr>\n')
        f_out.write('<td><figure>\n<figcaption>%s</figcaption>\n<img src="%s" width=450>\n</figure></td>\n'%(pngName,line.strip()))
    else:
        f_out.write('<td><figure>\n<figcaption>%s</figcaption>\n<img src="%s" width=450>\n</figure></td>\n'%(pngName,line.strip()))
        f_out.write('</tr>\n')

f_out.write('</table>\n<hr>\n')
if args.type == '':
    f_out.write('<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>\n')
if args.type == 'SW':
    f_out.write('<br><br>\n')
if args.type == 'BL':
    f_out.write('<br><br><br><br><br><br><br><br><br><br><br><br><br>\n')

count=0
for key in proteins.keys():
    count+=1
    f_out.write('<hr>\n')
    cntLogo=0
    for logo in proteins.get(key):
        cntLogo+=1
        f_out.write('%s\n'%(logo[logo.rfind('/')+1:]))
        f_out.write('<img src="%s" width=1100 height=500><br>\n'%logo)
        if not ((count == len(proteins.keys())) and (cntLogo == len(proteins.get(key)))):
            f_out.write('<br><br><br><br><br><br><br><br><br><br><br><br><br>\n')
    

f_out.write('</html>')

f_out.close()
IN.close()

#print proteins

print 'Generating .pdf report'

call('wkhtmltopdf %s %s' %(output, output[:-4]+'pdf'), shell = True)

print 'OK'
