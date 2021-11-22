import email
import smtplib
import argparse
from argparse import ArgumentParser
from email import MIMEMultipart
import mimetypes
import email.mime.application

arg_obj = ArgumentParser()
arg_obj.add_argument('-to', metavar = 'addressee', help='addressee')
arg_obj.add_argument('-pw', help = 'password', metavar = 'password')
arg_obj.add_argument('-attach', help = 'filename of attachment')
args = arg_obj.parse_args()   

sender = 'pbmanalysis@gmail.com'

msg = email.mime.Multipart.MIMEMultipart()
msg['Subject'] = 'PBM Analysis Pipeline Completed: %s'%args.attach
msg['From'] = sender
msg['To'] = args.to

filename = args.attach
fp=open(filename, 'rb')
att = email.mime.application.MIMEApplication(fp.read(), _subtype="pdf")
fp.close()
att.add_header('Content-Disposition','attachment', filename=filename)
msg.attach(att)

s=smtplib.SMTP('smtp.gmail.com',587)
s.starttls()
s.login(sender,args.pw)
s.sendmail(sender,args.to, msg.as_string())
s.quit()

print 'OK'
