#!/usr/bin/env python
#
# This file is part of Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-2007
# All rights reserved
#
# author: Damien Farrell, 2009

import urllib2
import smtplib
import StringIO

class URLChecker:
    """Basic checks for website in url
       cgitb should be disabled in the cgi script, or it won't work
    """
    def __init__(self):
        return

    url = "http://enzyme.ucd.ie/cgi-bin/titration_db/main.cgi"

    def check(self):
        try:
            stream = urllib2.urlopen(self.url)
            info = stream.info()
            print info
            status = info.status
            print status
            if status is not "":
                stream = None
                error = 'Error status %s' % str(status)
                print error
            else:
                lines = stream.readlines()
                string = "".join(lines)
                print string
                #stream = StringIO.StringIO(string)
                #print stream
        except urllib2.HTTPError, e:
            print e.code
            self.sendMail(e.code)
        return

    def sendMail(self, errcode=None):
        """send a mail"""
        to = 'damien.farrell@ucd.ie'
        mail_user = 'damienfa'
        sender = 'Titration_DB admin'
        smtpserver = smtplib.SMTP("mail.ucd.ie")
        smtpserver.ehlo()
        smtpserver.starttls()
        smtpserver.ehlo
        header = 'To:' + to + '\n' + 'From: ' + sender + '\n' + 'Subject: titdb problem \n'
        print header
        msg = header + '\n There is a problem with the titration_db website.. \n\n'
        smtpserver.sendmail(mail_user, to, msg)
        print 'done!'
        smtpserver.close()
        return

C=URLChecker()
C.check()
#C.sendMail()
