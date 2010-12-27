#!/usr/bin/env python
#
# pKa - various programs and scripts for pKa value analysis, calculation and redesign
# Copyright (C) 2010 Jens Erik Nielsen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact information: 
# Email: Jens.Nielsen_at_gmail.com
# Normal mail:
# Jens Nielsen
# SBBS, Conway Institute
# University College Dublin
# Dublin 4, Ireland

import sys, os

resources='/home/webserver/pKaWebServer/resources/'

def send_email(address,message):
    """Send an email to the address saying that prepration for the dir is done"""
    import smtplib
    from email.MIMEText import MIMEText
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMENonMultipart import MIMENonMultipart
    from email.MIMEBase import MIMEBase
    from email import Encoders
    import mimetypes
    #
    # Prepare the message
    #
    text=message
    text_mes=MIMEText(text)
    message=MIMEMultipart()
    message.preamble='pKD message'
    message.epilogue = ''
    message['Subject']='Message from pKD server'
    message['From']='pKD Server <pka@ucd.ie>' 
    message['Reply-To']='pka@ucd.ie'
    message['To']=address
    message.attach(text_mes)
    #
    # Send the email
    #
    #sys.stdout=X
    try:
        server=smtplib.SMTP('193.1.169.34') #mail.ucd.ie
        server.set_debuglevel(1)
        server.login('jensniel','alpha2')
        server.sendmail('pka@ucd.ie',address,message.as_string())
        server.quit()
    except:
        pass
    return

#
# ----
#

def do_pKa_calc(pdbfile):
    """Perform a pKa calculation on the PDB file"""
    #
    # Get the name of the PDB file
    #
    import os
    pdb_name=os.path.split(pdbfile)[1]
    pdb_dir=os.path.split(pdbfile)[0]
    topdir=os.getcwd()
    #
    # Prepare the PDB file
    #
    import pKarun.WI_tools
    logfile,files=pKarun.WI_tools.dellig_delhyd_corall(pdbfile=pdbfile,readall=None,setcha=None)
    #
    # Overwrite the file
    #
    for line in logfile:
        print line,
    print
    print files.keys()
    if len(files.keys())==0:
        #
        # Write failed file
        #
        fd=open(os.path.join(pdb_dir,'failed'),'w')
        return None
    lines=files[files.keys()[0]]
    
    fd=open(pdbfile,'w')
    for line in lines:
        fd.write(line)
    fd.close()
    #
    # copy the PDB file and the parameter files
    #
    import os
    source=resources
    files=['DELRAD.DAT','DELCRG.DAT','TOPOLOGY.H']
    for file in files:
        r_s=os.path.join(source,file)
        r_d=os.path.join(pdb_dir,file)
        if not os.path.isfile(r_d):
            os.link(r_s,r_d)
    os.system('chmod -R ugo+rwx %s' %pdb_dir)
    #
    # Start pKa calc
    #
    import pKarun.pKarun_base
    params={'dbcrit':1000,'lowph':0.1,'highph':20.0,'phstep':0.1,'pairene':1.0}
    #
    # Save a record of how the pKa calc was done
    #
    fdx=open(os.path.join(pdb_dir,'Invocation_dict'),'w')
    import pickle
    pickle.dump(params,fdx)
    fdx.close()
    #
    # Run the pKa calculation
    #
    print 'Starting to run the pKa calculation'
    import sys
    sys.stdout.flush()
    Y=pKarun.pKarun_base.pKarun(os.getcwd(),pdbfile,params)
    Y.runseq()
    #
    # Check if all went well
    #
    import pKaTool.pKaIO
    X=pKaTool.pKaIO.pKaIO(pdbfile)
    calculation_success=1
    if not X.calculation_completed:
        print 'Calculation failed in ',pdbdir
        calculation_success=None
    #
    # Write a flag file to alert the scheduler
    #
    if calculation_success:
        print 'Writing Ready file'
        fd=open(os.path.join(pdb_dir,'pKacalc_done'),'w')
    else:
        #
        # Calculation failed
        #
        send_email('Jens.Nielsen@ucd.ie','pKa calculation failed for '+pdb_dir)
        print 'Writing failed file'
        fd=open(os.path.join(pdb_dir,'failed'),'w')
    #
    # Write the file
    #
    fd.write('pKa calculations\n')
    fd.close()
    #
    # Running file
    #
    runningfile=os.path.join(pdb_dir,'running')
    if os.path.isfile(runningfile):
        os.unlink(runningfile)
        print 'Deleted runningfile',runningfile
    #
    # done
    #
    print 'Everything is done'
    return 1

#
# -----
#

def run_sugelm(pdbfile):
    """Run the SugELM command"""
    #
    # Get the name of the PDB file
    #
    import os
    pdb_name=os.path.split(pdbfile)[1]
    pdb_dir=os.path.split(pdbfile)[0]
    topdir=os.getcwd()
    #
    # Check if pKa values were calculated for this file
    #
    import pKaTool.pKaIO
    X=pKaTool.pKaIO.pKaIO(pdbfile)
    if not X.calculation_completed:
        print 'pKa calculation not completed for %s' %pdbfile
        return None
    #
    # Read the pKa values
    #
    pkas=X.readpka()
    target=pkas.keys()[2]
    #
    # Create the SUGELM file
    #
    import Design_pKa
    tparams=Design_pKa.get_defaults()
    params={}
    for key in tparams.keys():
        params[key]=tparams[key][0]
    #
    # Set the other parameters
    #
    params['pHstart']=0.1
    params['pHstop']=12.0
    params['pHstep']=0.05
    params['pKMCsteps']=200000,
    params['pKas']=target+'=+1.0'
    params['MC']=0
    params['TR']=0
    params['min_target_dist']=1.0
    params['save_solutions']=None
    params['silent']=0
    params['pdb']=pdbfile
    params['generate_mutations']=True
    #
    # Log message
    #
    import sys
    print 'Starting pKD mutation preparation'
    sys.stdout.flush()
    try:
        X=Design_pKa.Design_pKa(params)
        X.get_interaction_energies()
        print 'Calculated interaction energies'
    except:
        send_email('Jens.Nielsen@ucd.ie','SUGELM failed for '+pdbfile)
        return None
    #
    # Write a flag file to alert the scheduler
    #
    fd=open(os.path.join(pdb_dir,'ready'),'w')
    fd.write('SUGELM is done\n')
    fd.close()
    return 1

#
# -----
#
    
if __name__=='__main__':
    import sys
    do_pKa_calc(sys.argv[1])
    run_sugelm(sys.argv[1])


