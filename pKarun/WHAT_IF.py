#
# $Id: WHAT_IF.py 3060 2008-06-16 10:30:51Z nielsen $
#
# Run WHAT IF from a python program.
#
# (C) Rob W.W. Hooft, 1996
#
# Modified extensively
#
# Jens Erik Nielsen 1999
#
#---------------------------------------------------------------------
# WHAT_IF.py
#
# The WHAT_IF class allows WHAT IF to be run from PYTHON.
#
# The class needs a file in your homedir called ~/.WHAT_IF. This
# file should contain two lines. The first line specifies the path
# to the WHAT IF executable. The second line is either "None" or
# specifies the path to the WHATIF.FIG to be used. If you have no
# special preferences set in your WHATIF.FIG you should enter None
# in the second line.
# This class is also capable of creating the ~/.WHAT_IF file automatically.
#
# Functions are:
# 
#    def RunIfNew(self,command,logfile):
#    def Run(self,command,logfile,cleanup=1):
#    def RunNoClean(self,command,logfile):
#    def RedirectOutput(self,outfilnam):
# and should be self explanatory.
#
# address any questions to: Jens.Nielsen@UCD.IE
#
#Example usage:
#import WHAT_IF
# outlog='Whatif.out.log'
# logfile='Whatif.log'
# whatif=WHAT_IF.Whatif()
# whatif.RedirectOutput(outlog)
# commandline='getmol 1crn Y xxx \n grafic \n shoall 1 q \n center \n go \n '
#Call WHAT IF
# try:
#       whatif.Run(commandline,logfile,1)
#except 'WHAT_IFError',arg:
#	print 'WHAT IF failed for this molecule!!\n'
#
# ====================================================================
# WHAT IF notes:
#
# The WHAT IF 'go' can be executed via this module. Execution halts
# until 'CHAT' is clicked in the graphics window.
# Setting WIFPAR 710 to 1 ('setwif 710 1') will make WHAT IF continue
# execution even if 'go' is executed.
#
# The WHAT IF Website: http://www.cmbi.kun.nl/whatif
#

import os,string,stat,time
import re
Rcompress=re.compile('\\.gz$')
TrashFiles={'WHATIF.FIG':1,'STARTUP.FIL':1,'TEXTABLE.DAT':1,'ALTERR.LOG':1,'TEXSTORE.DAT':1}
Rtrash=re.compile('^\\('+string.joinfields(['FOR0[0-9][0-9]\\.DAT',
                                               'TAPE\\(IN\\|OUT\\)\\.DAT',
                                               'STARTUP\\.FIL',
                                               'pdbout\\.\\(tex\\|txt\\)',
                                               'check\\.db',
                                               'WHATIF\\.FIG',
                                               'TEX\\(STORE\\|TABLE\\)\\.DAT',
                                               'PICK\\.IDX',
                                               'SCATTER\\.\\(DAT\\|SCC\\|fig\\)',
                                               '\\(eps\\|sct\\)0...\\.eps'],'\\|')
                     +'\\)$')
              
class Whatif:
    def __init__(self,wiffile=None):
         """exceptions"""
         self.error=Exception
         """If the file ~/.WHAT_IF exists then we take the WHAT IF location from there"""
         if os.environ.has_key('WIF'):
             wiffile=os.environ['WIF']
         if not wiffile:
             homedir=os.environ['HOME']
             wiffile=homedir+'/.WHAT_IF'
         if os.path.isfile(wiffile):
              fd=open(wiffile)
              self.progrm=fd.readline()[:-1]
              self.wiffig=fd.readline()[:-1]
              fd.close()
         else:
              """Locate WHAT IF"""
              print
              print 'This is the first time you are using the WHAT_IF class.'
              print 'This class needs the location of the WHAT IF exeutable',
              print 'and the location of WHATIF.FIG',
              print 'specified in the file '+wiffile
              print
              answer=raw_input ('Should a search for WHAT IF executables be initiated? (y/n) ')
              if string.lower(answer)=='y':
                  path=raw_input("Please enter root of search path (default=/) ")
                  if path[-1]!='/':
                      path=path+'/'
                  answer=raw_input('Do you want to locate WHATIF.FIG too (y/n) (default=n) ')
                  if string.lower(answer)!='y':
                      self.restrict=0
                  else:
                      self.restrict=1
                  print 'Please wait.......'
                  X=WIfind()
                  locations=X.find('whatif',path)
                  if len(locations)==0:
                      raise self.error,'No executable found'
                  execlist=[]
                  for file in locations:
                      if os.path.isfile(file):
                          wifsrcdir,junk=os.path.split(file)
                          wifdir,junk=os.path.split(wifsrcdir)
                          dbdata=os.path.join(wifdir,"dbdata")
                          cconfi=os.path.join(dbdata,"CCONFI.FIG")
                          wiffig=os.path.join(wifdir,"WHATIF.FIG")
                          if self.restrict==1:
                              if (os.path.isfile(cconfi) and os.path.isfile(wiffig)):
                                  date=time.strftime('%a %d %b %Y %H:%M',time.localtime(os.stat(file)[stat.ST_MTIME]))
                                  execlist.append([file,date])
                          else:
                              date=time.strftime('%a %d %b %Y %H:%M',time.localtime(os.stat(file)[stat.ST_MTIME]))
                              execlist.append([file,date])
                  """Ask the user to identify the version"""
                  if len(execlist)==0:
                      raise self.error,'No executable found'
                  elif len(execlist)==1:
                      self.progrm=execlist[0][0]
                      
                  else:
                      x=0
                      print 'The following WHAT IF versions were found.'
                      for prog in execlist:
                          x=x+1
                          print '%2i %30s last modified %20s' %(x,str(prog[0]),str(prog[1]))
                      print
                      answer=raw_input ('Which version do you want to use? ')
                      if string.atoi(answer)<1 or string.atoi(answer)>len(execlist):
                          raise self.error,'Invalid number entered'
                      self.progrm=execlist[string.atoi(answer)-1][0]
                  """Write the path to the ~/.WHAT_IF file"""
                  print 'Creating ',wiffile
                  if self.restrict==1:
                      wifsrcdir,junk=os.path.split(self.progrm)
                      wifdir,junk=os.path.split(wifsrcdir)
                      dbdata=os.path.join(wifdir,"dbdata")
                      self.wiffig=os.path.join(wifdir,"WHATIF.FIG")
                  else:
                      self.wiffig='None'
                  fd=open(wiffile,'write')
                  fd.write(self.progrm+'\n')
                  fd.write(self.wiffig+'\n')
                  fd.close()
              else:
                  print
                  print 'Please enter the location of the WHAT IF executable in ',wiffile
                  print
                  return
         if not os.path.isfile(self.progrm):
             print self.progrm
             raise self.error,'WHAT IF executable not found'
         #print 'Using executable: ',self.progrm
         self.outfilnam=None
         return 

    def _prepare(self):
	if os.path.exists("STARTUP.FIL"):
	    os.unlink("STARTUP.FIL")
        if os.path.exists("WHATIF.FIG"):
            os.unlink("WHATIF.FIG")
        if self.wiffig!='None':
            os.symlink(self.wiffig,'WHATIF.FIG')
            print 'Using: ',self.wiffig
            
    def _cleanup(self):
	for file in os.listdir("."):
	    if Rtrash.match(file)>=0 or TrashFiles.has_key(file):
		os.unlink(file)

    def _dorun(self,command,logfile):
	commandfile=open("STARTUP.FIL","w")
	commandfile.writelines(
		["setwif 622 0\n",
		 "SETWIF 1194 1\n",
		 "dolog\n",
		 logfile+"\n",
		 "WHAT IF run from Python\n",
		 "\n",
                 'setwif 1816 30000\n',
                 'setwif 792 1\n',
                 'setwif 105 0\n',
                 'setwif 1194 1\n',
                 'setwif 1012 0\n',
		 command+"\n",
		 "nolog\n",
		 "fullst y\n"
		 ])
	commandfile.close()
	try:
	    if self.outfilnam:
		status=os.system(self.progrm+" > "+self.outfilnam+" 2>&1 </dev/null")
	    else:
                print self.progrm+" 2>&1 </dev/null"
		status=os.system(self.progrm+" 2>&1 </dev/null")
	except KeyboardInterrupt:
	    status=2
	if status==2:
	    if os.path.exists(logfile):
                print 'unlinking: ',logfile
		os.unlink(logfile)
	    self._cleanup()
	    raise KeyboardInterrupt,"during whatif run"
	elif status!=0:
	    print "Error code ",status
	    raise self.error,"nonzero exit code from whatif"
        return
    #
    # User interface
    #
    def RunIfNew(self,command,logfile):
	if os.path.exists(logfile):
	    return 1
	else:
	    self.Run(command,logfile)
	    return 0

    def Run(self,command,logfile,cleanup=1):
        """Run the WHAT IF command"""
        #
	if os.path.exists(logfile) and logfile[0:5]!="/dev/":
	    os.unlink(logfile)
	compress=0
	if Rcompress.search(logfile)>=0:
            compress=1
	    ilogfile=re.sub(Rcompress,'',logfile)
	else:
            ilogfile=logfile
	self._prepare()
	self._dorun(command,ilogfile);
	if cleanup:
	    self._cleanup();
	if compress:
	    status=os.system("gzip "+ilogfile)
	    if status==2:
		if os.path.exists(logfile):
		    os.unlink(logfile)
		raise KeyboardInterrupt,'during gzip call'
	    elif status>0:
		print "Error code ",status
		raise self.error,"nonzero exit code from gzip"
	return 0

    def RunNoClean(self,command,logfile):
	self.Run(command,logfile,0)

    def RedirectOutput(self,outfilnam):
	self.outfilnam=outfilnam

    def getHTOP(self,destination=os.getcwd()):
        """Symlink TOPOLOGY.H to the pwd"""
        wifsrcdir,junk=os.path.split(self.progrm)
        wifdir,junk=os.path.split(wifsrcdir)
        dbdata=os.path.join(wifdir,"dbdata")
        if os.path.isfile(os.path.join(dbdata,'TOPOLOGY.H')):
            if os.path.isfile(os.path.join(destination,'TOPOLOGY.H')):
                raise error,"TOPOLOGY.H already present"
            else:
                 os.symlink(os.path.join(dbdata,'TOPOLOGY.H'),os.path.join(destination,'TOPOLOGY.H'))
        else:
            raise error,"TOPOLOGY.H not found"
        return

#
# --------------------------------------
#

class WIfind:
     def __init__(self):
          self._debug = 0
          self._prune = ['(*)']
          return

     def find(self,pattern, dir = './'):
          import fnmatch
          import os
          list = []
          try:
               names = os.listdir(dir)
          except os.error,arg:
               return []
          names.sort()
          for name in names:
               if name in (os.curdir, os.pardir):
                    continue
               fullname = os.path.join(dir, name)
               if fnmatch.fnmatch(name, pattern):
                    list.append(fullname)
               if os.path.isdir(fullname) and not os.path.islink(fullname):
                    for p in self._prune:
                         if fnmatch.fnmatch(name, p):
                              if self._debug: print "skip", `fullname`
                              break
                         else:
                              if self._debug: print "descend into", `fullname`
                              list = list + self.find(pattern, fullname)
          return list
