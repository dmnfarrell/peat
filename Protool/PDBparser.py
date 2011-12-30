#--------------------------------------------------------------------------------------------
#--------------- This class can (sometimes) parse PDB files ---------------------------------
#--------------------------------------------------------------------------------------------
import string
from errors import *

class PDBparser:
    
    def __init__(self,lines,parse_terms):
        #
        # Initialise atom arrays
        #
        import string
        self.lines=[]
        for line in lines:
            self.lines.append(string.strip(line))
        self.sequence=[]
        self.residues={}
        self.atoms={}
        self.attribute={}
        #
        # The remaining arrays are for the header information
        #
        self.header=''
        self.compnd=''
        self.source=''
        self.author=''
        self.revdat=''
        self.remark=''
        self.seqres=[]
        self.helix=[]
        self.sheet=[]
        self.turn=[]
        self.ssbond=[]
        self.cryst=[]
        self.orig=[]
        self.scale=[]
        self.NMRmodel=False
        self.ignore_NMRmodel=False
        self.readmodel=True
        self.readmodels=1
        self.models_already_read=0
        #
        # Parse TER records and assign new chain IDs?
        #
        self.parse_terms=parse_terms
        #
        # Initialise some help variables
        #
        self.linenumber=0
        self.atomlinedefs={}
        self.chains={}
        self.useids=['J','K','L','M','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        for x in string.letters:
            self.useids.append(string.lower(x))
        for x in string.digits:
            self.useids.append(string.upper(x))
        self.oldchainid=None
        self.freeresnum=9000
        #
        # Defining the junk headers
        #
        self.junkheaders={'TITLE':'REMARK','KEYWDS':'REMARK','DBREF':'REMARK','SEQADV':'REMARK','HET':'REMARK','HETNAM':'REMARK','CISPEP':'REMARK','LINK':'REMARK','MTRIX1':'REMARK','MTRIX2':'REMARK','MTRIX3':'REMARK','SITE':'REMARK','HETSYN':'REMARK','ANISOU':'JUNK','SPRSDE':'REMARK','MODRES':'REMARK','HYDBND':'REMARK','TVECT':'REMARK','SLTBRG':'REMARK','SIGUIJ':'REMARK','SIGATM':'REMARK','CAVEAT':'REMARK','JNRL':'JRNL'}
        return

    def parse(self):
        #
        # Loop through all lines and extract all readable information
        # (for now only ATOM lines are parsed)
        #
        for line in self.lines:
            self.parseline(line)
        #
        # Check for name clashes and atoms with identical positons.
        #
        self.atoms=self.nameclash()
        #
        # Renumber
        #
        self.renumber()
        return

    def nameclash(self):
        """
        # This routine checks for atom/residue name clashes in self.atomlinedefs
        # and gives later entries the key 'ALT'
        # So 'A:231:CD' would become 'A:231:CD,ALT' if there was a 'A:231:CD' earlier
        # in the pdb file. Identical atom names in different molecules (separated by a 'TER')
        # are not affected since Protool does not allow identical chain identifiers.
        # If there is already an '*,ALT' then the name becomes '*,ALT:ALT' and so forth...
        """
        clashdict={}
        lines=self.atomlinedefs.keys()
        lines.sort()
        for line in lines:
            atom=self.atomlinedefs[line]
            uniqueid=atom['CHAINID']+':'+atom['RESNUM']+':'+str(atom['ATNAME'])
            clashdict,newuniqueid=self.enternewname(clashdict,uniqueid,atom)
        return clashdict

    def enternewname(self,clashdict,uniqueid,atom):
        #
        # inserts uniqueid:atom in clashdict. If clashdict already contains
        # uniqueid, then uniqueid gets ',ALT' appended
        #
        import string
        if clashdict.has_key(uniqueid):
            #
            # If the residue names are identical, or if the last three letters are identical
            # [like in AGLN], then the add the ALT tag
            # If the resnames are different, then we assume that it is an insertion, and add
            # an 'A' or a 'B' to the residue number
            # This should be written more robustly so that insertations that are identical in
            # type to the previous residue are handled correctly...
            #
            if atom['RESNAME']==clashdict[uniqueid]['RESNAME'] or (len(atom['RESNAME'])==4 and atom['RESNAME'][1:]==clashdict[uniqueid]['RESNAME'][1:]):
                if uniqueid[-4:]==',ALT':
                    uniqueid=uniqueid+':ALT'
                else:
                    uniqueid=uniqueid+',ALT'
                clashdict,uniqueid=self.enternewname(clashdict,uniqueid,atom)
                return clashdict,uniqueid
            else:
                resnum=string.split(uniqueid,':')[1]
                if resnum==resnum[:-1]+'A':
                    resnum=resnum[:-1]+'B'
                    uniqueid=string.split(uniqueid,':')[0]+':'+resnum+':'+string.join(string.split(uniqueid,':')[2:],':')
                else:
                    resnum=resnum+'A'
                    uniqueid=string.split(uniqueid,':')[0]+':'+resnum+':'+string.join(string.split(uniqueid,':')[2:],':')
                clashdict,uniqueid=self.enternewname(clashdict,uniqueid,atom)
                return clashdict,uniqueid
        else:
            clashdict[uniqueid]=atom
            return clashdict,uniqueid
    

    def parseline(self,line):
        """
        # This function calls other routines to parse the line
        """
        import string
        type=string.strip(string.upper(line[:6]))
        if len(type)==0:
            return
        #
        # There are so many junk headers, so we simplify a bit
        #
        if self.junkheaders.has_key(type):
            type=self.junkheaders[type]
        #
        # Should we ignore this model?
        #
        
        if self.readmodel is False and type!='ENDMDL':
            #print 'Skipping model no: %d' %self.NMRmodel
            return
        #
        # Now call the appropritate function
        #
        if type=='ATOM':
            self.parseatom(line)
        elif type=='HETATM':
            self.parsehet(line)
        elif type=='HEADER':
            self.parseheader(line)
        elif type=='COMPND':
            self.parsecompnd(line)
        elif type=='SOURCE':
            self.parsesource(line)
        elif type=='AUTHOR':
            self.parseauthor(line)
        elif type=='REVDAT':
            self.parserevdat(line)
        elif type=='REMARK':
            self.parseremark(line)
        elif type=='SEQRES':
            self.parseseqres(line)
        elif type=='HELIX':
            self.parsehelix(line)
        elif type=='SHEET':
            self.parsesheet(line)
        elif type=='TURN':
            self.parseturn(line)
        elif type=='MDLTYP':
            pass
        elif type=='SSBOND':
            self.parsessbond(line)
        elif type=='CRYST1' or type=='CRYSTL' or type=='CRYST':
            self.parsecryst(line)
        elif type=='ORIGX1' or type=='ORIGX2' or type=='ORIGX3':
            self.parseorig(line)
        elif type=='SCALE1' or type=='SCALE2' or type=='SCALE3':
            self.parsescale(line)
        elif type=='CONECT':
            self.parseconnect(line)
        elif type=='MASTER':
            self.parsemaster(line)
        elif type=='TER':
            self.parseter(line)
        elif type=='JRNL':
            self.parsejournal(line)
        elif type=='FTNOTE':
            self.parsefootnote(line)
        elif type=='HET':
            self.parsehetheader(line)
        elif type=='FORMUL':
            self.parseformul(line)
        elif type=='MODEL':
            self.parsemodel(line)
        elif type=='ENDMDL':
            self.parseendmodel(line)
        elif type=='END':
            self.parseend(line)
        elif type=='EXPDTA':
            self.parseexpdta(line)
        elif type=='JUNK':
            pass
        elif type=='NUMMDL':
            pass
        else:
            #print 'Ignored header: %s' %type
            #raise ParseLineError, 'Unknown Header type: %s' %type
            pass
        return

    #
    # ----
    #
        
    def parseatom(self,line,hetatm=None):
        """
        #
        # Parse a single ATOM line
        #"""
        import typecheck, string, structure
        G=structure.structure()
        #
        # This routine parses every atom line and stores it in self.atomlinedefs
        #
        self.linenumber=self.linenumber+1
        #
        # Check the atomnumber
        #
        atomnumber=line[6:11]
        if typecheck.isint(atomnumber):
            atomnumber=string.atoi(atomnumber)
        else:
            print '\nAtom number:',atomnumber
            raise ParseAtomError,'Non-int in atonnumber'
        #
        # All we now about the atomname is that it must contain a nondigit character
        #
        atomname=string.strip(line[12:16])
        if not typecheck.containsletter(atomname):
            raise ParseAtomError,'Atom name does not contain a character: %s,\nLINE: %s' %(atomname,line)
        #
        # Also the residuename must contain a letter
        #
        residuename=string.strip(line[17:21])
        if not typecheck.containsletter(residuename):
            print '\nResidue name:',residuename
            raise ParseAtomError,'Residuename does not contain a character'
        # ---------------------/------------------------------------------
        #
        # The ChainID must be a letter, a blank or a digit
        #
        if self.NMRmodel is False or self.readmodels==1 or self.ignore_NMRmodel:
            chainid=line[21]
        else:
            chainid='%d_%s' %(self.NMRmodel,line[21])
        #
        # Check the chainID
        #
        if not typecheck.containsletter(chainid) and chainid!=' ' and not typecheck.isint(chainid):
            print '\nChainID',chainid
            raise ParseAtomError,'Chain ID must be a letter,blank or digit'
        else:
            chainid=string.strip(chainid)
            #print 'Chainid: %s' %chainid
        newchainid=chainid
        #
        # Check that we are still building the same chain
        #
        if not self.chains.has_key(chainid):
            for cid in self.chains.keys():
                if self.chains[cid]=='building' or (self.oldchainid!=newchainid and self.chains[cid]=='building - dummy'):
                    self.chains[cid]='terminated'
            self.chains[chainid]='building'
        else:
            if self.chains[chainid]!='building':
                fail=1
                if chainid=='':
                    #
                    # Find a new chainid
                    #
                    for cid in self.chains.keys():
                        if self.chains[cid]=='building - dummy':
                            chainid=cid
                            fail=None
                            break
                    if chainid=='':
                        for cid in self.useids:
                            if not self.chains.has_key(cid):
                                self.chains[cid]='building - dummy'
                                chainid=cid
                                fail=None
                                break
                if fail:
                    for cid in self.chains.keys():
                        if self.chains[cid]=='building - dummy' and self.oldchainid==newchainid:
                            chainid=cid
                            fail=None
                            break
                    if fail:
                        for cid in self.useids:
                            if not self.chains.has_key(cid):
                                self.chains[cid]='building - dummy'
                                chainid=cid
                                fail=None
                                break
                    if fail:
                        print self.chains.keys()
                        raise ProtoolError, 'I ran out of alternative chain ID names'
        self.oldchainid=newchainid
        # -------------------------------------------------------------------------------------
        #
        # Check the residuenumber
        #
        residuenumber=line[22:27]
        done=False
        suffix=''
        while not done:
            if typecheck.isint(residuenumber):
                residuenumber=string.zfill(string.atoi(residuenumber),G.length_of_residue_numbers)
                residuenumber=residuenumber+suffix
                done=True
            else:
                suffix=residuenumber[-1]+suffix
                residuenumber=residuenumber[:-1]
        #
        # Get and check the coordinates
        #
        if len(line)<54:
            xcoord=None
            ycoord=None
            zcoord=None
        else:
            xcoord=line[30:38]
            ycoord=line[38:46]
            zcoord=line[46:54]
            if not typecheck.isnumber(xcoord) or not typecheck.isnumber(ycoord) or not typecheck.isnumber(zcoord):
                print '\nX,Y,Z',xcoord,ycoord,zcoord
                raise ParseAtomError,'All coordinates must be numbers'
            else:
                xcoord=string.atof(xcoord)
                ycoord=string.atof(ycoord)
                zcoord=string.atof(zcoord)
        #
        # Get the occupancy
        #
        if len(line)<60:
            print '\nNo occupancy for this line'
            occupancy=None
        else:
            occupancy=line[55:60]
            if not typecheck.isnumber(occupancy):
                print '\nOccupancy:',occupancy
                raise ParseAtomError,'Occupancy is not a number'
            else:
                occupancy=string.atof(occupancy)
        #
        # Get the B-factor
        #
        if len(line)<66:
            print '\nNo B-factor for this atom',line
            bfactor=None
        else:
            bfactor=line[60:66]
            if not typecheck.isnumber(bfactor):
                print 'B-factor:',bfactor
                bfactor=None
            else:
                bfactor=string.atof(bfactor)
        #
        # Put the rest of the line in a junk array
        #
        if len(line)>66:
            junk=line[66:]
        else:
            junk=None
        #
        # Finally put everything in self.atomlinedefs
        #
        self.atomlinedefs[self.linenumber]={'NUMBER':atomnumber}
        self.atomlinedefs[self.linenumber]['ATNAME']=atomname
        self.atomlinedefs[self.linenumber]['CHAINID']=chainid
        self.atomlinedefs[self.linenumber]['RESNAME']=residuename
        self.atomlinedefs[self.linenumber]['RESNUM']=residuenumber
        if xcoord!=None:
            self.atomlinedefs[self.linenumber]['X']=xcoord
            self.atomlinedefs[self.linenumber]['Y']=ycoord
            self.atomlinedefs[self.linenumber]['Z']=zcoord
        if occupancy or occupancy==0.0:
            self.atomlinedefs[self.linenumber]['OCCUPANCY']=occupancy
        if bfactor or bfactor==0.0:
            self.atomlinedefs[self.linenumber]['B-FACTOR']=bfactor
        if hetatm:
            self.atomlinedefs[self.linenumber]['HETATM']=1
        if junk:
            #
            # This is a tag
            #
            import string
            self.atomlinedefs[self.linenumber]['tag']=string.strip(junk)
                
        #
        # What is this atom bound to?
        #
        self.atomlinedefs[self.linenumber]['BOUND_TO']=[]
        return

    def parsehet(self,line):
        self.parseatom(line,1)
        return

    def parseheader(self,line):
        #if not hasattr(self,'header'):
        #    self.header=''
        self.header=self.header+line
        self.header=self.header.replace('HEADER','').strip()
        return

    def parsecompnd(self,line):
        self.compnd=self.compnd+line.replace('COMPND','').strip()
        return
    
    def parsesource(self,line):
        return

    def parseauthor(self,line):
        return

    def parserevdat(self,line):
        return

    def parseremark(self,line):
        return

    def parsescale(self,line):
        self.scale.append(line.strip())
        return

    def parsejournal(self,line):
        return

    def parseseqres(self,line):
        return

    def parsehelix(self,line):
        return
    
    def parsesheet(self,line):
        return

    def parseturn(self,line):
        return

    def parsessbond(self,line):
        return

    def parsecryst(self,line):
        self.cryst.append(line.strip())
        #
        # Find the space group
        #
        pos=self.cryst[0].find(' P')
        self.spacegroup=self.cryst[0][pos:].strip()
        if len(self.spacegroup.split())>4:
            sp=self.spacegroup.split()
            import string
            self.spacegroup=string.join(sp[:4],' ').strip()
        sg=self.spacegroup.strip()
        if sg[0]=='P' and sg[1]!=' ':
            sg=sg[0]+' '+sg[1:]
        conv={'1':'P 1','2':'P 2','P 1 21 1':'P 21','P21':'P 21'}
        if conv.has_key(sg):
            sg=conv[sg]
        self.spacegroup=sg
        return

    def parseorig(self,line):
        self.orig.append(line.strip())
        return

    def parseconnect(self,line):
        return

    def parsemaster(self,line):
        return

    def parseter(self,line):
        #
        # Molecule stops here
        #
        for cid in self.chains.keys() :
            if (self.chains[cid]=='building' or self.chains[cid]=='building - dummy') and self.parse_terms:
                self.chains[cid]='terminated'
                break
        return

    def parsejournal(self,line):
        return

    def parsefootnote(self,line):
        return

    def parsehetheader(self,line):
        return

    def parseformul(self,line):
        return

    def parseend(self,line):
        return

    def parseexpdta(self,line):
        return

    def parsemodel(self,line):
        self.NMRmodel=int(line.strip().split()[1])
        if self.models_already_read>=self.readmodels:
            self.readmodel=False
        else:
            self.readmodel=True
            self.models_already_read=self.models_already_read+1
        return

    def parseendmodel(self,line):
        self.NMRmodel=False
        self.readmodel=True
        return

    def renumber(self):
        #
        # Sort all atoms and renumber them
        #
        atoms=self.atoms.keys()
        atoms.sort()
        count=1
        for atom in atoms:
            self.atoms[atom]['NUMBER']=count
            count=count+1
        return
    
#
# ------------
#

class fast_PDBparser(PDBparser):
    """Fast class for parsing PDBlines"""

    def parseline(self,line):
        type=string.strip(string.upper(line[:6]))
        if type=='ATOM':
            self.parseatom(line)
        elif type=='HETATM':
            self.parsehet(line)
        else:
            pass
        return    

    #
    # ---------------------------------------
    #

    def parseatom(self,line,hetatm=None):
        self.linenumber=self.linenumber+1
        self.length_of_residue_numbers=4
        atomnumber=line[6:11]
        atomname=string.strip(line[12:16])
        residuename=string.strip(line[16:21])
        chainid=string.strip(line[21])
        newchainid=chainid
        if not self.chains.has_key(chainid):
            for cid in self.chains.keys():
                if self.chains[cid]=='building' or (self.oldchainid!=newchainid and self.chains[cid]=='building - dummy'):
                    self.chains[cid]='terminated'
            self.chains[chainid]='building'
        else:
            if self.chains[chainid]!='building':
                fail=1
                if chainid=='':
                    #
                    # Find a new chainid
                    #
                    for cid in self.chains.keys():
                        if self.chains[cid]=='building - dummy':
                            chainid=cid
                            fail=None
                            break
                    if chainid=='':
                        for cid in self.useids:
                            if not self.chains.has_key(cid):
                                self.chains[cid]='building - dummy'
                                chainid=cid
                                fail=None
                                break
                if fail:
                    for cid in self.chains.keys():
                        if self.chains[cid]=='building - dummy' and self.oldchainid==newchainid:
                            chainid=cid
                            fail=None
                            break
                    if fail:
                        for cid in self.useids:
                            if not self.chains.has_key(cid):
                                self.chains[cid]='building - dummy'
                                chainid=cid
                                fail=None
                                break
                    if fail:
                        raise ProtoolError, 'I ran out of alternative chain ID names'
        #
        # -----
        #
        self.oldchainid=newchainid
        residuenumber=string.zfill(string.atoi(line[22:27]),self.length_of_residue_numbers)
        xcoord=string.atof(line[30:38])
        ycoord=string.atof(line[38:46])
        zcoord=string.atof(line[46:54])
        self.atomlinedefs[self.linenumber]={'NUMBER':atomnumber}
        self.atomlinedefs[self.linenumber]['ATNAME']=atomname
        self.atomlinedefs[self.linenumber]['CHAINID']=chainid
        self.atomlinedefs[self.linenumber]['RESNAME']=residuename
        self.atomlinedefs[self.linenumber]['RESNUM']=residuenumber
        self.atomlinedefs[self.linenumber]['X']=xcoord
        self.atomlinedefs[self.linenumber]['Y']=ycoord
        self.atomlinedefs[self.linenumber]['Z']=zcoord
        if hetatm:
            self.atomlinedefs[self.linenumber]['HETATM']=1
        return

