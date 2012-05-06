#!/usr/bin/env python

import numpy 
import matplotlib.pyplot as plt
    
#
# Define the generic error for chemical shift data
#
errors={'N':0.1,'H':0.03,'HA':0.03}


#
# ----
#

def yasaraplot():
    """Make a plot of all ghosts on the structure"""
    ghosts=getdata()
    #
    rescolor={}
    tgcount=0
    tgcolstep=360.0/len(ghosts.keys())
    print 'tgcolstep: %5.2f' %tgcolstep
    #
    colors=[240,0,60,120,180,300,1030,1150]
    del ghosts[':0048:ASP']
    del ghosts[':0066:ASP']
    for tg in sorted(ghosts.keys()):
        for residue in ghosts[tg].keys():
            if not rescolor.has_key(residue):
                rescolor[residue]={}
            for atom in ['N','H']:
                if ghosts[tg][residue].has_key(atom):
                    if ghosts[tg][residue][atom]!='absent':
                        if not rescolor[residue].has_key(atom):
                            rescolor[residue][atom]=[colors[tgcount]]
                        else:
                            rescolor[residue][atom].append(colors[tgcount])
        print tg, colors[tgcount]
        tgcount=tgcount+1
    #
    import Yasara_handler
    Y=Yasara_handler.yasara_handler()
    obj=Y.load_mol('2LZT.pdb')[0]
    for plot in ['N','H']:
        for residue in sorted(rescolor.keys()):
            resnum=int(residue[1:])
            if rescolor[residue].has_key(plot):
                #if len(rescolor[residue][plot])==1:
                
                color=sum(rescolor[residue][plot])/len(rescolor[residue][plot])
                print resnum,rescolor[residue][plot],color
                #else:  
                #    color=1270
            else: 
                color=1059
            Y.yasara.ColorRes('%d' %(resnum),'%d' %color)
        raw_input('is this ok?')
    return

#
# -------------------------------------------
#
# Dummy instance to use scoring function
#    
import Create_epsmap
class dummy_epsmap(Create_epsmap.iterative_epsmap):
    def __init__(self):
        class X:
            def __init__(this):
                this.use_absent=True
                return
                
        self.options=X()
        return

#
# ----
#

class ghost_analyzer:

    def __init__(self,options):
        """Init and store options"""
        self.options=options
        #
        # Create dirs
        #
        import os
        if not os.path.isdir(options.plotdir):
            os.mkdir(options.plotdir)
        #
        # Load all data
        #
        exp_ghosts=self.get_expdata(options.expfile)
        del exp_ghosts[':0015:HIS'] # We trust none of the His15 data # This needs to be fixed
        #
        # filter my conformational change
        #
        self.s_change={}
        if options.filter_confchange:
            deleted=0
            self.read_confchange(exp_ghosts.keys())
            for tg in exp_ghosts.keys():
                if tg in self.conf_change.keys():
                    #
                    # Remove the ghosts of the residues that are moving
                    #
                    atom=options.atom
                    for residue in exp_ghosts[tg].keys():
                        if exp_ghosts[tg][residue].has_key(atom):
                            #print exp_ghosts[tg][residue][atom]
                            if self.conf_change[tg].has_key(atom):
                                intres=int(residue[1:])
                                if intres in self.conf_change[tg][atom]['res']:
                                    if not self.s_change.has_key(tg):
                                        self.s_change[tg]={}
                                    if not self.s_change[tg].has_key(residue):
                                        self.s_change[tg][residue]={}
                                    self.s_change[tg][residue][atom]=exp_ghosts[tg][residue][atom]
                                    del exp_ghosts[tg][residue][atom]
                                    deleted=deleted+1
                else:
                    print 'Deleting %s because we have no conformational restraints' %tg
                    del exp_data[tg]
            print 'Removed %d ghosts' %deleted
        #
        # Read and filter based on CPMG data
        #
        if options.CPMG_filter:
            fd=open('all_restraints/CPMG_data/CPMG_maxchange.csv')
            lines=fd.readlines()
            fd.close()
            for line in lines:
                if line[0]=='#':
                    continue
                line=line.strip().split(',')
                cpmgres=int(line[0])
                cpmgval=float(line[1])
                
                print line
                for tg in exp_ghosts.keys():
                    import string
                    residue=':%s' %(string.zfill(cpmgres,4))
                    if not exp_ghosts[tg][residue].has_key('N'):   
                        continue
                    Nghost=exp_ghosts[tg][residue]['N']
                    number=int(residue[1:])
                    if Nghost=='absent':
                        if cpmgres==number:
                            del exp_ghosts[tg][residue]
                    else:
                        if cpmgres==number:
                            if Nghost.find(';')!=-1:
                                Nghost=Nghost.split(';')[1]
                            import string
                            if not Nghost[0] in string.digits:
                                Nghost=float(Nghost[1:])
                            print type(cpmgval),type(Nghost)
                            if cpmgval>abs(0.5*float(Nghost)):
                                del exp_ghosts[tg][residue]
                            
        #
        # Display stats on experimental data
        #
        if options.stats:
            self.stats_forpaper(exp_ghosts)
            return

        if options.anatype=='eps':
            #
            # Get the calculated ghosts
            #
            if options.dofull:
                self.load_all_calculated(exp_ghosts)
                return
            #
            # Just load normally - i.e. just make correlation for a single epstype
            #
            calc_ghosts=self.load_calcghosts(options)
            #
            # Find optimal eps?
            #
            big_satisfied=self.find_opteps(exp_ghosts,calc_ghosts)
            #
            # Then do the plots/analyses that were requested
            #
            if options.restraints:
                self.restraint_plot(big_satisfied)
        elif options.anatype=='energies':
            self.energy_analysis(exp_ghosts)
 
        else:
            raise Exception('Unknown analysis type')


    def interpret_ranges(self,exp_value,atom):
        """Given a number of a range, get the true value + and error"""
        exp_error=errors[atom]
        if type(exp_value) is type('jkfjd'):
            if len(exp_value.split(';'))==2:
                s=[]
                for val in str(exp_value).split(';'):
                    if val[0] in ['<','>']:
                        val=val[1:]
                    s.append(float(val))
                exp_error=abs(s[0]-s[1])
                exp_value=float(sum(s))/len(s)
        return exp_value,exp_error

    #
    # ----------------
    #

    def energy_analysis(self,exp_ghosts):
        """
        # Calculate energies from ghost titrations
        """
        #
        # First instantiate the class for converting between chemical shift and energies
        #
        import get_dEF
        EF=get_dEF.map_struct(self.options.pdbfile)
        #EF.build_Hs()
        #
        # Now load all ghosts and calculate energies
        #
        excludes=self.find_excludes()
        Edict={}
        for titgroup in [':0035:GLU']: #sorted(exp_ghosts.keys()):
            Edict[titgroup]={}
            for residue in sorted(exp_ghosts[titgroup].keys()):
                if excludes.has_key(titgroup):
                    if residue in excludes[titgroup]:
                        continue
                vals=[]
                D=exp_ghosts[titgroup][residue]
                for nucleus in D.keys():
                    if nucleus=='HA':
                        continue
                    thisval=None
                    if D[nucleus]!='absent':
                        if D[nucleus][0]=='q':
                            thisval=D[nucleus][1:]
                        else:
                            thisval=D[nucleus]
                    else:
                        thisval=0.0
                    #
                    # Deal with range
                    #
                    exp_value,exp_error=self.interpret_ranges(thisval,nucleus)
                    # convert to kJ/mol
                    exp_value=EF.get_energy(exp_value,'%s:%s' %(residue,nucleus),titgroup)
                    vals.append([exp_value,exp_error,nucleus])
                Edict[titgroup][residue]=vals[:]
        #
        # Load PDB file and calculate distance from E35
        #
        import Protool
        X=Protool.structureIO()
        X.readpdb(self.options.pdbfile)
        #
        # Plot everything
        #
        tg=':0035:GLU'
        xs=[]
        ys=[]
        for residue in sorted(Edict[tg].keys()):
            if not X.residues.has_key(residue):
                continue
            dist=X.dist(':0035:CG','%s:N' %residue)
            for exp_value,exp_error,atom in Edict[tg][residue]:
                xs.append(dist)
                val=float(exp_value)
                ys.append(abs(val))
        import pylab
        pylab.plot(xs,ys,'ro',label='Ghosts')
        print 'Average Ghost energy: %5.2f kJ/mol' %(sum(ys)/len(ys))
        #
        # Read ddGcat values
        #
        fd=open('ddGcat.csv')
        lines=fd.readlines()
        fd.close()
        cat=[]
        for line in lines[1:]:
            sp=line.split(',')
            resnum=int(sp[0])
            mutation=sp[1]
            try:
                ddG=float(sp[2])
            except:
                continue
            cat.append([resnum,mutation,ddG])
        #
        # Plot it
        #
        xs=[]
        ys=[]
        for resnum,mutation,ddG in cat:
            import string
            residue=':%s' %(string.zfill(resnum,4))
            md=9999.9
            for atom in X.residues[residue]:
                md=min(md,X.dist(':0035:CG',atom))
            xs.append(md)
            ys.append(abs(ddG))
        pylab.plot(xs,ys,'bo',label='ddGcat')
        print 'Average ddGcat: %5.1f kJ/mol' %(sum(ys)/len(ys))
        pylab.xlabel('Distance from E35CG')
        pylab.ylabel('Energy change (kJ/mol)')
        pylab.legend()
        pylab.title('Correlation between ghosts and change in activity')
        pylab.show()
        return
                
            
 
    def _calc_Efield(self,dCS,angle,atom,charge):
        #
        # Maybe we want to do the sexy comparison with changes in activities?
        #
        """Given a change in chemical shift and an angle, calculate the electric field change"""
        cos_angle=math.cos(math.radians(angle))
        val=dCS/(NSP[atom]*e*ppm_au*cos_angle*1E6*-1)
        return val
        
	#
	# ------
	#
	
    def read_confchange(self,tgs):
        """Read the files with information on the conformational change"""
        import os
        self.conf_change={}
        for tg in tgs:
            self.conf_change[tg]={}
            for atom in ['N','H']:
                sp=tg.split(':')
                filename='%3s%d%1s_confchange.csv' %(sp[2],int(sp[1]),atom)
                filename=os.path.join('all_restraints','confchange',filename)
                if os.path.isfile(filename):
                    print 'Reading conformational change data from %s' %filename
                    data=self.read_csv(filename)
                    self.conf_change[tg][atom]=data.copy()
        return

    #
    # ----
    #

    def read_csv(self,csvfile):
        """Read a CSV file"""
        import csv
        f=open(csvfile,'r')
        cr = csv.reader(f)
        tmp_names = cr.next()
        names=[]
        for name in tmp_names:
            names.append(name.lower().strip())
        data={}
        #
        for n in names:
            n=n.lower()
            data[n]=[]
        for row in cr:
            for i in range(len(row)):
                n=names[i]
                if n=='res':
                    x=int(row[i])
                elif n=='shift':
                    x=float(row[i])
                else:
                    x=row[i]
                data[n].append(x)
        return data
		
		
        
    #
    # -----
    #
    
    def get_expdata(self,expfile):
        """Load the experimental data"""
        print 'Reading experimental data from: %s' %expfile
        import pickle
        fd=open(expfile)
        data=pickle.load(fd)
        fd.close()
        #
        Ndata=[]
        Hdata=[]
        #
        # Fill in the keys for the groups where we have no ghosts
        #
        for add in [':0048:ASP',':0066:ASP']:
            data[add]={}
            ress=data[':0035:GLU'].keys()
            for r in ress:
                data[add][r]={}
        return data
    
    #
    # -----
    #

    def RMSD(self,xs,ys,experrors,calcerrors):
        # Get the RMSD of the prediction
        import math
        diffs=[]
        for count in range(len(xs)):
            diff=abs(xs[count]-ys[count])
            if diff<abs(experrors[count])+abs(calcerrors[count]):
                diff=0.0
            else:
                diff=diff-abs(experrors[count])-abs(calcerrors[count])
            diffs.append(math.pow(diff,2))
        rmsd=math.sqrt(sum(diffs)/len(diffs))
        return rmsd
    
    #
    # --------
    #

    def find_excludes(self,options=None,exp_ghosts=None):
        """find residues to exclude based on an analysis of the structure"""
        excludes={':0007:GLU':[':0004:',':0005',':0006:'],
        ':0015:HIS':[':0011',':0012',':0013',':0014',':0016',':0088',':0090',':0089',':0092'],
        ':0018:ASP':[':0019',':0025',':0026',':0027',':0017',':0021',':0022',':0023',':0024'],
        ':0035:GLU':[':0110',':0036'],':0052:ASP':[':0053',':0045'],
        ':0052:ASP':[':0059','0058',':0045',':0044',':0047',':0053',':0051'],
        ':0101:ASP':[':0102',':0098'],
        ':0119:ASP':[':0121'],
        ':0129:CTERM':[':0128']}
        return excludes
    
    #
    # ------
    #

    def find_excludes_from_structure(self):
        """This is an old function for selecting residues to exclude. Does not work that great"""
        excludes={}
        #
        # Get the PDB file
        #
        import Protool
        X=Protool.structureIO()
        X.readpdb(options.pdbfile)
        #
        # Find all residues next to tg and next to Hbond-tg + residues with bb in VDW distance of sidechain
        #
        for tg in sorted(exp_ghosts.keys()):
            excludes[tg]=[]
            import string
            tg_res=string.join(tg.split(':')[:2],':')
            tg_num=int(tg.split(':')[1])
            CID=tg.split(':')[0]
            # Before - After
            excludes[tg].append('%s:%s' %(CID,string.zfill(tg_num-1,4)))
            excludes[tg].append('%s:%s' %(CID,string.zfill(tg_num+1,4)))
            #
            # Find Hbonding and VdW partners for side chain atoms
            #        
            for res2 in X.residues.keys():
                if res2==tg_res:
                    continue
                for atom in X.residues[tg_res]:
                    if X.is_backbone(atom) or X.is_hydrogen(atom):
                        continue
                    for atom2 in X.residues[res2]:
                        if X.dist(atom,atom2)<6.0:
                            if not res2 in excludes[tg]:
                                excludes[tg].append(res2)
                                break
        for tg in sorted(excludes.keys()):
            print tg
            print excludes[tg]
            print '-------'
        return excludes
        
    #
    # --------
    #    

    def calculate_average(self,data):
        """Calculate the average ghost observed and the standard deviation"""
        for datatype in data.keys():
            for diel in data[datatype].keys():
                for TG in data[datatype][diel].keys():
                    for residue in data[datatype][diel][TG].keys():
                        for nucleus in data[datatype][diel][TG][residue].keys():
                            try:
                                values=data[datatype][diel][TG][residue][nucleus]
                            except KeyError:
                                print 'Skipping %d for %s' %(diel)
                                continue
                            import stats
                            avg=stats.average(values)
                            SD=stats.stddev(values)
                            data[datatype][diel][TG][residue][nucleus]=[avg,SD]
        return data

    #
    # -----
    #
    
    def load_all_calculated(self,exp_ghosts):
        """Load all current predictions of ghosts"""
        predictions=[['Xray_avg.pickle','Xray'],['MD_calcghosts_nowat.pickle','MD_nowat'],['average_calcghosts_predrag_all.pickle','MD'],['Xrays/all_ghosts_2LZT_H.pdb.pickle','2LZT']] # ['average_calcghosts_predragGB.pickle','GB'],
        alldata={}
        for filename,ensemble_name in predictions:
            import os
            filename=os.path.join(os.getcwd(),self.options.calcghostdir,filename)
            calc_ghosts=self.load_calcghosts(self.options,filename)
            print 
            print '=================%s====================' %ensemble_name
            print
            big_satisfied=self.find_opteps(exp_ghosts,calc_ghosts)
            alldata[ensemble_name]=self.restraint_plot(big_satisfied)
        #
        # Do the full plot
        #
        import pylab
        colours=['k','r','b','g','y','c','m','grey','orange','pink']
        pylab.clf()
        count=0
        sumsum=[]
        for ensname in sorted(alldata.keys()):
            for xs,present,absent,method in alldata[ensname]:
                text='%s:%s' %(ensname,method)
                col=colours[count]
                marker=''
                if ensname=='GB':
                    marker='o'
                pylab.plot(xs,present,'-',marker=marker,color=col,label=text,linewidth=3)
                pylab.plot(xs,absent,'--',marker=marker,color=col,linewidth=3)
                count=count+1
                if count==len(colours):
                    count=0
                #
                # Find the best eps
                #
                ecount=0
                maxval=-999.9
                besteps=0.0
                for val in present:
                    if val>maxval:
                        maxval=val
                        besteps=ecount
                    ecount=ecount+1
                print 'Best eps for %s present: %5.1f with %5.1f %% satisfied, at this eps absent is %5.1f' %(text,xs[besteps],maxval,absent[besteps])
                sumsum.append([maxval+absent[besteps],xs[besteps],maxval,absent[besteps],text])
        sumsum.sort()
        print 'Sum, Eps, %present, %absent, method'
        for sum,eps,present,absent,text in sumsum:
            print '%5.1f, %5.1f, %5.1f, %5.1f %s' %(sum,eps,present,absent,text)
                    


        #
        # Finish the plot
        #
        pylab.ylim((0.0,100.0))
        pylab.xlim((0.0,30.0))
        pylab.legend()
        pylab.title('Satisfied restraints for %s' %self.options.atom)
        pylab.xlabel('Dielectric constant')
        pylab.ylabel('% restraints satisfied')
        pylab.savefig('bigplot.png',dpi=300)
        pylab.show()
        return
    
    #
    # -----
    #

    def load_calcghosts(self,options,loadfile=None): 
        """Load the calculated ghosts from a single file or a directory"""
        if options.loaddir:
            # Load all files in the directory specified
            import os
            if os.path.isdir(options.loaddir):
                files=os.listdir(options.loaddir)
                okfiles=[]
                for fn in files:
                    txt='_H.pdb.pickle'
                    size=len(txt)
                    if fn[-size:]!=txt:
                        continue
                    okfiles.append(fn)
                
                print 'Found %3d files with calculated ghosts' %len(okfiles)
                C=False
                count=1
                for fn in sorted(okfiles):
                    print 'Reading # %3d with name: %s' %(count,fn)
                    count=count+1
                    realfile=os.path.join(options.loaddir,fn)
                    fd=open(realfile)
                    import pickle
                    data=pickle.load(fd)
                    fd.close()
                    if not C:
                        C=data.copy()
                    #
                    for datatype in data.keys():
                        for diel in data[datatype].keys():
                            if not C[datatype].has_key(diel):
                                C[datatype][diel]={}
                            for TG in data[datatype][diel].keys():
                                if not C[datatype][diel].has_key(TG):
                                    C[datatype][diel][TG]={}
                                for residue in data[datatype][diel][TG].keys():
                                    if not C[datatype][diel][TG].has_key(residue):
                                        C[datatype][diel][TG][residue]={}
                                    for nucleus in data[datatype][diel][TG][residue].keys():
                                        if not C[datatype][diel][TG][residue].has_key(nucleus):
                                            C[datatype][diel][TG][residue][nucleus]=[]
                                        #
                                        value=data[datatype][diel][TG][residue][nucleus]
                                        
                                        if not type(C[datatype][diel][TG][residue][nucleus]) is type([]):
                                            C[datatype][diel][TG][residue][nucleus]=[]
                                        C[datatype][diel][TG][residue][nucleus].append(value)
                #
                # Get the averages and standard deviations    
                #
                calc_ghosts=self.calculate_average(C)
                avgfile=options.avgfile
                print 'Saving average ghosts to %s' %avgfile
                fd=open(avgfile,'w')
                import pickle
                pickle.dump(calc_ghosts,fd)
                fd.close()
            else:
                raise Exception('Not a directory: %s' %options.loaddir)
        else:
            # 
            # Load just a single file
            #
            filename=options.calcghosts
            if loadfile:
                filename=loadfile 
            fd=open(filename)
            import pickle
            print 'Reading calculated ghosts from %s' %filename
            calc_ghosts=pickle.load(fd)
            fd.close()
        return calc_ghosts
    
    #
    # -------
    #
    
    def stats_forpaper(self,exp_ghosts):
        """Get the stats for the paper"""
        # Get the excludes
        excludes=self.find_excludes(options,exp_ghosts)
        # Loop over all ghosts
        for atom in ['N','H','HA']:
            values=[]
            absent=[]
            for tg in sorted(exp_ghosts.keys()):
                for residue in sorted(exp_ghosts[tg].keys()):
                    if excludes.has_key(tg):
                        if residue in excludes[tg]:
                            continue
                    #
                    # Get the value
                    #
                    if exp_ghosts[tg][residue].has_key(atom):
                        exp_value=exp_ghosts[tg][residue][atom]
                        exp_error=errors[atom]
                        if exp_ghosts[tg][residue].has_key(atom+'_error'):
                            exp_error=exp_ghosts[tg][residue][atom+'_error']
                        # 
                        # Deal with ranges
                        #
                        if exp_value[0]=='q' and options.use_questionable:
                            exp_value=exp_value[1:]
                        # 
                        # Deal with ranges  - we need to do this better
                        #
                        if len(exp_value.split(';'))==2:
                            s=[]
                            for val in exp_value.split(';'):
                                if val[0] in ['<','>']:
                                    val=val[1:]
                                s.append(float(val))
                            exp_error=abs(s[0]-s[1])
                            exp_value=float(sum(s))/len(s)
                        if exp_value=='absent':
                            absent.append(residue)
                        else:
                            values.append(abs(float(exp_value)))
            #
            # Calculate average
            #
            import stats
            avg=0.0
            SD=0.0
            if len(values)>0:
                avg=stats.average(values)
                SD=stats.stddev(values)
            print '%2s ghost titrations: %3d, avg: %5.2f (%5.2f), %3d absent ghost titrations' %(atom,len(values),avg,SD,len(absent))
        return
        
    #
    # ------
    #
    
    def find_opteps(self,exp_ghosts,calc_ghosts):
        """Find the optimal eps value for each residue and overall"""
        #
        # Instantiate the dummy scoring function
        #
        S=dummy_epsmap()
        #
        # Find the ghosts to exclude
        #
        excludes=self.find_excludes(options,exp_ghosts)
        #
        # Initialize arrays
        #
        dep={}
        bestepses={}
        bigdata={}
        big_satisfied={}
        #
        # Start looping
        #
        for method in calc_ghosts.keys():
            #
            # Smoothed PBE calculations gave no improvement
            #
            if method=='sPBE':
                continue
            #
            # Method is the calculational method
            #
            dep[method]={}
            bestepses[method]={}
            big_satisfied[method]={}
            #
            bigdata[method]={}
            #
            for tg in sorted(exp_ghosts.keys()):
                tgcount=0
                wrong_sign=[]
                #
                # First find best eps
                #
                scores=[]
                data={}
                bigdata[method][tg]={}
                big_satisfied[method][tg]={}
                for eps in sorted(calc_ghosts[method]):
                    #
                    # Get the RMSD for this eps
                    #
                    xs=[]
                    ys=[]
                    experrors=[]
                    calcerrors=[]
                    residues=[]
                    big_satisfied[method][tg][eps]={}
                    #
                    # Loop over all residues 
                    #
                    atom=options.atom
                    if not calc_ghosts[method][eps].has_key(tg):
                        #print 'no TG in calc',tg
                        continue
                    #print 'CHECK',tg,eps,method
                    for residue in sorted(exp_ghosts[tg].keys()):
                        
                        #
                        # Catch for uncalculated ghosts from MD sims
                        #

                        if not calc_ghosts[method][eps][tg].has_key(residue):
                            #print 'missing residue',residue
                            
                            continue
                        #print calc_ghosts[method][eps][tg].keys()
                        #
                        #
                        #
                        if excludes.has_key(tg):
                            if residue in excludes[tg]:
                                continue
                        if exp_ghosts[tg][residue].has_key(atom):
                            exp_value=exp_ghosts[tg][residue][atom]
                            exp_error=errors[atom]
                            if exp_ghosts[tg][residue].has_key(atom+'_error'):
                                exp_error=exp_ghosts[tg][residue][atom+'_error']
                            #
                            # Get the calculated value - and the error if specified
                            #
                            if not calc_ghosts[method][eps][tg][residue].has_key(atom):
                                #print 'No calculated result for %s %d %s %s %s' %(method,eps,tg,residue,atom)
                                continue
                            #else:
                            #    print 'Found calc result for  %s %d %s %s %s' %(method,eps,tg,residue,atom)
                            calc_value=calc_ghosts[method][eps][tg][residue][atom]
                            calc_error=0.0
                            if type(calc_value) is type([]):
                                calc_error=calc_value[1]
                                calc_value=calc_value[0]
                            #
                            if options.no_calcerror:
                                calc_error=0.0
                            #
                            error,satisfied,abs_sat,tot_restraints,tot_abs,real_error=S.get_error_sub(exp_value,calc_value,exp_error,atom,calc_error=calc_error)
                            big_satisfied[method][tg][eps][residue]=[satisfied,abs_sat,tot_abs]
                            exp_error=errors[atom]
                            if exp_value[0]=='q':
                                exp_value=exp_value[1:]
                            if len(exp_value.split(';'))==2:
                                s=[]
                                for val in exp_value.split(';'):
                                    if val[0] in ['<','>']:
                                        val=val[1:]
                                    s.append(float(val))
                                exp_error=abs(s[0]-s[1])
                                exp_value=float(sum(s))/len(s)
                            #
                            # Append the values to xs and ys
                            #
                            if exp_value!='absent':
                                exp_value=float(exp_value)
                                xs.append(exp_value)
                                experrors.append(exp_error)
                                #
                                ys.append(calc_value)
                                calcerrors.append(calc_error)
                                residues.append(residue)
                                #
                                if exp_value*calc_value<0.0:
                                    wrong_sign.append(residue)
                    #
                    # Compared all ghosts for this tg at this eps
                    #
                    if len(xs)==0:
                        continue
                    scores.append([self.RMSD(xs,ys,experrors,calcerrors),eps])
                    bigdata[method][tg][eps]=[xs,ys,experrors,calcerrors,residues]
                    data[eps]=[xs,ys,residues,experrors,calcerrors]
                #
                # Looped over all eps
                #
                dep[method][tg]=None
                if len(scores)==0:
                    print 'Skipping %s because I have no points left' %tg
                    continue
                #
                # Format for dep
                #
                nscore=[]
                neps=[]
                for rmsd,eps in scores:
                    nscore.append(rmsd)
                    neps.append(eps/10.0)
                dep[method][tg]=[nscore[:],neps[:]]
                scores.sort()
                best_eps=scores[0][1]
                print 'Best eps for %12s using %12s is %5.1f, RMSD: %5.2f, numpoints: %3d' %(tg,method,best_eps/10.0,scores[0][0],len(data[best_eps][0]))
                #
                # Dig out the points we need
                #
                xs,ys,residues,experrors,calcerrors=data[best_eps]
                bestepses[method][tg]=data[best_eps][:]+[best_eps]
                import pylab
                xs1=[]
                ys1=[]
                err1=[]
                calc_err1=[]
                #
                xs2=[]
                ys2=[]
                err2=[]
                calc_err2=[]
                #
                
                for count in range(len(xs)):
                    res=residues[count].split(':')[-1]
                    if int(res)>40 and int(res)<=85:
                        xs1.append(xs[count])
                        ys1.append(ys[count])
                        err1.append(experrors[count])
                        calc_err1.append(calcerrors[count])
                    else:
                        xs2.append(xs[count])
                        ys2.append(ys[count])
                        err2.append(experrors[count])
                        calc_err2.append(calcerrors[count])
                #
                if len(xs1)>0:
                    pylab.errorbar(xs1,ys1,xerr=err1,yerr=calc_err1,fmt='ro',label='beta')
                if len(xs2)>0:
                    pylab.errorbar(xs2,ys2,xerr=err2,yerr=calc_err2,fmt='bo',label='alpha')
                #
                import stats
                if len(xs)==0:
                    continue
                #print xs
                #print ys
                #corr=stats.correlation(xs,ys)
                rmsd=self.RMSD(xs,ys,experrors,calcerrors)
                #
                pylab.plot(xs,xs,'y-')
                pylab.plot([min(xs),max(xs)],[0,0],'g-')
                pylab.xlabel('Experimental dCS')
                pylab.ylabel('Calculated dCS')
                pylab.title('%s %s, atom: %s, RMSD: %5.3f, eps_opt: %5.1f' %(method,tg,atom,rmsd,best_eps/10.0))
                pylab.legend(loc=8)
                pylab.savefig('%s/Opteps_%s_%s_%s.png' %(options.plotdir,method,tg,atom))
                #pylab.show()
                pylab.clf()
        # --------------------------
        #
        # Done with everything for this method
        # Find best overall RMSD
        #
        print 
        print 'Finding best overall eps for each method'
        print
        alldata={}
        for method in bigdata.keys():
            tgs=sorted(bigdata[method].keys())
            tg0=':0035:GLU'
            rmsds=[]
            alldata[method]={}
            for eps in bigdata[method][tg0].keys():
                xs=[]
                ys=[]
                experrors=[]
                residues=[]
                for tg in tgs:
                    if bigdata[method][tg].has_key(eps):
                        xs=xs+bigdata[method][tg][eps][0]
                        ys=ys+bigdata[method][tg][eps][1]
                        experrors=experrors+bigdata[method][tg][eps][2]
                        calcerrors=calcerrors+bigdata[method][tg][eps][3]
                        residues.append([tg,bigdata[method][tg][eps][4]])
                rmsds.append([self.RMSD(xs,ys,experrors,calcerrors),eps])
                alldata[method][eps]=[xs,ys,experrors,calcerrors,residues]
            rmsds.sort()
            #print rmsds
            #print tgs
            #print bigdata[method][tg0].keys()
            besteps=rmsds[0][1]
            print 'Best eps for %12s is %5.1f with RMSD: %6.3f' %(method,besteps,rmsds[0][0])
            import pylab
            xs,ys,experrors,calcerrors,residues=alldata[method][besteps]
        return big_satisfied
        
    #
    # -----
    #
    
    def structure_size_error(self):
        """Plot the size of the erors on the structure"""
        #
        # Let Yasara show the structure and color it according to the size of the error
        #
        raise Exception('function not completed')
        import Yasara
        Y=Yasara.yasara_handler()
        obj=Y.load_mol('2LZT.pdb')
        Y.yasara.run('HideAtom Sidechain')
        tgnum=int(tg.split(':')[1])
        Y.yasara.ColorRes('%d' %tgnum,'magenta')
        Y.yasara.run('ShowRes %s' %tgnum)
        for residue in sorted(rescolor.keys()):
            resnum=int(residue[1:])
            if rescolor.has_key(residue):
                color=rescolor[residue]
            else: 
                color=1059
            Y.yasara.ColorRes('%d' %(resnum),'%d' %color)
        print sorted(wrong_sign)
        raw_input('is this ok?')
        return
        
    #
    # -----
    #
        
    def restraint_plot(self,big_satisfied):
        """
        Make a plot of restraints satisfied vs. eps for each method
        """
        import pylab
        plotdata=[]
        for method in big_satisfied.keys():
            xs=[]
            present=[]
            absent=[]
            tgs=big_satisfied[method].keys()
            for eps in sorted(big_satisfied[method][tgs[0]].keys()):
                sat=[]
                abs_sat=[]
                for tg in tgs:
                    for residue in big_satisfied[method][tg][eps].keys():
                        data=big_satisfied[method][tg][eps][residue]
                        if data[2]==1:
                            abs_sat.append(data[1])
                        else:
                            sat.append(data[0])
                sat=sum(sat)/float(len(sat))*100.0
                abse=sum(abs_sat)/float(len(abs_sat))*100.0
                present.append(sat)
                xs.append(eps/10.0)
                absent.append(abse)
            #
            pylab.plot(xs,present,'-',label=method,linewidth=3)
            pylab.plot(xs,absent,'--',label='abs %s' %method,linewidth=3)
            plotdata.append([xs,present,absent,method])
        pylab.legend()
        pylab.ylim((0.0,100.0))
        pylab.xlabel('Dielectric constant')
        pylab.ylabel('% of restraints satisfied')
        pylab.savefig('%s/restraints.png' %self.options.plotdir,dpi=300)
        #pylab.show()
        return plotdata
        
        
        
#
# --------        
#
    
if __name__=='__main__':
    print
    print 'Analysis of ghost titrations'
    print 'Jens Erik Nielsen, 2011'
    print
    import optparse
    parser=optparse.OptionParser()
    #
    parser.add_option('-x','--experimentalghosts',dest='expfile',default='all_restraints/new_restraints.pickle',type='string',
        help='Pickle file containing the experimentally determined ghost restraints. Default: %default')
    #
    parser.add_option('--heatmap',dest='heatmap',action='store_true',default=False,
        help='Make heatmap plot of titrations. Default: %default')
    parser.add_option('-y','--yasarafig',dest='yasarafig',action='store_true',default=False,
        help='Show the ghosts on the structure. Default: %default')
    parser.add_option('-t','--titgroup',dest='titgroup',action='store',default=':0035:GLU',
        help='Titratable group selection. Default: %default')
    parser.add_option('-a','--atom',dest='atom',action='store',default='N',
        help='Atom to examine. Default: %default')
    parser.add_option('--pdb',dest='pdbfile',action='store',default='2LZT_H.pdb',
        help='PDB file. Default: %default')
    parser.add_option('--useabsent',dest='useabsent',action='store_true',default=False,
        help='Include absent ghosts in the plots and RMSDs (only enabled for cubescanplot. Default: %default')
    parser.add_option('--use_questionable',dest='use_questionable',default=True,action='store_true',
        help='Use questionable experimental ghosts in the plots and statistics. Default: %default')
    #
    # Options for specifying files
    #
    parser.add_option('--loaddir',dest='loaddir',default=None,type='string',
        help='Directory to load multiple calculated ghosts from. Used for loading calculated ghosts from MD snapshots. Default: %default')
    parser.add_option('-g','--calcghosts',dest='calcghosts',default='calc_ghosts/Xrays/all_ghosts_2LZT_H.pdb.pickle',help='File to load all calculated ghosts from. Default: %default')
    parser.add_option('--avgfile',dest='avgfile',default=None,help='File that averaged ghosts will be saved to. Default: %default')

    parser.add_option('--calcghostdir',dest='calcghostdir',default='calc_ghosts',help='Directory where all predictions reside. Deafult: %default')
    #
    # Which plots?
    #
    parser.add_option('--restraints',dest='restraints',action='store_true',default=False,
        help='Plot the number of restraints satisfied. Default: %default')
    parser.add_option('--stats',dest='stats',action='store_true',default=False,
        help='Print stats on experimental data. Default: %default')
    parser.add_option('-s','--cubescanplot',dest='cubescanfile',type='string',default='',
        help='If set to a cubescan output filename, then plot cube scores and make 3D epsmap. Default: "%default"')
    parser.add_option('--plotdir',dest='plotdir',type='string',default='plots',help='Dir used for saving all plots. Default; %default')
    parser.add_option('--dofull',dest='dofull',action='store_true',help='Make plot where predictions for all calculation types are loaded. Default: %default', default=False)

    parser.add_option('--filter_confchange',dest='filter_confchange',default=False,action='store_true',help='Filter by conformational change as determined by absolute chemical shift. Default: %default')
    parser.add_option('--filter_CPMG',dest='CPMG_filter',default=False,action='store_true',help='Filter residues where the CPMG signal is at least 50% of the ghost. Default: %default')
    
    parser.add_option('--no_calcerror',dest='no_calcerror',action='store_true',default=True,help='Disable the use of errors on calculated values. Default: %default')

    #
    # --------------
    #
    # Execution mode
    #
    parser.add_option('--anatype',dest='anatype',default='eps',help='Analysis mode. choices are: eps (compare calculated and experimental ghosts and find optimal eps), energies (convert experimental ghosts to energies). Default: %default')
    
    (options, args,) = parser.parse_args()
    #
    # Catch errors
    #
    if options.calcghosts and options.loaddir:
        options.calcghosts=False
    if options.avgfile and not options.loaddir:
        raise Exception('--avgfile has no meaning when --loaddir is not specified')
        
    #
    GA=ghost_analyzer(options)


    
