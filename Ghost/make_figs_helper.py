def error_size(options):
    """Unused code for looking at the size of errorrs"""
    #
    # Split into large and small errors
    #
    sxs=[]
    sys=[]
    serr=[]
    bxs=[]
    bys=[]
    berr=[]
    calc_bigerr=[]
    calc_small_err=[]
    for count in range(len(xs)):
        x=xs[count]
        y=ys[count]
        err=experrors[count]
        c_err=calcerrors[count]
        diff=abs(x-y)
        if diff<abs(err)+abs(c_err):
            diff=0.0
        else:
            diff=diff-abs(err)-abs(c_err)
        if diff>0.2:
            bxs.append(x)
            bys.append(y)
            berr.append(err)
            calc_bigerr.append(c_err)
        else:
            sxs.append(x)
            sys.append(y)
            serr.append(err)
            calc_smallerr.append(c_err)
    #
    # Plot
    #
    pylab.errorbar(sxs,sys,xerr=serr,fmt='go',label='Small errors')
    pylab.errorbar(bxs,bys,xerr=berr,fmt='ro',label='Large errors')
    pylab.plot(xs,xs,'y-')
    pylab.xlabel('Experimental chemical shift')
    pylab.ylabel('Calculated chemical shift')
    pylab.title('All titratable groups, N only: %s, eps = %5.1f, RMSD= %6.3f' %(method,besteps/10.0,rmsds[0][0]))
    pylab.savefig('%s.png' %method)
    pylab.show()
    pylab.clf()
    #
    # Reformat the residues info
    #
    longres=[]
    for tg,ress in residues:
        for residue in ress:
            longres.append([tg,residue])
    #
    # Find the worst predictions
    #
    nxs=[]
    nys=[]
    nerrs=[]
    import Yasara
    #Y=Yasara.yasara_handler()
    #obj=Y.load_mol('2LZT.pdb')
    #Y.yasara.run('HideAtom Sidechain')
    #Y.yasara.run('Color Green')
    print 'Worst predictions for %s' %method
    tgs=sorted(bigdata[method].keys())
    for count in range(len(xs)):
        x=xs[count]
        y=ys[count]
        experr=abs(experrors[count])
        tg,residue=longres[count]
        #
        color=tgs.index(tg)*30
        tgnum=int(tg.split(':')[1])
       
        #
        diff=abs(x-y)
        if diff<experr:
            diff=0.0
        else:
            diff=diff-experr
        #
        if diff>0.2:
            print '%6s residue: %6s, exp: %6.3f, calc: %6.3f, error: %6.3f' %(tg,residue,x,y,diff)
            #Y.yasara.run('ShowRes %s' %tgnum)
            #Y.yasara.ColorRes(tgnum,color)
            #Y.yasara.run('ColorAtom N Res %d, %d' %(int(residue[1:]),color))
        else:
            nxs.append(x)
            nys.append(y)
            nerrs.append(experr)
    import pylab
    pylab.errorbar(nxs,nys,xerr=nerrs,fmt='ro')
    pylab.plot(nxs,nxs,'g-')
    pylab.xlabel('Exp')
    pylab.ylabel('calc')
    rmsd=RMSD(nxs,nys,nerrs)
    pylab.title('Removing large errors for %s RMSD: %6.3f' %(method,rmsd))
    pylab.show()
    raw_input('is this ok?')
    return
    

#
# -----
#
class score_class:
    """Class for scoring calculated ghosts and compiling various statistics"""
    
    def __init__(self,options):
        self.options=options 
        self.exp_ghosts=getdata(options)
        del self.exp_ghosts[':0015:HIS']
        #
        # Instantiate the scoring function
        #
        self.S=dummy_epsmap()
        #
        # Get the residues that we should exclude
        #
        self.excludes=find_excludes(options,self.exp_ghosts)
        return
        
    #
    # -----
    #
        
    def score_ghosts(self,calc_ghosts):
        """Score all the calculated ghosts given"""
        big_satisfied={}
        atom=options.atom
        big_satisfied[atom]={}
        xs=[]
        ys=[]
        experrors=[]
        #
        for tg in sorted(self.exp_ghosts.keys()):
            big_satisfied[atom][tg]={}
            for residue in sorted(self.exp_ghosts[tg].keys()):
                #
                # Exclude this residue for this tg?
                #
                if self.excludes.has_key(tg):
                    if residue in self.excludes[tg]:
                        continue
                #
                # Get the experimental values and its error
                #
                if self.exp_ghosts[tg][residue].has_key(atom):
                    exp_value=self.exp_ghosts[tg][residue][atom]
                    exp_error=errors[atom]
                    if self.exp_ghosts[tg][residue].has_key(atom+'_error'):
                        exp_error=self.exp_ghosts[tg][residue][atom+'_error']
                    #
                    # Get calculated value and score it
                    #
                    calc_value=calc_ghosts[tg][residue][atom]
                    error,satisfied,abs_sat,tot_restraints,tot_abs,real_error=self.S.get_error_sub(exp_value,calc_value,exp_error,atom)
                    #
                    if exp_value!='absent':
                        exp_error=0.05
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
                    else:
                        if self.options.useabsent:
                            exp_value=0.0
                            exp_error=errors[options.atom]
                        else:
                            continue
                    #
                    xs.append(float(exp_value))
                    ys.append(calc_value)
                    experrors.append(exp_error)
                    big_satisfied[atom][tg][residue]={'error':error,'real_error':real_error,'calc_value':calc_value,'exp_value':exp_value,'exp_error':exp_error}
        return xs,ys,experrors,big_satisfied

#
# -----
#

def cubescanplot(options):
    """Read a cubescan output file and show the results"""
    fd=open(options.cubescanfile)
    import pickle
    data=pickle.load(fd)
    fd.close()
    #
    cubegrid=data[0]
    scanresults=data[1]
    #
    # Instantiate scoring class
    #
    SC=score_class(options)
    #
    # Score the ghosts from each cubescan
    #
    scores=[]
    for cube in sorted(scanresults.keys()):
        calc_ghosts=scanresults[cube]
        xs,ys,experrors,satisfied=SC.score_ghosts(calc_ghosts)
        rmsd=RMSD(xs,ys,experrors)
        scores.append([rmsd,cube])
        #import pylab
        #pylab.errorbar(xs,ys,xerr=experrors,fmt='ro')
        #pylab.plot(xs,xs,'g-')
        #pylab.xlabel('Experimental dCS')
        #pylab.ylabel('Calculated dCS')
        #pylab.title('Cubescan of cube %4d, atom: %s, RMSD: %5.3f' %(cube,options.atom,rmsd))
        #pylab.savefig('Cubescan_%d.png' %(cube))
        #pylab.clf()
    rmsds=[]
    scores.sort()
    import Protool
    P=Protool.structureIO()
    P.readpdb('2LZT_H.pdb')
    count=0
    for rmsd,cube in scores[:25]:
        print '%4d, rmsd: %5.2f' %(cube,rmsd)
        center=cubegrid[cube]['coord']
        P.add_atom('X:%4d:CS' %(count+1000),
                 atomnumber=0,atomname='CS',
                 chainid='X',residuename='CSS',residuenumber='999',
                 xcoord=center[0],ycoord=center[1],zcoord=center[2],update=1,BFACTOR=None,OCCUPANCY=None,CHARGE=None,RADIUS=None,tag=None,accept_duplicate=False)
        count=count+1
        rmsds.append(rmsd)
    P.writepdb('cubescan.pdb')
    import pylab
    pylab.hist(rmsds)
    pylab.savefig('Cubescanhist.png')
    return


def heatmap():
    """|Make a heatmap - not sure what this function does..."""
    data=getdata()
    TGs=sorted(data.keys())
    for titgroup in TGs:
        thisN=[]
        thisH=[]
        for residue in sorted(data[titgroup].keys()):
            Hval=0.0
            Nval=0.0
            for atom in ['H','N']:
                if data[titgroup][residue].has_key(atom):
                    x=data[titgroup][residue][atom]
                    #
                    # -----
                    #
                    if x=='absent':
                        print 'abs'
                        if atom=='H':
                            Hval=-20
                        else:
                            Nval=-20
                    else:
                        if x[0]=='q':
                            x=x[1:]
                        x=x.split(';')[0]
                        if x[0] in ['<','>']:
                            x=x[1:]
                        x=float(x)
                        if atom=='H':
                            Hval=20
                        else:
                            Nval=20
                    #
                    #
            thisN.append(Nval)
            thisH.append(Hval)
        for y in range(10):
            Ndata.append(thisN)
            Hdata.append(thisH)
    #
    #
    import TwoDplots
    blue_red1 = TwoDplots.LinearSegmentedColormap('BlueRed1', cdict1)
    #
    plt.imshow(Ndata,interpolation='nearest',cmap=blue_red1)
    cbar=plt.colorbar()
    cbar.set_label('score')
    #
    plt.yticks(range(5,105,10),TGs)
    plt.ylim([0,100])
    plt.show()
    return
    
