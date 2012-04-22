#!/usr/bin/env python
R=8.3145

def resact(T,time,Ea,A,startconc):
    """Calculate the residual activity at a given temp after a given time"""
    import math
    k=A*math.exp(-Ea/(R*T))
    Conc=startconc*math.exp(-k*time)
    return Conc-startconc
    
    #
    # -----
    #
    
def T50(Ea,A,t):
    """Calculate the T50"""
    import math
    tmp=math.log(2.0)/(t*A)
    result=-Ea/(math.log(tmp)*R)
    return result
    
    #
    # ------
    #
    
def eq_unfold(dH,dS,dCp,T1,T2):
    """Calculate the fraction of protein in the unfolded state"""
    import math
    dG=dH+dCp*(T2-T1)+dS+dCp*math.log(T2/T1)
    K=math.exp(-dG/(R*T2))
    return 1/(1+K)

    #
    # ------
    #

def unfold():
    t=25
    Ea=34000
    A=20000
    import numpy as np
    xs=[]
    ys=[]
    zs=[]
    #for T in np.arange(270,370,5):
    #    xs.append(T-273.15)
    #    ys.append(resact(T,t,Ea,A))
    wt=T50(Ea,A,t)
    print 'WT T50 is %5.2f' %(wt-273.15)
    for ddG in np.arange(-10000,10000,500):
        xs.append(ddG/1000.0)
        new_Ea=Ea+ddG
        Tm=T50(new_Ea,A,t)
        ys.append(Tm-wt)
        # 
        # Equilibrium unfolding
        #
        Tm=T50(Ea,A-ddG,t)
        zs.append(Tm-wt)
    
    import pylab
    pylab.plot(xs,ys,'ro',label='exp decay')
    pylab.plot(xs,zs,'bo',label='Eq. unfold')
    pylab.xlabel('ddG (kJ/mol)')
    pylab.ylabel('dTm (K)')
    pylab.legend()
    pylab.show()
    return
    
    #
    # ------
    #
    
def kinetic(options):
    """Simulate kinetic stability"""
    if options.reactions==[]:
        # Use the standard setup

        options.reactions=[1E10,100000.0,10000.0,1E10]
        options.reactions=options.kin_sites*options.reactions

    
    if options.experiment=='residual':
        #
        #
        import numpy as np
        xs=[]
        ys=[]
        Tm=None
        last_temp=False
        for temp in np.arange(options.starttemp,options.endtemp,1.0):
            conc=100.0
            for time in np.arange(0.0,options.endtime,options.dt):
                #
                # Calculate residual activity at the simulation time
                #
                dconc=[]
                for reaction in options.reactions:
                    #sp=reaction.split(',')
                    A_uf=float(reaction[0])
                    Ea_uf=float(reaction[1])
                    A_fold=float(reaction[2])
                    Ea_fold=float(reaction[3])
                    change=resact(temp,options.endtime,Ea_uf,A_uf,conc)
                    dconc.append(change)
                #
                # Apply the concentration changes
                #
                for change in dconc:
                    conc=conc+change
                # Find the Tm
                if last_temp:
                    if conc<50.0 and not Tm:
                        Tm=temp
                last_temp=temp
            ys.append(conc)
            xs.append(temp-273.15)
            #
        #
        # Graph it
        #
        if options.plot:
            import pylab
            pylab.plot(xs,ys,'ro-',label='T50: %5.1f C' %(Tm-273.15))
            pylab.xlabel('Temperature (C)')
            pylab.ylabel('Residual activity (%)')
            pylab.ylim([0.0,100.0])
            pylab.legend()
            pylab.show()
    return xs,ys,Tm-273.15
            
                    

    
    #
    # ------
    #
    
def main(options):
    """Parse options and find out what to do"""
    if options.graph=='numsites':
        tms=[]
        sites=[]
        import pylab
        for sitenum in range(1,10):
            options.kin_sites=sitenum
            options.plot=False
            xs,ys,Tm=kinetic(options)
            pylab.plot(xs,ys,'o-',label='%d sites. Tm: %5.1f' %(sitenum,Tm),linewidth=2)
            sites.append(sitenum)
            tms.append(Tm)
        pylab.xlabel('Temperature (C)')
        pylab.ylabel('Residual activity (%)')
        pylab.ylim([0.0,100.0])
        pylab.legend()
        pylab.title('Number of protease sites')
        pylab.show()
        #
        pylab.clf()
        pylab.plot(sites,tms,'ro-',linewidth=2)
        pylab.xlabel('Number of sites')
        pylab.ylabel('T50')
        pylab.show()
        return
    elif options.graph=='twosites':
        # fast / slow site, + stabilization
        import pylab
        std_reactions=[[1E10,1E5,5E4,1E10],
                       [1E10,9E4,5E4,1E10]]
        for stab1,stab2 in [[0,0],[options.stab,0],[0,options.stab],[options.stab,options.stab]]:
            import string,copy
            options.reactions=copy.deepcopy(std_reactions)
            options.reactions[0][1]=options.reactions[0][1]+stab1*1000.0
            options.reactions[1][1]=options.reactions[1][1]+stab2*1000.0
            
            #options.reactions[0]=string.join(options.reactions[0],',')
            #options.reactions[1]=string.join(options.reactions[1],',')
            #
            options.plot=False
            #print options.reactions
            #print std_reactions,'R'
            xs,ys,Tm=kinetic(options)
            pylab.plot(xs,ys,'o-',label='%5.1f,%5.1f, T50: %5.1f' %(stab1,stab2,Tm))
        pylab.legend()
        
        pylab.suptitle('Ea1: %3.2e kJ/mol, Ea2 %3.2e kJ/mol stab: %5.1f kJ/mol' %(std_reactions[0][1],std_reactions[1][1],options.stab))
        pylab.title('Two sites')
        pylab.show()
    
    
    if options.unfoldtype=='kinetic':
        kinetic(options)
    elif options.unfoldtype=='equilibrium':
        equilibrium(options)
    else:
        raise Exception('unknown unfolding type: %s' %options.unfoldtype)
    #
                        
            
    
    return
    
    
if __name__=='__main__':
    print
    print 'Simulation of different unfolding/stability scenarios'
    print
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
    parser.add_option('--unfoldtype',dest='unfoldtype',type='string',default='kinetic',help='Type of unfolding reaction. Options are: kinetic, equilibrium')
    parser.add_option('--experiment',dest='experiment',default='residual',help='Type of graph to show. Options: T50, residual activity (residual), ddG. Default: %default')
    parser.add_option('--reaction',dest='reactions',action='append',default=[],help='Specify an unfolding reaction. You can specify more than one using format A_unfold,Ea_unfold,A_fold,Ea_fold (in kJ/mol).')
    parser.add_option('-s','--sites',dest='kin_sites',default=1,type='int',help='Number of kinetic unfolding sites. Default: %default')
    #
    parser.add_option('--endtime',dest='endtime',type='float',default=100.0,help='End time of simulation (in seconds). Default: %default')
    parser.add_option('--dt',dest='dt',type='float',default=0.1,help='Time step in seconds. Default: %default')
    
    parser.add_option('--starttemp',dest='starttemp',default='20.0',type='float',help='Starting temperature (C). Default: %default')
    parser.add_option('--endtemp',dest='endtemp',default=100.0,type='float',help='Ending temperature (C). Default: %default')
    
    parser.add_option('-g','--graph',dest='graph',default=None)
    parser.add_option('--stab',dest='stab',default=10.0,type='float',help='Stabilization energy in kJ/mol. Default: %default')
    
    
    
    (options, args) = parser.parse_args()
    options.plot=True
    #
    # Convert temperatures to Kelvin
    #
    options.starttemp=options.starttemp+273.15
    options.endtemp=options.endtemp+273.15
    # Split reactions
    reacts=[]
    for reac in options.reactions:
        reacts.append(reac.split(','))
    options.reactions=reacts
    # Call main function
    main(options)