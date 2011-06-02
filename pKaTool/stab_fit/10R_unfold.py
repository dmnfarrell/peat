#!/usr/bin/env python
R=8.3145

def resact(T,time,Ea,A):
    """Calculate the residual activity at a given temp after a given time"""
    import math
    k=A*math.exp(-Ea/(R*T))
    #k=math.pow(k,1.5)
    Conc=100.0*math.exp(-k*time)
    return Conc
    
def T50(Ea,A,t):
    """Calculate the T50"""
    import math
    tmp=math.log(2.0)/(t*A)
    result=-Ea/(math.log(tmp)*R)
    return result
    
def eq_unfold(dH,dS,dCp,T1,T2):
    """Calculate the fraction of protein in the unfolded state"""
    import math
    dG=dH+dCp*(T2-T1)+dS+dCp*math.log(T2/T1)
    K=math.exp(-dG/(R*T2))
    return 1/(1+K)

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
    
    
    
if __name__=='__main__':
    unfold()