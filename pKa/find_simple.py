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

import Protool, Numeric
import dislin_driver


def plotit(xs,ys,title,legends):
    #
    # Do it 
    #
    num_points=0
    for x in xs:
        num_points=num_points+len(x)
    print num_points,'data points'
    #x=[1,2,3,4,5,6,7,8,9]
    #y=[2,4,6,8,10,12,14,16,18]
    mat_fix=[]
    vec_fix=[]
    for x in xs:
        for point in x:
            mat_fix.append([point,1.0])
    for y in ys:
        for point in y:
            vec_fix.append(point)
    import LinearAlgebra
    sols,rsq,rank,junk=LinearAlgebra.linear_least_squares(Numeric.array(mat_fix),
                                                          Numeric.array(vec_fix))
    slope=sols[0]
    intercept=sols[1]
    print rsq
    rsq=float(rsq[0])
    print 'Slope: %.2f, Intercept: %.2f, R^2: %.2f' %(slope,intercept,rsq)
    file=dislin_driver.graf_mult3(xs,ys,title,'Simple E','PBE_ene',legends)
    return

#
# -----
#

def matrix():
    import os
    dirs=os.listdir('data')
    xs=[]
    ys=[]
    for dir in dirs:
        print 'Processing %s' %dir
        #
        # find the PDB file
        #
        realdir=os.path.join(os.getcwd(),'data',dir)
        files=os.listdir(realdir)
        for file in files:
            realfile=os.path.join(realdir,file)
            if file[-4:]=='.pdb':
                import pKa.pKaTool.pKaIO
                X=pKa.pKaTool.pKaIO.pKaIO(realfile)
                X.assess_status()
                if X.calculation_completed==1:
                    #
                    # Hurra, the calc is complete. Load the matrix
                    #
                    PBEmatrix=X.read_matrix()
                    #
                    # Now calculate the same matrix with Protool
                    #
                    P=Protool.structureIO()
                    P.readpdb(realfile)
                    P.get_titratable_groups()
                    dist_matrix=P.Calculate_matrix(8)
                    #
                    # Plot it
                    #
                    x=[]
                    y=[]
                    for group1 in PBEmatrix.keys():
                        for group2 in PBEmatrix.keys():
                            PBE_ene=PBEmatrix[group1][group2][0]
                            try:
                                new_ene=dist_matrix[group1][group2]
                            except:
                                continue
                            #
                            # Load the values, distances in x, PBE_ene in y
                            #
                            if new_ene and PBE_ene:
                                x.append(abs(new_ene))
                                y.append(abs(PBE_ene))
                    #
                    # Append these result to the big arrays
                    #
                    ys.append(y)
                    xs.append(x)
    plotit(xs,ys,'Matrix',dirs)
    return

          
def desolvation():
    """Find correlation between desolvation and acc"""
    import os
    dirs=os.listdir('data')
    xs=[]
    ys=[]
    for dir in dirs:
        print 'Processing %s' %dir
        #
        # find the PDB file
        #
        realdir=os.path.join(os.getcwd(),'data',dir)
        files=os.listdir(realdir)
        for file in files:
            realfile=os.path.join(realdir,file)
            if file[-4:]=='.pdb':
                import pKa.pKaTool.pKaIO
                X=pKa.pKaTool.pKaIO.pKaIO(realfile)
                X.assess_status()
                if X.calculation_completed==1:
                    #
                    # Hurra, the calc is complete. Load the desolvation energies
                    #
                    PBEdesolv=X.read_desolv()
                    #
                    # Now calculate the desolvation energies with Protool
                    #
                    P=Protool.structureIO()
                    P.readpdb(realfile)
                    P.get_titratable_groups()
                    P.calculate_desolvation()
                    #
                    #
                    x=[]
                    y=[]
                    for residue in PBEdesolv.keys():
                        if P.desolv.has_key(residue):
                            x.append(P.desolv[residue])
                            y.append(PBEdesolv[residue])
                    xs.append(x)
                    ys.append(y)
    plotit(xs,ys,'Desolvation',dirs)
    return

def NN_desolvation():
    import os
    dirs=os.listdir('data')
    input=[]
    expect=[]
    groups={'ASP':1,'GLU':2,'ARG':3,'LYS':4,'HIS':5,'NTERM':6,'CTERM':7,'CYS':8,'TYR':9}
    #
    xs=[]
    ys=[]
    for dir in dirs:
        print 'Processing %s' %dir
        #
        # find the PDB file
        #
        realdir=os.path.join(os.getcwd(),'data',dir)
        files=os.listdir(realdir)
        for file in files:
            realfile=os.path.join(realdir,file)
            if file[-4:]=='.pdb':
                import pKaTool.pKaIO
                X=pKaTool.pKaIO.pKaIO(realfile)
                X.assess_status()
                if X.calculation_completed==1:
                    #
                    # Hurra, the calc is complete. Load the desolvation energies
                    #
                    PBEdesolv=X.read_desolv()
                    #
                    # Now calculate the desolvation energies with Protool
                    #
                    P=Protool.structureIO()
                    P.readpdb(realfile)
                    P.get_titratable_groups()
                    if not getattr(P,'atom_close',None):
                        P.get_atoms_around()
                    #
                    #
                    titgrps=P.titratable_groups
                    residues=titgrps.keys()
                    residues.sort()
                    for residue in residues:
                        for group in titgrps[residue]:
                            name1=residue+':'+P.resname(residue)+group['name']
                            if PBEdesolv.has_key(name1):
                                this_input=[0,0,0,0,0,0,0,0,0,0]
                                if group['name']!='':
                                    this_input[groups[group['name'][1:]]]=1
                                else:
                                    this_input[groups[P.resname(residue)]]=1
                            
                                #
                                # Add the info on how many atoms are close
                                #
                                this_input[0]=min(P.atom_close[name1]/20.0,1.0)
                                
                                this_expect=[min(max(15.0+PBEdesolv[name1],0.0)/30.0,1.0)]
                                print name1,this_input,this_expect
                                expect.append(this_expect)
                                input.append(this_input)
    #
    # Train the network
    #
    nettop=(10,5,1)
    import sys
    sys.path.append('/home/nielsen/lib/python/chresten/source2/methods')
    import NN
    network=NN.NN(nettop)
    print network.print_network()
    step=[]
    error=[]
    pred_value=[]
    real_value=[]
    for x in range(0,3000):
        res = network.train(input,expect)
    #print net for every 200 training rounds
        if(x%20 == 0 ):
            
            counter=0
            cum_error=0.0
            for inp in input:
                NNresult=network.test(inp)[0]
                real_result=expect[counter][0]
                counter=counter+1
                cum_error=cum_error+abs(NNresult-real_result)
            error.append(cum_error)
            step.append(x)
            print 'Step %d, cumulative error: %5.3f' %(x,cum_error)

    #make plot
    import dislin_driver
    file=dislin_driver.graf(step,error)
    import os
    os.system('eog %s' %file)
    return

#
# ----
#

def background():
    """Find correlation between background interaction energy and number of hbonds"""
    import os
    dirs=os.listdir('data')
    xs=[]
    ys=[]
    count=0
    #dirs=['1bli']
    for dir in dirs:
        print 'Processing %s' %dir
        #
        # find the PDB file
        #
        realdir=os.path.join(os.getcwd(),'data',dir)
        files=os.listdir(realdir)
        for file in files:
            realfile=os.path.join(realdir,file)
            if file[-4:]=='.pdb':
                import pKaTool.pKaIO
                X=pKaTool.pKaIO.pKaIO(realfile)
                X.assess_status()
                if X.calculation_completed==1:
                    #
                    # Hurra, the calc is complete. Load the desolvation energies
                    #
                    PBEbackground=X.read_backgr()
                    #
                    # Now calculate the desolvation energies with Protool
                    #
                    P=Protool.structureIO()
                    P.readpdb(realfile)
                    P.get_titratable_groups()
                    P.calculate_background()
                    #
                    #
                    x=[]
                    y=[]
                    residues= PBEbackground.keys()
                    residues.sort()
                    for residue in residues:
                        if P.background.has_key(residue) and PBEbackground[residue]<19.0:
                            x.append(P.background[residue])
                            #x.append(count)
                            count=count+1
                            y.append(PBEbackground[residue])
                            print '%12s %5.2f %5.2f' %(residue,P.background[residue],PBEbackground[residue])
                    xs.append(x)
                    ys.append(y)
    plotit(xs,ys,'Background',dirs)
    return
    

if __name__=="__main__":
    #matrix()
    NN_desolvation()
    #background()
