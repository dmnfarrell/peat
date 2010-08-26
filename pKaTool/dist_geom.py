#!/usr/bin/env python
#
# pKaTool - analysis of systems of titratable groups
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

from Tkinter import *
import numpy
import numpy.linalg
import sys
import math

def average(l):
    """Calculate the average value of the elements of a list"""
    sum=0.0
    for i in l:
        sum=sum+i
    return sum/float(len(l))

def avg_sq(l):
    sum=0.0
    for i in l:
        sum=sum+i
    return sum/float(len(l)*len(l))

#
# ----
#

class distance_optimisation:

    def __init__(self,distance_matrix,titration_curves):
        """Set up the control window and the PyOPenGL window"""
        #
        # Variables
        #
        self.P=None # Var for Protool instance
        #
        # Open window
        #
        self.Dcontrol=Toplevel()
        self.Dcontrol.geometry('+500+200')
        self.Dcontrol.title('Distance geometry control')
        #
        # Text window
        #
        self.coord_text=Text(self.Dcontrol,background='white',
                             foreground='black',
                             state=NORMAL,
                             exportselection=1,
                             width=60,
                             height=20)#,
        #yscrollcommand= detail_yscrollbar.set,
        #                  xscrollcommand=detail_xscrollbar.set)
        self.coord_text.grid(row=0,column=0,columnspan=3)
        self.specific_residue = StringVar()
        self.specific_residue.set('Chose residue')
        
        #
        self.eps=DoubleVar()
        self.eps.set(20)
        self.eps_scale=Scale(self.Dcontrol,from_=1,to=80,resolution=1,
                             orient='horizontal',relief='ridge',
                             variable=self.eps,
                             label='Effective dielectric constant')
        self.eps_scale.grid(row=8,column=0,columnspan=3)
        #
        Button(self.Dcontrol,text='Do new DG',command=self.update_all).grid(row=9,column=0)
        #
        Button(self.Dcontrol,text='Load PDB',command=self.load_PDB).grid(row=10,column=0)
        #
        Button(self.Dcontrol,text='Find match in PDB',command=self.find_match).grid(row=10,column=1)
        #
        Button(self.Dcontrol,text='Load pKa matrix file',command=self.load_matrix).grid(row=11,column=0)
        #
        Button(self.Dcontrol,text='Find match in matrix',command=self.find_matrix_match).grid(row=11,column=1)
        #
        Button(self.Dcontrol,text='Find match with specific residue',command=self.find_match_ineraction_with_specific_residue).grid(row=12,column=0)
        #
        Entry(self.Dcontrol,text='Residue',textvariable=self.specific_residue).grid(row=12,column=1)

#
#        mb = Menubutton(self.Dcontrol, textvariable=self.specific_residue)
#        menu = Menu(mb,tearoff=0)
#        mb['menu'] = menu
#        for i in range(len(titration_curves.keys())):
#            menu.add_radiobutton(label=titration_curves.keys()[i],variable=self.specific_residue,  value=titration_curves.keys()[i], indicatoron=1)
#        mb.grid(row=12,column=1)


        #
        #
        # Window for PDB search results
        #
        self.PDB_text=Text(self.Dcontrol,background='white',
                           foreground='black',
                           state=NORMAL,
                           exportselection=1,
                           width=60,
                           height=20)#,
        #yscrollcommand= detail_yscrollbar.set,
        #                  xscrollcommand=detail_xscrollbar.set)
        self.PDB_text.grid(row=13,column=0,columnspan=3)
        #
        # Do the first opt
        #
        self.distance_matrix=distance_matrix
        self.update_all()
        #
        # Start the OpenGL
        # 
        #X=OGL(self.coords,self)
        #glutMainLoop()

        
        return

    #
    # ----
    #

    def update_all(self):
        """Do new distance geometry"""
        self.coord_text.delete(1.0,END)
        #
        # Set the initial eps
        #
        matrix=self.set_eps(self.distance_matrix.copy(),self.eps.get())
        
        self.DG=distance_geometry(matrix)
        new_dist_matrix=self.DG.do_triangle_smoothing()
        self.coords=self.DG.construct_metric_matrix()
        # 
        # Insert the coordinates
        #
        groups=self.coords.keys()
        groups.sort()
        self.coord_text.insert(END,'%10s %7s %7s %7s\n' %('Group','X','Y','Z'))
        for group in groups:
            self.coord_text.insert(END,'%10s %7.2f %7.2f %7.2f\n' %(group,
                                                                self.coords[group]['X'],
                                                                self.coords[group]['Y'],
                                                                self.coords[group]['Z']))
        #
        # Insert distances
        #
        for g1 in range(len(groups)):
            for g2 in range(g1,len(groups)):
                group=groups[g1]
                group2=groups[g2]
                if g1==g2:
                    continue
                self.coord_text.insert(END,'%10s %10s %7.2f\n' %(group,group2,new_dist_matrix[group][group2]))
                                                             

    #
    # ----
    #

    def set_eps(self,matrix,eps):
        """Set a new eps for all interactions"""
        new_matrix={}
        for A in matrix.keys():
            new_matrix[A]={}
            for B in matrix[A].keys():
                if matrix[A][B]==0:
                    new_matrix[A][B]=0
                elif matrix[A][B]<999.9:
                    new_matrix[A][B]=matrix[A][B]/eps
                else:
                    # distances of 1000 are infinite
                    new_matrix[A][B]=matrix[A][B]
        return new_matrix

    #
    # ----
    #

    def load_PDB(self):
        """Load a PDB file"""
        import tkFileDialog, os
        filename=tkFileDialog.askopenfilename(defaultextension='.pdb',
                                              initialdir=os.getcwd(),
                                              parent=self.Dcontrol,
                                              filetypes=[("PDB file","*.pdb"),
                                                         ("All files","*.*")])
        if filename:
            import Protool
            self.P=Protool.structureIO()
            self.P.readpdb(filename)
        return

    #
    # -----
    #

    def find_match(self):
        """find match for the coordinates in teh PDB"""
        if not self.P:
            return
        #
        # Number of groups from distance geometry
        #
        groups=self.coords.keys()        
        num_groups=len(groups)
        #
        # Get the titratable groups in the protein
        #
        self.P.get_titratable_groups()
        tit_grps=self.P.titratable_groups
        tres=tit_grps.keys()
        tres.sort()
        #
        #
        #
        grp_centers={}
        for res in tres:
            name=self.P.resname(res)
            if not name in ['ARG','LYS','GLU','ASP','HIS']:
                continue
            for group in tit_grps[res]:
                atoms=[]
                for atom in group['atoms']:
                    atomname=res+':'+atom
                    if self.P.atoms.has_key(atomname):
                        atoms.append(self.P.GetPosition(atomname))
                #
                # Get avg coordinate
                #
                avg=numpy.zeros([3])
                for c in atoms:
                    avg=avg+c
                avg=avg/float(len(atoms))
                #
                # Construct group name
                #
                gname=res+':'+group['name']
                grp_centers[gname]=avg
        #
        # Construct all permutations and find best match
        #
        # Number of groups to fit is num_groups
        #
        prot_groups=grp_centers.keys()
        prot_groups.sort()
        print 'Protein groups'
        print prot_groups
        self.permutations=self.construct_permutations(num_groups,prot_groups)
        #
        # Load the coordinates from the distance geometry calc
        #
        fit_coords=numpy.zeros([num_groups,3])
        count=0
        for group in groups:
            fit_coords[count][0]=self.coords[group]['X']
            fit_coords[count][1]=self.coords[group]['Y']
            fit_coords[count][2]=self.coords[group]['Z']
            count=count+1
        #
        # Search for good matches in all group_centers
        #
        print 'Number of permutations to search: %d' %(len(self.permutations.keys()))
        fit_results={}
        big_count=0
        import sys
        for perm in self.permutations.keys():
            ref_coords=numpy.zeros([num_groups,3])
            count=0
            for group_name in eval(perm):
                ref_coords[count]=grp_centers[group_name]
                count=count+1
            #
            # Do the fit
            #
            rot_mat,trans_vec,rms,rtv,ftv=self.P.superpose(ref_coords,fit_coords)
            self.permutations[perm]=rms
            big_count=big_count+1
            streng='\b\b\b\b\b\b\b\b\b\b %8d' %big_count
            print streng,
            sys.stdout.flush()
        #
        # Get the best ten solutions
        #
        sols=[]
        rmsds=self.permutations.values()
        rmsds.sort()
        print rmsds[:10]
        for rms in rmsds:
            if rms:
                for perm in self.permutations.keys():
                    if not perm in sols:
                        if self.permutations[perm]==rms:
                            sols.append(perm)
            if len(sols)==10:
                break
        #
        # Print them
        #
        self.PDB_text.delete(1.0,END)
        for sol in sols:
            self.PDB_text.insert(END,'%30s %7.2f\n' %(sol,self.permutations[sol])) 

        x=[[':0035:',':0078:',':0172:'],[':0035:',':0172:',':0078:'],
           [':0078:',':0035:',':0172:'],[':0078:',':0172:',':0035:'],
           [':0172:',':0035:',':0078:'],[':0172:',':0078:',':0035:']]
        for p in x:
            self.PDB_text.insert(END,'\n%30s %7.2f\n' %(str(p),self.permutations[str(p)]))
            
        return

    #
    # -----
    #

    def load_matrix(self):
        """Load a matrix file"""
        import tkFileDialog, os
        filename=tkFileDialog.askopenfilename(defaultextension='.MATRIX.DAT',
                                              initialdir=os.getcwd(),
                                              parent=self.Dcontrol,
                                              filetypes=[("pKa MATRIX file","*.MATRIX.DAT"),
                                                         ("All files","*.*")])
        if filename:
            import pKaIO
            self.M=pKaIO.pKaIO()
            self.pkamatrix=self.M.read_matrix(filename)

            for res1 in self.pkamatrix.keys():
                for res2 in self.pkamatrix.keys():
                    print res1,res2,self.pkamatrix[res1][res2][0]-self.pkamatrix[res1][res2][1]-self.pkamatrix[res1][res2][2]+self.pkamatrix[res1][res2][3]
            
        return

    #
    # ----
    #

    def find_matrix_match(self):
        """Find a match in the energies in pkamatrix"""
        #
        # Complete distances
        #
        distmat=self.distance_matrix.copy()
        for group in distmat.keys():
            distmat[group][group]=0.0
        #
        # Number of groups in distance matrix
        #
        groups=distmat.keys()
        num_groups=len(groups)
        #
        # Make it into a Numeric matrix with energies
        #
        import math
        EM=numpy.zeros([num_groups,num_groups])
        for x in range(num_groups):
            for y in range(num_groups):
                if distmat[groups[x]][groups[y]]>0.0:
                    EM[x][y]=(243.4*math.log(10))/distmat[groups[x]][groups[y]]
                else:
                    EM[x][y]=0.0
                 #distance=243.3*math.log(10)/(eps*E)
        self.PDB_text.insert(END,'REF: %30s %7.2f\n' %(groups,0.0))

        EMT = EM.transpose() # not really necessary as EM is symmetric
        print 'Energy matrix',EM
        print 'Transposed energy matrix',EMT
        #
        # Now get all the groups in the protein
        #
        tmp_titgrps=self.pkamatrix.keys()
        tmp_titgrps.sort()
        titgrps=[]
        for grp in tmp_titgrps:
            if not 'TYR' in grp:
                titgrps.append(grp)
        print titgrps
        permutations=self.construct_permutations(num_groups,titgrps)
        #
        # Calculate the differences between the matrices
        #
        diffs=[]
        for perm in permutations.keys():
 #           print perm
            this_perm=eval(perm)
            this_matrix=numpy.zeros([len(this_perm),len(this_perm)])
            count1=0
            for res1 in this_perm:
                count2=0
                for res2 in this_perm:
                    if res1==res2:
                        this_matrix[count1][count2]=0.0
                    else:
                        # takes the average of the two interaction energies
                        this_matrix[count1][count2]=((self.pkamatrix[res1][res2][0]-self.pkamatrix[res1][res2][1]-self.pkamatrix[res1][res2][2]+self.pkamatrix[res1][res2][3])
                                                     +(self.pkamatrix[res2][res1][0]-self.pkamatrix[res2][res1][1]-self.pkamatrix[res2][res1][2]+self.pkamatrix[res2][res1][3]))/2

# this_matrix[count1][count2]=self.pkamatrix[res1][res2][0]-self.pkamatrix[res1][res2][1]-self.pkamatrix[res1][res2][2]+self.pkamatrix[res1][res2][3]
# this_matrix[count1][count2]=min((self.pkamatrix[res1][res2][0]-self.pkamatrix[res1][res2][1]-self.pkamatrix[res1][res2][2]+self.pkamatrix[res1][res2][3]),(self.pkamatrix[res2][res1][0]-self.pkamatrix[res2][res1][1]-self.pkamatrix[res2][res1][2]+self.pkamatrix[res2][res1][3]))



                    count2=count2+1
                count1=count1+1
            #
            # Get the difference
            #
            diff=EM-this_matrix
            sum=0.0
            for x in diff:
                for y in x:
                    sum=sum+abs(y)
            #
            # Calculate Frobenius inner product
            #
            product = numpy.dot(EMT,this_matrix)


#            normThis_matrix = math.sqrt(numpy.trace(numpy.dot(this_matrix.transpose(),this_matrix)))
                                      
            Frobenius  = math.sqrt(numpy.trace(numpy.dot(diff.transpose(),diff)))
            #Frobenius = numpy.trace(product)

            diffs.append([Frobenius,sum,perm])
        #
        # Sort and report
        #
        diffs.sort()
        #diffs.reverse()

        for fro,sum,sol in diffs:#[:10]:
            self.PDB_text.insert(END,'%30s %7.2f %7.2f\n' %(sol,sum,fro))
        self.PDB_text.insert(END,'Number of permutations: %d' %(len(permutations.keys())))
        return
    #
    # ----
    #

    def find_match_ineraction_with_specific_residue(self):

        if not len(self.distance_matrix.keys()) == 2:
            print 'Method can only be used when two titrating groups are fitted'
            return

        groups = self.distance_matrix.keys()
        print 'groups ',groups

        this_energy = (243.4*math.log(10))/self.distance_matrix[groups[0]][groups[1]]
        this_res = self.specific_residue.get()

        try:
            print 'all ',self.pkamatrix.keys()
            all_groups = self.pkamatrix.keys()
        except:
            return
        
        if not this_res in all_groups:
            print 'residue',this_res,' not found'
            return

        diffs = []
        for g in all_groups:
            E1 = self.pkamatrix[g][this_res][0]-self.pkamatrix[g][this_res][1]-self.pkamatrix[g][this_res][2]+self.pkamatrix[g][this_res][3]
            E2 = self.pkamatrix[this_res][g][0]-self.pkamatrix[this_res][g][1]-self.pkamatrix[this_res][g][2]+self.pkamatrix[this_res][g][3]
            E=(E1+E2)/2

            print 'E',E,'E1',E1,'E2',E2, 'this_energy',this_energy
            diff = abs(this_energy-E)
            diffs.append([diff,g,E])

        diffs.sort()
        print diffs
        self.PDB_text.insert(END,'Residue, Energy, Diff\n') 
        for diff,group,E in diffs[:100]:
            self.PDB_text.insert(END,'%10s %7.2f %7.2f\n' %(group,E,diff)) 
        return



    def construct_permutations(self,permutation_size,choices):
        """Construct all permutations of a combination of <permutation_size> items of <choices>"""
        print 'Construct all permutations of a combination of %d items of %d' %(permutation_size,len(choices))
        permutations={}
        x=[]
        count=[]
        for level in range(permutation_size):
            x.append(choices)
            count.append(0)
        done=None
        while not done:
            this_perm=[]
            #
            # Construct this permutation
            #
            for pos in range(len(count)):
                c_value=count[pos]
                if not x[pos][c_value] in this_perm:
                    this_perm.append(x[pos][c_value])
            #
            # Is this a valid permutation?
            #
            if len(this_perm)==permutation_size:
                this_perm=str(this_perm)
                if not permutations.has_key(this_perm):
                    permutations[this_perm]=None
            #
            # Increment count
            #
            count[0]=count[0]+1
            for pos in range(len(count)):
                if count[pos]==len(x[pos]):
                    if pos+1==len(count):
                        done=1
                        break
                    count[pos+1]=count[pos+1]+1
                    count[pos]=0
        return permutations
    
#
# ------
#

class distance_geometry:

    def __init__(self,D_matrix):
        """Store the matrix"""
        self.D_matrix=D_matrix
        self.complete_matrix()
        print self.D_matrix
        return

    #
    # ----
    #

    def complete_matrix(self):
        """Fill in distances of 0 with self"""
        for group in self.D_matrix.keys():
            self.D_matrix[group][group]=0.0

    #
    # ----
    #

    def do_triangle_smoothing(self):
        """Make sure that the distances are consistent"""
        anychange=1
        count=1
        while anychange:
            anychange=None
            for A in self.D_matrix.keys():
                for B in self.D_matrix.keys():
                    if B==A:
                        continue
                    for C in self.D_matrix.keys():
                        if C==B or C==A:
                            continue
                        #
                        # Flag for changes
                        #
                        changed=None
                        #
                        # Get the distances
                        #
                        AB=self.D_matrix[A][B]
                        AC=self.D_matrix[A][C]
                        BC=self.D_matrix[B][C]
                        #
                        # AC cannot be bigger than AB+BC
                        #
                        if AC>AB+BC:
                            AC=AB+BC
                        #
                        # AC dist must at least be AB-CB
                        #
                        if AC<AB-BC:
                            AC=AB-BC
                        #
                        # Put vals back into the matrix
                        #
                        if changed:
                            self.D_matrix[A][B]=AB
                            self.D_matrix[B][A]=AB
                            self.D_matrix[A][C]=AC
                            self.D_matrix[C][A]=AC
                            self.D_matrix[B][C]=BC
                            self.D_matrix[C][B]=BC
                            anychange=1
            print 'Triangle smoothing round %3d done' %count
            count=count+1
        print 'Triangle smoothing converged'
        print self.D_matrix
        print
        return self.D_matrix

    #
    # -----
    #

    def construct_metric_matrix(self):
        """Construct the metric matrix"""
        #
        # Set up groups
        #
        groups=self.D_matrix.keys()
        groups.sort()
        #
        # Distance from centre of coordinates
        #
        import math
        centre_dist_sq={}
        #
        # Avg of all distances
        #
        for A in groups:
            #
            # The avg sq dist for this group
            #
            sum=[]
            for B in groups:
                sum.append(math.pow(self.D_matrix[A][B],2))
            term1=average(sum)
            #
            # term2
            #
            term2=[]
            for x in range(1,len(groups)):
                for y in range(0,x):
                    term2.append(math.pow(self.D_matrix[groups[x]][groups[y]],2))
            print term2
            sum=0.0
            for i in term2:
                sum=sum+i
            term2=sum/float(math.pow(len(groups),2))
            print term2
            #
            # Get the dist
            #
            centre_dist_sq[A]=term1-term2
            print 'Distance to centre',A,centre_dist_sq[A]
                    
            
        #
        # Metric matrix
        #
        numgroups=len(groups)
        metric=numpy.zeros([numgroups,numgroups])
        for x in range(len(groups)):
            for y in range(len(groups)):
                A=groups[x]
                B=groups[y]
                metric[x][y]=(centre_dist_sq[A]+centre_dist_sq[B]-pow(self.D_matrix[A][B],2.0))/2.0
        #
        print 'Metric matrix'
        print metric
        print
        vals,vecs=numpy.linalg.eig(metric)
        print
        print 'Eigenvalues'
        print vals
        print 'Eigenvectors'
        print vecs
        l=[]
        for val in vals:
            l.append(val)
        l.sort()
        l.reverse()
        #
        done=[]
        count=0
        for best_val in l[:3]:
            for x in range(len(vals)):
                if best_val==vals[x] and not x in done:
                    done.append(x)
        print done
        
        count=0
        coords={}
        for group in groups:
            coords[group]={}
            coords[group]['X']=vecs[count][done[0]]*math.sqrt(abs(l[0]))
            coords[group]['Y']=vecs[count][done[1]]*math.sqrt(abs(l[1]))
            coords[group]['Z']=vecs[count][done[2]]*math.sqrt(abs(l[2]))
            print '%s %5.2f %5.2f %5.2f' %(group,coords[group]['X'],coords[group]['Y'],coords[group]['Z'])
            count=count+1
        return coords
            
        

# Some api in the chain is translating the keystrokes to this octal string
# so instead of saying: ESCAPE = 27, we use the following.
ESCAPE = '\033'

# Number of the glut window.
window = 0

# Rotation angle for the triangle. 
rtri = 0.0

# Rotation angle for the quadrilateral.
rquad = 0.0

class OGL(Frame):

    def __init__(self,coords,parent):
        self.parent=parent
        #
        # Vars
        #
        self.x_center=1
        self.y_center=1
        self.z_center=1
        self.count=0
        self.count_add=-1
        self.add=numpy.array([0,0,-1])
        #
        #
        self.main(coords) # draw the spheres
        
        #self.sphere_center=numpy.array([-1.5,0.0,-6.0])
        #self.add2=-numpy.array([0,0,0.1])
        #
        # Open a tkinter window
        #
        #self.geom_control=Toplevel()
        #self.geom_control.title('Tkinter')
        #self.geom_control.geometry('+700+100')
        #self.X=DoubleVar()
        #self.X.set(-1.5)
        #Scale(self.geom_control,variable=self.X,from_=-10,to=10,resolution=0.1).grid(row=0,column=0)
        #Button(self.geom_control,text='Quit',command=self.exit_application).grid(row=1,column=0)
        return

    #
    # ----
    #

    def exit_application(self,event=None):
        """Quit the application"""
        self.geom_control.destroy()
        import time
        time.sleep(1)
        glutDestroyWindow(window) #This apparently destroys everything...
        return

    #
    # ------
    # 

        
    def DrawGLScene(self):
        #
        # The main drawing function. 
        #
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)    # Clear The Screen And The Depth Buffer
        #
        # Colours
        #
        cols=[[1,0,0],
              [0,1,0],
              [0,0,1],
              [1,1,0],
              [1,0,1],
              [0,1,1]]
        #
        # Draw the spheres
        #
        count=0
        for group in self.coords.keys():
            self.coords[group]=self.coords[group]+self.add
            center=self.coords[group]
            self.count=self.count+self.count_add
            self.sphere(center,0.1,cols[count])
            count=count+1
        #
        # Update counter
        #
        if self.count<-70 or self.count>0:
            self.count_add=-self.count_add
            self.add=-1*self.add
        #
        #  since this is double buffered, swap the buffers to display what just got drawn. 
        #
        glutSwapBuffers()
        #
        # Sleep
        #
        import time
        time.sleep(0.1)
        #
        # Update Tkinter
        #
        self.parent.Dcontrol.update()
        return

    #
    # -----
    #

    # A general OpenGL initialization function.  Sets all of the initial parameters. 
    def InitGL(self,Width, Height):                # We call this right after our OpenGL window is created.
        glClearColor(0.0, 0.0, 0.0, 0.0)    # This Will Clear The Background Color To Black
        glClearDepth(1.0)                    # Enables Clearing Of The Depth Buffer
        glDepthFunc(GL_LESS)                # The Type Of Depth Test To Do
        glEnable(GL_DEPTH_TEST)                # Enables Depth Testing
        glShadeModel(GL_SMOOTH)                # Enables Smooth Color Shading

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()                    # Reset The Projection Matrix
                                            # Calculate The Aspect Ratio Of The Window
        gluPerspective(45.0, float(Width)/float(Height), 0.1, 100.0)

        glMatrixMode(GL_MODELVIEW)
        return


    # The function called when our window is resized (which shouldn't happen if you enable fullscreen, below)
    def ReSizeGLScene(self,Width, Height):
        if Height == 0:                        # Prevent A Divide By Zero If The Window Is Too Small 
            Height = 1

        glViewport(0, 0, Width, Height)        # Reset The Current Viewport And Perspective Transformation
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45.0, float(Width)/float(Height), 0.1, 100.0)
        glMatrixMode(GL_MODELVIEW)
        return

 

    #
    # ----
    #

    def sphere(self,position,radius,rgb):
        glLoadIdentity()                    # Reset The View
        glTranslatef(position[0],position[1],position[2])                # Move to the position
        quadratic=gluNewQuadric()
        glColor3f(rgb[0],rgb[1],rgb[2])            # Set The Color To Blue
        gluSphere(quadratic,radius,100,100)
        return




    # The function called whenever a key is pressed. Note the use of Python tuples to pass in: (key, x, y)  
    def keyPressed(self,*args):
        # If escape is pressed, kill everything.
        if args[0] == ESCAPE:
            sys.exit()

    def main(self,coords):
        """Set up everything"""
        self.coords={}
        for group in coords.keys():
            self.coords[group]=numpy.array([coords[group]['X'],coords[group]['Y'],coords[group]['Z']])
        #
        #
        #
        global window
        glutInit(sys.argv)

        # Select type of Display mode:   
        #  Double buffer 
        #  RGBA color
        # Alpha components supported 
        # Depth buffer
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)

        # get a 640 x 480 window 
        glutInitWindowSize(640, 480)

        # the window starts at the upper left corner of the screen 
        glutInitWindowPosition(0, 0)

        # Okay, like the C version we retain the window id to use when closing, but for those of you new
        # to Python (like myself), remember this assignment would make the variable local and not global
        # if it weren't for the global declaration at the start of main.
        window = glutCreateWindow("3D representation")

        # Register the drawing function with glut, BUT in Python land, at least using PyOpenGL, we need to
        # set the function pointer and invoke a function to actually register the callback, otherwise it
        # would be very much like the C version of the code.    
        glutDisplayFunc(self.DrawGLScene)

        # Uncomment this line to get full screen.
        # glutFullScreen()

        # When we are doing nothing, redraw the scene.
        glutIdleFunc(self.DrawGLScene)

        # Register the function called when our window is resized.
        glutReshapeFunc(self.ReSizeGLScene)

        # Register the function called when the keyboard is pressed.  
        glutKeyboardFunc(self.keyPressed)

        # Initialize our window. 
        self.InitGL(640, 480)
        return


#
# ---------
#

