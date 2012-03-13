from numpy import *
from numpy.linalg import * 
inverse=inv
Float=float

class LM_functions:

    def fit_LM_ghost(self):
        """Do Levenberg-Marquardt fitting"""
        J,E =self.get_jacobian_ghost()
        JT = transpose(J)
        JTE = dot(JT,E)
        JTJ = dot(JT,J)
        JTJd = JTJ + self.LM_damper*identity(shape(JTJ)[0])
        invJTJd = inv(JTJd)
        q = -dot(JTE,invJTJd)
        
        self.cubes=sorted(self.cube_eps.keys())
        for cube in self.cubes:
            print 'Cube %4d q: %5.2f, actual: %5.2f' %(cube,q[cube],self.cube_eps[cube]+q[cube])
            self.cube_eps[cube]=min(100.0,max(1.0,self.cube_eps[cube]+q[cube])) # We don't allow any dielectric constant to drop below 1.0 or get above 100.0
        return
    
    #
    # ----------
    #

    def get_jacobian_ghost(self):
        """Get the Jacobian matrix and errors of the data points"""
        #
        # Get the number of data points
        #
        no_data_points=0
        for titgroup in sorted(self.exp_ghosts.keys()):
            for residue in sorted(self.exp_ghosts[titgroup].keys()):
                for atom in ['H','N']:
                    no_data_points=no_data_points+1
        #
        #no_data_points = len(self.exp_data)
        errors = resize(array(0,float),[no_data_points])        
        jacobian = resize(array(0,float),[no_data_points,len(self.cube_eps.keys())])
        #
        # Precalculate the variation of all parameters
        #
        now=self.current_ghosts #self.get_spans()
        variations=[]
        step = 1
        #
        self.cubes=sorted(self.cube_eps.keys())
        self.cube_grad=self.get_cube_scores()
        for cube in self.cubes:
            variations.append(self.cube_grad[cube])
        #
        # construct jacobian
        #
        data_id=0
        #
        # Experimental data points
        #
        for titgroup in sorted(self.exp_ghosts.keys()):
            self.errors[titgroup]={}
            for residue in sorted(self.exp_ghosts[titgroup].keys()):
                for atom in ['H','N']:
                    if self.exp_ghosts[titgroup][residue].has_key(atom):
                        errors[data_id] = self.get_error(titgroup,residue,atom,now)[0]
                        #
                        # Find the derivative for this data point for this cube
                        #
                        diff=resize(array(0,float),[len(self.cubes)])
                        count=0
                        for variation in variations:
                            diff[count]=(now[titgroup][residue][atom]-variation[titgroup][residue][atom])/step
                            count=count+1
                        jacobian[data_id]=diff
                        data_id=data_id+1
        return jacobian,errors

        
           
