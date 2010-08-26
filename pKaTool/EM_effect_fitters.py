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

import LM_Fitter
class flexible_fitter(LM_Fitter.LM_Fitter):
    """This class fits anything depending on fit_variables: pKa values, positions and effective dielectric constants for any number of groups

        fit_variables is a dictionary: {<group>:{<var1>:startval,<var2>:startval2}, <group2>:{<var1>:startval1}}

    """

    def __init__(self,exp_data,fit_variables,nonfit_variables,EM_effect_instance,callback=None,calc_type='Coulomb'):
        """parse fit_variables - construct two lists:
            start_values - just holds the start values
            var_name: defines the type of each variable - order corresponds to start_values
        """
        def_startvals={'pka':4.0,'deff':8.0,'X':0.0,'Y':0.0,'Z':0.0}
        start_values=[]
        self.value_position={}
        self.groups=fit_variables.keys()
        self.groups.sort()
        self.EI=EM_effect_instance
        count=0
        for group in fit_variables.keys():
            self.value_position[group]={}
            for var in fit_variables[group].keys():
                self.value_position[group][var]=count
                if not fit_variables[group][var]:
                    start_values.append(def_startvals[var])
                else:
                    start_values.append(fit_variables[group][var])
                count=count+1
        #
        self.calc_type=calc_type
        self.nonfit_variables=nonfit_variables
        LM_Fitter.LM_Fitter.__init__(self,start_values,exp_data,callback)
        return
        
    #
    # ----
    #
    
    def get_value(self,function_variables,data_point):
        """Calculate the function value at the data point"""
        #
        # Data point
        #
        residue=data_point[0]
        atom_type=data_point[1]
        pH=data_point[2]
        charge=data_point[3]
        #
        # If we work with spans, then pH must be set to -99
        #
        if pH<-90:
            mode='span'
        else:
            mode='titration'
        #
        # Function variables
        #
        sum_value=0.0
        for group in self.value_position.keys():
            import numpy, math
            if self.value_position[group].has_key('X'):
                x_index=self.value_position[group]['X']
                xpos=function_variables[x_index]
                
                y_index=self.value_position[group]['Y']
                ypos=function_variables[y_index]
                
                z_index=self.value_position[group]['Y']
                zpos=function_variables[z_index]
            else:
                xpos=self.nonfit_variables[group]['X']
                ypos=self.nonfit_variables[group]['Y']
                zpos=self.nonfit_variables[group]['Z']
            titpos=numpy.array([xpos,ypos,zpos])
            #
            # Get the charge
            #
            if mode=='span':
                charge=charge
            elif mode=='titration':
                if self.value_position[group].has_key('pKa'):
                    pka_index=self.value_position[group]['pKa']
                    pKa=function_variables[pka_index]
                    charge=-1.0/(1.0+math.pow(10.0,pKa-pH))
                else:
                    print 'No pKa variable found?'
                    charge=charge
            else:
                print 'Something is really wrong'
                raise Exception
            #
            # Get the dielectric constant
            #
            # Write code for global eps, and for 3D epsmap
            #   
            if self.value_position[group].has_key('deff'):
                deff_index=self.value_position[group]['deff']
                deff=function_variables[deff_index]
            else:
                deff=self.nonfit_variables[group]['deff']
            #
            # Calculation type - Coulomb or PBE?
            #
            if self.calc_type=='Coulomb':
                value,degrees=self.EI.get_dEF(residue,titpos,charge,deff,atom_type=atom_type)
            elif self.calc_type=='PBE':
                value=self.EI.get_span_from_PBE(residue,titpos,charge,deff,atom_type=atom_type)
                angle=0.0
            else:
                raise Exception,'Calculation method not known'
            sum_value=sum_value+value
        return sum_value



import LM_Fitter
class pKa_pos_deff(LM_Fitter.LM_Fitter):
    """This class fits the pKa values, positions and effective dielectric constants for any number of groups"""

    def __init__(self,exp_data,start_values,pka_groups,callback=None,calc_type='Coulomb'):
        self.EI=EM_instance
        self.pka_groups=pka_groups
        self.calc_type=calc_type
        LM_Fitter.LM_Fitter.__init__(self,start_values,exp_data,callback)
        return
    
    def get_value(self,function_variables,data_point):
        #
        # Format for variables: [pka1, xpos1, ypos1, zpos1,pka2, xpos2, ypos2, zpos2,...., deff]
        # Data points [residue, pH, chemshift]
        #
        # Data point
        #
        residue=data_point[0]
        pH=data_point[1]
        #
        # Function variables
        #
        sum_value=0.0
        for pkagroup in range(self.pka_groups):
            import numpy, math
            pKa=function_variables[0+pkagroup*4]
            xpos=function_variables[pkagroup*4+1]
            ypos=function_variables[pkagroup*4+2]
            zpos=function_variables[pkagroup*4+3]
            deff=function_variables[-1]
            titpos=numpy.array([xpos,ypos,zpos])
            charge=-1.0/(1.0+math.pow(10.0,pKa-pH))
            if self.calc_type=='Coulomb':
                value,degrees=self.EI.get_dEF(residue,titpos,charge,deff,atom_type='N')
            elif self.calc_type=='PBE':
                span=self.EI.get_span_from_PBE(residue,titpos,charge,deff,atom_type='N')
                angle=0.0
            else:
                raise Exception,'Calc method not known'
            sum_value=sum_value+value
        return sum_value
                    

#
# ----
#

import LM_Fitter
class pos_deff(LM_Fitter.LM_Fitter):
    """This class fits the position and effective dielectric constant for a single group
    from the residue number and the span

    We always assume full charge"""

    def __init__(self,exp_data,start_values,callback=None,calc_type='Coulomb'):
        self.EI=EM_instance
        self.calc_type=calc_type
        LM_Fitter.LM_Fitter.__init__(self,start_values,exp_data,callback)
        return
        
    #
    # -----
    #
    
    def get_value(self,function_variables,data_point):
        #
        # Format for variables: [xpos1, ypos1, zpos1, deff1]
        # Data points [residue, span]
        #
        # Data point
        #
        residue=data_point[0]
        span=data_point[1]
        #
        # Function variables
        #
        import numpy, math
        xpos=function_variables[0]
        ypos=function_variables[1]
        zpos=function_variables[2]
        deff=function_variables[3]
        titpos=numpy.array([xpos,ypos,zpos])
        charge=-1.0
        if self.calc_type=='Coulomb':
            value,degrees=self.EI.get_dEF(residue,titpos,charge,deff,atom_type='N')
        elif self.calc_type=='PBE':
            value=self.EI.get_span_from_PBE(residue,titpos,charge,deff,atom_type='N')
        else:
            raise Exception,'Calc method not known'
        return value
        
        
#
# -------
#

class deff(LM_Fitter.LM_Fitter):
    """This class fits effective dielectric constant for a single group
    from the residue number and the span

    We always assume full charge"""

    def __init__(self,exp_data,start_values,callback=None,calc_type='Coulomb'):
        self.callback=callback
        self.calc_type=calc_type
        LM_Fitter.LM_Fitter.__init__(self,start_values,exp_data,callback)
        return
        
    #
    # -----
    #
    
    def get_value(self,function_variables,data_point):
        #
        # Format for variables: [deff1]
        # Data points [titposx, titposy, titposz, charge, atom_type, residue, span]
        #
        # Data point
        #
        xpos=data_point[0]
        ypos=data_point[1]
        zpos=data_point[2]
        charge=data_point[3]
        atom_type=data_point[4]
        residue=data_point[5]
        span=data_point[6]
        #
        # Function variables
        #
        import numpy, math
        deff=function_variables[0]
        titpos=numpy.array([xpos,ypos,zpos])
        #
        # Calculate the span using all the variables. Use either Coulomb's law or PBE
        #
        if self.calc_type=='Coulomb':
            value,degrees=self.EI.get_dEF(residue,titpos,charge,deff,atom_type=atom_type)
        elif self.calc_type=='PBE':
            value=self.EI.get_span_from_PBE(residue,titpos,charge,deff,atom_type=atom_type)
        else:
            raise Exception,'Calc method not known'
        return value

    #
    # -----
    #
   
