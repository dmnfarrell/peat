#!/usr/bin/env python
#
# Constants
#
ppm_au=1.9447E-18
NSP={'N':977,'H':90}
e=1.602E-19 # electron charge
e0=8.85E-12 # vacuum permittivity
kB=1.3806503E-23 #Boltzmann's constant

import EM_effect, numpy

class map_struct(EM_effect.EM_effect):

    def __init__(self,pdbfile):
        """Load the PDB file"""
        import Protool
        self.PI=Protool.structureIO()
        self.PI.readpdb(pdbfile)
        return

    #
    # -----
    #
    
    def _calc_eff_eps(self,dCS,dist,angle,atom,charge):
        import math
        #E=1.0/(4*math.pi*e0*deff)*e/(r*1E-10)**2*ppm_au*NSP[atom_type]*1E6*cos_angle*charge
        cos_angle=math.cos(math.radians(angle))
        val=dCS*4*math.pi*e0*(dist*1E-10)**2*-charge
        val=val/(NSP[atom]*e*ppm_au*cos_angle*1E6*-1)
        return 1.0/val

    #
    # -----
    #

    def get_eff_eps(self,dCS,titgroup,atom,charge):
        """Get the effective dielectric constant if given a change in chemical shift, a
        titratable group and an atom where the change in chemical shift is measured"""
        titpos=self.get_charge_center(titgroup)
        atom_type=atom.split(':')[-1]
        resid=self.PI.resnum(atom)
        cos_angle,distance=self.get_angle(residue=resid,titpos=titpos,atom_type=atom_type)
        import math
        angle=math.degrees(math.acos(cos_angle))
        #print 'Distance: %5.3fA, angle: %5.3f' %(distance,angle)
        #
        # Get the atom type
        #
        atomtype=atom[-1]
        if atomtype not in ['N','H']:
            raise Exception('I only know of N and H')
        return self._calc_eff_eps(dCS,distance,angle,atomtype,charge)

    #
    # ----
    #

    def get_Eint(self,deff,distance,T=298.15):
        """if given a distance and an effective dielectric constant, return
        interaction energy in kT/e"""
        import math
        conv=kB*T/e #Conversion from volts to kT/e
        return 1.0/(4*math.pi*e0*deff)*e/(distance*1E-10)/conv

    #
    # ----
    #

    def get_dCS(self,atom,titgroup,charge=1.0,method='PBE',diel=8.0):
        """Calculate the change in the chemical shift of atom due
        to the charge present on the titgroup"""
        #
        # Get the residue ID
        #
        if not self.PI.atoms.has_key(atom):
            return None
        residue=self.PI.resnum(atom)
        #
        # Center of titgroup
        #
        if not hasattr(self,'titpos_centers'):
            self.titpos_centers={}
        if self.titpos_centers.has_key(titgroup):
            titpos=self.titpos_centers[titgroup]
        else:
            titpos=self.get_charge_center(titgroup)
            self.titpos_centers[titgroup]=titpos
        # Atom type
        atom_type=atom.split(':')[-1]
        #
        # Get the value of the span
        #
        if method=='Coulomb':
            dCS,angle=self.get_dEF(residue,titpos,charge,diel,atom_type=atom_type)
        elif method=='PBE':
            dCS=self.get_span_from_PBE(residue,titpos,charge,diel,atom_type=atom_type)
            angle=0.0
        elif method=='APBS':
            dCS=self.get_span_from_PBE(residue,titpos,charge,diel,atom_type=atom_type,APBS=1)
            angle=0.0
        else:
            raise Exception,'Unknown method for calculating electric field'
        return dCS

    #
    # ----
    #

    def build_Hs(self):
        """Build all HNs on the protein"""
        import Protool.hydrogens
        H=Protool.hydrogens.Hydrogens(self.PI)
        H.build_all_HNs()
        return
    

def main():
    """Calculate effective dielectric constants"""
    X=map_struct('2lzt.pka.pdb')
    eff_eps=X.get_eff_eps(0.45,titgroup=':0052:ASP',atom=':0035:GLU:N')
    print 'Effective eps',eff_eps
    # Do the reverse calc
    span,angle=X.get_dEF(residue=':0035',
                         titpos=X.get_charge_center(':0052:ASP'),
                         charge=-1,deff=eff_eps,atom_type='N')
    print 'Backcalc, span: %5.3f, angle: %5.3f' %(span,angle)
    print
    eff_eps=X.get_eff_eps(0.2,titgroup=':0035:GLU',atom=':0052:ASP:N')
    print 'ghost in D52',eff_eps
    print 'Eint',X.get_Eint(21,7.3)


if __name__=='__main__':
    main()
