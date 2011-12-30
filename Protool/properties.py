#
# # Protool - Python class for manipulating protein structuress
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

def calculate_mass(PI):
    """Calculate the mass of all protein atoms in a Protool instance"""
    weight=0.0
    for residue in PI.residues:
        if PI.isaa(residue):
            weight=weight+PI.all_attribute[PI.three_to_one[PI.resname(residue)]]['Mw']
    return weight-16.0
            
def calculate_unit_cell_volume(PI):
    """calculate the volume of the unit cell"""
    cr=PI.cryst[0].split()
    a=float(cr[1])
    b=float(cr[2])
    c=float(cr[3])
    from math import *
    alpha=radians(float(cr[4]))
    beta=radians(float(cr[5]))
    gamma=radians(float(cr[6]))
    Vuc= a*b*c*sqrt(1.0-pow(cos(alpha),2)-pow(cos(beta),2)-pow(cos(gamma),2)+2.0*cos(alpha)*cos(beta)*cos(gamma))
    return Vuc

SGdata={'P 42cm': 8, 'P 6(3)/mmc': 24, 'I4(1)/amd': 32, 'P 4 3 2': 24, 'P 42/mcm': 16, 'I4/mcm': 32, 'P-42(1)m': 8, 'Pn-3': 24, 'P4/nnc': 16, 'P-3c1': 12, 'Abm2': 8, 'P 21': 2, 'P-42(1)c': 8, 'P4/mbm': 16, 'P4/ncc': 16, 'P 41 21 2': 8, 'P3(1)12': 6, 'Pmc2(1)': 4, 'I-42d': 16, 'I4(1)22': 16, 'Cmm2': 8, 'I-42m': 16, 'P6(5)22': 12, 'R32': 18, 'C222(1)': 8, 'R-3c': 36, 'P321': 6, 'R-3': 18, 'R3': 9, 'P-6c2': 12, 'I 4 3 2': 48, 'R-3m': 36, 'Pbcn': 8, 'Pbcm': 8, 'Pbca': 8, 'R3m': 18, 'P6/mcc': 24, 'Pmn2(1)': 4, 'Cmma': 16, 'Cmmm': 16, 'P6(3)22': 12, 'R3c': 18, 'Iba2': 8, 'Pm-3m': 48, 'F432': 96, 'Ia-3': 48, 'Ama2': 8, 'P 43 2 2': 8, 'P 21 3': 12, 'F222': 16, 'P 4 21 2': 8, 'P3(2)21': 6, 'P 4/n': 8, 'Pmmn': 8, 'Pmmm': 8, 'Imma': 16, 'Pmma': 8, 'P 21/c': 4, 'P 2/m': 4, 'I4(1)/acd': 32, 'P6mm': 12, 'I-4': 8, 'P3m1': 6, 'P 2/c': 4, 'P 21/m': 4, 'Pccn': 8, 'I2(1)3': 24, 'Pccm': 8, 'P-4c2': 8, 'Pnn2': 4, 'Pcca': 8, 'Pbam': 8, 'P 41 2 2': 8, 'Pban': 8, 'P 42 /ncm': 16, 'P-3m1': 12, 'P4(2)mc': 8, 'P 62 2 2': 12, 'P 41': 4, 'Ibca': 16, 'P 2': 2, 'P3': 3, 'P 1': 1, 'P6': 6, 'P4(3)32': 24, 'P 4': 4, 'P 42 /mbc': 16, 'P4(2)/m': 8, 'P4(2)/n': 8, 'Pmna': 8, 'P4(2)bc': 8, 'P-4b2': 8, 'Pmm2': 4, 'Fd-3': 96, 'Imm2': 8, 'P3(2)': 3, 'P6(3)': 6, 'P-6': 6, 'P4/mcc': 16, 'P 41 32': 24, 'P6(1)': 6, 'I-43d': 48, 'Pba2': 4, 'Pnnm': 8, 'P 42': 4, 'Pcc2': 4, 'Pnna': 8, 'P6(1)22': 12, 'P 43 21 2': 8, 'C2/c': 8, 'P-62c': 12, 'Pnma': 8, 'P4/m': 8, 'P6/mmm': 24, 'Amm2': 8, 'P3c1': 6, 'P 21 21 21': 4, 'Fmm2': 16, 'I 2 3': 24, 'P-6m2': 12, 'P4/nbm': 16, 'I4/mmm': 32, 'F4(1)32': 96, 'Fm-3m': 192, 'Immm': 16, 'P3(1)': 3, 'P4/mmm': 16, 'P4(2)/nmc': 16, 'P4cc': 8, 'Pa-3': 24, 'P6(3)cm': 12, 'Aba2': 8, 'P 43': 4, 'Fddd': 32, 'P6cc': 12, 'I4(1)': 8, 'Ibam': 16, 'C2': 4, 'P-4m2': 8, 'I-43m': 48, 'I422': 16, 'Pca2(1)': 4, 'P3(2)12': 6, 'P 42 32': 24, 'I2(1)2(1)2(1)': 8, 'I4(1)cd': 16, 'I4/m': 16, 'P6(2)': 6, 'P 2 3': 12, 'Fmmm': 32, 'F-43c': 96, 'P6/m': 12, 'P-43n': 24, 'P6(4)22': 12, 'P4bm': 8, 'Fdd2': 16, 'Cm': 4, 'Pc': 2, 'Cc': 4, 'I222': 8, 'F23': 48, 'P-31c': 12, 'P31c': 6, 'I-4m2': 16, 'P3(1)21': 6, 'P-31m': 12, 'P4nc': 8, 'I4mm': 16, 'I4(1)32': 48, 'I4(1)/a': 16, 'P4/mnc': 16, 'I-4c2': 16, 'P622': 12, 'Pna2(1)': 4, 'Im-3': 48, 'P6(3)/mcm': 24, 'P 42 21 2': 8, 'P312': 6, 'P6(4)': 6, 'F4-3m': 96, 'Pn-3n': 48, 'Pn-3m': 48, 'P 2 2 21': 4, 'I4': 8, 'Ccca': 16, 'C2/m': 8, 'P-42c': 8, 'P-42m': 8, 'Pm': 2, 'Cccm': 16, 'C222': 8, 'P31m': 6, 'P4/nmm': 16, 'Cmcm': 16, 'Pnnn': 8, 'Fm-3': 96, 'Cmca': 16, 'P4mm': 8, 'P422': 8, 'P 42 nm': 8, 'P4(2)/nbc': 16, 'Pma2': 4, 'P6(3)/m': 12, 'P6(3)mc': 12, 'P-62m': 12, 'Pnc2': 4, 'P-4': 4, 'Ccc2': 8, 'P-4n2': 8, 'I4cm': 16, 'Cmc2(1)': 8, 'P-1': 2, 'P-43m': 24, 'P-3': 6, 'P 2 2 2': 4, 'P 42/mnm': 16, 'P 65': 6, 'Pm-3': 24, 'P 42 2 2': 8, 'P 21 21 2': 4, 'P 42/nnm': 16, 'I4(1)md': 16, 'Pm-3n': 48, 'P 42/mmc': 16, 'Ima2': 8}

def get_asym_units(PI):
    SG=PI.spacegroup
    print 'Looking for spacegroup: %s' %SG
    if SGdata.has_key(SG):
        print 'Found!'
        return SGdata[SG]
    # No exact match, then we search for the closest match
    for cut in range(len(SG),3,-1):
        nSG=SG[:cut]
        if SGdata.has_key(nSG):
            print 'Found %s' %nSG
            return SGdata[nSG]
    print 'Spacegroup not found: %s' %SG
    return 1

def calculate_Vm(PI):
    """Calculate Vm, Protein volume, Solvent volume and crystal density"""
    # Mass
    Mw=calculate_mass(PI)
    # Number of assymmetric units in unit cell
    z=get_asym_units(PI)
    # Unit cell volume
    Vuc=calculate_unit_cell_volume(PI)
    # Calculate Vm
    Vm=Vuc/(Mw*z)
    # Calculate protein volume
    Na=6.02E23
    densprot=1.35E-24 # g / A**3
    protein_vol=100.0*Mw*z/(densprot*Na)/Vuc
    solvent_vol=100.0-protein_vol
    cryst_dens=(protein_vol*densprot*1E24+(100.0-protein_vol)*1.00)/100.0
    return Vm,protein_vol,solvent_vol,cryst_dens,Vuc,Mw*z
