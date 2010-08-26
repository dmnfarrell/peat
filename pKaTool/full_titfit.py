# -*- coding: iso-8859-15 -*-
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


class titfitter:

    def __init__(self,filename=None,pH_start=None,pH_end=None):
        """Load a full set of titration data from an Ekin project file"""
        import tkFileDialog, os
        savedir=os.getcwd()
        if not filename:
            filename=tkFileDialog.askopenfilename(defaultextension="*.Ekinprj",
                                                  initialdir=savedir,
                                                  filetypes=[("Ekin project file","*.Ekinprj"),("All files","*.*")])
        if filename:
            fd=open(filename)
            import pickle
            titdata=pickle.load(fd)
            fd.close()
            #
            # Extract data into a more reasonable format
            #
            self.alltitdata={}
            groups=titdata.keys()
            groups.sort()
            for group in groups:
                
                if group[:2]=='__' or group=='mode':
                    pass
                else:
                    #
                    # Get the atom to map data to
                    #
                    # In the future we will get this from the Data source mapping
                    #
                    if titdata.has_key('__structure_mappings__'):
                        if titdata['__structure_mappings__'].has_key(group):
                            for mapping in titdata['__structure_mappings__'][group].keys():
                                thismapping=titdata['__structure_mappings__'][group][mapping]
                                if 'Data source' in thismapping['Data property']:
                                    print group,thismapping
                                    print 'and then write the code for getting the info'
                    #
                    # For now we just grab the residue number from the key
                    #
                    resnumber=group[1:]
                    Npos=resnumber.find('N')
                    import string
                    resnumber=resnumber[:Npos]
                    try:
                        int(resnumber[1:])
                    except:
                        continue
                    resnumber=':'+string.zfill(resnumber,4)
                    import types
                    self.alltitdata[resnumber]=[]
                    if type(titdata[group][0]) is types.DictType:
                        for datapoint in titdata[group][0].keys():
                            import types
                            if type(titdata[group][0][datapoint]) is types.DictType:
                                if titdata[group][0][datapoint].has_key('var'):
                                    pH=titdata[group][0][datapoint]['var']
                                    chemshift=titdata[group][1][datapoint]['var']
                                    #
                                    # Convert to floats and filter out emtry strings
                                    #
                                    import string
                                    if string.strip(pH)=='':
                                        pH=None
                                    if string.strip(chemshift)=='':
                                        chemshift=None
                                    if chemshift and pH:
                                        #
                                        # pH range criteria
                                        #
                                        if float(pH)>=pH_start and float(pH)<=pH_end:
                                            self.alltitdata[resnumber].append([float(pH),float(chemshift)])
                    self.alltitdata[resnumber].sort()
            #
            # Draw the buggers
            #
            self.draw_titrations()
        return

    #
    # ----
    #

    def draw_titrations(self):
        """Draw all the titrations"""
        import EM_effect_struct_win
        self.EStruct=EM_effect_struct_win.EM_effect_struct(self)
        self.EStruct.draw_titrations(self.alltitdata)
        return

    #
    # ----
    #

    def put_data_in_fitter_format(self):
        """Rearrange all the data in the fitter format"""
        #
        # Subtract the first point to get rid of the abs value
        #
        exp_data=[]
        for residue in self.alltitdata.keys():
            if len(self.alltitdata[residue])<2:
                continue
            first_chemshift=self.alltitdata[residue][0][1]
            for ph,chemshift in self.alltitdata[residue]:
                exp_data.append([residue,ph,chemshift-first_chemshift])
        return exp_data

    
