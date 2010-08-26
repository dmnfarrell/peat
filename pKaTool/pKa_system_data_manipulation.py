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
import ftir_data

class data_manipulation:

    def pack_all_data(self):
        """Pack all data in a dictionary"""
        #
        # Dictionary for storing everything
        #
        pKaTool_data = {}
        #
        # Store groups
        #
        pKaTool_data['groups'] = {}
        no_groups = len(self.groups)
        pKaTool_data['no_groups'] = no_groups
        
        for i in range(no_groups):
            pKaTool_data['groups'][i] = {}
            pKaTool_data['groups'][i]['pka'] = self.groups[i].intpka.get()
            pKaTool_data['groups'][i]['name'] = self.groups[i].name.get()
            pKaTool_data['groups'][i]['ab'] = self.groups[i].acid_base.get()
            for j in self.groups[i].intenes.keys():
                pKaTool_data['groups'][i][j] = self.groups[i].intenes[j].get()
        #
        # Store experimental data
        #
        import copy
        if getattr(self,'titration_data',None):
            pKaTool_data['titration_data'] = copy.deepcopy(self.titration_data)
        #
        # store pH-activity profile
        #
        if getattr(self,'activity_data',None):
            pKaTool_data['pH-activity profile'] = copy.deepcopy(self.activity_data)
            pKaTool_data['CCPS'] = {}
            for state in self.act_state.keys():
                pKaTool_data['CCPS'][state]=self.act_state[state].get()

        #
        # store FTIR data
        #
        if getattr(self,'FTIR_win',None):
            if len(self.FTIR_win.ftir_data)>0:
                pKaTool_data['FTIR data'] = copy.deepcopy(self.FTIR_win.ftir_data)
                pKaTool_data['FTIR data offset'] = self.FTIR_win.offset.get()
                pKaTool_data['FTIR data c'] = self.FTIR_win.c.get()
                
                
        return pKaTool_data

    #
    # -----
    #

    def unpack_all_data(self,pKaTool_data):
        """Unpack all data"""
        #
        # load groups
        #
        print 'DATA'
        print pKaTool_data
        no_groups = pKaTool_data['no_groups']
        for i in range(no_groups):
            self.groups[i].intpka.set(pKaTool_data['groups'][i]['pka'])
            if pKaTool_data['groups'][i].has_key('name'):
                self.groups[i].name.set(pKaTool_data['groups'][i]['name'])
            else:
                self.groups[i].name.set('Group %d' %i)
            print i
            print pKaTool_data['groups'][i]
            self.groups[i].acid_base.set(pKaTool_data['groups'][i]['ab'])
            for j in range(i):
                self.groups[i].intenes[j].set(pKaTool_data['groups'][i][j]) 
        #
        # load experimental data
        #
        if pKaTool_data.has_key('titration_data'):
            self.titration_data = pKaTool_data['titration_data'].copy()
            self.display_loaded_curves.set(1)
        #
        # load pH-activity profile
        #
        if pKaTool_data.has_key('pH-activity profile'):
            self.micro_var.set(1)    
            self.activity_data = pKaTool_data['pH-activity profile'].copy()
            self.update_pkasystem_curves()
            for state in pKaTool_data['CCPS'].keys():
                self.act_state[state].set(pKaTool_data['CCPS'][state])
        #
        # load FTIR data
        #
        if pKaTool_data.has_key('FTIR data'):
            if not getattr(self,'FTIR_win',None):
                self.FTIR_win = ftir_data.FTIR_data(self)
                self.show_ftir.set(1)
                
            self.FTIR_win.ftir_data=pKaTool_data['FTIR data'].copy()
            self.FTIR_win.offset.set(pKaTool_data['FTIR data offset'])
            self.FTIR_win.c.set(pKaTool_data['FTIR data c'])

        #
        # Update display
        #
        self.update_scales_from_fit()
        self.update_pkasystem_curves()
        return

    #
    # -----
    #

    def send_system_to_EAT(self):
        """Get all data and send it to EAT"""
        #
        # Do we have an EAT instance?
        #
        EAT_instance=None
        if self.parent_application:
            if self.parent_application.ID=='EAT':
                EAT_instance=self.parent_application
            elif self.parent_application.ID=='Ekin':
                if self.parent_application.parent:
                    if self.parent_application.parent.ID=='EAT':
                        EAT_instance=self.parent_application.parent
        if not EAT_instance:
            import tkMessageBox
            tkMessageBox.showerror('Not called from EAT',
                                   'You have to call pKa_system from EAT or from EAT via Ekin to use this function',
                                   parent=self.window)
            return
        #
        # Get the data for the system
        #
        PSdata=self.pack_all_data()
        #
        # Add info for EAT if we have it
        #
        import copy
        if self.protein or self.field_name:
            PSdata['__EATinfo__']={'record':copy.deepcopy(self.protein),'column':copy.deepcopy(self.field_name)}
        #
        # Send the whole thing to EAT
        #
        import copy
        EAT_instance.pKa_system_data=copy.deepcopy(PSdata)
        #
        # Destroy window
        #
        if self.parent_application.ID=='EAT':
            self.window.destroy()
            self.win.destroy()
        else:
            #
            # We have to manually load data into EAT
            #
            EAT_instance.load_pKa_system_data(self.window)
        return
