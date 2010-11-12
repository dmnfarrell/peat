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
from titration_class import *
import numpy
import numpy.linalg
import random
import math
import Pmw
import pickle
import ftir_data


class file_menu:

    def load_system(self):
        """ loads parameters and experimental curves from a file """
        import tkFileDialog, os
        savedir=os.getcwd()
        filename=tkFileDialog.askopenfilename(defaultextension="*.pks",
                                              initialdir=savedir,
                                              filetypes=[("pKaTool system file","*.pks"),("All files","*.*")])
        if filename:
            if os.path.isfile(filename):
                file = open(filename,'r')
                pKaTool_data = pickle.load(file)
                self.unpack_all_data(pKaTool_data)
        return

    #
    # -----
    #

    def save_system(self):
        """ saves parameters and experimental curves to a file """
        import tkFileDialog, os
        savedir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension="*.pks",
                                                initialdir=savedir,
                                                filetypes=[("pKaTool system file","*.pks"),("All files","*.*")],
                                                parent=self.window)
        if filename:
            pKaTool_data=self.pack_all_data()
            fd=open(filename,'w')
            pickle.dump(pKaTool_data,fd)
            fd.close()

        return

    #
    # -----
    #

    def load_curves(self, filename = '', acid_bases = None):
        """Read titration curves from a text file"""
        import tkFileDialog, os
        savedir=os.getcwd()
        if filename == '':
            filename=tkFileDialog.askopenfilename(defaultextension="*.tcrv",
                                                  initialdir=savedir,
                                                  filetypes=[("pKaTool titration curve file","*.tcrv"),
                                                             ("Text file","*.txt"),
                                                             ("CSV file",".csv"),
                                                             ("All files","*.*")])
        if filename:
            if os.path.isfile(filename):
                #
                # Get the type
                #
                if filename[-4:]=='.csv':
                    filetype='csv'
                    split_char=','
                else:
                    filetype='spacesep'
                    split_char=None
                #
                # Read the file
                #
                fd=open(filename)
                local_titration_data={}
                header=fd.readline()
                header=header.split(split_char)
                exp_groups=header[1:]
                #
                # store experimental data in local_titration_data
                #
                keys = exp_groups
                for key in keys:
                    local_titration_data[key]={}

                line=fd.readline()
                while line.find('END')==-1 and len(line)>2:
                    sp=line.split(split_char)
                    pH=float(sp[0])

                    count=1
                    for key in keys:
                        try:
                            local_titration_data[key][pH]=float(sp[count])
                        except:
                            print 'Could not load data point for',key,'at pH',pH
                        count=count+1
                    #
                    # Read the next line
                    #
                    line=fd.readline()
                fd.close()

                self.convert_titration_data(local_titration_data, acid_bases=acid_bases)

            else:
                import tkMessageBox
                tkMessageBox.showwarning('File does not exist',
                                         'Could not find %s.' %filename)
        return

    def load_titdb(self, filename = '', acid_bases = None):
        """Read titration curves from a titdb exported file - added by damien"""
        import tkFileDialog, os
        savedir=os.getcwd()
        if filename == '':
            filename=tkFileDialog.askopenfilename(defaultextension="*.tcrv",
                                                  initialdir=savedir,
                                                  filetypes=[("CSV file",".csv"),
                                                             ("All files","*.*")])
        if filename:
            import csv
            if os.path.isfile(filename):
                # Read the file
                csvfile = open(filename)
                reader = csv.reader(csvfile)
                local_titration_data={}
                for row in reader:
                    if len(row)==1: #dataset name
                        current = row[0]
                        local_titration_data[current]={}
                    else:
                        ph=float(row[0]); val=float(row[1])
                        local_titration_data[current][ph]=val

                print local_titration_data

                self.convert_titration_data(local_titration_data, acid_bases=acid_bases)
            else:
                import tkMessageBox
                tkMessageBox.showwarning('File does not exist',
                                         'Could not find %s.' %filename)
        return

    #
    # ---
    #

    def convert_titration_data(self,exp_tit_data, acid_bases = None):
        """Scale the loaded curves and fix the experimental -> theoretical curve mapping"""
        if not getattr(self,'titration_data',None):
            self.titration_data={}
        #
        # associate experimental curves with groups
        #
        do_scaling = IntVar()
        keys = exp_tit_data.keys()

        # experimental -> theoretical group mapping, if not provided
        if acid_bases == None:
            acid_bases = []
            for key in keys:
                # exp residue, acid/base,Group name
                temp = [key, StringVar(),StringVar()]
                acid_bases.append(temp)
            popup = Dialog(self.window, acid_bases, do_scaling, self.titration_curves,self.groups)
            self.window.wait_window(popup.top)

        else:
            do_scaling.set(1)

        for i in range(len(acid_bases)):
            if not acid_bases[i][2].get() == 'None':
                #print 'Matching',acid_bases[i][2].get(),'with',acid_bases[i][0]
                #
                # Name groups in gui
                #
                j = int(acid_bases[i][2].get()[1:5]) #This is probably a stupid way to find the group index
                self.titration_data[acid_bases[i][2].get()]= exp_tit_data[acid_bases[i][0]]

                self.groups[j].name.set(keys[i])
        #
        # scale data to fit titration curves
        #
        if do_scaling.get() == 1:
            for key in self.titration_data.keys():
                min = 1e6
                max = -1e6
                pHatMax = 0
                pHatMin = 0
                for ph in self.titration_data[key].keys():
                    cur = self.titration_data[key][ph]
                    if cur > max:
                        max = cur
                        pHatMax = ph
                    if cur < min:
                        min = cur
                        pHatMin = ph

                type = 1
                if pHatMax < pHatMin:
                    type = 2

                a_b = ''
                for i in range(len(acid_bases)):
                    if key == acid_bases[i][2].get():
                        a_b = acid_bases[i][1].get()

                for ph in self.titration_data[key].keys():
                    if type == 1:
                        if a_b == 'acid':
                            self.titration_data[key][ph] = -(self.titration_data[key][ph]-min)/(max - min)
                        elif a_b == 'base':
                            self.titration_data[key][ph] = -(self.titration_data[key][ph]-max)/(max - min)
                    elif type == 2:
                        if a_b == 'acid':
                            self.titration_data[key][ph] = (self.titration_data[key][ph]-max)/(max - min)
                        elif a_b == 'base':
                            self.titration_data[key][ph] = (self.titration_data[key][ph]-min)/(max - min)

            #
            # Set the flag for displaying loaded curves
            #
            self.display_loaded_curves.set(1)
            self.update_pkasystem_curves()
        return

    #
    # ----
    #


    def save_curves(self):
        """Save titration curves in a file"""
        import tkFileDialog, os
        savedir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension='.tcrv',
                                                initialdir=savedir,
                                                filetypes=[("pKaTool titration curve file","*.tcrv"),("All files","*.*")])
        if filename:
                #
                # If the file is already there then delete it
                #
                if os.path.isfile(filename):
                    os.unlink(filename)
                #
                # Write everything nicely
                #
                fd=open(filename,'w')
                groups=self.titration_curves.keys()
                groups.sort()
                pHs=self.titration_curves[groups[0]].keys()
                pHs.sort()
                text='pH\t'
                for group in groups:
                    text=text+'%s\t' %group
                fd.write(text+'\n')
                for pH in pHs:
                    if pH=='pKa':
                        continue
                    text='%5.2f\t' %pH
                    for group in groups:
                        text=text+'%9.2f\t' %self.titration_curves[group][pH]
                    fd.write(text+'\n')
                fd.write('END\n')
                fd.close()
        return

    #
    # -----
    #

    def load_pH_activity_profile(self, filename=None,parent=None):
        """Read pH activity profile from a text file"""
        if not parent:
            parent=self.win
        import tkFileDialog, os
        if not filename:
            savedir=os.getcwd()
            filename=tkFileDialog.askopenfilename(defaultextension="*.txt",
                                                  initialdir=savedir,
                                                  filetypes=[("Text file","*.txt"),
                                                             ("Comma separated","*.csv"),
                                                             ("All files","*.*")],
                                                  parent=parent)

        if filename:
            if os.path.isfile(filename):
                if filename[-4:]=='.csv':
                    splitter=','
                else:
                    splitter=None
                #
                # Read the file
                #
                fd=open(filename)
                self.activity_data={}
                header=fd.readline()
                header=header.split()

                line = fd.readline()
                while len(line)>2 and not line =='END':
                    words = line.split(splitter)
                    print words
                    self.activity_data[float(words[0])] = float(words[1])
                    print 'loaded',words[0],':  ',words[1]
                    line = fd.readline()
                fd.close()
                #
                # normalise profile
                #
                max = 0
                for key in self.activity_data.keys():
                    if self.activity_data[key] > max:
                        max = self.activity_data[key]
                #print 'max',max
                for key in self.activity_data.keys():
                    self.activity_data[key] = self.activity_data[key]/max

                #print self.activity_data

                self.micro_var.set(1.0)
                self.update_pkasystem_curves()


            else:
                import tkMessageBox
                tkMessageBox.showwarning('File does not exist',
                                         'Could not find %s.' %filename,
                                         parent=parent)

        return

    #
    # ----
    #

    def load_2nd_pH_activity_profile(self,filename=None,parent=None):
        """Load a 2nd pH-activity profile"""
        if getattr(self,'activity_data',None):
            self.backup=self.activity_data.copy()
        else:
            self.backup={}
        self.load_pH_activity_profile(filename,parent=parent)
        self.secondary_activity_data=self.activity_data.copy()
        #
        # copy the old one back
        #
        self.activity_data=self.backup.copy()
        self.update_pkasystem_curves()
        return

    #
    # --------
    #

    def add_group(self):
        return

    def remove_exp_curve(self):
        win = Toplevel()
        win.title('Remove experimental curves')

        row = 1
        col = 1

        # make window
        if getattr(self,'titration_data',None):
            removes = []
            for i in self.titration_data.keys():
                v = IntVar()
                removes.append(v)

            Label(win,text='Remove curves:').grid(row=row, column=col)
            for i in range(len(self.titration_data.keys())):
                j = int(self.titration_data.keys()[i][1:5]) #This is probably a stupid way to find the group index
                row+=1
                Checkbutton(win,text='%s ' %self.groups[j].name.get(),variable=removes[i]).grid(row=row, column=col)

        else:
            Label(win,text='No titration curves loaded').grid(row=row, column=col)
        row += 1
        Button(win,text='Ok',command=win.destroy).grid(row=row, column=col)
        self.window.wait_window(win)

        # remove selected curve(s)
        if getattr(self,'titration_data',None):
            i = 0
            total = len(self.titration_data.keys())
            while i < total:
                if removes[i].get() == 1:
                    j = int(self.titration_data.keys()[i][1:5])
                    self.groups[j].name.set('Group %d'%j)
                    del self.titration_data[self.titration_data.keys()[i]]
                    del removes[i]
                    total-=1
                else:
                    i+=1
            # delete the self.titration_data dictionary, if no more exp titration curves are present
            if len(self.titration_data.keys()) == 0:
                del self.titration_data
                self.display_loaded_curves.set(0)

        self.titwin.update()
        self.update_pkasystem_curves()
        return

    def print_table(self):
        return


    def load_FTIR_data(self):
        """Method for loading FTIR data"""

        if not getattr(self, 'FTIR_win',None):
            self.FTIR_win = ftir_data.FTIR_data(self)
            self.show_ftir.set(1)


        self.FTIR_win.load_ftir_data()

        return

#
# -----
#

class fit_menu:

    def fit_ftir(self):
        """Method for fitting FTIR data"""
        if getattr(self,'FTIR_win',None):
            self.FTIR_win.fit()

        return



    def fit_system_to_curves(self):
        """Given loaded titration curves, attempt to fit the system characteristics
        (intrinsic pKa value and interaction energies) so that we get the best match.
        We map the titration curves to groups (in the system) based on the charges"""
        #
        # Can we fit?
        #
        if not getattr(self,'titration_data',None):
            import tkMessageBox
            tkMessageBox.showwarning('No titration curves loaded',
                                     'Load titration curves first')
            return

        self.fit()
        return

    #
    # -----
    #

    def uniqueness_scan(self):
        """Attempts to find critiria for when a unique solution can be found. Combinatorial scans are carried out on generated titration curves based on varius pKa values and interaction energies """
        self.display_loaded_curves.set(1)
#        pka_dists = range(0,10)
#        for i in range(len(pka_dists)):
#            pka_dists[i]=pka_dists[i]/1.0
        pka_dists = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

#        int_dists = range(0,10)
#        for i in range(len(int_dists)):
#            int_dists[i]=int_dists[i]/1.0
        int_dists = [0.0, 0.5, 1.0, 1.5, 2.0]
        out = 'pka_dist   int_dist   no_close  max_para_dist'
        print out
        for pka_dist in pka_dists:
            for int_dist in int_dists:
                #
                # move to new position
                #
                self.groups[0].intpka.set(4)
                self.groups[0].intpka.set(4+pka_dist)
                self.groups[1].intenes[0].set(int_dist)

                now = self.eval(silent=1)
                #
                # copy this conformation
                #
                self.titration_data = {}
                for group in now.curves.keys():
                    self.titration_data[group] = {}
                    for ph in now.curves[group].keys():
                        if not ph == 'pKa':
                            if not ph == 'pka':
                                self.titration_data[group][ph] = float('%4.1f'%now.curves[group][ph])################ one decimal point
                #file name
                fn = 'close_'+str(int(pka_dist*10))+'_'+str(int(int_dist*10))+'.txt'

                #max_diff = self.test_uniqueness(silent=1)
                no_close,max_para_dist = self.combinatorial_scan(silent=1,filename=fn)
                out = '%3.2f      %3.2f      %5d      %3.2f' %(pka_dist,int_dist,no_close,max_para_dist)
                print out

        return
    #
    # ---
    #



    def test_uniqueness(self,silent=0):
        """Tests uniqueness by performing a large number of fits using different start conditions"""
        #
        # make a reference fit
        #
        out_file = open('test_uniqueness.txt','w')

        import random
        self.fit(silent=1)
        first_vars=[]
        out = ''

        for key in self.groups.keys():
            out += '%8s '%self.groups[key].name.get()
            for key2 in self.groups[key].intenes.keys():
                out += '%4s-%4s '%(self.groups[key].name.get(),self.groups[key2].name.get())
        out += '\n'

        for group in self.groups:
            first_vars.append(self.groups[group].intpka.get())
            out += '%8.4f '%self.groups[group].intpka.get()
            for group2 in self.groups[group].intenes.keys():
                first_vars.append(self.groups[group].intenes[group2].get())
                out += '%8.4f '%self.groups[group].intenes[group2].get()
        first_error = self.evaluate_fit(make_win = 0)[1]

        out += '%8.4f '%0.0 #difference in parameters from reference fit
        error = self.evaluate_fit(make_win=0)[1]
        out += '%8.4f \n'%error

        print out
        out_file.write(out)
        out_file.write('#---------------------------------------------------------\n')

        #
        # Try out a lot of fits with different starting parameters
        #
        max_max_diff = 0
        for i in range(1000):
            print '+++++++++++ %d +++++++++++++++++++++++++++++\n' %i
            #
            # Randomize start conditions
            #
            for group in self.groups:
                self.groups[group].intpka.set(random.uniform(2,12))
                for group2 in self.groups[group].intenes.keys():
                    self.groups[group].intenes[group2].set(random.uniform(0,10))
            self.update_scales_from_fit()
            self.titwin.update()
            #
            # do fit
            #
            status = self.fit(silent=1)
            #
            # compare results to reference fit
            #
            cur_error = self.evaluate_fit(make_win = 0)
            #            if cur_error < 0.1:
            cur_vars = []
            out = ''
            for group in self.groups:
                cur_vars.append(self.groups[group].intpka.get())
                out += '%8.4f '%self.groups[group].intpka.get()
                for group2 in self.groups[group].intenes.keys():
                    cur_vars.append(self.groups[group].intenes[group2].get())
                    out += '%8.4f '%self.groups[group].intenes[group2].get()

            #
            # calculate difference in parameters from reference fit
            #
            diff_vars = 0
            for i in range(len(cur_vars)):
                diff_vars += (cur_vars[i]-first_vars[i])*(cur_vars[i]-first_vars[i])
            diff_vars = math.sqrt(diff_vars)

            out += '%8.4f '%diff_vars

            #
            # calculate error of this fit
            #
            error = self.evaluate_fit(make_win=0)[1]
            out += '%8.4f %3s\n'%(error,status)

            print out
            out_file.write(out)

        out_file.close()
        return max_max_diff
    #
    # ---
    #

    def fit(self, silent=0):
        """Do 10000 rounds of fitting, or until converged"""
        #
        # Set damper
        #
        self.LM_damper = DoubleVar()
        self.LM_damper.set(0.1)
        #
        # Open the window
        #
        if not silent == 1:
            self.count_win=Toplevel()
            self.count_win.geometry('100x100+200+200')
            self.count_win.title('Fit in progress')
            self.count=IntVar()
            self.count.set(0)
            Label(self.count_win,text='Iteration #').grid(row=0,column=0)
            Label(self.count_win,textvariable=self.count).grid(row=0,column=1)
            #
            # Keep track of number of pKa calcs
            #
            self.pka_count=IntVar()
            self.pka_count.set(0)
            Label(self.count_win,text='pKa calc #').grid(row=1,column=0)
            Label(self.count_win,textvariable=self.pka_count).grid(row=1,column=1)
            #
            Button(self.count_win,text='Stop fit',command=self.stop_fit).grid(row=2,column=0,columnspan=2)
            self.count_win.update_idletasks()
        self.keep_running=1
        #
        # Make copies of all the variables
        #
        self.vars=[]
        for group in self.groups:
            self.vars.append(self.groups[group].intpka)
            for group2 in self.groups[group].intenes.keys():
                self.vars.append(self.groups[group].intenes[group2])
        #
        # Set up the loaded curves
        #
        self.exp_data=titration_curve(self.titration_data)
        #
        # Start iterating
        #
        old_diff=0.0
        status = 'MI'
        if not silent == 1:
            print 'Step    Diff   LM_damper'
        for x in range(1,1000):#1000

            if not self.keep_running:
                break
            self.fit_LM(silent=silent)
            #
            # Update more stuff
            #
            if not silent==1:
                self.count.set(self.count.get()+1)
                self.count_win.update()
            self.update_scales_from_fit()
            self.titwin.update()
            #self.master.update_idletasks()
            now_diff=self.eval(silent=silent)-self.exp_data
            #
            # Check convergence
            #
            if not silent == 1:
                print '%5d  %6.4f  %6.4f' %(x, now_diff, self.LM_damper.get())

            if abs(now_diff-old_diff)<0.00005:
                status = 'C'
                if not silent == 1:
                    print 'Converged',now_diff
                break
            else:
                old_diff=now_diff


        if not silent==1:
            self.count_win.destroy()
        self.update_pkasystem_curves()
        return status

    #
    # -----
    #



    def fit_LM(self,silent=0):
        """Do Levenberg-Marquardt fitting"""

        J,E =self.get_jacobian(silent=silent)
        JT = numpy.transpose(J)
        JTE = numpy.dot(JT,E)
        JTJ = numpy.dot(JT,J)
        JTJd = JTJ + self.LM_damper.get()*numpy.identity(numpy.shape(JTJ)[0])
        invJTJd = numpy.linalg.inv(JTJd)
        q = -numpy.dot(JTE,invJTJd)
        out1 ='bv '
        out2 ='av '
        for var in range(len(self.vars)):
            out1 += '%4.6f '%self.vars[var].get()
            self.vars[var].set(self.vars[var].get()+q[var])
            out2 += '%4.6f '%self.vars[var].get()
        return

    #
    # -----
    #


    def stop_fit(self,event=None):
        """Stop the fit"""
        self.keep_running=None
        return

    #
    # ------
    #

    def get_jacobian(self,silent=0):
        """Get the Jacobian matrix and errors of the data points"""
        #
        # Get the number of data points
        #
        count=0
        for group in self.exp_data.curves.keys():
            for ph in self.exp_data.curves[group].keys():
                count=count+1
        no_data_points = count
        #
        #
        #
        errors = numpy.resize(numpy.array(0,float),[no_data_points])
        jacobian = numpy.resize(numpy.array(0,float),[no_data_points,len(self.vars)])
        #temp = numpy.resize(numpy.array(0,float),[len(self.vars)])
        #
        # Precalculate the variation of all parameters
        #
        now=self.eval(silent=silent)
        variations=[]

        step = 1e-8

        for var in range(len(self.vars)):
            self.vars[var].set(self.vars[var].get()+step)
            variations.append(self.eval(silent=silent))
            self.vars[var].set(self.vars[var].get()-step)
        #
        # construct jacobian
        #
        data_id=0

        for group in self.exp_data.curves.keys():
            for ph in self.exp_data.curves[group].keys():
                data_point = self.exp_data.curves[group][ph]
                x=float(ph)
                y=float(data_point)

                errors[data_id] = y-self.titration_curves[group][ph]

                #
                # Find the derivative at this ph (x) for this data point
                #
                diff=numpy.resize(numpy.array(0,float),[len(self.vars)])
                count=0
                for variation in variations:
                    diff[count]=(now.curves[group][ph]-variation.curves[group][ph])/step
                    count=count+1
                jacobian[data_id]=diff
                data_id=data_id+1

        return jacobian,errors

    #
    # ------
    #

    def combinatorial_scan(self,silent=0,filename='close.txt'):
        """Do a full combinatorial scan"""
        if not getattr(self,'titration_data',None):
            import tkMessageBox
            tkMessageBox.showwarning('No titration curves loaded',
                                     'Load titration curves first')
            return


        #
        # Open config window
        #
        self.pka_min = StringVar()
        self.pka_min.set(4.0)
        self.pka_max = StringVar()
        self.pka_max.set(9.0)
        self.pka_step = StringVar()
        self.pka_step.set(0.1)

        self.int_min = StringVar()
        self.int_min.set(0.0)
        self.int_max = StringVar()
        self.int_max.set(4.0)
        self.int_step = StringVar()
        self.int_step.set(0.1)

        self.close_threshold = StringVar()
        self.close_threshold.set(10.0)

        if silent == 0:
            CS = CS_dialog(self)
            self.window.wait_window(CS.CS)

        #
        # Make copies of all the variables
        #
        self.pkas=[]
        self.interactions = []
        for group in self.groups:
            self.pkas.append(self.groups[group].intpka)
            for group2 in self.groups[group].intenes.keys():
                self.interactions.append(self.groups[group].intenes[group2])
        #
        # Set up the loaded curves
        #
        self.exp_data=titration_curve(self.titration_data)
        #
        # more setup
        #
        self.pka_values = range(int(float(self.pka_min.get())*10),
                                int(float(self.pka_max.get())*10),
                                int(float(self.pka_step.get())*10))
        for i in range(len(self.pka_values)):
            self.pka_values[i] = self.pka_values[i]/10.0

        self.interaction_values = range(int(float(self.int_min.get())*10),
                                        int(float(self.int_max.get())*10),
                                        int(float(self.int_step.get())*10))
        for i in range(len(self.interaction_values)):
            self.interaction_values[i] = self.interaction_values[i]/10.0

        #
        # calculate no. of combinations
        #
        self.no_comb=1
        for group in self.groups:
            self.no_comb=self.no_comb*len(self.pka_values)
            for group2 in self.groups[group].intenes.keys():
                self.no_comb=self.no_comb*len(self.interaction_values)


        self.record = 1e9
        self.pka_count=IntVar()
        self.pka_count.set(0)

        self.close_parameters = []



        #
        # Open counter window
        #
        self.stop_CS = IntVar()
        self.stop_CS.set(0)
        self.record_gui=StringVar()
        self.record_gui.set('None')
        self.count=IntVar()
        self.count.set(0)
        if not silent==1:
            self.count_win=Toplevel()
            self.count_win.geometry('400x200+200+200')
            self.count_win.title('Combinatorial scan')

            Label(self.count_win,text='Combination ').grid(row=0,column=0,sticky='w')
            Label(self.count_win,textvariable=self.count).grid(row=0,column=1,sticky='w')
            Label(self.count_win,text=' of %d'%(self.no_comb)).grid(row=0,column=2)
            Label(self.count_win,text='Smallest AUE so far:').grid(row=1,column=0,sticky='w')
            Label(self.count_win,textvariable=self.record_gui).grid(row=1,column=1,sticky='w')
            Label(self.count_win,text='Best configuration found so far').grid(row=2,column=1,sticky='w')

        #
        # count no of points
        #
        self.no_points = 0
        for t in self.exp_data.curves.keys():
            self.no_points = self.no_points + len(self.exp_data.curves[t].keys())
        #
        # make record vairables
        #

        self.record_vars=[]
        count = 0
        r = 3
        c = 0
        for group in self.groups:
            self.record_vars.append(StringVar())
            self.record_vars[count].set(self.groups[group].intpka.get())
            if not silent==1:
                Label(self.count_win,text='%s'%self.groups[group].name.get()).grid(row=r,column=c,sticky='w')
                c+=1
                Label(self.count_win,textvariable=self.record_vars[count]).grid(row=r,column=c,sticky='w')
            count +=1
            for group2 in self.groups[group].intenes.keys():
                self.record_vars.append(StringVar())
                self.record_vars[count].set(self.groups[group].intenes[group2].get())
                if not silent==1:
                    c+=1
                    Label(self.count_win,textvariable=self.record_vars[count]).grid(row=r,column=c,sticky='w')
                    count +=1
            c=0
            r+=1

        if not silent==1:
            Button(self.count_win,text='Stop',command=self.cs_stopper).grid(row=r,column=1)

        self.scan_pka(0,silent)

        if not silent==1:
            self.count_win.destroy()
        #
        # set parameters to optimal configuration
        #
        count=0
        for group in self.groups:
            self.groups[group].intpka.set(self.record_vars[count].get())
            count=count+1
            for group2 in self.groups[group].intenes.keys():
                self.groups[group].intenes[group2].set(self.record_vars[count].get())
                count=count+1
        self.update_scales_from_fit()
        self.titwin.update()
        #
        # write close configs to file
        #

        cf = open(filename,'w')

        cf.write('#A record of %f was found with the parameter set:\n' %self.record)
        for r in self.record_vars:
            cf.write('%15s'%r.get())
        cf.write('\n\n')

        cf.write('#Parameter sets with errors < %s * record = %f\n' %(self.close_threshold.get(), float(self.close_threshold.get())*self.record))
        cf.write('#[parameters], error, parameter_dist\n')
        max_para_dist=0
        for cp in self.close_parameters:
            line = ''
            para_dist = 0
            for i in range(len(cp)):
                line+='%15s'%cp[i]
                if not i == len(cp)-1:
                    para_dist += (float(cp[i])-float(self.record_vars[i].get()))*(float(cp[i])-float(self.record_vars[i].get()))
            para_dist = math.sqrt(para_dist)
            if para_dist>max_para_dist:
                max_para_dist=para_dist
            line += '%15.2f\n'%para_dist
            cf.write(line)
        cf.write('\n\n')
        cf.close()


        return len(self.close_parameters),max_para_dist
    #
    # ---
    #

    def scan_pka(self,i,silent=0):
        for pv in self.pka_values:
            self.pkas[i].set(pv)
            if not i == len(self.pkas)-1:
                self.scan_pka(i+1,silent)
            else:
                self.scan_int(0,silent)

            if(self.stop_CS.get()==1):
                return
        return

    def scan_int(self,i,silent=0):
        for iv in self.interaction_values:
            self.interactions[i].set(iv)
            if not i == len(self.interactions)-1:
                self.scan_int(i+1,silent)
            else:
                self.check_record(silent)
            if(self.stop_CS.get()==1):
                return
        return


    def check_record(self,silent=0):
        now_diff=(self.eval(silent)-self.exp_data)/(float(self.no_points))

        # Is this parameter set close to being a record?
        if now_diff <= float(self.close_threshold.get())*self.record:
            close_vars = []
            for group in self.groups:
                close_vars.append(self.groups[group].intpka.get())
                for group2 in self.groups[group].intenes.keys():
                    close_vars.append(self.groups[group].intenes[group2].get())

            close_vars.append('%2.5f'%now_diff)

#            abs_sum_error,scaled_abs_error = self.evaluate_fit(make_win=0)

#            close_vars.append('%2.5f'%abs_sum_error)
#            close_vars.append('%2.5f'%scaled_abs_error)

            self.close_parameters.append(close_vars[:])
            # is it also a new record?
            if now_diff <= self.record:
                count = 0
                for group in self.groups:
                    self.record_vars[count].set(self.groups[group].intpka.get())
                    count+=1
                    for group2 in self.groups[group].intenes.keys():
                        self.record_vars[count].set(self.groups[group].intenes[group2].get())
                        count+=1
                self.record = now_diff

                # given the new record, should we remove some of the old close-to-be-records?
                remove_these_close =[]
                for c in self.close_parameters:
                    if float(c[len(c)-1]) >= float(self.close_threshold.get())*float(self.record):
                        remove_these_close.append(c[:])

                for c in remove_these_close:
                    self.close_parameters.remove(c)

                self.record_gui.set('%4.3f'%self.record)
                self.update_scales_from_fit()
                self.titwin.update()
        if not silent ==1:
            self.count.set(self.count.get()+1)
            self.count_win.update()
        return

    def cs_stopper(self):
        self.stop_CS.set(1)
        return


    def show_close(self, filename='close.txt'):
        """Displays close parameter sets found in a combinatorial scan"""
        import tkMessageBox
        if not getattr(self,'titration_data',None):
            tkMessageBox.showwarning('No titration curves loaded',
                                     'Make sure to load the titration curves used in the last combinatorial scan first')
            return
        try:
            lines = open(filename,'r')
        except:
            tkMessageBox.showwarning('Could not open file',
                                     'Could not open file %s'%filename)
            return

        record_vars = []
        self.close_vars = []

        first = 1
        for line in lines:
            if not line[0:1] == '#' and len(line) > 4:
                words = line.split()

                for i in range(len(words)):
                    words[i]=float(words[i])

                # record vars
                if first == 1:
                    record_vars = words
                    first = 0
                # close vars
                else:
                    self.close_vars.append(words)


        #print 'record_vars',record_vars
        #        print 'close_vars',close_vars

        self.close_vars.sort(lambda x,y: cmp(x[len(x)-2],y[len(y)-2]))

        #        for i in range(len(old)):
        #            print '--------------------------------------------'
        #            print close_vars[i]



        scd = Show_close_dialog(self)
        return

    def identify_no_groups(self):
        return



    def evaluate_fit(self, make_win=1):
        import tkMessageBox
        if not getattr(self,'titration_data',None):
            tkMessageBox.showwarning('No titration curves loaded',
                                     'Load titration curves first')
            return

        self.exp_data=titration_curve(self.titration_data)

        no_points = 0
        for t in self.exp_data.curves.keys():
            no_points = no_points + len(self.exp_data.curves[t].keys())

        abs_sum_error=self.eval(silent=1)-self.exp_data
        avr_abs_error=abs_sum_error/(float(no_points))
        #
        # Do scaled MUE
        #
        scaled_abs_error=self.eval(silent=1).sub_scaled(self.exp_data)
        #
        # Do HHd scaled MUE
        #
        calc_pkas = {}
        for group in self.groups.keys():
            calc_pkas[self.names[group]] = self.groups[group].pkavalue.get()
        HHd_scaled_abs_error=self.eval(silent=1).sub_HHd_scaled(self.exp_data,calc_pkas)
        avr_HHd_scaled_abs_error = HHd_scaled_abs_error/no_points


        if make_win == 1:
            text = """
Sum of absolute errors: %3.3f
Number of data points:  %d
Average absolute error: %3.5f
Scaled absolute error:  %3.5f
Sum of H-H dev. scaled error: %3.5f
Avr. H-H dev. scaled error: %3.5f
------------------------------------------------
Avr. absolute errors of individual curves:
""" %(abs_sum_error,
      no_points,
      avr_abs_error,
      scaled_abs_error,
      HHd_scaled_abs_error,
      avr_HHd_scaled_abs_error)


            individual_diffs = self.exp_data.subtract_individually(self.eval(silent=1))
            print 'individual_diffs',individual_diffs
            for t in range(len(self.exp_data.curves.keys())):
                j = int(self.exp_data.curves.keys()[t][1:5]) #This is probably a stupid way to find the group index
                text += 'Error of %-5s: %f\n' %(self.groups[j].name.get(),individual_diffs[t]/len(self.exp_data.curves[self.exp_data.curves.keys()[t]]) )


            tkMessageBox.showinfo('Evaluation of fit',text)


        return abs_sum_error,avr_abs_error,scaled_abs_error,avr_HHd_scaled_abs_error


    def eval(self,silent=0):
        """Evaluate the function"""
        X,pKa_values,prot_states=self.calc_pKas_from_scales(self.groups)
        curve=titration_curve(X.prot_states)
        if not silent == 1:
            self.pka_count.set(self.pka_count.get()+1)
            self.count_win.update()
        return curve



    def X(self):

        files = ['Data_asp_low_salt.txt','Data_glu_low_salt.txt']
        groups_in_files = {}
        all_groups = []



        #
        # Get groups from file headers
        #
        for file in files:
            f = open(file,'r')
            heads = f.readline()
            f.close()
            heads = heads.split()
            heads = heads[1:]
            for h in heads:
                all_groups.append(h)
            groups_in_files[file] = heads

        print 'groups_in_files',groups_in_files
        print 'all_groups',all_groups

        out = open('X.txt','w')

        out.write('Doing cross-fitting of the groups:\n')
        for file in files:
            out.write('%s: '%file)
            for g in groups_in_files[file]:
                out.write('%s '%g)
            out.write('\n')
        out.write('g1, g2, [parameters],abs_sum_error, avr_abs_error, scaled_abs_error, avr_HHd_scaled_abs_error\n\n')
        #
        # Cycle through all combinations
        #

        for g1 in all_groups:
            filename = ''
            for file in files:
                if g1 in groups_in_files[file]:
                    filename = file

            acid_bases = self.get_import_scheme(g1,groups_in_files[filename],0)
            print 'filename',filename

            for a in acid_bases:
                print '1:::', a[0],a[1].get(),a[2].get()
            self.load_curves(filename=filename,acid_bases=acid_bases)

            for g2 in all_groups:
                if not g1 == g2:
                    filename = ''
                    for file in files:
                        if g2 in groups_in_files[file]:
                            filename = file

                    acid_bases = self.get_import_scheme(g2,groups_in_files[filename],1)

                    print 'filename - 2',filename
                    self.load_curves(filename=filename,acid_bases=acid_bases)
                    self.update()

                    print 'Doing',g1,'and',g2

                    self.fit(silent=1)
                    res = '%5s %5s' %(g1,g2)
                    for group in self.groups:
                        res += '%6.2f '%self.groups[group].intpka.get()
                        for group2 in self.groups[group].intenes.keys():
                            res += '%6.2f '%self.groups[group].intenes[group2].get()



                    res += '%10f %10f %10f %10f\n' %(self.evaluate_fit(make_win=0))
                    print res
                    out.write(res)


        return

    def get_import_scheme(self, group, groups, index):
        """Makes import scheme for automatic load of curves - used by X """
        acid_bases = []
        for i in range(len(groups)):
            temp_ab = ['', StringVar(),StringVar()]
            acid_bases.append(temp_ab)

        for i in range(len(groups)):
            if group == groups[i]:
                acid_bases[i][0] = group
                acid_bases[i][1].set('acid')
                acid_bases[i][2].set(self.titration_curves.keys()[index])
            else:
                acid_bases[i][0] = ''
                acid_bases[i][1].set('acid')
                acid_bases[i][2].set('None')

        return acid_bases


    def estimate_experimental_uncertainty(self):
        """Estimates uncertainty of fits using experimental uncertainties """


        if not getattr(self,'titration_data',None):
            import tkMessageBox
            tkMessageBox.showwarning('No titration data',
                                     'You need to load experimental titration curves first!')
            return
        #
        # get parameters from user
        #
        titration_uncertainty = DoubleVar(value=0.01)
        pH_uncertainty = DoubleVar(value=0.1)
        activity_uncertainty = DoubleVar(value=0.05)
        no_runs = IntVar(value = 10)

        win = Toplevel(self)
        win.title('Experimental uncertainty')

        row = 1
        col = 1

        Label(win,text='Uncertainty of NMR data [charge units]:').grid(row=row,column=col)
        col+=1
        Entry(win,textvariable=titration_uncertainty).grid(row=row,column=col)
        col=1
        row+=1

        Label(win,text='Uncertainty of Activity data [normalised units]:').grid(row=row,column=col)
        col+=1
        Entry(win,textvariable=activity_uncertainty).grid(row=row,column=col)
        col=1
        row+=1

        Label(win,text='Uncertainty of pH values:').grid(row=row,column=col)
        col+=1
        Entry(win,textvariable=pH_uncertainty).grid(row=row,column=col)
        col=1
        row+=1

        Label(win,text='Number of runs:').grid(row=row,column=col)
        col+=1
        Entry(win,textvariable=no_runs).grid(row=row,column=col)
        col=1
        row+=1

        Button(win, text='Go',command=win.destroy).grid(row=row,column=col)

        self.window.wait_window(win)

        #
        # backup the experimental data
        #
        backup_titration_data = self.titration_data.copy()
        if getattr(self,'activity_data',None):
            backup_activity_data = self.activity_data.copy()

        #
        # Make arrays for parameters
        #
        parameters=[]
        gui_parameters=[]
        gui_parameters_standard_deviations=[]
        for group in self.groups:
            parameters.append(numpy.resize(numpy.array(0,float),[no_runs.get()]))
            gui_parameters.append(DoubleVar(value=0.0))
            gui_parameters_standard_deviations.append(DoubleVar(value=0.0))
            for group2 in self.groups[group].intenes.keys():
                parameters.append(numpy.resize(numpy.array(0,float),[no_runs.get()]))
                gui_parameters.append(DoubleVar(value=0.0))
                gui_parameters_standard_deviations.append(DoubleVar(value=0.0))
        #
        # Make arrays for errors
        #
        errors=[]
        gui_errors=[]
        gui_standard_deviations=[]
        for i in range(4):
            errors.append(numpy.resize(numpy.array(0,float),[no_runs.get()]))
            gui_errors.append(DoubleVar(value=0.0))
            gui_standard_deviations.append(DoubleVar(value=0.0))
        #
        # Make counter win
        #
        current_run = IntVar(value=0)

        self.eeu_win = Toplevel(self)
        self.eeu_win.title('Estimating experimental uncertianty')
        row=1
        col=1
        Label(self.eeu_win,text='Run number ').grid(row=row,column=col)
        col+=1
        Label(self.eeu_win,textvariable=current_run).grid(row=row,column=col)
        col+=1
        Label(self.eeu_win,text=' of ').grid(row=row,column=col)
        col+=1
        Label(self.eeu_win,textvariable=no_runs).grid(row=row,column=col)
        col=1
        row+=1
        #
        # headings
        #
        col=4
        Label(self.eeu_win,text='Average values').grid(row=row,column=col,columnspan=2)
        col+=2
        Label(self.eeu_win,text='Standard deviations').grid(row=row,column=col,columnspan=2)
        col=1
        row+=1
        #
        # average parameters and standard deviations
        #
        count_p=0
        count_s=0
        line=0
        for group in self.groups:
            Label(self.eeu_win,text='%s'%self.groups[group].name.get()).grid(row=row,column=col,sticky='w',columnspan=3)
            col+=3
            Label(self.eeu_win,textvariable=gui_parameters[count_p]).grid(row=row,column=col,sticky='w')
            count_p+=1
            for group2 in self.groups[group].intenes.keys():
                col+=1
                Label(self.eeu_win,textvariable=gui_parameters[count_p]).grid(row=row,column=col,sticky='w')
                count_p +=1
            col+=len(self.groups)-line
            Label(self.eeu_win,textvariable=gui_parameters_standard_deviations[count_s]).grid(row=row,column=col,sticky='w')
            count_s+=1
            for group2 in self.groups[group].intenes.keys():
                col+=1
                Label(self.eeu_win,textvariable=gui_parameters_standard_deviations[count_s]).grid(row=row,column=col,sticky='w')
                count_s +=1
            col=1
            row+=1
            line+=1
        #
        # Errors and standard deviations
        #
        Label(self.eeu_win,text='Absolut sum of errors:').grid(row=row,column=col,sticky='w',columnspan=3)
        col+=3
        Label(self.eeu_win,textvariable=gui_errors[0]).grid(row=row,column=col,sticky='w',columnspan=2)
        col+=2
        Label(self.eeu_win,textvariable=gui_standard_deviations[0]).grid(row=row,column=col,sticky='w',columnspan=2)
        row+=1
        col=1
        Label(self.eeu_win,text='Average absolut sum of errors:').grid(row=row,column=col,sticky='w',columnspan=3)
        col+=3
        Label(self.eeu_win,textvariable=gui_errors[1]).grid(row=row,column=col,sticky='w',columnspan=2)
        col+=2
        Label(self.eeu_win,textvariable=gui_standard_deviations[1]).grid(row=row,column=col,sticky='w',columnspan=2)
        row+=1
        col=1
        Label(self.eeu_win,text='Scaled absolut sum of errors:').grid(row=row,column=col,sticky='w',columnspan=3)
        col+=3
        Label(self.eeu_win,textvariable=gui_errors[2]).grid(row=row,column=col,sticky='w',columnspan=2)
        col+=2
        Label(self.eeu_win,textvariable=gui_standard_deviations[2]).grid(row=row,column=col,sticky='w',columnspan=2)
        row+=1
        col=1
        Label(self.eeu_win,text='Avr. HHd scaled error:').grid(row=row,column=col,sticky='w',columnspan=3)
        col+=3
        Label(self.eeu_win,textvariable=gui_errors[3]).grid(row=row,column=col,sticky='w',columnspan=2)
        col+=2
        Label(self.eeu_win,textvariable=gui_standard_deviations[3]).grid(row=row,column=col,sticky='w',columnspan=2)
        row+=1
        col=1
        #
        # stop and quit buttons
        #
        Button(self.eeu_win,text='Stop',command=self.stop_estimate_uncertainty).grid(row=row,column=col)
        col+=1
        Button(self.eeu_win,text='Quit',command=self.quit_estimate_uncertainty).grid(row=row,column=col)
        #
        # do the runs
        #
        import random
        self.stop_uncertainty_estimation = 0
        for i in range(no_runs.get()):
            if self.stop_uncertainty_estimation == 1:
                return
            #
            # update run no
            #
            current_run.set(i+1)

            #
            # mutate titration data
            #
            print 'mutating'
            self.titration_data.clear()

            for group in backup_titration_data:
                self.titration_data[group] = {}
                for ph in backup_titration_data[group]:
                    mut_ph = ph + random.uniform(-pH_uncertainty.get(),pH_uncertainty.get())
                    mut_charge = backup_titration_data[group][ph] + random.uniform(-titration_uncertainty.get(),titration_uncertainty.get())
                    self.titration_data[group][mut_ph] = mut_charge
            #
            # mutate activity data (if present)
            #
            if getattr(self,'activity_data',None):
                self.activity_data.clear()

                for ph in backup_activity_data:
                    mut_ph = ph + random.uniform(-pH_uncertainty.get(),pH_uncertainty.get())
                    mut_activity = backup_activity_data[ph] + random.uniform(-activity_uncertainty.get(),activity_uncertainty.get())
                    self.activity_data[mut_ph] = mut_activity
            #
            # Update windows
            #
            self.update_pkasystem_curves()
            self.win.update()
            #
            # Do the fit
            #
            print 'fitting'
            if getattr(self,'activity_data',None):
                self.include_tcs_in_fit = 1
                self.fit_to_ph_activity_profile(silent=1)
                self.include_tcs_in_fit = 0
            else:
                self.fit(silent=1)

            #
            # Save parameters and errors
            #
            count=0
            for group in self.groups:
                parameters[count][i] = self.groups[group].intpka.get()
                count +=1
                for group2 in self.groups[group].intenes.keys():
                    parameters[count][i] = self.groups[group].intenes[group2].get()
                    count +=1
            #
            # calculate average parameters and standard deviations
            #
            print 'calc parameters'
            for k in range(len(parameters)):
                avr = 0.0
                for j in range(i+1):
                    avr += parameters[k][j]
                avr=avr/(i+1)
                gui_parameters[k].set('%5.4f'%avr)
            for k in range(len(parameters)):
                temp = parameters[k][:i+1]
                print 'std for parameter',temp,'is', numpy.std(temp)
                gui_parameters_standard_deviations[k].set('%5.4f'%numpy.std(temp))


            abs_sum_error,avr_abs_error,scaled_abs_error,avr_HHd_scaled_abs_error = self.evaluate_fit(make_win=0)
            errors[0][i] = abs_sum_error
            errors[1][i] = avr_abs_error
            errors[2][i] = scaled_abs_error
            errors[3][i] = avr_HHd_scaled_abs_error
            #
            # Calculate average errors
            #
            print 'calc errors'
            for k in range(4):
                temp = 0.0
                for j in range(i+1):
                    temp+=errors[k][j]

                temp=temp/(i+1)
                gui_errors[k].set('%5.4f'%temp)
            #
            # Calculate standard deviations of errors
            #
            for k in range(4):
                temp = errors[k][:i+1]
                print 'std for errors of ',temp,' is ', numpy.std(temp)
                gui_standard_deviations[k].set('%5.4f'%numpy.std(temp))
            #
            # Reset parameters
            #
            count=0
            for group in self.groups:
                self.groups[group].intpka.set(4.0)
                for group2 in self.groups[group].intenes.keys():
                    self.groups[group].intenes[group2].set(0.0)



            print 'parameters',parameters
            print 'errors',errors

        #
        # recover original exp curves
        #
        self.titration_data = backup_titration_data.copy()
        self.update_pkasystem_curves()
        self.win.update()

        return

    def stop_estimate_uncertainty(self):
        self.stop_uncertainty_estimation = 1
        self.keep_running_activity = 0
        self.stop_fit()
        return

    def quit_estimate_uncertainty(self):
        self.stop_estimate_uncertainty()
        self.eeu_win.destroy()
        return

#
# -----
#

class Dialog:

    def __init__(self, parentwindow, acid_bases, do_scaling, theoretical_curves, groups):

        #
        # Make window
        #
        self.top = Toplevel(parentwindow)
        import Window_geom
        Window_geom.set_geometry(parentwindow,self.top)
        self.do_scaling = do_scaling
        self.acid_bases = acid_bases
        self.groups = groups
        r = 0

        self.headline = StringVar()
        self.headline.set('Chose mapping of experimental curves to theoretical curves')
        Label(self.top, textvariable=self.headline).grid(row=r,columnspan=4,sticky='we')
        r+=1
        Label(self.top, text='Groups in file').grid(row=r,column=0)
        Label(self.top, text='Acid/base').grid(row=r,column=1,columnspan=2)
        Label(self.top, text='Theoretical group').grid(row=r,column=3)
        #
        # Append gui name mapping
        #
        for ab in acid_bases:
            var=StringVar()
            var.set('None')
            ab.append(var)
        #
        for i in range(len(acid_bases)):
            r+=1
            ab = acid_bases[i]
            Label(self.top, text=ab[0]).grid(row=r,column=0)
            a = Radiobutton(self.top,text="Acid", variable=ab[1], value='acid')
            a.select()
            a.grid(row=r,column=1)
            Radiobutton(self.top, text="Base", variable=ab[1], value='base').grid(row=r,column=2)

            ab[2].set('None')

            #
            # exp to theoretical curve mapping menu
            #
            mb = Menubutton(self.top, textvariable=ab[3])
            menu = Menu(mb,tearoff=0,postcommand=self.update_menu)
            mb['menu'] = menu
            #
            #
            #theoret_keys=theoretical_curves.keys()
            #theoret_keys.sort()
            for i in range(len(theoretical_curves.keys())):
                j = int(theoretical_curves.keys()[i][1:5]) #This is probably a stupid way to find the group index
                menu.add_radiobutton(label=groups[j].name.get(),
                                     variable=ab[2],
                                     value=theoretical_curves.keys()[i],
                                     indicatoron=1,
                                     command=self.update_menu)

            menu.add_radiobutton(label='None', variable=ab[2], value='None', indicatoron=1,command=self.update_menu)
            mb.grid(row=r,column=3)



        c = Checkbutton(self.top, text='Scale loaded titration data between 0 and 1/-1?', variable=do_scaling)
        c.select()
        c.grid(row=r+1,column=0,columnspan=4)

        ok = Button(self.top, text="Load curves", command=self.ok).grid(row=r+2,column=0)

        ok = Button(self.top, text="Cancel", command=self.top.destroy).grid(row=r+2,column=1)
        return

    #
    # ----
    #

    def update_menu(self):
        """Updates gui names in menus"""
        for ab in self.acid_bases:
            if not ab[2].get() == 'None':
                j = int(ab[2].get()[1:5])
                ab[3].set(self.groups[j].name.get())
            else:
                ab[3].set('None')
        return

    #
    # -----
    #

    def ok(self):
        #
        # Check that each group has been assigned one or no curve
        #
        for i in range(len(self.acid_bases)):
            test = self.acid_bases[i][2].get()
            count = 0
            for i in range(len(self.acid_bases)):
                if self.acid_bases[i][2].get() == test:
                    if not self.acid_bases[i][2].get() == 'None':
                        count = count+1

            if not count < 2:
                self.headline.set('Import data - choose one group for each exp group')
                return
        #
        # Check that acids and bases match
        #
#        for i in range(len(self.acid_bases)):
#            if self.acid_bases[i][1].get() == 'acid' and self.acid_bases[i][2].get()[-3:] == 'ARG':
#                self.headline.set('Import data - make sure that acids and bases match')
#                return
#            if self.acid_bases[i][1].get() == 'base' and self.acid_bases[i][2].get()[-3:] == 'ASP':
#                self.headline.set('Import data - make sure that acids and bases match')
#                return


        #print '0 1 2'
        #for ab in self.acid_bases:
        #    print ab[0],' ',ab[1].get(),' ',ab[2].get()

        self.top.destroy()

class CS_dialog:

    def __init__(self,parent):
        self.CS=Toplevel()
        self.CS.title('Combinatorial scan')
        self.CS.geometry('+200+200')

##         include_pkas = []
##         include_ints = []
##         for group in parent.groups:
##             temp = IntVar()
##             temp.set(1)
##             include_pkas.append(temp)
##             for group2 in parent.groups[group].intenes.keys():
##                 temp = IntVar()
##                 temp.set(1)
##                 include_ints.append(temp)


        row = 0
        col = 0
        pka_count = 0
        int_count = 0
       ##  L = Label(self.CS,text='Which parameters do you want to include in the combinatorial scan?')
##         L.grid(row=row,column=col,columnspan=2)
##         row=row+1
##         for group in parent.groups:
##             c = Checkbutton(self.CS, text='pka',onvalue=1, offvalue =0, variable=include_pkas[pka_count])
##             c.grid(row=row,column=col)
##             col = col+1
##             pka_count = pka_count+1
##             for group2 in parent.groups[group].intenes.keys():
##                 c = Checkbutton(self.CS, text='int',onvalue=1, offvalue =0, variable=include_ints[int_count])
##                 c.grid(row=row,column=col)
##                 col = col+1
##                 int_count=int_count+1
##             row=row+1
##             col=0

        L=Label(self.CS,text='Chose boundaries and steps for the scan')
        L.grid(row=row,column=col,columnspan=2)
        row=row+1

        L = Label(self.CS, text='Min. pKa')
        L.grid(row=row,column=col)
        col=col+1
        E = Entry(self.CS, textvariable=parent.pka_min)
        E.grid(row=row,column=col)
        col=col+1
        L = Label(self.CS, text='Max. pKa')
        L.grid(row=row,column=col)
        col=col+1
        E = Entry(self.CS, textvariable=parent.pka_max)
        E.grid(row=row,column=col)
        col=col+1
        L = Label(self.CS, text='pKa step')
        L.grid(row=row,column=col)
        col=col+1
        E = Entry(self.CS, textvariable=parent.pka_step)
        E.grid(row=row,column=col)
        col=0
        row=row+1
        L = Label(self.CS, text='Min. int.')
        L.grid(row=row,column=col)
        col=col+1
        E = Entry(self.CS, textvariable=parent.int_min)
        E.grid(row=row,column=col)
        col=col+1
        L = Label(self.CS, text='Max. int.')
        L.grid(row=row,column=col)
        col=col+1
        E = Entry(self.CS, textvariable=parent.int_max)
        E.grid(row=row,column=col)
        col=col+1
        L = Label(self.CS, text='Int step')
        L.grid(row=row,column=col)
        col=col+1
        E = Entry(self.CS, textvariable=parent.int_step)
        E.grid(row=row,column=col)
        col=col=0
        row=row+1

        L = Label(self.CS, text='Include parameter sets in \'close.txt\' based on threshold:')
        L.grid(row=row,column=col)
        col=col+1
        E = Entry(self.CS, textvariable=parent.close_threshold)
        E.grid(row=row,column=col)
        col=col=0
        row=row+1


        B = Button(self.CS,text='OK',command=self.ok)
        B.grid(row=row,column=col)

        return

    def ok(self):
        self.CS.destroy()


class Show_close_dialog:

    def __init__(self, parent):
        self.parent=parent
        self.win=Toplevel()
        self.win.title('Close parameter sets')
        self.win.geometry('+400+200')


        row = 1
        col = 1

        #
        # set up counter
        #
        self.counter = Pmw.Counter(self.win,
                labelpos = 'w',
                label_text = 'Close parameter set #:',
                orient = 'horizontal',
                entry_width = 10,
                entryfield_value = 1,
                entryfield_command = self.update,
                entryfield_validate = self._counter_validate)
        self.counter.grid(row=row,column=col)
        Button(self.win, text='Update',command=self.update).grid(row=row,column=col+1)
        Button(self.win, text='Print to txt file',command=self.print_close_configs).grid(row=row,column=col+2)
        #
        # set up parameters
        #
        count=0

        self.gui_vars = []
        self.gui_errors = []

        for group in self.parent.groups:
            row+=1
            Label(self.win,text=self.parent.groups[group].name.get()).grid(row=row,column=col)
            col+=1
            pka = StringVar()
            Label(self.win,textvariable=pka).grid(row=row,column=col,padx=10)
            self.gui_vars.append(pka)
            for group2 in self.parent.groups[group].intenes.keys():
                col+=1
                intes = StringVar()
                Label(self.win,textvariable=intes).grid(row=row,column=col,padx=10)
                self.gui_vars.append(intes)
            col = 1
        #
        # set error display
        #
        row+=1
        Label(self.win, text='Error: ').grid(row=row,column=col)
        col += 1
        error = StringVar()
        Label(self.win, textvariable=error).grid(row=row,column=col)
        self.gui_vars.append(error)

        col =1
        row+=1
        Label(self.win, text='Dist. in parameter space: ').grid(row=row,column=col)
        col += 1
        pd = StringVar()
        Label(self.win, textvariable=pd).grid(row=row,column=col)
        self.gui_vars.append(pd)

        col =1
        row+=1
        Label(self.win, text='Error diff.: ').grid(row=row,column=col)
        col += 1
        ed = StringVar()
        Label(self.win, textvariable=ed).grid(row=row,column=col)
        self.gui_errors.append(ed)

        col =1
        row+=1
        Label(self.win, text='Rel.error diff.: ').grid(row=row,column=col)
        col += 1
        red = StringVar()
        Label(self.win, textvariable=red).grid(row=row,column=col)
        self.gui_errors.append(red)
        col += 1
        Label(self.win, text='%').grid(row=row,column=col)

        col =1
        row+=1
        Label(self.win, text='Rel.scaled error diff. (cmp to #1): ').grid(row=row,column=col)
        col += 1
        rsed = StringVar()
        Label(self.win, textvariable=rsed).grid(row=row,column=col)
        self.gui_errors.append(rsed)
        col += 1
        Label(self.win, text='%').grid(row=row,column=col)


        col =1
        row+=1
        Label(self.win, text='Avr. HH dev.scaled error: ').grid(row=row,column=col)
        col += 1
        aHHd = StringVar()
        Label(self.win, textvariable=aHHd).grid(row=row,column=col)
        self.gui_errors.append(aHHd)

        col =1
        row+=1
        Label(self.win, text='Rel. diff. avr. HH dev.scaled error: ').grid(row=row,column=col)
        col += 1
        rdaHHd = StringVar()
        Label(self.win, textvariable=rdaHHd).grid(row=row,column=col)
        self.gui_errors.append(rdaHHd)
        col += 1
        Label(self.win, text='%').grid(row=row,column=col)

        exp_uncertainty, avr_exp_uncertainty =titration_curve(self.parent.titration_data).experimental_uncertainty()

        col =1
        row+=1
        Label(self.win, text='Estimated exp. uncertainty of curves: ').grid(row=row,column=col)
        col += 1
        Label(self.win, text='%3.5f'%exp_uncertainty).grid(row=row,column=col)

        col =1
        row+=1
        Label(self.win, text='Estimated exp. uncertainty per point: ').grid(row=row,column=col)
        col += 1
        Label(self.win, text='%3.5f'%avr_exp_uncertainty).grid(row=row,column=col)

        self.lowest_error = self.parent.close_vars[0][len(self.parent.close_vars[0])-2]
        self.lowest_scaled_error = -1

        self.update()
        return


    def update(self):
        """updates scd and parent window """

        no = int(self.counter.get())-1
        i = 0

        #
        # update scales
        #

        for group in self.parent.groups:
            self.parent.groups[group].intpka.set(self.parent.close_vars[no][i])
            i+=1
            for group2 in self.parent.groups[group].intenes.keys():
                self.parent.groups[group].intenes[group2].set(self.parent.close_vars[no][i])
                i+=1


        #
        # update windows
        #
        self.parent.update_scales_from_fit()
        self.parent.titwin.update()

        #
        # find lowest scaled error the first time this method is called
        #
        if self.lowest_scaled_error == -1:
            self.lowest_scaled_error = self.parent.evaluate_fit(make_win=0)[2]
            self.lowest_HHd_scaled_error = self.parent.evaluate_fit(make_win=0)[3]


        #
        # calc scaled error
        #
        abs_sum_error,avr_abs_error,scaled_abs_error,avr_HHd_scaled_abs_error = self.parent.evaluate_fit(make_win=0)
        rel_diff_scaled_error = (scaled_abs_error- self.lowest_scaled_error)/self.lowest_scaled_error
        rel_diff_HHd_scaled_error = (avr_HHd_scaled_abs_error- self.lowest_HHd_scaled_error)/self.lowest_HHd_scaled_error

        #
        # update gui variables
        #
        i=0
        while i < len(self.parent.close_vars[no]):
            self.gui_vars[i].set(self.parent.close_vars[no][i])
            i+=1

        self.gui_errors[0].set(self.parent.close_vars[no][i-2] - self.lowest_error)
        self.gui_errors[1].set('%2.2f'%((self.parent.close_vars[no][i-2] - self.lowest_error)/self.lowest_error*100))
        self.gui_errors[2].set('%2.2f'%(rel_diff_scaled_error*100))
        self.gui_errors[3].set('%2.6f'%(avr_HHd_scaled_abs_error))
        self.gui_errors[4].set('%2.2f'%(rel_diff_HHd_scaled_error*100))

        return

    def print_close_configs(self):
        no = 1

        out = open('configs.gnu','w')
        count = 1

        out.write('# count    [parameters]        Error    para_d  E_diff  rE_dif rsE_dif  HHd_s_e   r_HHd\n')

        while self._counter_validate(self.counter.get()) == 1:
            self.update()
            print self.counter.get()
            line = '%7d'%count
            for g in self.gui_vars:
                line += '%7s '%g.get()

            for g in self.gui_errors:
                line += '%7s '%g.get()

            line += '\n'

            self.counter.increment()
            count += 1

            out.write(line)

        return



    def _counter_validate(self,integer):
        res = -1
        if int(integer) >= 1 and int(integer) <= len(self.parent.close_vars):
            res = 1

        return res

