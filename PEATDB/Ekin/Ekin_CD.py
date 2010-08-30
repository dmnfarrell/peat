# $Id: Ekin_CD.py 3535 2008-10-11 12:04:31Z farrell $
#
# This file is part of Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
from Tkinter import *
import Pmw

class CD_handler:

    def import_CD_tempscan(self):
        """Import a temperature scan from a Jasco CD spec"""
        import tkFileDialog, os
        filename=tkFileDialog.askopenfilename(defaultextension='.txt',initialdir=os.getcwd(),
                                              filetypes=[("Jasco txt","*.txt"),
                                                         ("All files","*.*")],
                                              parent=self.ekin_win)
        if not filename:
            return
        #
        # If we got a filename then read the stuff
        #
        import os
        if not os.path.isfile(filename):
            import tkMessageBox
            tkMessageBox.showwarning('File not found',
                                     'I could not find %s' %filename,
                                     parent=self.ekin_win)
            return
        #
        # Open and read file
        #
        fd=open(filename)
        lines=fd.readlines()
        fd.close()
        #
        # Parse file
        #
        # Format is:
        # 00 TITLE	0.6mg/ml hewl sigma
        # 01 DATA TYPE
        # 02 ORIGIN	JASCO
        # 03 OWNER
        # 04 DATE	06/10/17
        # 05 TIME	17:59:52
        # 06 SPECTROMETER/DATA SYSTEM	J-810
        # 07 DELTAT	0
        # 08 TUNITS	Temperature [C]
        # 09 FIRSTT	25
        # 10 LASTT	35
        # 11 NPOINTST	3
        # 12 DELTAX	-0.5
        # 13 XUNITS	Wavelength [nm]
        # 14 YUNITS	CD[mdeg]
        # 15 YUNITS2	HT[V]
        # 16 FIRSTX	350
        # 17 LASTX	250
        # 18 NPOINTS	201
        # 19 XYDATA1	25	30	35
	#
        # For a file with scans of three temperatures
        #
        data={}
        data['TITLE']=lines[0]
        data['ORIGIN']=lines[1]
        data['DATE']=lines[4]+lines[5]
        data['SPEC']=lines[6]
        #
        # Get the range of temps
        #
        t_temps=lines[19].split()[1:]
        temps=[]
        for t in t_temps:
            temps.append(float(t))
        #
        # Get the number of wavelengths
        #
        lambda_points=int(lines[18].split()[-1])
        #
        # Read the data
        #
        raw_data={}
        for x in range(lambda_points):
            line=lines[x+20]
            line=line.split()
            wavlen=float(line[0])
            count=0
            for temp in temps:
                count=count+1
                mdeg=float(line[count])
                if not raw_data.has_key(temp):
                    raw_data[temp]=[]
                raw_data[temp].append([wavlen,mdeg])
        #
        # Insert the tabs
        #
        temp1=temps[0]
        dps=len(raw_data[temp1])
        for temp in temps:
            newtab_num=self.add_datalist('CD(T'+str(temp)+')',data_points=dps)
            count=0
            for wavlen,mdeg in raw_data[temp]:
                self.data[newtab_num][0][count]['var'].set(wavlen)
                self.data[newtab_num][1][count]['var'].set(mdeg)
                count=count+1
        self.mode_var.set(5)
        self.update_mode()
        self.redraw_graph(fitting=1)

        return

    #
    # -----
    #

    def insert_CD_temp_datatab(self):
        """Insert datatab for fitting the temperature dependence of CD data"""
        #
        # Get the wavelengths
        #
        thistab = self.nb.index(self.nb.getcurselection())
        wavlens_dps=self.data[thistab][0].keys()
        wavlens=[]
        for dp in wavlens_dps:
            if dp!='label' and dp!='label_widget':
                wavlens.append(self.data[thistab][0][dp]['var'].get())
        wavlens.sort()
        self.wavlen=StringVar()
        #
        # Open window for selecting the wavelength
        #
        # first check if there is any imported data - bad hack
        if self.tabnames[thistab]=='data' or self.tabnames[thistab]=='Temp-dependence':
            import tkMessageBox
            tkMessageBox.showwarning("Import Data First",
                                     "Please import data or choose the first tab.")
            return
        self.wav_win=Toplevel()
        self.wav_win.title('Select wavelength')
        #
        # Set the geometry

        #PEATDB.PEAT_window.set_geometry(self.ekin_win,self.wav_win)
        self.wav_win.protocol("WM_DELETE_WINDOW",self.wav_win_cancel)

        #
        # Entry field for wavelength instead - asked for by Una
        #
        wav_entry=Entry(self.wav_win,textvariable=self.wavlen,width=10,font='Courier 14',justify='center')
        wav_entry.grid(row=1,column=0,sticky='news',padx=3,pady=3,columnspan=2)
        # Label
        lb=Label(self.wav_win,text='Choose Wavelength',font='Arial 16',justify='center')
        lb.grid(row=0,column=0,sticky='news',padx=3,pady=3,columnspan=2)
        #
        increment = float(wavlens[1])-float(wavlens[0])
        self.wavscale=Scale(self.wav_win,resolution=increment,variable=self.wavlen,orient='horizontal',showvalue=0,
                            command=self.insert_wavdata,from_=wavlens[0],to=wavlens[-1])
        self.wavscale.grid(row=2,column=0,sticky='news',padx=3,pady=3,columnspan=2)
        self.wavlen.set(wavlens[0])
        #
        # Manually specify min and max mdeg
        #
        Label(self.wav_win,text='min mdeg').grid(row=3,column=0)
        self.min_mdeg=DoubleVar()
        self.min_entry=Entry(self.wav_win,textvariable=self.min_mdeg,width=10,font='Courier 14',justify='center')
        self.min_entry.grid(row=4,column=0,sticky='news',padx=3,pady=3)
        self.min_entry.bind('<KeyPress>',self.update_CDmin)
        #
        Label(self.wav_win,text='max mdeg').grid(row=3,column=1,padx=3,pady=3)
        self.max_mdeg=DoubleVar()
        self.max_entry=Entry(self.wav_win,textvariable=self.max_mdeg,width=10,font='Courier 14',justify='center')
        self.max_entry.grid(row=4,column=1,sticky='news',padx=3,pady=3)
        self.max_entry.bind('<KeyPress>',self.update_CDmax)
        #
        Button(self.wav_win,text='Cancel',command=self.wav_win_cancel).grid(row=5,column=0,sticky='news',padx=3,pady=4,columnspan=1)
        Button(self.wav_win,text='Done',command=self.wav_win_close).grid(row=5,column=1,sticky='news',padx=3,pady=4,columnspan=1)
        #
        # Add the new datatab
        #
        dps=len(self.data.keys())
        names=self.tabnames.values()
        #if already created temp dep. tab do not do it again
        if not 'Temp-dependence' in names:
            self.temp_dep_tab=self.add_datalist('Temp-dependence',data_points=dps)
        else:
            thistab_num=len(self.datatabs)-1
            self.nb.selectpage(self.tabnames[thistab_num])
        #
        # Plot the data
        #
        self.CD_min='auto'
        self.CD_max='auto'
        self.insert_wavdata()
        #
        # Update the fitting model
        #
        self.FIT.update_model('CD - Thermal denaturation')
        self.reprint_parameters()
        self.redraw_graph()
        return

    #
    # ----
    #

    def wav_win_close(self,event=None):
        """
        Scale the CD data between 0 and 1, and destroy wavelength selection window
        """
        #
        # Get data
        #
        chosen_wavelength=self.wavlen.get()
        #
        data=[]
        for tabnum in self.tabnames:
            if self.tabnames[tabnum][:4]=='CD(T':
                #
                # Grab data from here
                #
                temp=float(self.tabnames[tabnum][4:-1])
                for dp in self.data[tabnum][0].keys():
                    if dp!='label' and dp!='label_widget':
                        wavlen=self.data[tabnum][0][dp]['var'].get()
                        if wavlen==chosen_wavelength:
                            data.append([temp,wavlen,self.data[tabnum][1][dp]['var'].get()])

        #
        # Normalise
        #
        i=0
        for temp,wavlen,mdeg in data:
            mdeg=float(mdeg)+abs(self.min_mdeg.get())
            mdeg=float(mdeg)/(self.max_mdeg.get()-self.min_mdeg.get())
            data[i][2]=mdeg
            i=i+1
        #
        # Change fitting model
        #
        self.FIT.update_model('Thermal denaturation')
        self.reprint_parameters()
        self.redraw_graph()
        #
        # Insert normalised data
        #
        count=0
        for temp,wavlen,mdeg in data:
            self.data[self.temp_dep_tab][0][count]['var'].set(temp)
            self.data[self.temp_dep_tab][1][count]['var'].set(mdeg)
            count=count+1
        self.redraw_graph(fitting=1)
        #
        # Close window
        #
        if self.wav_win:
            self.wav_win.destroy()
        return

    #
    # ----
    #

    def wav_win_cancel(self,event=None):
        """
        Destroy wavelength selection window and delete the new tab
        """
        if self.wav_win:
            self.delete_datatab()
            self.wav_win.destroy()
        return

    #
    # -----
    #

    def update_CDmin(self,event=None):
        """Update the min value for CD signal"""
        self.CD_min='manual'
        self.insert_wavdata()
        return

    #
    # -----
    #

    def update_CDmax(self,event=None):
        """Update the max value for CD signal"""
        self.CD_max='manual'
        self.insert_wavdata()
        return

    #
    # ------
    #

    def insert_wavdata(self,junk=None):
        """Get the data for the wavelength that was selected and insert it in the sheet"""
        #
        # Loop over all datatabs and get the ones that are CD(T) data
        #
        chosen_wavelength=self.wavlen.get()
        #
        data=[]
        for tabnum in self.tabnames:
            if self.tabnames[tabnum][:4]=='CD(T':
                #
                # Grab data from here
                #
                temp=float(self.tabnames[tabnum][4:-1])
                for dp in self.data[tabnum][0].keys():
                    if dp!='label' and dp!='label_widget':
                        wavlen=self.data[tabnum][0][dp]['var'].get()
                        if wavlen==chosen_wavelength:
                            data.append([temp,wavlen,self.data[tabnum][1][dp]['var'].get()])

        #
        # find min and max of CD signal
        #
        min_val=0
        max_val=-1000
        i=0
        for temp,wavlen,mdeg in data:
            mdeg=float(mdeg)
            if mdeg>max_val:
                max_val=mdeg
            if mdeg<min_val:
                min_val=mdeg
            i=i+1
        #
        # Did the user specify one value or the other?
        #
        if self.CD_min=='auto':
            self.min_mdeg.set(min_val)
        #
        if self.CD_max=='auto':
            self.max_mdeg.set(max_val)
        #
        # Insert data from the wavelength
        #
        count=0
        for temp,wavlen,mdeg in data:
            self.data[self.temp_dep_tab][0][count]['var'].set(temp)
            self.data[self.temp_dep_tab][1][count]['var'].set(mdeg)
            count=count+1
        self.redraw_graph(fitting=1)
        return






