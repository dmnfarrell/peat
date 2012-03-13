#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
# $Id: EM_effect.py 4884 2009-06-09 06:29:41Z nielsen $
#
# This file is part of the pKaTool package
# Copyright Jens Erik Nielsen University College Dublin 2003-
# All rights reserved
#

def length(vector):
    # This function returns the length of vector
    import math
    sum=0.0
    for value in vector:
        sum=sum+math.pow(value,2)
    return math.sqrt(sum)


from Tkinter import *
import Pmw
import pKaTool.pKa_base

#
# Constants
#
ppm_au=1.9447E-18
NSP={'N':977,'H':90}

class EM_effect(Frame,pKaTool.pKa_base.pKa_base):

    def __init__(self,parent=None,Yasara=None,display=1,PI=None):
        """Init the class"""
        self.parent_application=None
        self.ID='EM_effect'
        if parent:
            self.parent_application=parent
        if PI:
            self.PI=PI
        if not display:
            return
        #
        # yasara?
        #
        self.yasara=None
        if Yasara:
            self.yasara=Yasara
        #
        # Savedir?
        #
        import os
        if os.environ.has_key('HOME'):
            self.savedir=os.getenv('HOME')
        #
        # Open the window
        #
        if not self.parent_application:
            Frame.__init__(self)
            self.EM_window=self.master
        else:
            self.EM_window=Toplevel()
        self.EM_window.title('Titration curve simulation')
        self.EM_window.geometry('+200+50')
        #
        # Draw the canvas
        #
        self.draw_canvas()
        self.balloon=Pmw.Balloon(self.master)
        #
        # Pulldown menu
        #
        self.menu=Menu(self.EM_window)
        #
        # File menu
        #
        self.file_menu=Menu(self.menu,tearoff=0)
        self.file_menu.add_command(label='Load PDB file',command=self.load_pdb)
        self.file_menu.add_command(label='Load span data for one group',command=self.load_csv)
        self.file_menu.add_command(label='Load full titration data',command=self.load_full_titdata)
        self.file_menu.add_command(label='Exit',command=self.quit)
        self.menu.add_cascade(label='File',menu=self.file_menu)
        #
        # Command menu
        #
        self.com_menu=Menu(self.menu,tearoff=0)
        self.com_menu.add_command(label='Select/Define titratable group',command=self.def_titgrp)
        self.com_menu.add_command(label='Print graph to file',command=self.print_chemshift_to_file)
        self.menu.add_cascade(label='Command',menu=self.com_menu)
        #
        # Configure the menu
        #
        self.EM_window.config(menu=self.menu)
        #
        # Add listbox to the bottom of graphic window
        #
        self.LB_added_groups=Listbox(self.EM_window,width=40,bg='lightblue')
        self.LB_added_groups.grid(row=1,column=0,columnspan=3,stick='news')
        #
        # Listbox for experimental data
        #
        self.LB_exp_data=Listbox(self.EM_window,width=40,bg='lightgreen')
        self.LB_exp_data.grid(row=1,column=10,columnspan=3,sticky='news')
        #
        # Buttons for controlling the calc for each group
        #
        self.colorbutton=Button(self.EM_window,text='Colour',state='disabled')
        self.colorbutton.grid(row=1,column=4)
        self.diel_var=DoubleVar()
        self.diel_var.set(8)
        self.eff_diel_scale=Scale(self.EM_window,
                                  from_=1.0,to=30.0,
                                  resolution=0.1,
                                  variable=self.diel_var,
                                  orient='horizontal',
                                  state='disabled',
                                  label='Effective eps')
        self.eff_diel_scale.grid(row=1,column=5)
        self.removebutton=Button(self.EM_window,text='Remove',state='disabled')
        self.removebutton.grid(row=1,column=6)
        #
        # Display checkbutton
        #
        self.display_var=IntVar()
        self.display_var.set(1)
        self.display_button=Checkbutton(self.EM_window,onvalue=1,offvalue=0,variable=self.display_var,state='disabled')
        self.display_button.grid(row=1,column=7)
        self.group_properties={}
        #
        # Fit button
        #
        self.ele_calc=StringVar()
        self.ele_calc.set('Coulomb')
        self.elecalc_button=Menubutton(self.EM_window,textvariable=self.ele_calc,relief=RAISED,state=DISABLED)
        self.elecalc_menu=Menu(self.elecalc_button,tearoff=0)
        self.elecalc_button['menu']=self.elecalc_menu
        #
        # Add the two options
        #
        for method in ['Coulomb','PBE']:
            self.elecalc_menu.add_radiobutton(label=method,
                                              variable=self.ele_calc,
                                              value=method,
                                              indicatoron=1,
                                              command=self.update_EM_effect_prediction)
                                 
        self.elecalc_button.grid(row=1,column=8,sticky='ew')
    
        Button(self.EM_window,command=self.fit,text='Fit').grid(row=1,column=9)
        #
        # Button for colouring structure in yasara
        #
        def update_EM_with_yasara():
            self.update_EM_effect_prediction(do_yasara=True)
            return
        Button(self.EM_window,command=update_EM_with_yasara,text='Colour structure').grid(row=2,column=5)
        #
        # Binding for when a group is chosen
        #
        self.LB_added_groups.bind('<ButtonRelease-1>',self.enable_group_editing)
        return

    #
    # ----
    #

    def draw_canvas(self,protein_length=129):
        """(Re)Draw the canvas"""
        if hasattr(self,'tc'):
            self.tc.destroy()
        #
        # Draw the charge vs. chem shift graph
        #
        self.linewidth=2
        self.pHstart=0.00001
        self.pHend=20.0
        self.maxcrg=10.0
        self.mincrg=0.0
        self.titwidth=1200
        self.titheight=450
        self.tc=Canvas(self.EM_window,bd=5,bg='white',width=self.titwidth,
                       height=self.titheight,
                       scrollregion=(0,0,self.titwidth,self.titheight))
        self.tc.xview("moveto", 0)
        self.tc.yview("moveto", 0)
        self.tc.grid(row=0,column=0,columnspan=20)
        #
        # Axes
        #
        self.set_graph_size(x_max=protein_length,
                            x_min=0.0,
                            y_max=1.5,
                            y_min=-1.5,
                            y_tick_level=1.0,
                            x_tick_level=int(protein_length/25))
        self.draw_ordinates(self.tc,y_label='Chemical shift change',x_label='Residue')
        return

    #
    # -----
    #

    def quit(self):
        self.EM_window.destroy()
        return

    #
    # ----
    #
    def get_geometry(self,widget):
        """Get the geometry of a widget
        Return width,height,xorg,yorg"""
        widget.update_idletasks()
        txt=widget.winfo_geometry()
        width=int(txt.split('x')[0])
        rest=txt.split('x')[1]
        height=int(rest.split('+')[0])
        xorg=int(rest.split('+')[1])
        yorg=int(rest.split('+')[2])
        return width,height,xorg,yorg
        
    #
    # ----
    #

    def set_geometry(self,pwidget,widget):
        """Set the position of widget in the middle of pwidget"""
        w,h,x,y=self.get_geometry(pwidget)
        pwidget.update()
        widget.update()
        sw,sh,dummy,dummy2=self.get_geometry(widget)
        xoffset=int((w-sw)/2)
        yoffset=int((h-sh)/2)
        widget.geometry('+%d+%d' %(x+xoffset,y+yoffset))
        return

    #
    # -----
    #

    def enable_group_editing(self,event):
        """Enable editing of the properties of a group"""
        group_num=int(self.LB_added_groups.curselection()[0])
        selected_group=self.selected_groups[group_num][0]
        #
        #
        def get_new_colour():
            import tkColorChooser
            rgb,colorstring=tkColorChooser.askcolor(self.group_properties[selected_group]['colour'])
            self.group_properties[selected_group]['colour']=colorstring
            self.update_EM_effect_prediction()
            return
        self.colorbutton.configure(state=ACTIVE)
        self.colorbutton.configure(fg=self.group_properties[selected_group]['colour'])
        self.colorbutton.configure(command=get_new_colour)
        
        self.diel_var.set(self.group_properties[selected_group]['diel'])
        self.eff_diel_scale.configure(state=ACTIVE)
        def update_diel(event=None):
            """Update teh dielectric constant for a calc"""
            diel=self.diel_var.get()
            self.group_properties[selected_group]['diel']=diel
            self.update_EM_effect_prediction()
            return
        self.eff_diel_scale.configure(command=update_diel)
        self.removebutton.configure(state=ACTIVE)
        #
        # Display checkbutton
        #
        def update_on_off():
            self.group_properties[selected_group]['display']=self.display_var.get()
            self.update_EM_effect_prediction()
            return
        self.display_var.set(self.group_properties[selected_group]['display'])
        self.display_button.configure(state=ACTIVE)
        self.display_button.configure(command=update_on_off)
        #
        # Elecalc button
        #
        self.elecalc_button.configure(state=ACTIVE)
        return


    #
    # ----
    #

    def load_pdb(self,filename=None):
        """Load a PDB file and build all the HNs"""
        import tkFileDialog, os
        if hasattr(self,'savedir'):
            savedir=self.savedir
        else:
            savedir=os.getcwd()
        if not filename:
            filename=tkFileDialog.askopenfilename(defaultextension="*.pdb",
                                                  initialdir=savedir,
                                                  filetypes=[("PDB file","*.pdb"),("All files","*.*")])
        if filename:
            if os.path.isfile(filename):
                wait_win=Toplevel()
                self.set_geometry(self.master,wait_win)
                wait_win.title('Please wait')
                Label(wait_win,text='Preparing PDB file').grid(row=0,column=0)
                statusvar=StringVar()
                statusvar.set('Reading PDB file')
                Label(wait_win,textvariable=statusvar).grid(row=1,column=0,sticky='news')
                self.master.update_idletasks()
                wait_win.update_idletasks()
                import Protool
                self.PI=Protool.structureIO()
                self.PI.readpdb(filename)
                self.PI.RemoveALT()
                #
                # Build all HNs
                #
                import Protool.hydrogens
                H=Protool.hydrogens.Hydrogens(self.PI)
                residues=self.PI.residues.keys()
                residues.sort()
                count=0
                for residue in residues:
                    statusvar.set('Building hydrogen #%d of %d' %(count,len(residues)))
                    wait_win.update_idletasks()
                    if self.PI.isaa(residue):
                        H.build_HN(residue)
                    count=count+1
                self.PI.Update()
                #
                # Count protein residues
                #
                count=0
                for residue in residues:
                    if self.PI.isaa(residue):
                        count=count+1
                wait_win.destroy()
                #
                # Redraw the canvas
                #
                self.draw_canvas(protein_length=count)
                #
                # Get the titratable groups
                #
                self.PI.get_titratable_groups()
                #
                # If we are started from yasara then load the pdb file
                #
                if self.yasara:
                    self.yasara.run('DelAll')
                    self.yasara.run('obj=LoadPDB %s' %filename)
                    self.yasara.run('Style Stick')
                    self.yasara.run('HideRes Hoh')
                    self.yasara.run('ColorAll Green')
                    self.yasara.run('HUD Off')
                    self.yasara.Finish()
            else:
                raise Exception,'PDB file not found: %s' %(filename)
        return

    #
    # ----
    #

    def def_titgrp(self):
        """Define a titratable group
        There are two steps
        1. Define the position
        2. Define the titrational behaviour
        
        """
        self.seldef_win=Toplevel()
        self.set_geometry(self.master,self.seldef_win)
        self.seldef_win.title('Select/Define titratable group')
        Label(self.seldef_win,text='Titratable group').grid(row=1,column=0)
        Button(self.seldef_win,text='Select from PDB file',command=self.select_titgrp_frompdb).grid(row=1,column=2)
        Button(self.seldef_win,text='Define XYZ manually',state='disabled').grid(row=1,column=3)
        
        #
        # Titration behaviour
        #
        Label(self.seldef_win,text='Titration curve').grid(row=10,column=0)
        Checkbutton(self.seldef_win,text='Copy from pKa_system',command=self.copy_from_pkasystem).grid(row=10,column=1)
        def model_unit_charge():
            self.charge=1.0
            self.sel_titgrp=StringVar()
            self.sel_titgrp.set('Unit charge')
            return
        self.unit_charge=IntVar()
        Checkbutton(self.seldef_win,text='Model as unit charge',variable=self.unit_charge,onvalue=1,offvalue=0,command=model_unit_charge).grid(row=10,column=2)
        model_unit_charge()
        #
        # Nitrogen or Hydrogen chemical shift
        #
        self.chemshift_type=IntVar()
        Label(self.seldef_win,text='Calculate chemical shift for').grid(row=11,column=0)
        Radiobutton(self.seldef_win,
                    text='Backbone amide nitrogen',
                    value=0,
                    variable=self.chemshift_type).grid(row=11,column=1)
        Radiobutton(self.seldef_win,
                    text='Backbone amide hydrogen',
                    value=1,
                    variable=self.chemshift_type).grid(row=11,column=2)
        self.chemshift_type.set(0)
        #
        # Add group button
        #
        def add_group():
            #
            # Add to listbox
            #
            if not hasattr(self,'selected_groups'):
                self.selected_groups=[]
            titgroups=self.PI.titratable_groups.keys()
            titgroups.sort()
            selection=self.LB_groups.curselection()
            selection=int(selection[0])
            
            #
            # Get the type (N or H), and construct the name from it
            #
            if self.chemshift_type.get()==0:
                name='%s - %s' %(self.titgroups_inPDB[selection],'N')
                self.group_properties[name]={'colour':'red','diel':8,'display':1,'atom':'N'}
            else:
                name='%s - %s' %(self.titgroups_inPDB[selection],'H')
                self.group_properties[name]={'colour':'red','diel':8,'display':1,'atom':'H'}
            #
            # Add to listbox
            #
            self.selected_groups.append([name,self.titgroups_inPDB[selection],self.sel_titgrp.get()])
            self.LB_added_groups.insert(END,str(self.selected_groups[-1]))
            #
            # Activate the callback from pKa_system 
            #
            if hasattr(self,'PS'):
                self.PS.add_callback(self.update_EM_effect_prediction)
            else:
                self.update_EM_effect_prediction()
            return
        Button(self.seldef_win,text='Add group',command=add_group).grid(row=13,column=0)
        return

    #
    # ----
    #

    def add_titgroup(self,name,groupname,model='Unit charge'):
        """Add a titgroup"""
        if not hasattr(self,'selected_groups'):
            self.selected_groups=[]
        self.selected_groups.append([name,groupname,model])
        self.LB_added_groups.insert(END,str(self.selected_groups[-1]))
        self.group_properties[name]={'colour':'red','diel':8,'display':1,'atom':'N'}
        self.update_EM_effect_prediction()
        return
        
        
    #
    # ------
    #
    
    def select_titgrp_frompdb(self,window=None):
        """Select a titratable group from the PDB file"""
        #
        # We can do this in any window
        #
        if not window:
            window=self.seldef_win
        #
        # Make the dialog
        #
 
        titgrps=self.PI.titratable_groups.keys()
        titgrps.sort()
        self.titgroups_inPDB=[]
        for res in titgrps:
            for titgrp in self.PI.titratable_groups[res]:
                name='%s:%s' %(res,titgrp['name'])
                self.titgroups_inPDB.append(name)
        #
        yscrollbar=Scrollbar(window,orient='vertical',width=10)
        yscrollbar.grid(row=2,column=1,rowspan=5,sticky='nws')
        self.LB_groups=Listbox(window,bg='white',selectmode=BROWSE,yscrollcommand= yscrollbar.set)
        self.LB_groups.grid(row=2,column=0,rowspan=5)
        yscrollbar.config(command=self.LB_groups.yview)
        for element in self.titgroups_inPDB:
            self.LB_groups.insert(END,element)
        return
    #
    # -----
    #
    
    def copy_from_pkasystem(self):
        if self.parent_application==None:
            import pKaTool.pKa_system
            self.PS=pKa_system.pKa_system(parent_application=self)
        elif self.parent_application.ID=='pKaTool':
            #
            # I must still start pKa_system
            #
            pass
        elif self.parent_application.ID=='pKa_system':
            #
            # Hurray - no action necessary
            #
            pass 
        #
        #
        #
        groups=self.PS.groups.keys()
        groups.sort()
        if not hasattr(self,'sel_titgrp'):
            self.sel_titgrp=StringVar()
            self.sel_titgrp.set(self.PS.groups[groups[0]].name.get())
            self.selgroup_button=Menubutton(self.seldef_win,textvariable=self.sel_titgrp,relief=RAISED,bg='white')
            self.selgroup_menu=Menu(self.selgroup_button,tearoff=0)
            self.selgroup_button['menu']=self.selgroup_menu
            #
            # Add the groups in pKa_system
            #
            for group in groups:
                self.selgroup_menu.add_radiobutton(label=self.PS.groups[group].name.get(),
                                                    variable=self.sel_titgrp,
                                                    value=self.PS.groups[group].name.get(),
                                                    indicatoron=1,
                                                    command=self.update_EM_effect_prediction)
            self.selgroup_button.grid(row=10,column=2)
        return
  
    #
    # ----
    #

    def update_EM_effect_prediction(self,do_yasara=None):
        """Recalculate all EM effects and redraw all"""
        #
        # Should we recolour yasara?
        #
        if self.yasara and do_yasara:
            do_yasara=True
        #
        # First delete the bars from last time
        #
        if not hasattr(self,'graph_bars'):
            self.graph_bars={}
        for bar in self.graph_bars.keys():
            self.tc.delete(bar)
        #
        # Loop over all groups we're looking at
        #
        yascom=[]
        for group_name,group_position,group_titration in self.selected_groups:
            #
            # Get the coordinate center of the charge that it titrating
            #
            print group_position
            resnumber=self.PI.resnum(group_position)
            group_type=self.PI.resname(group_position)
            group=self.PI.titgroups[group_type]
            xpos=0.0
            ypos=0.0
            zpos=0.0
            for atom in group['atoms']:
                xpos=xpos+self.PI.atoms[resnumber+':'+atom]['X']/float(len(group['atoms']))
                ypos=ypos+self.PI.atoms[resnumber+':'+atom]['Y']/float(len(group['atoms']))
                zpos=zpos+self.PI.atoms[resnumber+':'+atom]['Z']/float(len(group['atoms']))
            #
            # Get the properties
            #
            colour=self.group_properties[group_name]['colour']
            eff_diel=self.group_properties[group_name]['diel']
            charge=group['charge']
            display=self.group_properties[group_name]['display']
            if not display:
                continue

            atom_type=self.group_properties[group_name]['atom']
            #
            # Loop over all residues
            #
            residues=self.PI.residues.keys()
            residues.sort()
            import numpy
            #
            # Draw the new bars
            #
            for residue in residues:
                #
                # If not an AA, then skip
                #
                if not self.PI.isaa(residue):
                    continue
                #
                # Get the residue number
                #
                resnum=int(residue.split(':')[1])
                #
                # Get the value of the span
                #
                if self.ele_calc.get()=='Coulomb':
                    span,angle=self.get_dEF(residue,numpy.array([xpos,ypos,zpos]),charge,eff_diel,atom_type=atom_type)
                elif self.ele_calc.get()=='PBE':
                    span=self.get_span_from_PBE(residue,numpy.array([xpos,ypos,zpos]),charge,eff_diel,atom_type)
                    angle=0.0
                elif self.ele_calc.get()=='APBS':
                    span=self.get_span_from_PBE(residue,numpy.array([xpos,ypos,zpos]),charge,eff_diel,atom_type,APBS=True)
                else:
                    raise Exception,'Unknown method for calculating electric field'
                #
                # If we get None back then we could not calculate the value
                #
                if not span and not angle:
                    #
                    # Draw black circle
                    #
                    x,y=self.get_xy(resnum,0.2)
                    x1,y0=self.get_xy(resnum+1,-0.2)
                    handle=self.tc.create_oval(x,y,x1,y0,fill='black',stipple='warning')
                    self.balloon.tagbind(self.tc,handle,'Could not calculate span for %s' %residue)
                    #
                    # colour residue green in yasara
                    #
                    if do_yasara:
                        yascom.append('ColorRes %s,green' %(str(resnum)))
                else:
                    #
                    # value ok, plot it
                    #
                    x,y=self.get_xy(resnum,span)
                    x1,y0=self.get_xy(resnum+1,0)
                    handle=self.tc.create_rectangle(x,y,x1,y0,fill=colour)
                    self.balloon.tagbind(self.tc,handle,'%s, span: %5.3f, angle: %5.2f' %(residue,span,angle))
                    #
                    # colour residue in yasara
                    #
                    if do_yasara:
                        if span>0.0:
                            color=min(50*span,120)
                        else:
                            color=max(50*span,-120)
                        if color<0.0:
                            color=360+color
                        yascom.append('ColorRes(%s,%d)' %(str(resnum),int(color)))
                self.graph_bars[handle]=1
        if do_yasara:
            for com in yascom:
                self.yasara.run(com)
                print com
            self.yasara.Finish()
        return

    #
    # ----
    #

    def show_fit_prediction(self,xyz,deff):
        """Show the current fit prediction"""
        stop
        #
        # First delete the bars from last time
        #
        if not hasattr(self,'fitpred_bars'):
            self.fitpred_bars={}
        for bar in self.fitpred_bars.keys():
            self.tc.delete(bar)
        residues=self.PI.residues.keys()
        residues.sort()
        import numpy
        #
        # Draw the new bars
        #
        for residue in residues:
            if not self.PI.isaa(residue):
                continue
            #
            # Get the residue number
            #
            resnum=int(residue.split(':')[1])
            #
            # Get the value of the span
            #
            span,angle=self.get_dEF(residue,xyz,-1,deff,atom_type='N')
            #
            # If we get None back then we could not calculate the value
            #
            if not span and not angle:
                #
                # Draw black circle
                #
                x,y=self.get_xy(resnum,0.2)
                x1,y0=self.get_xy(resnum+1,-0.2)
                handle=self.tc.create_oval(x,y,x1,y0,fill='black',stipple='warning')
                self.balloon.tagbind(self.tc,handle,'Could not calculate span for %s' %residue)
            else:
                #
                # value ok, plot it
                #
                x,y=self.get_xy(resnum,span)
                x1,y0=self.get_xy(resnum+1,0)
                handle=self.tc.create_rectangle(x,y,x1-3,y0,fill='yellow',stipple='gray50')
                #self.balloon.tagbind(self.tc,handle,'%s, span: %5.3f, angle: %5.2f' %(residue,span,angle))
            self.fitpred_bars[handle]=1
        self.EM_window.update_idletasks()
        return
        
    #
    # ---- 
    #
    
    def show_fit_prediction2(self,values):
        """Show the current fit prediction"""
        #
        # First delete the bars from last time
        #
        if not hasattr(self,'fitpred_bars'):
            self.fitpred_bars={}
        for bar in self.fitpred_bars.keys():
            self.tc.delete(bar)

        import numpy
        #
        # Draw the new bars
        #
        for residue,value in values:
            if not self.PI.isaa(residue):
                continue
            #
            # Get the residue number
            #
            resnum=int(residue.split(':')[1])
            #
            # If we get None back then we could not calculate the value
            #
            if not value:
                #
                # Draw black circle
                #
                x,y=self.get_xy(resnum,0.2)
                x1,y0=self.get_xy(resnum+1,-0.2)
                handle=self.tc.create_oval(x,y,x1,y0,fill='black',stipple='warning')
                self.balloon.tagbind(self.tc,handle,'Could not calculate span for %s' %residue)
            else:
                #
                # value ok, plot it
                #
                x,y=self.get_xy(resnum,value)
                x1,y0=self.get_xy(resnum+1,0)
                handle=self.tc.create_rectangle(x,y,x1-3,y0,fill='yellow',stipple='gray50')
                #self.balloon.tagbind(self.tc,handle,'%s, span: %5.3f, angle: %5.2f' %(residue,span,angle))
            self.fitpred_bars[handle]=1
        self.EM_window.update_idletasks()
        return

    #
    # ----
    #

    def get_angle(self,residue,titpos,atom_type):
        """Get the cos(angle)and distance between a bond and a charge"""
        if atom_type=='N':
            try:
                prev_res=self.PI.PreviousResidue(residue)
            except:
                return None,None
            atoms=[residue+':N',prev_res+':C']
        elif atom_type=='H':
            atoms=[residue+':H',residue+':N']
        else:
            print 'Atom type %s unknown' %str(atom_type)
            raise Exception
        #
        # Check that both atoms are present
        #
        if not self.PI.atoms.has_key(atoms[0]) or not self.PI.atoms.has_key(atoms[1]):
            return None,None
        #
        # Get the bond vector
        #
        bond_vector=self.PI.GetPosition(atoms[1])-self.PI.GetPosition(atoms[0])
        #
        # Vector to charge
        #
        charge_vector=titpos-self.PI.GetPosition(atoms[0])
        #
        # Get angle
        #
        import numpy
        dp=numpy.dot(bond_vector,charge_vector)
        cos_angle=dp/(length(bond_vector)*length(charge_vector))
        #
        # Get the distance
        #
        dist=length(charge_vector)
        return cos_angle,dist

    #
    # ----
    #
        
    def get_dEF(self,residue,titpos,charge,deff,atom_type='N'):
        """When given a residue, calculate the backbone and possible side chain polarisation due to
        the titratable charge at position titpos. The titratable group has a charge = charge and
        interacts with the residue with an effective dielectric constant deff        """
        import numpy
        import os
        if not os.path.isdir(os.path.join(os.getcwd(),'Coulomb_potmap/')):
		    os.system('mkdir Coulomb_potmap')
		
        cos_angle,dist=self.get_angle(residue=residue,titpos=titpos,atom_type=atom_type)
        if not cos_angle or not dist:
            return None,None
        #
        # Calculate electric field
        #
        import math
        e=1.602E-19
        e0=8.85E-12
        #print deff,dist,cos_angle,charge
        #print residue,titpos
        E=1.0/(4*math.pi*e0*deff)*e/(dist*1E-10)**2*ppm_au*NSP[atom_type]*1E6*cos_angle*charge
        return E,math.degrees(math.acos(cos_angle))

    #
    # -----
    #

    def get_span_from_PBE(self,residue,titpos,charge,eff_diel,atom_type='N',APBS=None):
        """Get the span from Chresten's PBE solver or APBS"""
        #
        # Define the name of the run and see if we have it already
        #
        import os
        pot=None
        #potname=os.path.join(os.getcwd(),'%4.2f_%5.1f_%5.1f_%5.1f_eps_%4.1f.potmap' %(charge,titpos[0],titpos[1],titpos[2],eff_diel))
        if APBS:
			if not os.path.isdir(os.path.join(os.getcwd(),'APBS_potmap/')):
                            os.system('mkdir APBS_potmap')
			potname=os.path.join(os.getcwd(),'APBS_potmap/','%4.2f_%5.1f_%5.1f_%5.1f_eps_%4.1f.potmap' %(charge,titpos[0],titpos[1],titpos[2],eff_diel))
			potname=potname+'_APBS'
        else:
            self.potmaps={}
        #if not pot:
        #    #
        #    # Is the file present?
        #    #
        #    if os.path.isfile(potname):
        #        fd=open(potname#)
        #        import pickle
        #        pot=pickle.load(fd)
        #        fd.close()
        #        self.potmaps[potname]=pot.copy()
        #        print 'Loaded map from pickle file',potname
        #
        # If no potmap, then calculate it
        #
        if not pot:
            #
            # Get the pqr file
            #
            print 'Recalculating: %s' %potname
            import Protool.assign_parameters
            self.PQR=Protool.assign_parameters.assign_parameters(self.PI)
            self.PQR.clear_all_charges_and_radii()
            self.PQR.set_all_radii()
        if APBS:
            print 'Using APBS'
            self.PI. add_atom(uniqueid='Z:999:CHA',
                atomnumber=99999,atomname='CHA',    
                chainid='Z',residuename='DUM',residuenumber='999',
                xcoord=titpos[0],ycoord=titpos[1],zcoord=titpos[2],update=1,BFACTOR=None,OCCUPANCY=None,CHARGE=charge,RADIUS=1.55,tag=None)
            import Protool.APBS_interface as APBS
            A=APBS.APBS_calc(self.PI,pdie=eff_diel,sdie=80.0)
            pot=A.potentials.copy()
            self.PI.remove_atom('Z:999:CHA')
            #
            # Save the map
            #
            #import pickle
            #fd=open(potname,'w')
            #pickle.dump(pot,fd)
            #fd.close()
            #self.potmaps[potname]=pot.copy()
        #
        # Calculate difference in potential for the two atoms in question
        #
        import numpy
        if atom_type=='N':
            try:
                prev_res=self.PI.PreviousResidue(residue)
            except:
                return None
            atoms=[residue+':N',prev_res+':C']
        elif atom_type=='H':
            atoms=[residue+':H',residue+':N']
        #
        # Get potentials
        #
        if not APBS:
            value1=pbs.get_potential_at(self.PI.atoms[atoms[0]]['X'],
                                        self.PI.atoms[atoms[0]]['Y'],
                                        self.PI.atoms[atoms[0]]['Z'],
                                        method='linear_interpolation',potmap=pot)
            value2=pbs.get_potential_at(self.PI.atoms[atoms[1]]['X'],
                                        self.PI.atoms[atoms[1]]['Y'],
                                        self.PI.atoms[atoms[1]]['Z'],
                                        method='linear_interpolation',potmap=pot)
        else:
            value1=pot[atoms[0]]
            value2=pot[atoms[1]]
        #
        value=value2-value1
        value=value*0.025692*ppm_au*NSP[atom_type]*1E6
        return value/(self.PI.dist(atoms[0],atoms[1])*1E-10)
        '''cos_angle,dist=self.get_angle(residue=residue,titpos=titpos,atom_type=atom_type)
        value=value1*0.025692*ppm_au*NSP[atom_type]*1E6*cos_angle/(dist*1E-10)
        return value'''
    #
    # -----
    #

    def print_chemshift_to_file(self):
        """Print the chemical shift canvas to file"""
        import sys, os
        import tkFileDialog, os
        if not hasattr(self,'lastdir'):
            self.lastdir=os.getcwd()
        filename=tkFileDialog.asksaveasfilename(defaultextension='.ps',
                                                initialdir=self.lastdir,
                                                filetypes=[("Postscript files","*.ps"),("All files","*.*")])
        if filename:
            import os
            self.lastdir=os.path.split(filename)[0]
            if filename[-3:]!='.ps':
                filename=filename+'.ps'
            self.tc.postscript(colormode='color',file=filename)
        return

    #
    # ----
    #

    def load_csv(self,filename=None,group=None):
        """Open a file with chem shift data"""
        self.exp_group_selected=group
        #
        # Make sure we have the expdata array
        #
        if not hasattr(self,'expdata'):
            self.expdata={}
        
        #
        # Load data
        #
        import tkFileDialog, os
        savedir=os.getcwd()
        if not filename:
            filename=tkFileDialog.askopenfilename(defaultextension="*.csv",
                                                  initialdir=savedir,
                                                  filetypes=[("Comma separated","*.csv"),("All files","*.*")])
        if filename:
            if os.path.isfile(filename):
                if not self.exp_group_selected:
                    #
                    # Ask which titratable group the data should be associated with
                    #
                    exp_titgroup=Toplevel()
                    exp_titgroup.transient()
                    exp_titgroup.title('Select titratable group for experimental data')
                    self.select_titgrp_frompdb(window=exp_titgroup)
                    chemshift_type=IntVar()
                    Label(exp_titgroup,text='Chemical shift for').grid(row=11,column=0)
                    Radiobutton(exp_titgroup,
                                text='Backbone amide nitrogen',
                                value=0,
                                variable=chemshift_type).grid(row=11,column=1)
                    Radiobutton(exp_titgroup,
                                text='Backbone amide hydrogen',
                                value=1,
                                variable=chemshift_type).grid(row=11,column=2)
                    chemshift_type.set(0)
                    def select_groups():
                        self.exp_group_selected=self.titgroups_inPDB[int(self.LB_groups.curselection()[0])]
                        if chemshift_type.get()==1:
                            self.exp_group_selected=self.exp_group_selected+' - H'
                        else:
                            self.exp_group_selected=self.exp_group_selected+' - N'
                        exp_titgroup.destroy()
                        return
                    Button(exp_titgroup,text='Select Group(s)',command=select_groups).grid(row=10,column=0)
                    self.EM_window.wait_window(exp_titgroup)
                #
                # Read the data
                #
                fd = open(filename,'r')
                lines=fd.readlines()
                fd.close()
                self.expdata[self.exp_group_selected]=[]
                for line in lines[1:]:
                    import string
                    line=string.strip(line)
                    sp=line.split(',')
                    resnum=int(sp[0][1:-4])
                    span=float(sp[1])
                    self.expdata[self.exp_group_selected].append([resnum,span])
                    x,y=self.get_xy(resnum,span)
                    x1,y0=self.get_xy(resnum+1,0)
                    handle=self.tc.create_rectangle(x,y,x1-2,y0,fill='green',stipple='gray50')
                #
                # Add it to the listbox
                #
                self.LB_exp_data.insert(END,self.exp_group_selected)
            else:
                import tkMessageBox
                tkMessageBox.showwarning('Cannot load experimental data file',
                                         'File: %s not found' %filename, parent=self.EM_window)
                return
        return

    #
    # ----
    #

    def load_full_titdata(self,filename=None,pH_start=5.2,pH_end=7.2):
        """Load a full set of titration data from en Ekin project file"""
        import full_titfit
        self.fullfit=full_titfit.titfitter(filename,pH_start,pH_end)
        return

    #
    # ----
    #

    def fit(self):
        """Fit data"""
        expgroups=self.expdata.keys()
        expgroups.sort()
        if len(expgroups)==0:
            import tkMessageBox
            tkMessageBox.showwarning('Cannot fit data',
                                         "Cannot fit because you haven't loaded any experimental data./nLoad experimental data first.", parent=self.EM_window)
            return
        self.fit_window=Toplevel()
        self.set_geometry(self.EM_window,self.fit_window)
        headings=['Exp group','start X','start Y','start Z','X','Y','Z','start eps','eps']
        count=0
        for heading in headings:
            if heading in ['Z','eps']:
                Label(self.fit_window,text=heading).grid(row=0,column=count)
                count=count+1
            else:
                Label(self.fit_window,text=heading).grid(row=0,column=count)
            count=count+1
        # 
        # List all variables and checkboxes for selecting what to optimise
        #
        row=1
        self.show_fit={}
        self.include_var={}
        self.start_val={}
        for group in expgroups:
            self.show_fit[group]={}
            self.start_val[group]={}
            Label(self.fit_window,text=group).grid(row=row,column=0)
            #
            #
            # XYZ start values
            #
            self.start_val[group]['X']=DoubleVar()
            self.start_val[group]['Y']=DoubleVar()
            self.start_val[group]['Z']=DoubleVar()
            Entry(self.fit_window,textvariable=self.start_val[group]['X'],width=4,justify=CENTER,bd=0).grid(row=row,column=1,sticky='news')
            Entry(self.fit_window,textvariable=self.start_val[group]['Y'],width=4,justify=CENTER,bd=0).grid(row=row,column=2,sticky='news')
            Entry(self.fit_window,textvariable=self.start_val[group]['Z'],width=4,justify=CENTER,bd=0).grid(row=row,column=3,sticky='news')
            #
            # Use the coordinates for the group charge as start values
            #
            position=self.get_charge_center(group)
            self.start_val[group]['X'].set(position[0])
            self.start_val[group]['Y'].set(position[1])
            self.start_val[group]['Z'].set(position[2])
            
            self.include_var[group]={}
            self.show_fit[group]['X']=StringVar()
            self.show_fit[group]['X'].set('0.0')
            Label(self.fit_window,textvariable=self.show_fit[group]['X']).grid(row=row,column=4,sticky='news')
            #
            self.show_fit[group]['Y']=StringVar()
            self.show_fit[group]['Y'].set('0.0')
            Label(self.fit_window,textvariable=self.show_fit[group]['Y']).grid(row=row,column=5,sticky='news')
            #
            self.show_fit[group]['Z']=StringVar()
            self.show_fit[group]['Z'].set('0.0')
            Label(self.fit_window,textvariable=self.show_fit[group]['Z']).grid(row=row,column=6,sticky='news')  
            #
            # Checkbutton for including position in fit
            #
            self.include_var[group]['position']=IntVar()
            self.include_var[group]['position'].set(0)
            Checkbutton(self.fit_window,variable=self.include_var[group]['position']).grid(row=row,column=7,sticky='news')

            #
            self.start_val[group]['deff']=DoubleVar()
            self.start_val[group]['deff'].set(8.0)
            Entry(self.fit_window,textvariable=self.start_val[group]['deff'],width=4,justify=CENTER,bd=0).grid(row=row,column=8,sticky='news')
            self.show_fit[group]['deff']=StringVar()
            self.show_fit[group]['Z'].set('0.0')
            Label(self.fit_window,textvariable=self.show_fit[group]['deff']).grid(row=row,column=9,sticky='news')
            #
            # Checkbutton for including deff in fit
            #
            self.include_var[group]['deff']=IntVar()
            self.include_var[group]['deff'].set(0)
            Checkbutton(self.fit_window,variable=self.include_var[group]['deff']).grid(row=row,column=10)
            #
            # Update row
            #
            row=row+1
        #
        # Array holding the names of the parameters we optimise
        #
        self.optimise_vars=[]
        #
        # Define function for updating variables
        #
        def update_vars(difference,variables,value_position,fit_values):
            """Update all the Tkinter vars"""
            for group_name in self.optimise_vars.keys():
                for var_name in self.optimise_vars[group_name].keys():
                    position=value_position[group_name][var_name]
                    self.show_fit[group_name][var_name].set('%5.2e' %variables[position])

            #
            # Redraw the fit
            #
            values=[]
            for datapoint,fit_span in fit_values:
                residue=datapoint[0]
                values.append([residue,fit_span])
            self.show_fit_prediction2(values)
            #
            # Update the sum of differences
            #
            self.difference.set('%5.3e' %difference)
            #
            # Tell Tkinter to update the windows
            #
            self.fit_window.update_idletasks()
            self.EM_window.update_idletasks()
            return
        #
        # Fill in all the experimental data - that doesn't change even when we choose to 
        # fit a different number of variables
        #
        import string
        exp_data=[]
        for group in expgroups:
            # Get the group type
            titgroup=group.split('-')[0]
            import string
            titgroup=string.strip(titgroup)
            titgroup_type=titgroup.split(':')[-1]
            charge=self.PI.titgroups[titgroup_type]['charge']
            for resnum,span in self.expdata[group]:
                pH=-99.9 #Span mode
                atom_type=group.split('-')[-1]
                import string
                atom_type=string.strip(atom_type)
                # Data point format: [residue number, atom_type, pH, charge, span]
                exp_data.append([':'+string.zfill(resnum,4),atom_type,pH,charge,span])
        #
        # Define functions for fitting
        #
        def do_parm_fit():
            """Function for starting the fit"""
            #
            # Find the variables that we should fit
            #
            self.optimise_vars={}
            nonfit_vars={}
            for grp in self.include_var.keys():
                for var in self.include_var[group].keys():
                    if self.include_var[group][var].get()==1:
                        if not self.optimise_vars.has_key(group):
                            self.optimise_vars[group]={}
                        if var=='position':
                            self.optimise_vars[group]['X']=self.start_val[group]['X'].get()
                            self.optimise_vars[group]['Y']=self.start_val[group]['Y'].get()
                            self.optimise_vars[group]['Z']=self.start_val[group]['Z'].get()
                        else:
                            self.optimise_vars[group][var]=self.start_val[group][var].get()
                        #
                    else:    
                        if not nonfit_vars.has_key(group):
                            nonfit_vars[group]={}
                        if var=='position':
                            nonfit_vars[group]['X']=self.start_val[group]['X'].get()
                            nonfit_vars[group]['Y']=self.start_val[group]['Y'].get()
                            nonfit_vars[group]['Z']=self.start_val[group]['Z'].get()
                        else:
                            nonfit_vars[group][var]=self.start_val[group][var].get()
            #
            # Do we have any variables to fit?
            #
            #print 'Fit vars'
            #print self.optimise_vars
            #print
            #print 'nonfit vars'
            #print nonfit_vars
            #print 
            if self.optimise_vars=={}:
                import tkMessageBox
                tkMessageBox.showwarning('No variables to optimise','You must select at least one variable to fit',parent=self.fit_window)
                return
            import EM_effect_fitters
            FIT=EM_effect_fitters.flexible_fitter(exp_data=exp_data,
                                                    fit_variables=self.optimise_vars,
                                                    nonfit_variables=nonfit_vars,
                                                    EM_effect_instance=self,
                                                    callback=update_vars,
                                                    calc_type=self.ele_calc.get())
            FIT.fit(step=0.2,convergence_crit=1E-8)
            return
        # End of do_parm_fit function
        # --->
        #
        def ID_titgroup():
            return
        #
        # General info and convergence criteria
        #
        Label(self.fit_window,text='R^2').grid(row=row,column=0,sticky='news')
        self.difference=StringVar()
        Label(self.fit_window,textvariable=self.difference).grid(row=row,column=1,sticky='news')
        row=row+1
        #
        # Add buttons for doing different kinds of fits
        #
        Button(self.fit_window,text='Fit parameters',command=do_parm_fit).grid(row=row,column=0,columnspan=2)
        Button(self.fit_window,text='ID titgroup',command=ID_titgroup,state=DISABLED).grid(row=row,column=2,columnspan=3)
        return
    
    #
    # ------
    #
    
    def get_charge_center(self,group):
        """Get the center of the charge for a group"""
        import string
        if group.split('-')>1:
            titgroup=string.strip(group.split('-')[0])
        else:
            titgroup=group
        tit_type=titgroup.split(':')[-1]
        residue=self.PI.resnum(titgroup)
        centers=self.PI.calculate_titgroup_center(residue)
        if centers.has_key(residue):
            if centers[residue].has_key(tit_type):
                return centers[residue][tit_type]
        raise Exception,'Could not calculate titgroup center'




if __name__=="__main__":
    X=EM_effect()
    X.update_idletasks()
    #X.load_pdb('2lzt.pka.pdb')
    #
    #
    #X.load_csv(filename='0035GLU.csv',group=':0035:GLU - N')
    #X.add_titgroup(':0035:GLU - N',':0035 GLU',model='Unit charge')
    #X.ele_calc.set('Coulomb')
    #X.update_EM_effect_prediction()
    #X.fit()
    #X.load_full_titdata('HEWL_wt.Ekinprj',pH_start=5.0,pH_end=7.2)
    #X.fit()
    #num_groups=1
    #exp_data=X.fullfit.put_data_in_fitter_format()
    #start_values=[]
    #import random
    #for dummy in range(num_groups):
    #    start_values.append(4.0)
    #    start_values=start_values+[1,2,3]
    #
    # Dielectric constant
    #
    #start_values.append(8.0)
    #FIT=pKa_pos_deff(exp_data=exp_data,start_values=start_values,EM_instance=X,pka_groups=num_groups)
    #print 'Fitting'
    #status,variables=FIT.fit()
    #print status,variables
    #x=[]
    #y=[]
    #for tpH in range(20,100,5):
    #    ph=tpH/10.0
    #    x.append(ph)
    #    y.append(FIT.get_value(start_values,[':0036',ph,0.0]))
    #import python.dislin_driver
    #python.dislin_driver.graf(x,y)
    

    mainloop()
                
