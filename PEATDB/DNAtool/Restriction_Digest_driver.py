#!/usr/bin/env python
#
# Protein Engineering Analysis Tool DataBase (PEATDB)
# Copyright (C) 2010 Damien Farrell & Jens Erik Nielsen
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
# 

"""Driver for the Restriction digest class in mutation.py"""

from Tkinter import *
import tkFont

def alphabetical(text1,text2):
	import string
	text1=string.lower(text1)
	text2=string.lower(text2)
	if text1==text2:
		return 0
	if text1<text2:
		return -1
	else:
		return 1

class Restriction_Digest:

        def __init__(self):

                return		
        
	def init_vars(self):
		
		# Restriction digest instance  
		import mutation
		self.RS=mutation.restriction_digest()
		
		# Variables for the restriction digest		
		self.show_only_unique_sites=IntVar()
		
		# Min length of recognition sequence		
		self.min_restr_enz_seq=IntVar()
		self.min_restr_enz_seq.set(6)
		
		# Ignore the cut position when defining isoschiziomers		
		self.ignore_cut_position=IntVar()
		self.ignore_cut_position.set(1)
		
		# Exclude promiscuous enzymes		
		self.exclude_promiscuous_enzymes=IntVar()
		self.exclude_promiscuous_enzymes.set(0)
		self.update_promenz()
		
		# Maximum number of sites		
		self.max_num_sites=IntVar()
		self.max_num_sites.set(100)

		#colors for line
		if 'Darwin' in self.currplatform:
			self.linecol='lightblue'
		else:
			self.linecol='blue'
		return

	def restriction_digest(self):
		
		#Get the enzymes to use		
		self.select_enzymes()		
		# Do the restriction digest		
		self.data['cut_pos']=self.RS.get_restriction_sites(self.data['DNAseq'],self.data['used_enzymes'])
		enzymes=self.data['cut_pos'].keys()
		self.plot_restriction_sites(self.data['cut_pos'])
		return

	def plot_restriction_sites(self,cut_pos,colour='black',direction='up',add_text='',
			temporary=None,delete_temporary=1,underline_unique=1):
		"""
		Plot all restriction sites in enzymes
		"""
		enzymes=cut_pos.keys()  #enzyme names, for string label of each site

		#
		# Control vars
		#
		if not getattr(self,'new_seq_win_objs',None):
			self.new_seq_win_objs={}
		#
		step=0
		colour_cycle=[colour]#,'green','magenta','cyan','blue']
		
		# Keep track of all the sites we plot
		
		if not getattr(self,'donepos',None):
			self.donepos={}
		if not getattr(self,'temp_sites',None):
			self.temp_sites={}
			self.temp_objs={}
		
		# Delete all temporary sites from last time
		
		if delete_temporary:
			for obj in self.temp_objs.keys():
				self.seqframe.delete(obj)
			self.temp_objs={}
			
			# Remove the influence of the tempsites on self.donepos			
			for site in self.temp_sites.keys():
				self.donepos[site]=self.donepos[site]-self.temp_sites[site]
			#
			self.temp_sites={}
		
		# Plot up or down?
		
		if direction=='up':
			operator=-1
			base_shift=0
			y_org_shift=0
		else:
			# Down
			operator=1
			base_shift=60
			y_org_shift=10
		#
		#List for storing restr. label objects
		#
		if not getattr(self,'sites',None):
			self.sites=list()
		if not getattr(self,'tempsites',None):
			self.tempsites=list()

		i=0
		#
		# Main Loop
		#
		for enzyme in enzymes:
			positions=cut_pos[enzyme]
			for site in positions:
				#print site
				for subsite in site:
					#print subsite
					x,y=self.get_base_pos_on_screen(subsite+1)
					yorg=y # Store org y position as end-point for line
					yorg=yorg+y_org_shift
					#Apply the base shift for drawing upside or downside
					y=y+base_shift

					# If a comparison sequence is loaded raise more
					if self.show_comp_sequence.get()==1:
						if direction=='up':
							if self.maxseqlevel==0:
								y=(y-15)/self.y_scale
							else:
								y=y-(self.maxseqlevel*15)/self.y_scale
							if self.primer_displayed==1 and self.maxseqlevel>0:
								y=(y-15)/self.y_scale
						else:
							y=(y-15)*self.y_scale

					maxlevel=self.canvas_height-50
					#determine if more than one cut at this site
					l=15
					level=1
					if self.donepos.has_key(subsite):
					   if self.donepos[subsite]>1:
							y=y+l*operator
							level=level+1

					#iterate over donepos to see if nearby sites, and if so shift up
					for donesite in self.donepos.keys():
							if abs(donesite-subsite)<15:
								y=y+l*operator
								level=level+1
					#
					# Should we only show unique sites?
					#
					if self.show_only_unique_sites.get()==1 and len(positions)>1:
						continue
					#
					# Maximum number of sites (other than the unique site button)
					#
					if len(positions)<=self.max_num_sites.get() or self.max_num_sites.get()==100:
						#
						# Select colour
						#
						col=colour_cycle[step]
						step=step+1
						if step==len(colour_cycle):
							step=0

						# Underline unique sites
						font = self.restr_font.get() + " 10"
						uniquetag='nu'
						if underline_unique:
							if len(positions)==1:
								uniquetag='u'
								font=font+' underline'
						
						# Draw site/enzyme label here
						
						obj=self.seqframe.create_text(x,(y-45)/self.y_scale,text=enzyme+add_text,
										activefill='red',font=font,anchor='sw',fill=col,
										tags=('textlabel',uniquetag,level,enzyme+add_text,subsite))

						if temporary:
							self.temp_objs[obj]=1
						else:
							self.new_seq_win_objs[obj]=1
						
						# Keep track of where we've plotted everything						
						if not self.donepos.has_key(subsite):
							self.donepos[subsite]=0
						self.donepos[subsite]=self.donepos[subsite]+1
						#
						if temporary:
							if not self.temp_sites.has_key(subsite):
								self.temp_sites[subsite]=0
							self.temp_sites[subsite]=self.temp_sites[subsite]+1

						#add site text item to list of those currently on canvas
						#if not temporary:
						#	self.sites.append(obj)


		#try to tidy up labels on canvas - testing
		self.tidy_restriction_labels(direction,temporary,underline_unique)
		
		# Restriction details?		
		self.restriction_details()
		return


        def move_restriction_label(self, obj, y):
                """Handle restriction label movement"""                
                c=self.seqframe
                oldy = c.coords(obj)[1]
                enzyme = c.gettags(obj)[3]
                site = c.gettags(obj)[4]
        
                #print enzyme, site, obj                   
                rects=c.find_withtag('labelrect')+c.find_withtag('templabelrect')
                lines=c.find_withtag('line')+c.find_withtag('templine')
                
                for item in rects:
                        if site in c.gettags(item) and enzyme in c.gettags(item):
                                rect = item                                
                for item in lines:
                        tags=c.gettags(item)
                        if site in tags and enzyme in tags:
                                line = item
                                x1,y1,x2,y2 = c.coords(line)
                                c.delete(line)
                                line=c.create_line(x1, y1,x1,y, fill=self.linecol,
                                                width=2,stipple='gray25',
                                                tag=tags)
                                c.tag_lower(line)
                
                c.move(obj, 0, y-oldy)
                c.move(rect, 0, y-oldy)
                c.tag_raise(rect)
                c.tag_raise(obj)                
                return
        
	def tidy_restriction_labels(self,direction='up',temporary=None,underline_unique=None):
		"""This function is for tidying up the restr. site placement after initial drawing
                   It also draws the boxes and lines for the sites"""
		c = self.seqframe
		if direction=='up':
		   l=-15
		   base_shift=0
		else:
			l=15
			base_shift=30

		#clear the temporary rects and lines from previous time
		self.seqframe.delete('templabelrect')
		self.seqframe.delete('templine')
		currentsites=[]
		#use seperate lists depending on whether temp sites are being drawn
		if temporary:
			for obj in self.temp_objs.keys():
				currentsites.append(obj)
		else:
			for obj in self.new_seq_win_objs.keys():
				currentsites.append(obj)
		#for each obj, see if it overlaps nearest, and if so - move it up/down
		for obj in currentsites:
			box = c.bbox(obj)
			x1=box[0]
			y1=box[1]
			x2=box[2]
			y2=box[3]
			overlap=[]
			overlap=c.find_overlapping(x1,y1,x2,y2)
			if overlap:
				for overobj in overlap:
					if overobj != obj:
					#if overlap is a text label, move it upwards until
					#it no longer overlaps another label
						if 'textlabel' in c.gettags(overobj):
							c.move(overobj, 0, l)

		#Now plot rectangle, lift text and then remove sites from list for next time
		#first determine if unique site from tag
		for obj in currentsites:
			#first draw the lines..
			coords = c.coords(obj)
			enzyme = c.gettags(obj)[3]
			sitetag = c.gettags(obj)[4]
			x=coords[0]
			y1=coords[1]
			y2=self.seq_row+base_shift
			yorg_here=y2
			# If we have a protein sequence then lift the line a bit
			if self.data.has_key('ORF_selected'):
				if direction=='up':
					yorg_here=(y2-15)/self.y_scale
				else:
					yorg_here=(y2-15)*self.y_scale
                        # If a comparison sequence is loaded raise line more
			if self.show_comp_sequence.get()==1:
				if direction=='up':
					if self.maxseqlevel==0:
						yorg_here=(yorg_here-15)/self.y_scale
					else:
						yorg_here=yorg_here-(self.maxseqlevel*15)/self.y_scale
  					if self.primer_displayed==1 and self.maxseqlevel>0:
						yorg_here=(yorg_here-15)/self.y_scale
				else:
				    yorg_here=(yorg_here-15)*self.y_scale

			if temporary:
				line=self.seqframe.create_line(x,yorg_here-5,x,y1,fill=self.linecol,
						 width=2,stipple='gray25',tag=('templine',enzyme,sitetag,direction))
			else:
				line=self.seqframe.create_line(x,yorg_here-5,x,y1,fill=self.linecol,
						width=2,stipple='gray25',tag=('line',enzyme,sitetag,direction))
			#lower the lines behind other objects
			c.tag_lower(line)
			if 'u' in c.gettags(obj):
				fillcolour='pink'
			else:
				fillcolour='lightyellow'
			box = c.bbox(obj)			
			if temporary:
				rect = c.create_rectangle(box,tag=('templabelrect',enzyme,sitetag),fill=fillcolour)
			else:
				rect = c.create_rectangle(box,tag=('labelrect',enzyme,sitetag),fill=fillcolour)				
			c.tag_raise(obj)

	
	def show_item(self, event):
		"""Raise items(s) to top if selected with right click
                   and show recognition seq if it's a site label, called from DNAtool main"""
		c=self.seqframe
		box = c.bbox(CURRENT)
		x1=box[0]
		y1=box[1]
		x2=box[2]
		y2=box[3]
		items=[]
		#make selection rectangle one pixel larger to include rect and text
		items=c.find_enclosed(x1-1,y1-1,x2+1,y2+1)
		#get this for recog sequence
		enzymes=self.RS.enzymes_regexs
		
		sfont = tkFont.Font (family='Arial', size=12,weight='bold')
		for obj in items:
			c.tag_raise(obj)
			#if item is text, get recog sequence and display
			if 'textlabel' in c.gettags(obj):
				name=c.itemcget(obj, 'text')
				name=name.rstrip('(+-)')
				seq=self.get_sequence(enzymes[name]['realseq'])
				obj=c.create_text(x2+2,y1-2,text=seq,tags='recogseqlabel',
							font=sfont,width=120,anchor='nw')
				box = c.bbox(obj)
				rect = c.create_rectangle(box,tag='recogseqlabel',fill='yellow')
				c.lift(obj)


	def remove_recog_label(self, event):
		"""Remove recognition sequence label"""
		c=self.seqframe
		c.delete('recogseqlabel')
		return

	def clear_restriction_details(self):
		"""Remove all restriction sites from the sequence window"""
	
		if getattr(self,'new_seq_win_objs',None):
			for obj in self.new_seq_win_objs.keys():
				self.seqframe.delete(obj)
			self.new_seq_win_objs={}
			self.donepos={}
			#
			for obj in self.temp_objs.keys():
				self.seqframe.delete(obj)
			self.temp_objs={}
			self.temp_sites={}
			self.seqframe.delete('labelrect')
			self.seqframe.delete('line')
			self.seqframe.delete('templabelrect')
			self.seqframe.delete('templine')
			#also clear the sites list - this is used in tidying and rendering lines/rects

		return


	def restriction_details(self,junk=None):
		
		# Should we do this?		
		if self.restr_detail.get()==0:
			return
		
		# do we have a digest?		
		if not self.data['cut_pos']:
			self.warning('No restriction digest','Load a DNA sequence first')
			self.restr_detail.set(0)
			return
		
		# Open window		
		if not self.digest_window:
			self.digest_window=Toplevel()
			self.digest_window.transient(self.master)
			self.digest_window.geometry('+100+350')
			self.digest_window.title('Restriction digest summary')
			
			# Sort-by label and button			
			l=Label(self.digest_window,text='Sort enzymes')
			l.grid(row=0,column=0,sticky='nes')
			
			self.sort_by=StringVar()
			self.sort_by.set('by # of cuts')
			self.sortby_button=Menubutton(self.digest_window,textvariable=self.sort_by,relief=RAISED)
			self.sortby_menu=Menu(self.sortby_button,tearoff=0)
			self.sortby_button['menu']=self.sortby_menu
			
			# Set the choices			
			choices=['by # of cuts','alphabetically']
			for choice in choices:
				self.sortby_menu.add_radiobutton(label=choice,variable=self.sort_by,value=choice,indicatoron=0,command=self.restriction_details)
			self.sortby_button.grid(row=0,column=1)
			
			# Close button			
			cls=Button(self.digest_window,text='Close window',command=self.restr_detail_close)
			cls.grid(row=0,column=2)
			
			# Scrollbar			
			scrollbar=Scrollbar(self.digest_window,orient='vertical')
			scrollbar.grid(row=1,column=3,sticky='nws')
			
			# Textbox			
			self.textbox=Text(self.digest_window,background='white',foreground='black',height=30,width=60,state=NORMAL,exportselection=1,yscrollcommand=scrollbar.set)
			self.textbox.grid(row=1,column=0,columnspan=3,sticky='NEWS')
			scrollbar.config(command=self.textbox.yview)
		
		# Clear the textbox
		self.textbox.config(state=NORMAL)
		self.textbox.delete(1.0,END)
		
		# List all restriction enzymes and their cutting positions		
		enzymes=self.data['cut_pos'].keys()

		import string
		self.textbox.insert(END,'Enzymes that cut the sequence\n')
		self.textbox.insert(END,'\n')
		self.textbox.insert(END,'# of sites   enzyme    positions cut (base pos))\n')
		self.textbox.insert(END,'------------------------------------------\n')
		#
		line=2
		if self.sort_by.get()=='by # of cuts':
			#
			# Sort the enzymes according to # of cuts
			#
			cuts={}
			for enz in enzymes:
				poss=[]
				for pos in self.data['cut_pos'][enz]:
					for subsite in pos:
						poss.append(subsite)
				poss.sort()
				
				# Put the enzyme in the correct bin				
				ncuts=len(poss)
				if not cuts.has_key(ncuts):
					cuts[ncuts]={}
		
				postxt=[]
				for pos in poss:
					postxt.append(str(pos))
				cuts[ncuts][enz]=string.join(postxt,', ')
		
			allcuts=cuts.keys()
			allcuts.sort()
			for cut in allcuts:
				these_enzs=cuts[cut].keys()
				these_enzs.sort()
				for enz in these_enzs:
					txt='%2d sites: %12s: %s\n' %(cut,enz,cuts[cut][enz])
					index='%d.%d' %(line,0)
					line=line+1
					#self.textbox.insert(index,txt)
					self.textbox.insert(END,txt)
		else:
                        # Sort the enzymes alphabetically			
			enzymes.sort(alphabetical)
			for enz in enzymes:
				poss=[]
				for pos in self.data['cut_pos'][enz]:
					for subsite in pos:
						poss.append(subsite)
				poss.sort()
				postxt=[]
				for pos in poss:
					postxt.append(str(pos))
				txt='%2d sites %12s: %s\n' %(len(poss),enz,string.join(postxt,', '))
				index='%d.%d' %(line,0)
				line=line+1
				#print index,txt
				#self.textbox.insert(index,txt)
				self.textbox.insert(END,txt)

		# List enzymes that did not cut
		
		import mutation
		all_enzymes=self.data['used_enzymes']
		all_enzymes.sort(alphabetical)
		self.textbox.insert(END,'\n------------------\nEnzymes in the set that do not cut\n')
		count=0
		for enzyme in all_enzymes:
			if not enzyme in enzymes:
				self.textbox.insert(END,'%10s, ' %enzyme)
				count=count+1
				if count==5:
					self.textbox.insert(END,'\n')
					count=0
		
		# Lock the textbox		
		self.textbox.config(state=DISABLED)

		return


	def win_select_enzymes(self):
		#
		# Open a new window for setting the restriction enzymes
		#
		self.set_digest_win=Toplevel()
		self.set_digest_win.geometry('+300+450')
		self.set_digest_win.title('Restriction digest controls')
		#
		self.lbl1=Label(self.set_digest_win,text='Show only unique sites',relief='ridge')
		self.lbl1.grid(row=0,column=0)
		self.uni_sit=Checkbutton(self.set_digest_win,
									  var=self.show_only_unique_sites,onvalue=1,offvalue=0)
		self.uni_sit.grid(row=0,column=1, sticky='wens')
		#
		# Number of cuts to show
		#
		row=1
		l=Label(self.set_digest_win,text='Show only enzymes that cut max:',relief='ridge')
		l.grid(row=row,column=0)
		scl=Scale(self.set_digest_win,from_=1,to=100,resolution=1,orient='horizontal',relief='ridge',variable=self.max_num_sites,
							   label='times')
		scl.grid(row=row,column=1, sticky='wens')
		#
		# Ignore cut position
		#
		row=2
		self.lbl2=Label(self.set_digest_win,text='Ignore cut position when finding isoschiziomers',relief='ridge')
		self.lbl2.grid(row=row,column=0)
		self.uni_sit=Checkbutton(self.set_digest_win,
									  var=self.ignore_cut_position,onvalue=1,offvalue=0)
		self.uni_sit.grid(row=row,column=1, sticky='wens')
		#
		# Length of recognition sequence
		#
		row=3
		self.lbl=Label(self.set_digest_win,text='Min length of recognition sequence',relief='ridge')
		self.lbl.grid(row=row,column=0,sticky='wens')
		self.restr_scale=Scale(self.set_digest_win,from_=2,to=15,resolution=1,orient='horizontal',relief='ridge',variable=self.min_restr_enz_seq,
							   label='bases')
		self.restr_scale.grid(row=row,column=1, sticky='wens')
		#
		# Exclude promiscuous enzymes (enzymes that don't have a single recognition sequence)
		#
		row=4
		self.lbl4=Label(self.set_digest_win,text='Exclude promiscuous enzymes')
		self.lbl4.grid(row=row,column=0)
		self.promiscuous_enz=Checkbutton(self.set_digest_win,
									  var=self.exclude_promiscuous_enzymes,onvalue=1,offvalue=0,
										 command=self.update_promenz)
		self.promiscuous_enz.grid(row=row,column=1, sticky='wens')
		# ------------
		row=5
		#
		# Apply Button
		#
		b = Button(self.set_digest_win, text="Apply new settings",command=self.update_sequence_window)
		b.grid(row=row,column=1,sticky='wens')
		#
		# Close button
		#
		c=Button(self.set_digest_win,text='Close window',command=self.close_win_select_enzymes)
		c.grid(row=row,column=0,sticky='wens')
		return

	#
	# ---------
	#

	def update_promenz(self):
		#
		# Update the promiscuous enzymes flag
		#
		self.RS.exclude_promiscuous_enzymes=self.exclude_promiscuous_enzymes.get()
		self.RS.compile_regexs()
		return

	#
	# -------------------
	#

	def close_win_select_enzymes(self):
		self.set_digest_win.destroy()
		return

	#
	# ----------------
	#

	def select_enzymes(self,include=[],exclude=[]):
		#
		# Select enzymes based on length and mandatory inclusions and exclusions
		#
		enzymes=self.RS.enzymes_regexs
		#
		# Calculate length and define isoschiziomers
		#
		preferred_enzymes=['Aat II','Acc I','Acc III','Alu I',
						   'Bbu I','Bcl I','Bgl II','Bst98 I','BstE II','BstZ I',
						   'Cfo I','Sac I','EcoR I','EcoR V','Fok I','Hae II',
						   'Hae III','Hinc II','Hind III','Hinf I','Hpa I','Hpa II',
						   'Kpn I','Mbo I','Mbo II','Mlu I','Msp I','Nae I','Nar I',
						   'Nci I','Nco I','Nde I','Nde II','Nhe I','Not I','Nru I',
						   'Nsi I','Pst I','Sac II','Sal I','Sca I','Sfi I','Sgf I',
						   'Sin I','Sma I','SnaB I','Spe I','Sph I','Ssp I','Stu I',
						   'Sty I','Taq I','Tru9 I','Tth111 I','Vsp I','Xba I','Xho I',
						   'Xho II','Xmn I']
		#
		# Check that we have all enzymes
		#
		#for enzyme in preferred_enzymes:
		#    if not enzyme in enzymes:
		#         print 'Missing',enzyme
		iso_s={}
		rseq_len={}
		for enzyme in enzymes.keys():
			import string
			seq=self.get_sequence(enzymes[enzyme]['realseq'])
			#
			# Do we have another enzyme that cuts like this?
			#
			if iso_s.has_key(seq):
				iso_s[seq].append(enzyme)
			else:
				iso_s[seq]=[enzyme]
			#
			# Clean the recognition sequence
			#
			import string
			#seq=string.replace(seq,'^','')
			c_seq=string.replace(seq,'(','')
			seq=string.replace(seq,')','')
			#
			# Do we have numbers?
			#
			numbers=[]
			#last_was_num=None
			thisnum=''
			for char in c_seq:
				if char in string.digits:
					thisnum=thisnum+char
				else:
					if thisnum!='':
						numbers.append(thisnum)
						thisnum=''
			if thisnum!='':
				numbers.append(thisnum)
			#
			# Calculate length
			#
			rseq_len[enzyme]=len(c_seq)
			for number in numbers:
				#print 'Old: %3d, -: %2d, +:%2d' %(rseq_len[enzyme],len(number),int(number))
				rseq_len[enzyme]=rseq_len[enzyme]-len(number)+int(number)-1
			#print numbers
			#print enzyme,enzymes[enzyme].keys()[0],rseq_len[enzyme]
			#print
		#
		# ----------------------------------------------------------
		#
		# Apply the selection criteria
		# iso_s_added holds all the enzymes that should be used
		#
		iso_s_added={}
		self.data['used_enzymes']=[]
		for enz in rseq_len.keys():
			seq=self.get_sequence(enzymes[enz].keys()[0])
			add=None
			if rseq_len[enz]>=self.min_restr_enz_seq.get():
				add=1
			if enz in include:
				add=1
			if enz in exclude:
				add=None
			#
			# Add the enzyme?
			#
			if add:
				self.data['used_enzymes'].append(enz)
				if not iso_s_added.has_key(seq):
					iso_s_added[seq]=[enz]
				else:
					iso_s_added[seq].append(enz)
		#
		# ------------------------------------------------------------
		#
		# Check that we use the correct and preferred iso_schiziomers
		#
		new_selection=[]
		iso_s_added_new={}
		for enzyme in self.data['used_enzymes']:
			#
			# Is this one of the redundant enzymes?
			#
			seq=self.get_sequence(enzymes[enzyme]['realseq'])
			#
			# Check if we should add it.
			#
			add=1
			if iso_s_added.has_key(seq):
				#
				# Was another isoschiziomer already added?
				#
				if iso_s_added_new.has_key(seq):
					add=None
				elif len(iso_s[seq])==1:
					#
					# If there is only one enz then just add it
					#
					iso_s_added_new[seq]=1
				else:
					#
					# If not, then is this our preferred name?
					#
					name_control=None
					import string
					for enz_name in iso_s[seq]:
						if enz_name in preferred_enzymes:
								name_control=enz_name
					#
					# If name control then check this is the correct one
					# (and check that the preferred name isn't exluded)
					#
					if name_control:
						#
						# Is the preferred enzyme included in the set?
						#
						if name_control in iso_s_added[seq]:
							if not enzyme in preferred_enzymes:
								add=None
								#print 'Skipping %s because of name control' %enzyme
							else:
								#
								# Mark that we added the preferred enz
								#
								iso_s_added_new[seq]=1
								#print 'Adding because it is preferred',enzyme
						else:
							#
							# OK, the preferred enzyme is not added
							# Then we exclude this enzyme as well
							#
							add=None
					else:
						#
						# No name control just add this one
						#
						iso_s_added_new[seq]=1
			if add:
				import string
				new_selection.append(enzyme)
		#
		# Overwrite the old selection
		#
		self.data['used_enzymes']=new_selection
		return

	#
	# ---------------
	#

	def get_sequence(self,seq):
		#
		# Get the restriction enzyme recognition sequence and
		# see if we should ingore the position of the cut
		#
		import string
		if self.ignore_cut_position.get()==1:
			seq=string.upper(string.strip(string.replace(seq,'^','')))
		return seq

	#
	# ---------------
	#

	def restr_on_off(self,junk=None):
		newstate=self.restr_detail.get()
		if newstate==1:
			if self.digest_window:
				#
				# Already open
				#
				return
			#
			# open win
			#
			self.restriction_details()
		else:
			#
			# Kill?
			#
			if not self.digest_window:
				#
				# Already dead
				#
				return
			else:
				self.digest_window.destroy()
				self.digest_window=None
				return
		return

	def restr_detail_close(self,junk=None):
		#
		# Close the restriction detail window
		#
		self.restr_detail.set(0)
		self.restr_on_off()
		return

#
# -------------------------
#

