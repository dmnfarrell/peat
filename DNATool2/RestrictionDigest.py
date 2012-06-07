#!/usr/bin/env python
#
# DNATool - A program for DNA sequence manipulation
# Copyright (C) 2012- Damien Farrell & Jens Erik Nielsen
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
# Email: farrell.damien_at_gmail.com

"""Restriction digest methods for DNATool"""

import string
#import Mutation

#RS = mutation.restriction_digest()

def selectEnzymes(self,include=[],exclude=[]):
    """Select enzymes based on length and mandatory
    inclusions and exclusions"""

    enzymes = RS
    
    # Calculate length and define isoschiziomers
    
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
    
    # Check that we have all enzymes
    
    #for enzyme in preferred_enzymes:
    #    if not enzyme in enzymes:
    #         print 'Missing',enzyme
    iso_s={}
    rseq_len={}
    for enzyme in enzymes.keys():        
        seq=self.get_sequence(enzymes[enzyme]['realseq'])
        
        # Do we have another enzyme that cuts like this?
        
        if iso_s.has_key(seq):
            iso_s[seq].append(enzyme)
        else:
            iso_s[seq]=[enzyme]
        
        # Clean the recognition sequence
                
        #seq=string.replace(seq,'^','')
        c_seq=string.replace(seq,'(','')
        seq=string.replace(seq,')','')
        
        # Do we have numbers?
        
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
        
        # Calculate length
        
        rseq_len[enzyme]=len(c_seq)
        for number in numbers:
            #print 'Old: %3d, -: %2d, +:%2d' %(rseq_len[enzyme],len(number),int(number))
            rseq_len[enzyme]=rseq_len[enzyme]-len(number)+int(number)-1
        #print numbers
        #print enzyme,enzymes[enzyme].keys()[0],rseq_len[enzyme]
        #print

    # Apply the selection criteria
    # iso_s_added holds all the enzymes that should be used
    
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
        
        # Add the enzyme?        
        if add:
            self.data['used_enzymes'].append(enz)
            if not iso_s_added.has_key(seq):
                iso_s_added[seq]=[enz]
            else:
                iso_s_added[seq].append(enz)    
   
    # Check that we use the correct and preferred iso_schiziomers    
    new_selection=[]
    iso_s_added_new={}
    for enzyme in self.data['used_enzymes']:        
        # Is this one of the redundant enzymes?        
        seq=self.get_sequence(enzymes[enzyme]['realseq'])
        
        # Check if we should add it
        
        add=1
        if iso_s_added.has_key(seq):
            
            # Was another isoschiziomer already added?
            
            if iso_s_added_new.has_key(seq):
                add=None
            elif len(iso_s[seq])==1:                
                # If there is only one enz then just add it                
                iso_s_added_new[seq]=1
            else:                
                # If not, then is this our preferred name?                
                name_control=None                
                for enz_name in iso_s[seq]:
                    if enz_name in preferred_enzymes:
                            name_control=enz_name
                
                # If name control then check this is the correct one
                # (and check that the preferred name isn't exluded)
                
                if name_control:                    
                    # Is the preferred enzyme included in the set?                    
                    if name_control in iso_s_added[seq]:
                        if not enzyme in preferred_enzymes:
                            add=None
                            #print 'Skipping %s because of name control' %enzyme
                        else:                            
                            # Mark that we added the preferred enz                            
                            iso_s_added_new[seq]=1
                            #print 'Adding because it is preferred',enzyme
                    else:                        
                        # OK, the preferred enzyme is not added
                        # Then we exclude this enzyme as well                        
                        add=None
                else:                    
                    # No name control just add this one                    
                    iso_s_added_new[seq]=1
        if add:        
            new_selection.append(enzyme)
    
    # Overwrite the old selection    
    self.data['used_enzymes']=new_selection
    return



