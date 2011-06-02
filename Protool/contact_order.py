#!/usr/bin/env python
#
# # Protool - Python class for manipulating protein structures
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

class contact_order:

    def calculate_RCO(self,ss=False):
        """Calculate the relative contact order for the soup. Setting ss to True will lead to the sequence distance being determine using SS bridges"""
        cutoff=6.0
        contacts={}
        done={}
        if hasattr(self,'sep_results'):
            delattr(self,'sep_results')
        for residue1 in sorted(self.residues.keys()):
            print residue1, 
            import sys
            contacts[residue1]=self.calculate_residue_CO(residue1,ss=ss)
            sys.stdout.flush()
        #
        # Calculate the contact order
        #
        #print contacts
        sum_cont=[]
        RCOs={}
        for res in contacts.keys():
            RCO=sum(contacts[res])/float(len(contacts[res]))
            sum_cont=sum_cont+contacts[res]
            RCOs[res]=RCO
            
        RCO=sum(sum_cont)/(float(len(self.residues.keys()))*len(sum_cont))
        RCOs['all']=RCO
        return RCOs
        
    #
    # ------
    #
    
    def res_separation(self,ID1,ID2,ss=False):
        """Get the number of residues between the two IDs"""
        #print 'Finding shortest route between',ID1,ID2
        res1=self.resid(ID1)
        res2=self.resid(ID2)
        #
        # Keep a dictionary of results
        #
        if not hasattr(self,'sep_results'):
            self.sep_results={}
            for r in self.residues.keys():
                self.sep_results[r]={}
        if self.sep_results[res1].has_key(res2):
            return self.sep_results[res1][res2]
        #
        # Use Dijkstra's algorithm to find the shortest path between the two
        #
        G={}
        for residue in self.residues.keys():
            G[residue]={}
            neighbours=self.get_neighbours(residue,ss=ss)
            G[residue]['NB']={}
            for r2 in neighbours:
                G[residue]['NB'][r2]=1 # 1 is the distance between neighbours
            G[residue]['dist']=9999999
            G[residue]['visited']=False
        #
        # We have the graph, now search it
        #
        current=res1
        G[res1]['dist']=0.0
        done=False
        while not done:
            for nb in G[current]['NB']:
                if not G[nb]['visited']:
                    thisdist=G[current]['dist']+G[current]['NB'][nb]
                    if thisdist<G[nb]['dist'] and thisdist<G[res2]['dist']:
                        G[nb]['dist']=thisdist
            G[current]['visited']=True
            #
            # Find new current node
            #
            mindist=99999999999999
            bestnode=None
            for node in G.keys():
                if not G[node]['visited']:
                    if G[node]['dist']<mindist:
                        mindist=G[node]['dist']
                        bestnode=node
            current=bestnode
            if current is None or G[current]['dist']>G[res2]['dist']:
                done=True
        #
        # Shortest distance
        #
        self.sep_results[res1][res2]=G[res2]['dist']
        self.sep_results[res2][res1]=G[res2]['dist']
        #print 'Shortest distance between %7s and %7s is %5.1f' %(res1,res2,G[res2]['dist'])
        return G[res2]['dist']
        
    #
    # ------
    #
        
    def get_neighbours(self,residue,ss=False):
        """Get the neighbours in sequence for this residue"""
        NBs=[]
        try:
            next=self.NextResidue(residue)
            NBs.append(next)
        except:
            pass
        try:
            prev=self.PreviousResidue(residue)
            NBs.append(prev)
        except:
            pass
        #
        # Use SS bonds?
        #
        if ss:
            partner=self.is_SSbonded(residue)
            if partner:
                NBs.append(partner)
        return NBs

    #
    # ------
    #
        
    def calculate_Elec_CO(self,residue1):
        """Calculate the Electrostatic contact order"""
        for residue2 in sorted(self.residues.keys()):
            if residue1==residue2:
                continue
            if self.dist('%s:CA' %(residue1),'%s:CA' %(residue2))>25.0:
                continue
            for atom1 in self.residues[residue1]:
                for atom2 in self.residues[residue2]:
                    if self.dist(atom1,atom2)<=cutoff:
                        contacts.append(self.res_separation(atom1,atom2,ss=ss))    
                        

        return

    #
    # -----
    #

        
    def calculate_residue_CO(self,residue1,ss=False):
        """Calculate the contact order for a single residue"""
        cutoff=6.0
        contacts=[]
        for residue2 in sorted(self.residues.keys()):
            if residue1==residue2:
                continue
            if self.dist('%s:CA' %(residue1),'%s:CA' %(residue2))>25.0:
                continue
            for atom1 in self.residues[residue1]:
                for atom2 in self.residues[residue2]:
                    if self.dist(atom1,atom2)<=cutoff:
                        contacts.append(self.res_separation(atom1,atom2,ss=ss))                       
        #
        # Calculate the contact order
        #
        #print contacts
        #RCO=sum(contacts)/float(len(contacts))
        #print 'RCO for residue: %5.3f' %RCO
        return contacts
        
        
        
        
 
        