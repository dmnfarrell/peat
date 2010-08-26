#!/bin/env python
#
# UFFBAPS - Unified Force Field for Binding And Protein Stability
# Copyright (C) 2010 Jens Erik Nielsen & Chresten Soendergaard
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


import cgi,random,cgitb,os,time
cgitb.enable()


class BACI:
    # Binding Affinity Calculation Initiator
    def __init__(self):
        self.calculation_dir = ''
        self.base_dir = '/var/www/cgi-bin/uffbaps/calculations/'
        return

    def start_binding_affinity_calculation(self):
        self.out = """Content-Type: text/html

"""
        res = ''
        self.form = cgi.FieldStorage()
        if self.form.has_key("protein") and self.form.has_key("ligand"):
            res = self.do_run()
        else:
            self.out += """You will need to upload structure files -
go <a href=\"http://peat.ucd.ie/cgi-bin/uffbaps/index.cgi\">back</a> and try again
"""
        print self.out
        return res

    def do_run(self):

        while (os.path.isdir(self.calculation_dir) or self.calculation_dir==''):
            self.generate_calculation_dir()

        self.make_files(self.form["protein"].value, self.form["ligand"].value)
        self.make_run_file(self.form["ligand"].value)
        self.run()
        self.get_result()
        res= self.make_xml_file()

        
        self.out += "<br>"
        self.out += self.form["protein"].filename
        self.out += "<br>"
        self.out += self.form["ligand"].filename
        self.out += "<br>"
        self.out += self.calculation_dir
        
        return res

    def generate_calculation_dir(self):
        d = time.strftime("%Y-%m-%d-",time.gmtime())
        ##
        i = random.randint(1000,9999)
        self.calculation_dir = self.base_dir+d+str(i)+'/'
        return



    def make_files(self, protein_file, ligand_file):

        os.system('cp -r /var/www/cgi-bin/uffbaps/calculations/new_binding_affinity_run '+self.calculation_dir)
        # write protein
        u = open(self.calculation_dir+'protein.pdb','w')
        for line in protein_file:
            u.write(line)
        u.close()

        # write ligand
        u = open(self.calculation_dir+'ligand.mol2','w')
        for line in ligand_file:
            u.write(line)
        u.close()
        

        return
    
    def make_run_file(self,ligand_file):

        # count number of models in mol2 file
        no_model = ligand_file.count('@<TRIPOS>MOLECULE')

        # more than one model => we need to overwrite the standard run.txt
        if no_model > 1:
            u = open(self.calculation_dir+'run.txt','w')
            u.write('\nload protein.pdb protein\nload ligand.mol2 ligand\n\n')
            u.write('task energy {')            
            for i in range(1,no_model+1):
                u.write('protein ')
            u.write('}{')
            for i in range(1,no_model+1):
                u.write('ligand_%d '%i)
            u.write('}\n')

        self.out += 'Number of models:%d'%no_model

        return
        
    def run(self):
        os.chdir(self.calculation_dir)
        os.system('/usr/local/bin/genialtNavn '+self.calculation_dir+'run.txt > '+self.calculation_dir+'output.txt')
        os.chdir(self.base_dir)
        return

    def get_result(self):
        lines = open(self.calculation_dir+'output.txt','r').readlines()

        no = 1
        
        self.out+="""<br><pre>
          |                  COMPONENTS                       |  DG total       
----------+---------------------------------------------------+--------------
Structure |     vdW     ES     HB     HP     PE     LE     WE |  Prediction 
----------+---------------------------------------------------+--------------
"""
        self.components = {}
        self.predictions = {}
        
        for line in lines:
            if 'protein   |' in line:
                self.out+='Ligand %3d|'%no + line[11:-22]+line[-12:]

                # store components
                temp_components = line[11:-26].split(' ')
                while temp_components.count('')>0:
                    temp_components.remove('')

                self.components[no] = temp_components

                # store energy
                temp_energy = line[-13:-1]
                while ' ' in temp_energy:
                    temp_energy = temp_energy.replace(' ','')

                self.predictions[no] = temp_energy

                no+=1

        self.out+='</pre>'
        return
        

    def make_xml_file(self):

        res = '<UFFBAPSresult type=\"BA\">'

        for i in range(len(self.predictions.keys())):
            res +="""
  <Ligand>
    <vdW>%s</vdW>
    <ES>%s</ES>
    <HB>%s</HB>
    <HP>%s</HP>
    <PE>%s</PE>
    <LE>%s</LE>
    <WE>%s</WE>
    <ENERGY>%s</ENERGY>
    <number>%d</number>
  </Ligand>
"""%(
    self.components[i+1][0],
    self.components[i+1][1],
    self.components[i+1][2],
    self.components[i+1][3],
    self.components[i+1][4],
    self.components[i+1][5],
    self.components[i+1][6],
    self.predictions[i+1],
    i+1
    )

        res +='</UFFBAPSresult>\n'
        
        fd = '/uffbaps/'+self.calculation_dir[-15:-1]+'_results.xml'
        u = open('/var/www/html'+fd,'w')
        u.write(res)
        u.close()

        self.out+='<br><a href=%s>XML</a>' %fd

        return res



if (__name__ == "__main__"):
    b = BACI()
    b.start_binding_affinity_calculation()
