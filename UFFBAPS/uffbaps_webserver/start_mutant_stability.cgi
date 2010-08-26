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


class MSCCI:
    # Mutant Stability Change Calculation Initiator
    def __init__(self):
        self.calculation_dir = ''
        self.base_dir = '/var/www/cgi-bin/uffbaps_webserver/calculations/'
        return

    def start_mutant_stability_change_calculation(self):
        self.out = """Content-Type: text/html

"""
        res = ''
        self.form = cgi.FieldStorage()
        if self.form["wt"].value!='' and self.form["mt"].value!='':
            res = self.do_run()
        else:
            self.out += """You will need to upload structure files -
go <a href=\"http://peat.ucd.ie/cgi-bin/uffbaps_webserver/index.cgi\">back</a> and try again
"""
        print self.out
        return res

    def do_run(self):

        ## store name
        self.name = self.form["msc_name"].value
        
        while (os.path.isdir(self.calculation_dir) or self.calculation_dir==''):
            self.generate_calculation_dir()

            
        self.make_files(self.form["wt"].value, self.form["mt"].value)
        self.run()
        self.get_result()

        res = self.make_xml_file()
        self.clean_up()

        return res

    def generate_calculation_dir(self):
        d = time.strftime("%Y-%m-%d-",time.gmtime())
        ## generate random dir
        i = random.randint(1000,9999)
        self.calculation_dir = self.base_dir+d+str(i)+'/'

        ## check that it does not already exist
        if os.path.isdir(self.calculation_dir):
            self.calculation_dir = ''

        return



    def make_files(self, protein_file, ligand_file):

        os.system('cp -r /var/www/cgi-bin/uffbaps_webserver/calculations/new_run_dir '+self.calculation_dir)
        # write protein
        u = open(self.calculation_dir+'wt.pdb','w')
        for line in protein_file:
            u.write(line)
        u.close()

        # write ligand
        u = open(self.calculation_dir+'mt.pdb','w')
        for line in ligand_file:
            u.write(line)
        u.close()

        return
    
    def run(self):
        os.chdir(self.calculation_dir)
        os.system('/local/chresten/UFFBAPS/trunk/Uffbaps '+self.calculation_dir+'run_msc.txt > '+self.calculation_dir+'output.txt')
        os.chdir(self.base_dir)
        return

    def get_result(self):

        lines = open(self.calculation_dir+'stability.m','r').readlines()
        self.dd_energies = lines[1].split()

        self.energy = 0.0
        for i in range(len(self.dd_energies)):
            self.energy += float(self.dd_energies[i])



        self.out+='<h1>The results for your calculation %s</h1>'%self.name

        self.out+= """
<table rules=all cellspacing=30 cellpadding=10>
  <tr>
    <th>Van der Waals</th>
    <th>Electrostatic</th>
    <th>Hydrogen Bonds</th>
    <th>Desolvation</th>
    <th>Backbone Conformation</th>
    <th>Protein Entropy</th>
    <th>Ligand Entropy</th>
    <th>Water Bridges</th>
    <th><b>Total Stability Change</b></th>
  </tr>
  <tr>
    <td align=right>%s</td>
    <td align=right>%s</td>
    <td align=right>%s</td>
    <td align=right>%s</td>
    <td align=right>%s</td>
    <td align=right>%s</td>
    <td align=right>%s</td>
    <td align=right>%s</td>
    <td align=right><b>%s</b></td>
  </tr>
</table>
<br>
All energies are in kJ/mol. The calculation was done using the files %s (wild type) and %s (mutant).
"""%(
    self.dd_energies[0],
    self.dd_energies[1],
    self.dd_energies[2],
    self.dd_energies[3],
    self.dd_energies[4],
    self.dd_energies[5],
    self.dd_energies[6],
    self.dd_energies[7],
    self.energy,
    self.form["wt"].filename,
    self.form["mt"].filename
    )


        

        return
        

    def make_xml_file(self):

        res = """
<UFFBAPSresult type=\"MSC\">
  <Name>%s</Name>
  <D_vdW>%s</D_vdW>
  <D_ES>%s</D_ES>
  <D_HB>%s</D_HB>
  <D_DE>%s</D_DE>
  <D_BB>%s</D_BB>
  <D_PE>%s</D_PE>
  <D_LE>%s</D_LE>
  <D_WB>%s</D_WB>
  <ENERGY>%s</ENERGY>
</UFFBAPSresult>
"""%(
    self.name,
    self.dd_energies[0],
    self.dd_energies[1],
    self.dd_energies[2],
    self.dd_energies[3],
    self.dd_energies[4],
    self.dd_energies[5],
    self.dd_energies[6],
    self.dd_energies[7],
    self.energy
    )

        fd = '/uffbaps/calculations/'+self.calculation_dir[-14:-1]+'_results.xml'
        u = open('/var/www/html'+fd,'w')
        u.write(res)
        u.close()

        self.out+='<br>Here is a permanent link to your results in a <a href=%s>XML</a> file.' %fd
        self.out += '<br><a href=\"http://peat.ucd.ie/cgi-bin/uffbaps_webserver/index.cgi\">New calculation</a><br>'

        return res


    def clean_up(self):
        c = 'rm -rf '+self.calculation_dir
        os.system(c)
        return


if (__name__ == "__main__"):
    b = MSCCI()
    b.start_mutant_stability_change_calculation()

