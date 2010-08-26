#!/bin/env python
#
# Written by Chresten Sondergaard (2008)
#
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


def make_front_page():
    out  = """Content-Type: text/html

<html>
<body>
<center>

<table border="5" frame="void" rules="cols" cellspacing="20%">
  <tr>
    <td colspan=2 align=center>
      <h1>UFFBAPS Server</h1>
    </td>
  </tr>

  <tr>
    <td width=40%><h2>Binding affinity prediction</h2></td>
    <td width=40%><h2>Mutant stability change prediction</h2></td>
  </tr>
  <tr>
    <td>
      <form action=\"start_binding_affinity.cgi\" method=post  enctype=\"multipart/form-data\">
        <table>
          <tr>
            <td>
              Protein structure (PDB file) 
            </td>
            <td>
              <input type=\"file\" name =\"protein\" disabled>
            </td>
          </tr>
          <tr>
            <td>
              Ligand structure (mol2 file) 
            </td>
            <td>
              <input type=\"file\" name =\"ligand\" disabled>
            </td>
          </tr>
          <tr>
            <td>
              <input type=\"submit\" value=\"Submit\" disabled>
            </td>
          </tr>
        </table>
      </form>
    </td>
    <td>

      <form action=\"start_mutant_stability.cgi\" method=post  enctype=\"multipart/form-data\">
        <table>
          <tr>
            <td>
              Name*
            </td>
            <td>
              <input type=\"text\" name =\"msc_name\">
            </td>
          </tr>
          <tr>
            <td>
              Wild type structure (PDB file) 
            </td>
            <td>
              <input type=\"file\" name =\"wt\">
            </td>
          </tr>
          <tr>
            <td>
              Mutant structure (PDB file) 
            </td>
            <td>
              <input type=\"file\" name =\"mt\">
            </td>
          </tr>
          <tr>
            <td>
              <input type=\"submit\" value=\"Submit\">
            </td>
          </tr>
          <tr>
            <td colspan=2>
              * For future reference give your calculation a name that makes sense to you 
            </td>
          </tr>
        </table>
      </form>
      


    </td>
  </tr>
</table>

</center>
</body>
</html>"""
    
    print out
    return

if (__name__ == "__main__"):
    make_front_page()
