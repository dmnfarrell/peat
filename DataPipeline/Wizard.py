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

from Base import BaseImporter

"""Implements a simple wizard in the application to allows users to
select and appropriate Importer and some other base settings"""

from Tkinter import *

class WizardDialog(Frame):
    
    def __init__(self, parent=None):
        self.parent=parent         
        if not self.parent:
            Frame.__init__(self)
            self.main=self.master
        else:
            self.main=Toplevel()
            self.master=self.main
        self.main.title('Wizard')
        self.main.geometry('500x500+250+200')
        self.spanel = Frame(self.main)
        self.spanel.pack(fill=BOTH)
        self.step1()
        
    def step1(self):
        b = self.getButton(self.spanel)
        b.grid(row=0,column=0)
        
    def getButton(self, parent):
        self.img = case1()
        button = Button(parent,image=self.img)
        return button       
    
def case1():
    img = PhotoImage(format='gif',data=
              'R0lGODlhgAAuAKU4AAAAAAQEBA4ODhMTEzg4OD8/P0BAQEJCQkRERIs3NYw4'
            +'NjdbhzdbiDdciThdijlgjjxklLNKSD1mlz5ml7VLSUBom8BQTcNRTkVxqEl4'
            +'sUp4sbJqL7NrMLVsMLZtMLZtMU+BvXCIPnGJP3OLQHSMQeGIP+eMQeqOQu+R'
            +'Q5CuUpGvU5KwU/eWRvuYR6ioqJq6WJu7WZy8Wb29vcDAwMbGxszMzODg4Pv7'
            +'+////////////////////////////////yH+EUNyZWF0ZWQgd2l0aCBHSU1Q'
            +'ACH5BAEKAD8ALAAAAACAAC4AAAb+wB9jQSwaj0gk48dsJJ9QY4MpjFqVVOd1'
            +'S2w0MJmweEwulzHTH8TMbo8hTbB7TkYz1/R8GOLQgP6AgYKDgxoOTBOEiouB'
            +'E0x9jJGFhz+Jkpd/E5CYkoaInJeOP5ugi56VpZGafqmKp5athKKksYGvtbK0'
            +'uCC3u42PrL68lLDCq8J/vcizwb7Kxrq4z77MyMOf1sfI07vV28TWIJoSFeXm'
            +'5+jp6RKUD+rv8OcPj+Tx9ujsTO73/OXzVAADChxIsKDBgwgTKlzIsKHDhxAj'
            +'SgyYoKLFixgzaqyogAmAjyBDihxJ8qPHkihTgjypsqVIBRQuyJxJs6ZNmxQ6'
            +'/kCAo6f+z59AgwpF4LGG0KNIfdYAwIRn0qdBESiIYKGq1atYs2aNoNMA1K8+'
            +'DXicAfbrDKY/vJaFamCq1rdwrXJlonZtUrE/AJC1i/QsXb533cYdjHVuWsBI'
            +'8epFLNTvYcZB21IlTNmC4bqQw47N/NMxZs6SK1O+zPmn4r2lPZcOK1h0XNKr'
            +'cZyOrTp2aNevu8aWvXl17dVtNwofzpGly+PGj7tMrjzlxOfQo0ufTr269YQh'
            +'smvfzr27d+3Mm5cML54k+fIhQ6R4wb69+/fw4acIwaTAjfv48+vfz7+ARxn8'
            +'BSggfjKgZd+ACO5XgAgqwODggxBGKKGEKojwl229pYbWZ5n+GcDghCCG+GCF'
            +'FwKXIWe/leZhgyK2GCGJj5mYF2oobrjbii7m6CCMHEI2m2822vahji7yeOOJ'
            +'maUI2pBEimgkhjPSFiRwIqwQw5VYZqnllluuYOEPB+x2gEc0xEYDWmHGdkAI'
            +'I5Dg5ptwximnnCPQ90MAA+Sp55589ulnAB4J4OeghOopAFp4Fqpon4Be5+ij'
            +'kEYq6aQGdcDBpZhmqummnHLgwXnomZRXqM6NSipJH5SAwqqsturqq6+W8AET'
            +'BNhg66245qrrrgR45MKuwAZ7qwto1SrssboS8IEJLDTr7LPQRhutCbPGqCKS'
            +'kCnZ4bLSduuts9SWeG2UQIoLGrdi36YLbbjWgoYtY9r6iK669LLbI2M/amju'
            +'tszSW2+19yKWb437ytuvv+naeyS5+ra77QktRCzxxBRXXPEJ1Tq1GlF5GbXa'
            +'Uk3thsAHHGxg8skop6yyyhxUe+p4pr4sEqihBgEAOw==')

    return img
    
