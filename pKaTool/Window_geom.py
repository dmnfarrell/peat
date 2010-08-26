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


def get_geometry(widget):
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

def set_geometry(pwidget,widget):
    """Set the position of widget in the middle of pwidget"""
    w,h,x,y=get_geometry(pwidget)
    pwidget.update()
    widget.update()
    sw,sh,dummy,dummy2=get_geometry(widget)
    xoffset=int((w-sw)/2)
    yoffset=int((h-sh)/2)
    widget.geometry('+%d+%d' %(x+xoffset,y+yoffset))
    return

