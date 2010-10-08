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

acidbase={'ARG':1,'HIS':1,'LYS':1,'TYR':-1,'ASP':-1,'GLU':-1,
          'CYS':-1,'CTERM':-1,'NTERM':1,'SER':-1,'THR':-1,
          'TPO':-1,'SEP':-1,'ACID':-1,'BASE':1,'ATP':-1}
modelpKas={'NTERM':8.00,'LYS':10.40,'GLU':4.40,'HIS':6.30,
           'ASP':4.00,'TYR':9.6,'ARG':13.0,'CTERM':3.80,
           'CYS':8.7,'SER':15.0,'THR':15.0,
           'TPO':7.0,'SEP':7.0,'ATP':6.8}
