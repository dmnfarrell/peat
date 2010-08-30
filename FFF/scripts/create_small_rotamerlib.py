#!/usr/bin/env python
#
# FFF - Flexible Force Field
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
import Protool.rotamer_lib as RL
aas=RL.rots.keys()
aas.sort()
import random
for aa in aas:
    for x in range(1000):
        rot=random.choice(RL.rots[aa])
        #for rot in RL.rots[aa]:
        #print rot
        txt='%3s %4d %4d %d %d %d %d %d %f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f' %(aa,rot['phi'],rot['psi'],
            0,0,0,0,0,rot['p'],rot['chi1'],rot['chi2'],rot['chi3'],rot['chi4'],0,0,0,0)
        print txt
        
