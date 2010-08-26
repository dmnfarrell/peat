#!/usr/bin/env python
#
# pKa - various programs and scripts for pKa value analysis, calculation and redesign
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

def main():
    import sys, types
    for file in sys.argv[1:]:
        fd=open(file)
        import pickle
        d=pickle.load(fd)
        fd.close()
        #
        #
        for target in d.keys():
            if target=='wt_full' or target=='pdbfile':
                continue
            for num in d[target].keys():
                for dist in d[target][num].keys():
                    if not type(d[target][num][dist]) is types.DictType:
                        d[target][num][dist]={'MC':d[target][num][dist]}
                    if not d[target][num][dist].has_key('MC'):
                        d[target][num][dist]={'MC':d[target][num][dist]}
        #
        fd=open(file,'w')
        pickle.dump(d,fd)
        fd.close()

if __name__=="__main__":
    main()
