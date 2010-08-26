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

def get_real_mutations(mutations):
    number=0
    if not eval(mutations):
        return 0
    for mut in eval(mutations):
        if mut:
            number=number+1
    return number

#
# ---------
# 

def main():
    #
    # Get the file
    #
    import sys, types
    file=sys.argv[1]
    #
    # Open
    #
    fd=open(file)
    import pickle
    dict=pickle.load(fd)
    fd.close()
    #
    # Get the correct residues
    #
    if file.find('2lzt')!=-1:
        print 'This is lysozyme'
        donor=':0035:GLU'
        nuc=':0052:ASP'
    elif file.find('1xnb')!=-1:
        print 'This is Bacillus circulans xylanase'
        donor=':0172:GLU'
        nuc=':0078:GLU'
    else:
        raise 'Unknown protein'
    #
    # Construct the result dictionary
    #
    max_dist=15.0
    max_muts=10
    matrix={}
    method='MC'
    for design in dict.keys():
        if design=='wt_full' or design=='pdbfile':
            continue
        this_design=dict[design]
        for num_muts in this_design.keys():
            if num_muts>max_muts:
                continue
            matrix[num_muts]={}
            for dist in this_design[num_muts].keys():
                if dist>max_dist:
                    continue
                matrix[num_muts][dist]=[]

    #
    # Now count for real
    #
    for design in dict.keys():
        if design=='wt_full' or design=='pdbfile':
            continue
        this_design=dict[design]
        for num_muts in this_design.keys():
            for dist in this_design[num_muts].keys():
                if dist>max_dist:
                    continue
                if type(this_design[num_muts][dist][method]) is types.DictType:
                    sol_dict=this_design[num_muts][dist][method]
                    #
                    # Loop over all the solutions
                    #
                    for solution in sol_dict.keys():
                        if not eval(solution):
                            matrix[real_muts][dist].append(0.0)
                            continue
                        real_muts=get_real_mutations(solution)
                        targets=sol_dict[solution].keys()
                        if len(targets)>2:
                            raise 'Too many targets'
                        if len(targets)==2:
                            diff=sol_dict[solution][donor]-sol_dict[solution][nuc]
                            if sol_dict[solution][nuc]<0.3 and real_muts>0:
                                matrix[real_muts][dist].append(diff)
                        elif len(targets)==1:
                            dpka=sol_dict[solution][donor]
                            matrix[real_muts][dist].append(dpka)
    #
    # Save only the best solution
    #
    for num_muts in matrix.keys():
        for dist in matrix[num_muts].keys():
            if matrix[num_muts][dist]==[]:
                matrix[num_muts][dist]=0.0
            else:
                matrix[num_muts][dist].sort()
                matrix[num_muts][dist]=matrix[num_muts][dist][-1]
    #
    # Plot it
    #
    import dislin_driver
    dislin_driver.colour_2D(matrix,title='Max deltapKa',
                            title2='#muts dist',
                            xlabel='# of muts',
                            ylabel='dist from target',
                            zlabel='abs(dpKa shift)',
                            file='two_dpka.tif')
                        

    return


if __name__=="__main__":
    main()
