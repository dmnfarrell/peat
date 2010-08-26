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

import os

convert={'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE','K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER','T':'THR','V':'VAL','W':'TRP','Y':'TYR','TPO':'Z','SEP':'X','-':'---','X':'XXX'}


def main(options,args):
    """Extract mutations from a multiple sequence alignment"""
    import PEATDB.sequence_alignment as SA
    alignment=SA.read_clustal_aln_file(args[0])

    HEWL=options.wt
    HEWL_seq=alignment[HEWL]


    fd=open('mutations.txt','w')
    fd2=open('frequencies.csv','w')

    aas=convert.keys()
    aas.sort()
    import string
    fd2.write('WT Residue number, %s\n' %(string.join(aas,',')))
    #
    real_pos=0
    lines = []
    for position in range(0,len(HEWL_seq)):
        res_observed={}
        if HEWL_seq[position]=='-':
            continue
        real_pos=real_pos+1
        #print 'Now looking at position',real_pos
        for seq in alignment.keys():
            res_pos=alignment[seq][position]
            if not res_observed.has_key(res_pos):
                res_observed[res_pos]=0
            res_observed[res_pos]=res_observed[res_pos]+1
        #
        # Calculate frequencies of observation
        #
        total=sum(res_observed.values())
        text='%3d' %real_pos
        #print res_observed.keys()
        for aa in aas:
            if res_observed.has_key(aa):
                text=text+',%5.1f' %(float(res_observed[aa])/total*100.0)
            else:
                text=text+',%5.1f' %(0)
        fd2.write(text+'\n')

        #
        #print '%3d   %d' %(real_pos,len(res_observed.keys()))
        lines += ['%3d   %d\n' %(real_pos,len(res_observed.keys()))]
        for mut in res_observed.keys():
            if mut=='-':
                continue
            if mut==HEWL_seq[position]:
                continue
            #if mut=='Y':
            #    continue
            org_res=HEWL_seq[position]
            new_res=mut
            import string
            if org_res=='X' or new_res=='X':
                continue
            #
            # Within the PDB file?
            #
            if real_pos<options.start_aa or real_pos>options.end_aa:
                pass
            else:
                muttext='%s:%s:%s:%s' %(options.CID,string.zfill(real_pos+options.offset-options.noffset,4),convert[org_res],convert[new_res])
                fd.write('%s,%s\n' %(muttext,muttext))
                print muttext
    fd.close()
    fd2.close()
    return

if __name__=='__main__':
    print
    print 'Extract mutations from a multiple sequence alignment'
    print 'Copyright (c) Jens Erik Nielsen, 2009-2010 All rights reserved'
    print
    import sys, os
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] <file>',version='%prog 1.0')
    parser.add_option('-w',"--wt",type='string',dest='wt',action='store',default='HEWL',
                      help='The name of the wild type sequence')
    parser.add_option('-c','--chainID',type='string',dest='CID',action='store',
                      help='The chain ID to mutate in the PDB file. Default: %s',default='')
    parser.add_option('-b','--begin',type='int',dest='start_aa',action='store',
                      help='First residue in PDB file. Default: %default',default=1)
    parser.add_option('-e','--end',type='int',dest='end_aa',action='store',
                      help='Last residue in PDB file. Default: %default',default=9999)
    parser.add_option('-o','--offset',type='int',dest='offset',action='store',
                      help='Offset. This number will be added to the alignment residue number. Default: %default',default=0)
    parser.add_option('-n','--negative_offset',type='int',dest='noffset',action='store',
                      help='Offset. This number will be subtracted from the alignment residue number. Default: %default',default=0)
    (options, args) = parser.parse_args()
    if len(args)!=1:
        parser.error('You must specify the name of a clustal .aln file')
    #
    # Call main
    #
    main(options,args)

