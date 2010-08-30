
# This file is part of Protein Engineering Analysis Tool (PEAT)
# (C) Copyright Jens Erik Nielsen, University College Dublin 2003-
# All rights reserved
#
#
import UserDict
import PEAT_dict

class Align_dict(PEAT_dict.PEAT_dict):
    """Keep track of the alignments"""

    def __init__(self,changelist,protein_name,parentapp):
        self.data={}
        self.parentapp=parentapp
        self.changelist=PEAT_dict.Changelist()
        self.protein_name=protein_name
        return

    def __getitem__(self,parent):
        self.changelist=PEAT_dict.Changelist()
        if not self.data.has_key(parent):
            self.data[parent]=self.align(self.protein_name,parent)
        return self.data[parent]

    def align(self,protein,parent):
        print 'Aligning %s %s' %(protein,parent)
        return self.parentapp.isparent(protein,parent)


class Ancestry(dict):
    """On-the-fly alignment of sequences"""

    def __init__(self,parentapp):
        """Set up the container"""
        self.alignments={}
        self.parentapp=parentapp
        return

    def __getitem__(self,protein):
        if not self.alignments.has_key(protein):
            self.alignments[protein]=Align_dict(None,protein_name=protein,parentapp=self.parentapp)
        #
        # Get the alignment
        #
        return self.alignments[protein]

    def keys(self):
        return self.alignments.keys()

    def __repr__(self):
        return str(self.alignments)



