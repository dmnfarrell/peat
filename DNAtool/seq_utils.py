#! /bin/env python

import sys, string, sys

#
# ---------
#

def get_base_count (seq, base):
    """Get the number of occurrences of the given base(s) in the sequence."""
    count=0 

    if base=='G':
        for b in seq:
            if b=='G':
                count=count+1
    elif base=='C':
        for b in seq:
            if b=='C':
                count=count+1
    elif base=='T':
        for b in seq:
            if b=='T':
                count=count+1
    elif base=='A':
        for b in seq:
            if b=='A':
                count=count+1
    elif base=='GC':
        for b in seq:
            if b=='G' or b=='C':
                count=count+1    
    elif base=='AT':
        for b in seq:
            if b=='A' or b=='T':
                count=count+1
                
    return count  
    

