#
# $Id: large_dataset_manager.py 165 2005-11-23 18:47:46Z nielsen $
#
"""A class for managing very large amount of data"""
#
# Copyright (C) 2005 Jens Erik Nielsen, University College Dublin, All Rigths Reserved
#

class data_manager:

    def __init__(self,dir):
        """Initialise by setting up the dir that holds the data"""
        #
        # We split data into 10000 files per dir
        #
        self.basedir=dir
        self.dir_count=100000
        #
        # We store the data as 'name':filename
        #
        # self.dirs_index holds count of the number of files in each dir
        # Dirs are numbered by 7-digits
        #
        self.dir_digits=7
        return
    
 
    def __call__(self,name):
        return self.get(name)

    #
    # -----
    #

    def get_dir_name(self,name):
        import os
        dir_num=int(float(name)/float(self.dir_count))
        newdir=str(dir_num).zfill(self.dir_digits)
        new_dir=os.path.join(self.basedir,newdir)
        #
        # If the dir does not exist then create it
        #
        if not os.path.isdir(new_dir):
            os.mkdir(new_dir)
        return new_dir

    #
    # -----
    #

    def find_file(self,name,check_it=1):
        """Give the filename, find it in a dir"""
        import os
        newdir=self.get_dir_name(name)
        filename=os.path.join(newdir,str(name))
        if check_it:
            if os.path.isfile(filename):
                return filename
            print 'Could not find',filename
            raise Exception('Crap routine')
        else:
            return filename

    #
    # -----
    #
    
    def get(self,name):
        """Given its name get a value"""
        import os, pickle
        #
        # Can we find this file
        #
        try:
            filename=self.find_file(name)
        except:
            return None
        fd=open(filename)
        value=pickle.load(fd)
        fd.close()
        return value

    #
    # ----
    #

    def put(self,name,value):
        """Store a value"""
        #
        # Do we already have this name?
        #
        import os
        filename=None
        filename=self.find_file(name,None)
        #
        # Open the file and store the data
        #
        fd=open(filename,'w')
        import pickle
        pickle.dump(value,fd)
        fd.close()
        return



