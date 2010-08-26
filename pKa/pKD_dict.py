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

import UserDict

class sub_dict(UserDict.DictMixin):
    
    def __init__(self,data,parent,parentkey):
        self.data=data
        self.parent=parent
        self.parentkey=parentkey
        return
        
    def __getitem__(self,key):
        if self.data.has_key(key):
            if type(self.data[key]) is type({}):
                return self.__class__(self.data[key],self,key)
            else:
                return self.data[key]
        raise KeyError(key)
        
    def __setitem__(self,key,value):
        self.data[key]=value
        self.parent.__setitem__(self.parentkey,self.data)
        return
        
    def keys(self):
        return self.data.keys()
        
    def has_key(self,key):
        return self.data.has_key(key)
        
    def __delitem__(self,key):
        if self.data.has_key(key):
            del self.data[key]
            self.parent.__setitem__(self.parentkey,self.data)
        else:
            raise KeyError(key)
        return
        

class pKD_dict(UserDict.DictMixin):
    """This class is used to emulate a dictionary for all the potential maps used
    in pKD"""
    
    def __init__(self):
        """Init the pKD_dict"""
        self.data={} 
        import os
        self.datadir=os.path.join(os.getcwd(),'.pKDtmp')
        if not os.path.isdir(self.datadir):
            os.mkdir(self.datadir)
        return

    def __getitem__(self,key):
        """Get an item from disk"""
        if self.data.has_key(key):
            filename=self.data[key]['filename']
            fd=open(filename)
            import pickle
            data=pickle.load(fd)
            fd.close()
        else:
            raise KeyError(key)
        #
        # If a dictionary then return an instannce that will call this __setitem__
        #
        if type(data) is type({}):
            return sub_dict(data,parent=self,parentkey=key)
        return data

    def __setitem__(self,key,value):
        """Save the item to disk"""
        #
        # save as file
        #
        filename=False
        if self.data.has_key(key):
            if self.data[key].has_key('filename'):
                filename=self.data[key]['filename']
        #
        if not filename:
            import tempfile
            oshandle,filename=tempfile.mkstemp(dir=self.datadir)
            import os
            os.close(oshandle)
            self.data[key]={'filename':filename}
        #
        fd=open(filename,'w')
        import cPickle
        cPickle.dump(value,fd)
        fd.close()
        return

    def __delitem__(self,key):
        if self.data.has_key(key):
            filename=self.data[key]['filename']
            import os
            if os.path.isfile(filename):
                os.unlink(filename)
            del self.data[key]
        else:
            raise KeyError('No such key: %s' %key)
        return
    
    def has_key(self,key):
        """Do we have this key"""
        if self.data.has_key(key):
            return True
        return False
    
    def keys(self):
        return self.data.keys()

    def __del__(self):
        """Destructor"""
        for key in self.data.keys():
            self.__delitem__(key)


#
# Test nested dicts
#
if __name__=='__main__':
    D=pKD_dict()
    D['a']={'value':'original','d3':{'haha':'does not work'}}
    D['a']['value']='new'
    print D['a']['value']
    D['a']['d3']['haha']='oh yes'
    print D['a']['d3']['haha']
    print 'keys',D['a']['d3'].keys()
    del D['a']['d3']['haha']
    print D
