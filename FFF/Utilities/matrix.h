/*
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
*/

#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

template<class Type>
class matrix_class {
 public:
  //
  // Default constructor
  //
  matrix_class(): _size1(0),_size2(0),_val() {};

  //
  // Constructor, args: int,int
  //
  matrix_class(const int s1,const int s2) :
    _size1(s1),_size2(s2),_val(s1*s2) {}
  //
  // Constructor: int, int, t
  // the matrix will be filled with the value t
  //
  matrix_class(const int s1,const int s2,const Type t) :
    _size1(s1),_size2(s2),_val(s1*s2,t) {}

  //
  // Using default destructor
  //
  //~matrix_class()

  //
  // Allocate space in val for s1*s2 elements
  //
  void Alloc(const int s1,const int s2) {
    if (s1!=_size1 || s2!=_size2) _val.resize(s1*s2);
    _size1=s1;_size2=s2;
  }

  //
  // Set all elements to t. if _val needs to be resized then do it
  //
  void Alloc(const int s1,const int s2,const Type t) {
    if (s1==_size1 && s2==_size2) SetAll(t);
    else {
      _val.resize(s1*s2);_size1=s1;_size2=s2;
      SetAll(t);
    }
  }

  //
  // set all elements of _val to t
  //
  void SetAll(const Type t) {for (int i=0;i<_size1*_size2;i++) _val[i]=t;}

  //
  // Operators
  //
  Type& operator()(const int s1,const int s2) {
    if (s1>=_size1 ||s2>=_size2) {
      printf ("Error: out of matrix dimension. x: %d y: %d\n",s1,s2);
      exit(0);
    }
    return _val[s1*_size2+s2];
  }

  Type operator()(const int s1,const int s2) const {
    return _val[s1*_size2+s2];
  }

  //
  // Transpose the matrix
  //
  void Transpose(const matrix_class<Type> &rhs) {
    Alloc(rhs.Size2(),rhs.Size1());
    for (int i=0;i<_size1;i++) {
      for (int j=0;j<_size2;j++) {
	(*this)(i,j)=rhs(j,i);
      }
    }
    return;
  }


  //
  // Declaration of the Print function
  //
  void print() {
    for (int i=0;i<_size1;i++) {
      for (int j=0;j<_size2;j++) cout<<"["<<i<<"]["<<j<<"]="<<_val[i*_size2+j]<<" ";     cout<<endl;
    }
    return;
  };


  //
  // Return the sizes
  //
  int Size1() const {return _size1;}
  int Size2() const {return _size2;}

 protected:
  int _size1,_size2;
    std::vector<Type*> _val;
};

#endif
