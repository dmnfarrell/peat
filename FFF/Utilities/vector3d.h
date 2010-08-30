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
#ifndef VECTOR3D_H
#define VECTOR3D_H
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string>

#ifndef EPS
#define EPS 1.0E-7
#endif

class vector3d {
 public:
  double x,y,z;
  vector3d(): x(0.0),y(0.0),z(0.0) {};
  //~vector3d() use default
  vector3d(double xx,double yy,double zz): x(xx),y(yy),z(zz) {};
  vector3d(const vector3d &rhs): x(rhs.x), y(rhs.y), z(rhs.z) {}
  void operator=(const vector3d &rhs) { x=rhs.x;y=rhs.y;z=rhs.z;}
  void Normalize() {
    double l=sqrt(x*x+y*y+z*z);
    if (l>EPS) {x/=l;y/=l;z/=l;}
  }
  //
  // A print function
  //
  void print_xyz() {
    cout << x << " " << y << " " << z << endl;
  }
};

//
// Define inline operators for some common vector operations
//
inline vector3d operator+(const vector3d &r1,const vector3d &r2) {
  return vector3d(r1.x+r2.x,r1.y+r2.y,r1.z+r2.z);}
inline vector3d operator-(const vector3d &r1,const vector3d &r2) {
  return vector3d(r1.x-r2.x,r1.y-r2.y,r1.z-r2.z);}

inline vector3d operator%(const vector3d &r1,const vector3d &r2)
 { //cross product
  return vector3d(r1.y*r2.z-r1.z*r2.y,r1.z*r2.x-r1.x*r2.z,r1.x*r2.y-r1.y*r2.x);
  }
inline double operator*(const vector3d &r1,const vector3d &r2) {
  return r1.x*r2.x+r1.y*r2.y+r1.z*r2.z;
  }

inline double Length2(const vector3d &r) {return r.x*r.x+r.y*r.y+r.z*r.z;}
inline double Length(const vector3d &r) {return sqrt(r.x*r.x+r.y*r.y+r.z*r.z);}
inline double Dist2(const vector3d &r1,const vector3d &r2) {return Length2(r2-r1);}
inline double Dist(const vector3d &r1,const vector3d &r2) {return Length(r2-r1);}

inline double Dihedral(const vector3d &i,const vector3d &j,const vector3d &k,const vector3d &l) {
  double dih;
  vector3d jk=k-j;
  vector3d c=(i-j)%jk; //cross product
  c.Normalize();
  vector3d d=(l-k)%jk;
  d.Normalize();
  double scal=c*d;
  if (fabs(scal+1.0)<EPS) dih=180.0;
  else if (fabs(scal-1.0)<EPS) dih=0.0;
  else dih=57.2958*acos(scal);
  double chiral=(c%d)*jk;
  if (chiral<0) dih=-dih;
  return dih;
}

inline void print(vector3d &r) {
  printf ("X: %6.3f, Y: %6.3f, Z: %6.3f\n",r.x,r.y,r.z);
  return;
}

//
// Convert degrees to radians
//
template<class T>
inline T deg2rad(T angle) {
  T rad=angle*atan(1.0)/45.0;
  return rad;
}


#endif
