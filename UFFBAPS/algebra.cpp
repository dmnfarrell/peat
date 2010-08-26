/*
 #
 # UFFBAPS - Unified Force Field for Binding And Protein Stability
 # Copyright (C) 2010 Jens Erik Nielsen & Chresten Soendergaard
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

// 3D vectors and matrices
//
//

#include "algebra.h"
#define PI 3.14159265


Vector::Vector()
  {
    x=0.0;
    y=0.0;
    z=0.0;
  }
  
Vector::Vector(float xi, float yi, float zi)
{
  x=xi;
  y=yi;
  z=zi;
}

Vector Vector::operator-(Vector other)
{
  Vector res;
  res.x = x - other.x;
  res.y = y - other.y;
  res.z = z - other.z;
  
  return res;
}

Vector Vector::operator+(Vector other)
{
  Vector res;
  res.x = x + other.x;
  res.y = y + other.y;
  res.z = z + other.z;
  
  return res;
}

void Vector::operator=(Vector other)
{
  x = other.x;
  y = other.y;
  z = other.z;
}

Vector Vector::operator*(float constant)
{
  return Vector(x*constant,y*constant,z*constant);
}


void Vector::normalise()
{
  float len = length();
  x = x/len;
  y = y/len;
  z = z/len;
 
}

         
float Vector::length(){return sqrt(x*x + y*y + z*z);}

float Vector::operator*(Vector other)
{
  return x*other.x + y*other.y + z*other.z; 
}


Vector Vector::operator^(Vector other)
{
  return Vector(y*other.z - z*other.y,
		z*other.x - x*other.z,
		x*other.y - y*other.x);
}



Vector Vector::get_perpendicular()
// return a perpendicular vector
{

  Vector random(1.0, 1.0, 1.0);
  Vector res = (*this)^random;

  if(res.length() < 0.000001)  //random and vector were parallel - try another random
    { 
      Vector random(1.0, 1.0, 0.0);
      res = (*this)^random;
    }

  res.normalise();

  return res;
}

















void Vector::rotate(float angle, Vector p1,Vector p2)
{
  //http://local.wasp.uwa.edu.au/~pbourke/geometry/rotate/

   Vector u,q1,q2;
   float d;

   // translate space so that the rotation axis passes through the origin
   q1.x = x - p1.x;
   q1.y = y - p1.y;
   q1.z = z - p1.z;

   // norm vector
   u.x = p2.x - p1.x;
   u.y = p2.y - p1.y;
   u.z = p2.z - p1.z;
   u.normalise();

   d = sqrt(u.y*u.y + u.z*u.z);

   // rotate space about the x axis so that the rotation axis lies in the xz plane
   if (d != 0) {
      q2.x = q1.x;
      q2.y = q1.y * u.z / d - q1.z * u.y / d;
      q2.z = q1.y * u.y / d + q1.z * u.z / d;
   } else {
      q2 = q1;
   }

   //rotate space about the y axis so that the rotation axis lies along the z axis
   q1.x = q2.x * d - q2.z * u.x;
   q1.y = q2.y;
   q1.z = q2.x * u.x + q2.z * d;

   //perform the rotation by theta about the z axis 
   q2.x = q1.x * cos(angle * PI/180) - q1.y * sin(angle * PI/180);
   q2.y = q1.x * sin(angle * PI/180) + q1.y * cos(angle * PI/180);
   q2.z = q1.z;

   /* Inverse of step 3 */
   q1.x =   q2.x * d + q2.z * u.x;
   q1.y =   q2.y;
   q1.z = - q2.x * u.x + q2.z * d;

   /* Inverse of step 2 */
   if (d != 0) {
      q2.x =   q1.x;
      q2.y =   q1.y * u.z / d + q1.z * u.y / d;
      q2.z = - q1.y * u.y / d + q1.z * u.z / d;
   } else {
      q2 = q1;
   }

   /* Inverse of step 1 */
   x = q2.x + p1.x;
   y = q2.y + p1.y;
   z = q2.z + p1.z;

   if (isnan(x))
     {
       printf( "NAN in algebra:");
       p1.display("p1");
       p2.display("p2");
     }

   //   printf("                                hopsan (%8.2f, %8.2F %8.2f)\n",x,y,z);

   
}







void Vector::rotate(Vector axis, float angle)
{
  /*

tX^2 + c       tXY + sZ          tXZ - sY      0
tXY-sZ         tY^2 + c          tYZ + sX      0
tXY + sY       tYZ - sX          tZ^2 + c       0
0              0                 0             1
Where c = cos (theta), s = sin (theta), t = 1-cos (theta), and <X,Y,Z> is the unit vector representing the arbitary axis

http://www.cprogramming.com/tutorial/3d/rotation.html

  */


  float 
    c = cos(angle * PI/180),
    s = sin(angle * PI/180),
    t = 1-c;

  Matrix4x4 rot(t*x*x + c   , t*x*y + s*z , t*x*z - s*y , (float) 0.0,
		t*x*y - s*z , t*y*y + c   , t*y*z + s*x , (float) 0.0,
		t*x*y + s*y , t*y*z - s*x , t*z*z + c   , (float) 0.0,
		(float) 0.0 , (float) 0.0 , (float) 0.0 , (float) 1.0);



  
  //*this =  mult4(*this, rot);


  Vector res = mult4(*this, rot);
  res.display("In rot");
  

}

void Vector::display(string name)
{
  //  printf("%10s  x:%8.2e  y:%8.2e  z:%8.2e\n",name.c_str(), x, y, z);
  printf("%10s  %8.2f  %8.2f  %8.2f\n",name.c_str(), x, y, z);
  
}






void Vector::test()
{
  Vector v1(1.0, 0.0, 0.0), v2(0.0, 1.0, 0.0);

  v1.display("v1");
  v2.display("v2");

  Vector c = v1^v2;
  
  c.display("c = v1 X v2");
  


}



//////////////////// Matrix implementations ////////////////////


Matrix4x4::Matrix4x4(float a11, float a12, float a13, float a14,
		     float a21, float a22, float a23, float a24,
		     float a31, float a32, float a33, float a34,
		     float a41, float a42, float a43, float a44)
{
  a[0][0] = a11;
  a[0][1] = a12;
  a[0][2] = a13;
  a[0][3] = a14;

  a[1][0] = a21;
  a[1][1] = a22;
  a[1][2] = a23;
  a[1][3] = a24;

  a[2][0] = a31;
  a[2][1] = a32;
  a[2][2] = a33;
  a[2][3] = a34;

  a[3][0] = a41;
  a[3][1] = a42;
  a[3][2] = a43;
  a[3][3] = a44;
  
}




////////////// Methods using both matrices and vectors ///////////////////

Vector mult4(Vector v, Matrix4x4 m)
{

  float 
    x = m.a[0][0]*v.x + m.a[0][1]*v.y + m.a[0][2]* v.z + m.a[0][3]*1.0,
    y = m.a[1][0]*v.x + m.a[1][1]*v.y + m.a[1][2]* v.z + m.a[1][3]*1.0,
    z = m.a[2][0]*v.x + m.a[2][1]*v.y + m.a[2][2]* v.z + m.a[2][3]*1.0;
  //a = m.a[3][0]*v.x + m.a[3][1]*v.y + m.a[3][2]* v.z + m.a[3][3]*1.0;


  return Vector(x,y,z);




}


/////////// Sigmoidal function for desolvation factors ////////////////

Sigmoidal::Sigmoidal()
{

  offset = 0.5;
  factor = 10;
  step = 0.01;

  // precalculate the sigmoidal function
  float x = 0;
  while (x<=1.0)
    {
      precalcs.push_back(calculate(x));
      x += step;
    }
}


float Sigmoidal::value(float desolvation_factor)
{
  if (desolvation_factor<0)
    desolvation_factor=0;
  if (desolvation_factor>1)
    desolvation_factor=1;

  index = (int) (desolvation_factor/step);

  return precalcs[index];
}


float Sigmoidal::calculate(float x)
{  return ((x-offset)*factor/sqrt(1+(x-offset)*factor*(x-offset)*factor)+1)/2;}
