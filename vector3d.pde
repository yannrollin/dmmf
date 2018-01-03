/**
Copyright 2017-2018 Yann Rollin
Contact yann.rollin@univ-nantes.fr

  This file is part of DMMF (Discrete Moment Map Flow program).

    DMMF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DMMF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DMMF.  If not, see <http://www.gnu.org/licenses/>.
**/    

// The class V3d defines a vector in R^3

class V3d {


 float x ;
 float y ;
 float z ;


  V3d(){
    x=0.0;
     y=0.0;
     z=0.0;
  }

/*
This function rotates the vertex vert using the rotation matrix. 
The vertex is then projected on the 3-sphere then on R^3 via stereographic projection
*/
void compute(V4d vert, float[][] matrix) {
  double ts, tt, tu, tv,r;
  //apply the rotation matrix and store the result in (ts,tt,tu,tv)
  ts = matrix[0][0] * vert.s + matrix[0][1] * vert.t + matrix[0][2] * vert.u + matrix[0][3] * vert.v ;
  tt = matrix[1][0] * vert.s + matrix[1][1] * vert.t + matrix[1][2] * vert.u + matrix[1][3] * vert.v ;
  tu = matrix[2][0] * vert.s + matrix[2][1] * vert.t + matrix[2][2] * vert.u + matrix[2][3] * vert.v ;
  tv = matrix[3][0] * vert.s + matrix[3][1] * vert.t + matrix[3][2] * vert.u + matrix[3][3] * vert.v ;

  //Project the point onto the unit 3-sphere
  r=vert.norm();
  ts = ts / r;
  tt = tt / r ; 
  tu = tu / r ; 
  tv = tv / r ; 

  //perform a stereographic projection
  
  x= (float)(ts *2 /(1-tv));
  y= (float)(tt *2 /(1-tv));
  z= (float)(tu *2 /(1-tv));
  
}

}