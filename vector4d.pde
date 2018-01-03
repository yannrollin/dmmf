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

// The class V4d defines a vector in R^4 identified to C^2
// real numbers are represented with double type for more accuracy

class V4d {
// s+it and u+iv are the complex coordinates
// for the identification between C^2 and R^4

 double s ;
 double t ;
 double u ;
 double v ;

V4d(){
  s=0.0;
   t=0.0;
   u=0.0;
   v=0.0;
}

 //this function stores the vector A to the current vector
V4d(V4d A){
 s = A.s; 
 t = A.t; 
 u = A.u; 
 v = A.v; 
}

//this function returns the magnitude of the vector
 double norm(){
   return(sqrt((float)normsq()));
}

//this function returns the square of the magnitude
double normsq(){
  return(s*s+t*t+u*u+v*v);
}

//this function computes the vector AB from the coordinates of A and B
//and stores the resut in the vector
void bipoint(V4d A, V4d B){
 s = B.s- A.s; 
 t = B.t- A.t;
 u = B.u- A.u;
 v = B.v- A.v;
}


//this function adds the vector A to the current vector
void add(V4d A){
 s += A.s; 
 t += A.t; 
 u += A.u; 
 v += A.v; 

}

//this function substract the vector A to the current vector
void sub(V4d A){
 s -= A.s; 
 t -= A.t; 
 u -= A.u; 
 v -= A.v; 

}
//this function stores the vector A to the current vector
void store(V4d A){
 s = A.s; 
 t = A.t; 
 u = A.u; 
 v = A.v; 

}

//this function rotates the vector using the complex structure J
void J(){
 double ss=s,tt=t,uu=u,vv=v;
 s= - tt;
 t= ss;
 u= -vv;
 v= uu;
}

// computes the Euclidean inner product with another vector
double eucl_inner_product(V4d A){
 return(s*A.s + t* A.t + u * A.u + v* A.v); 
}

// computes the symplecti inner product with another vector
// Beware that the function modifies the vector of the class
double sympl_inner_product(V4d A){
  J();
 return(eucl_inner_product(A)); 
}

//multiply the vector by a scalar c
void scal_mult(double c){
 s *= c;
 t *= c;
 u *= c;
 v *= c;
 
}

}