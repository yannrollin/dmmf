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

/**
* Discrete Moment Map Flow program
* This program provides a numerical simulation of 
* an evolution equation for surfaces in the four
* dimensional space.
*
* For more details on the mathematical aspects, see our
* Research paper : Discrete Geometry and Isotropic Surfaces
* by F. Jauberteau, Y. Rollin and S. Tapie.
*/


//Global variables
PShape graph, vertices, group; //for opengl
float xmag, ymag = 0; // Mouse motions
float newXmag, newYmag = 0; //Mouse motions
int framerate = 60; 
//Initialization of the main class
Torus t=new Torus(25);

//This flag is used to trigger intro sequence
Boolean intro_flag=true;

void setup(){
//size(1280,640,P3D);
 fullScreen(P3D);
  frameRate(framerate);
//default text size and color
  textSize(16);
  fill(0); 
 
 }



void keyReleased() {
  switch(key){
      
    case 'm':
    t.movie_stop();
    break;
  }
}

void keyPressed() {
  switch(key){
      
    case ' ':
    t.play_switch();
    break;
  
    case 'p':
    saveFrame("torus-######.tiff");
    break;
 
   case 'm':
   t.movie_rec();
   break;
    
    case 'i':
    t.showinfo_switch();
    break;
  
    case 'r':
    t.restart();
    break;
  
    case 'f':
    t.fill_switch();
    break;

    case 'b':
    t.balls_switch();
    break;

    case 'h':
    t.help_switch();
    break;

    case 't':
    t.type_prev();
    break;
 
    case 'T':
    t.type_next();
    break;

    case 'n':
    t.less_noise();
    break;

    case 'N': t.more_noise();
    break;
    
    case 'q': t.less_quads();
    break;
 
    case 'Q': t.more_quads();
    break;
 
    case 'd': t.decrease_dt();
    break;
  
    case 'D':
    t.increase_dt();
    break;

    case 's': t.scale_down();
    break;

    case 'S':
    t.scale_up();
    break;

    case 'c': t.compute_down();
    break;

    case 'C':
    t.compute_up();
    break;
  }
}

void draw(){  
 
  //detect mouse moves 
  newXmag = mouseX/float(width) * TWO_PI;
  newYmag = mouseY/float(height) * TWO_PI;
  float diff = xmag-newXmag;
  if (abs(diff) >  0.01)  
    xmag -= diff/4.0; 
  diff = ymag-newYmag;
  if (abs(diff) >  0.01)  
    ymag -= diff/4.0; 
    
  // precompose 3d picture with rotations of R4
  // angles depend on mouse motions
  t.resetMatrix();
  t.rotate(1,2,-ymag); 
  t.rotate(1,3,-xmag); 

  //compute the 3 dimensional positions of points via radial projection, then stereographic projection
  t.compute3d();

  //Draw the opengl picture
  
  t.render();
  if(intro_flag==true) {
    intro(); 
  }
  
  //Apply discrete evolution of the flow
  t.evolve();
}

void intro() {
  fill(255,0,0);
  if(millis()<10000){
    text("Discrete Moment Map Flow", 10, height -10 );
  } else if(millis()<15000){
    text("...a program after the work of...", 10, height -10 );    
  } else if(millis()<25000){
    text("François Jauberteau, Yann Rollin and Samuel Tapie", 10, height - 10 );    
  } else if(millis()<30000){
    text("Built with Processing ", 10, height - 10);    
  } else
  intro_flag=false;
}