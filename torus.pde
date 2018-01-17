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

//This is the main class for 
//storing a quadrangular mesh
//apply the evolution equation
//rendre the opengl picture
//and interact with parameters.

class Torus {
 V4d[][]  vertex ; // position of the vertices of the mesh
 V4d[][]  motionvector ; // this vector for each vertex is used to move the quadrangulation
 V3d[][]  vertex3d ; // coordinates projected in R3
 double[][] phi; // a function on faces that contains the symplectic area of each face
 float[][] matrix; // 4d orthogonal matrix used to precompose with 3d projections
 float timer=0.0;//timeline parameter
 double mu_min=0.0;//minimum symplectic density
 double mu_max =0.0; //maximum symplectic density

 //most important parameters
 int N =20; //size of the quadrangulation with N^2 quads
 int N_MAX=100; // for safety N must be smaller than N_MAX

 int Delta_pow = -2;
 double Delta = pow(10.0,Delta_pow); // dt time step for the flow
 float noise = 0.0; // random noise parameter for initial torus
 int torus_type = 0; // Clifford = 0 ; Chekanov = 1, etc..
 int general_scale = 50 ; // overall scaling factor

 //Accelaration factor number of streps of 
 //computation the evolution flow
 //before drawing a picture
 int  Boost_comp= 1;
 
 //Overall number of initial tori implemented here 
 int number_torus_types = 9; 
 
 //some colors presets to make nice gradients as a function of 
 //symplectic density
 color isotropic = color(0,0,255);
 color isointerm = color(255,255,0);
 color notiso = color (255,0,0);
 float thr_max=1.0;
 float thr_interm = .5;
 float thr_min=.001;

 //vertices
 int ball_size =10;
 color black=color(0,0,0);
 color far = color(255,0,0);
 color close = color(0,0,255);
 color interm = color(0,255,0);
 float thr_balls_max=4.0;
 float thr_balls_interm = 2.0;
 float thr_balls_min=.1;



 //some flags for the running status
 boolean run = true; //running state (run if true or stop if false)
 boolean showinfo = false; //toggle show information
 boolean restart =true; // if true, generate a new surface
 boolean fill =true; // if true, fill the triangles
 boolean balls = false; // if true draw colored balls at the vertices
 boolean help = false; // if true show help (and override info)
 boolean movie_flag = false ; //if true save all pictures as a continuous movie
 //constructor
 Torus(){
   matrix = new float[4][4]; //rotation matrix to apply to each point   
   resetMatrix();
   resetvertex();
 }
 //another constructor 
 Torus(int k){
   N=k;
   matrix = new float[4][4]; //rotation matrix to apply to each point   
   resetMatrix();
   resetvertex();
 }
 
//Setup the mesh and allocate memory
 void resetvertex(){
   vertex = new V4d[N][N];
   motionvector = new V4d[N][N];
   vertex3d = new V3d[N][N];
   phi = new double[N][N];

   for(int n=0 ; n<N;n++) {
   for(int m=0; m<N;m++)  {
    vertex[n][m]= new V4d();   
    motionvector[n][m]= new V4d();   
    vertex3d[n][m]= new V3d();   
     }
   }
 }
 
 //this function set the rotation matrix to the identity
 void resetMatrix(){
 for(int i = 0;i<4;i++){ 
   for (int j = 0; j<4;j++){
   if(i==j) 
     matrix[i][j]=1.0;
     else
      matrix[i][j]=0.0;
     }
   }
 }
 
//mutliply the rotation matrix on the left by a rotation in the (i,j)-plane
 void rotate(int i, int j, float phi){
  float [][] rmatrix = new float[4][4];
  float [][] oldmatrix = new float[4][4];
   
  //initialize the temporary matrix as the indentity
  // copy matrix in oldmatrix
  for(int ii=0;ii<4;ii++) {
     for(int jj=0;jj<4;jj++){
       oldmatrix[ii][jj]= matrix[ii][jj];
       if(ii==jj)
         rmatrix[ii][jj]=1.0;
       else
         rmatrix[ii][jj]=0.0;
       }
    }

  //modify the temporary matrix in a rotation matrix of angle phi
  rmatrix[i][i]=cos(phi);
  rmatrix[i][j]=sin(phi);
  rmatrix[j][i]=-sin(phi);
  rmatrix[j][j]=cos(phi);

  // perform matrix multiplication
  for(int ii=0;ii<4;ii++){
    for(int jj=0;jj<4;jj++){
      matrix[ii][jj]=rmatrix[ii][0]* oldmatrix[0][jj] 
        + rmatrix[ii][1]* oldmatrix[1][jj] 
        +rmatrix[ii][2]* oldmatrix[2][jj] 
        +rmatrix[ii][3]* oldmatrix[3][jj];
       }
    }
 }
 

//Implementation of several parametrizations of tori

void create_clifford_Torus(){ 
    for (int n=0; n<N; n++){
     for (int m=0; m<N; m++){
       vertex[n][m].s=cos(2*n*PI/N)+ (random(noise)-noise/2)/N;
       vertex[n][m].t=sin(2*n*PI/N)+ (random(noise)-noise/2)/N;
       vertex[n][m].u=cos(2*m*PI/N)+ (random(noise)-noise/2)/N;
       vertex[n][m].v=sin(2*m*PI/N)+ (random(noise)-noise/2)/N;
       }
    }
 }
 
void create_knot_Torus(int p, int q,float r)
  { 
    for (int n=0; n<N; n++){
     for (int m=0; m<N; m++){
       vertex[n][m].s=cos(2*n*PI/N*p)+ (random(noise)-noise/2)/N;
       vertex[n][m].t=sin(2*n*PI/N*p)+ (random(noise)-noise/2)/N;
       vertex[n][m].u=cos(2*n*PI/N*q)+cos(2*m*PI/N)*r+ (random(noise)-noise/2)/N;
       vertex[n][m].v=sin(2*n*PI/N*q)+sin(2*m*PI/N)*r+ (random(noise)-noise/2)/N;
       }
    }
  }
 
 //This knotted torus is a spin knot on q (p,q)-torus knot
 //This is one of the luttinger examples which does not
 //admit a lagrangian type
 void create_spin_knot_Torus(int p, int q,float r)
  { 
    float x,y,z;
    for (int n=0; n<N; n++){
     for (int m=0; m<N; m++){
       x= cos(2*m*p*PI/N) * (cos(2*m*q*PI/N)+ 2) ;
       y= sin(2*m*p*PI/N) * (cos(2*m*q*PI/N)+ 2 ) ;
       z= -sin(2*m*PI/N*q) ;
       vertex[n][m].s= cos(2*n*PI/N)+ (random(noise)-noise/2)/N;
       vertex[n][m].t= sin(2*n*PI/N)* (x+5) + (random(noise)-noise/2)/N;
       vertex[n][m].u=  y + (random(noise)-noise/2)/N;
       vertex[n][m].v= z +  (random(noise)-noise/2)/N;
       }
    }
  }
 
 
 //this function creates a flat torus contained in a complex line
void create_complex_Torus(){ //noise is the level of noise of the standard clifford torus
      for (int n=0; n<N; n++){
         for (int m=0; m<N; m++){
           vertex[n][m].s= (n*2.0-N)/N+ (random(noise)-noise/2)/N;
           vertex[n][m].t= (m*2.0-N)/N+ (random(noise)-noise/2)/N;
           vertex[n][m].u= (random(noise)-noise/2)/N;
           vertex[n][m].v= (random(noise)-noise/2)/N;
           }
        }
   }
 
void create_complex_Torus_B(){   
    for (int n=0; n<N; n++){
     for (int m=0; m<N; m++){
       vertex[n][m].s= (n*2.0-N)/N+ (random(noise)-noise/2)/N;
       vertex[n][m].t= (m*2.0-N)/N+ (random(noise)-noise/2)/N;
       vertex[n][m].u= 1.0+(random(noise)-noise/2)/N;
       vertex[n][m].v= (random(noise)-noise/2)/N;//random(noise);
       }
      }
   }
 
void create_double_circle() {
  float a= 2.0;
  float b= 2.0;
  for (int n=0; n<N; n++){
   for (int m=0; m<N; m++){
     vertex[n][m].s= - a +a * cos(PI+2*n*PI/N)+ random(noise)/N;
     vertex[n][m].t= b *sin(PI+2*n*PI/N)+ random(noise)/N;
     vertex[n][m].u=- a +a * cos(2*m*PI/N)+ random(noise)/N;
     vertex[n][m].v=  b *sin(2*m*PI/N)+ random(noise)/N;
   }
  }
}

void create_Checkanov_Torus() {
  float xc= 1.2;
  float yc= 1.2;
  for (int n=0; n<N; n++){
   for (int m=0; m<N; m++){
       vertex[n][m].s=cos(2*PI*m/N) * (cos(2*n*PI/N)+xc)
       -sin(2*PI*m/N) * (sin(2*n*PI/N)+yc)+ (random(noise) - noise/2)/N +1.5;
       vertex[n][m].t =cos(2*PI*m/N) * (sin(2*n*PI/N)+yc)
       +sin(2*PI*m/N) * (cos(2*n*PI/N)+xc)+ (random(noise)-noise/2)/N;
       vertex[n][m].u=cos(-2*PI*m/N) * (cos(2*n*PI/N)+xc)
       -sin(-2*PI*m/N) * (sin(2*n*PI/N)+yc)+ (random(noise)-noise/2)/N;
       vertex[n][m].v=cos(-2*PI*m/N) * (sin(2*n*PI/N)+yc)
       +sin(-2*PI*m/N) * (cos(2*n*PI/N)+xc)+ (random(noise)-noise/2)/N;
     }
  }
}



// This function computes the 3-dimensional representation of the torus
void compute3d() {
    for( int n =0;n<N;n++) {
      for (int m=0;m<N;m++) {
        vertex3d[n][m].compute(vertex[n][m], matrix);
        }
      }
  }
 
 // this funtion computes the symplectic area of each face
 // and stores the result in the member phi of the class
void computephi(){

  boolean flag_init = false;
  double mu;

  mu_min=0.0;
  mu_max=0.0;

   for( int n =0;n<N;n++) {
      for (int m=0;m<N;m++) {
       V4d ac= new V4d();
       V4d bd= new V4d();
       int nn = next_n(n);
       int nm = next_m(m);

       ac.bipoint(vertex[n][m],vertex[nn][nm]);
       bd.bipoint(vertex[nn][m],vertex[n][nm]);
       //   compute the symplectic area
       phi[n][m] = -.5 * ac.sympl_inner_product(bd);
       //compute the symplectic density
       mu=phi[n][m]*N*N;

      if(flag_init) {
        if(mu<mu_min)
          mu_min=mu;
        else if (mu>mu_max)
          mu_max=mu;
        } else
      {
        flag_init=true;
        mu_min = mu;
        mu_max = mu;
      }
    }
   }
 }

//This function computes the evolution vector of the flow
void computemotionvector(){
    int nn,nm, pn,pm;
    V4d[] D = new V4d[4];
    D[0]= new V4d();
    D[1]= new V4d();
    D[2]= new V4d();
    D[3]= new V4d();

    //We compute first the symplectic area of each face
    computephi();
    //Then we compute the vector
    for( int n =0;n<N;n++) {
      for (int m=0;m<N;m++) {
        pn=prev_n(n);
        pm=prev_m(m);
        nn=next_n(n);
        nm=next_m(m);
  
        D[0].bipoint(vertex[nn][m], vertex[n][nm]);
        D[1].bipoint(vertex[n][nm], vertex[pn][m]);
        D[2].bipoint(vertex[pn][m], vertex[n][pm]);
        D[3].bipoint(vertex[n][pm], vertex[nn][m]);

        D[0].scal_mult(phi[n][m]);
        D[1].scal_mult(phi[pn][m]);
        D[2].scal_mult(phi[pn][pm]);
        D[3].scal_mult(phi[n][pm]);

        D[0].add(D[1]);
        D[0].add(D[2]);
        D[0].add(D[3]);

        D[0].J();

        motionvector[n][m].store(D[0]);
      }
    }
 }
  
  
// Some color gradient functions for faces and vertices

color color_gradient_balls(float f){
  if (f>thr_balls_max) 
    return(far);
  if (f<thr_balls_min)
    return(black);
  if (f <thr_balls_interm)
    return(lerpColor(close,interm,f/thr_balls_interm));
  else
    return(lerpColor(interm,far,(f-thr_balls_interm)/(thr_balls_max-thr_balls_interm)));
}

color color_gradient(float f){
  if (f>thr_max) 
    return(notiso);
  if (f<thr_min)
    return(black);
  if (f <thr_interm)
    return(lerpColor(isotropic,isointerm,f*2));
  else
   return(lerpColor(isointerm,notiso,(f-thr_interm)*2));
}



  
//this function evolves the quadrangulation along the motion vector
  
void evolve(){  
  //check if a restart was requested
  if(restart) {
   restart=false;
   switch(torus_type) {
    case 0: create_clifford_Torus();
    break;

    case 1: create_Checkanov_Torus();
    break;

    case 2: create_complex_Torus();
    break;

    case 3: create_complex_Torus_B();
    break;
    
    case 4: create_knot_Torus(2,3,.7);
    break;
    
    case 5: create_knot_Torus(2,5,.6);
    break;
    
    case 6: create_knot_Torus(3,5,.5);
    break; 
      
    case 7: create_knot_Torus(5,7,.3);
    break;
 
    case 8: create_spin_knot_Torus(2,3,2.0);
    break;
 
   }
   computephi();
   timer=0;
  }

  for(int q =0; q<Boost_comp ;q++) {
    if(run){
      computemotionvector();
      for( int n =0;n<N;n++) {
        for (int m=0;m<N;m++) {
          motionvector[n][m].scal_mult(Delta);
          vertex[n][m].add(motionvector[n][m]);
          }
        }
      timer += Delta/(N*N*N*N);
        }
    }
}
  
  
 //the functions next_n and next_m return the correct index modulo N and modulo M
int next_n(int i){
  if (i==(N-1)) return (0);
    else return(i+1);
}

int next_m(int i){
  if (i==(N-1)) return (0);
    else return(i+1);
}
  
  //the functions prev_n and prev_m return the correct index modulo N and modulo M
int prev_n(int i){
  if (i==0) return (N-1);
    else return(i-1);
}
  
int prev_m(int i){
  if (i== 0) return (N-1);
    else return(i-1);
}
  
//rendering the opengl surface
void  render() {
 //graphic opengl rendering of the 3dimensional picture  
 graph =createShape();
 if(fill) {
 graph.beginShape(QUADS);
 graph.stroke(200);
 
 for (int n = 0; n<N ; n++) {
  for (int m = 0 ; m<N ; m++) {
    int n_next,m_next;
    if (n==(N-1)) 
      n_next = 0;
    else 
      n_next = n+1;
    if (m==(N-1)) 
      m_next = 0;
    else 
      m_next = m+1;
 
 
 graph.fill(color_gradient(N*N*abs((float)t.phi[n][m])));
  
 graph.vertex(t.vertex3d[n][m].x, t.vertex3d[n][m].y, t.vertex3d[n][m].z);
 graph.vertex(t.vertex3d[n_next][m].x, t.vertex3d[n_next][m].y, t.vertex3d[n_next][m].z);
 graph.vertex(t.vertex3d[n_next][m_next].x, t.vertex3d[n_next][m_next].y, t.vertex3d[n_next][m_next].z);
 graph.vertex(t.vertex3d[n][m_next].x, t.vertex3d[n][m_next].y, t.vertex3d[n][m_next].z);
  }
 }
 
 graph.endShape();
  
}
 else {
   graph.beginShape(QUADS);
   graph.noFill();
   graph.stroke(0);
   for (int n = 0; n<N ; n++) {
     for (int m = 0 ; m<N ; m++) {
      int n_next,m_next;
      if (n==(N-1)) 
       n_next = 0;
      else 
       n_next = n+1;
      if (m==(N-1)) 
        m_next = 0;
      else 
        m_next = m+1;
 
       graph.vertex(t.vertex3d[n][m].x, t.vertex3d[n][m].y, t.vertex3d[n][m].z);
       graph.vertex(t.vertex3d[n_next][m].x, t.vertex3d[n_next][m].y, t.vertex3d[n_next][m].z);
       graph.vertex(t.vertex3d[n_next][m_next].x, t.vertex3d[n_next][m_next].y, t.vertex3d[n_next][m_next].z);
       graph.vertex(t.vertex3d[n][m_next].x, t.vertex3d[n][m_next].y, t.vertex3d[n][m_next].z);
      }
   }
   graph.endShape();
 }
   
  
 if(balls) //if balls=true draw colored balls at each vertex
   {
    vertices = createShape();
    vertices.beginShape(POINTS);
    vertices.strokeWeight(ball_size);
    for (int n = 0; n<N ; n++) {
      for (int m = 0 ; m<N ; m++) {
        //define the color of verticces accordint to their distance to the origin
        float r = (float)t.vertex[n][m].normsq();
        vertices.stroke(color_gradient_balls(r));
        vertices.vertex(t.vertex3d[n][m].x,t.vertex3d[n][m].y,t.vertex3d[n][m].z);
        }
     }
    vertices.endShape();
  }

 group = createShape(GROUP);
 if(balls)  
   group.addChild(vertices);
 group.addChild(graph);
 
 background(255);
 pushMatrix(); 
 translate(width/2, height/2, -30);  
 scale(general_scale);
 shape(group);
 popMatrix();

 if(t.help){
   fill(0);
   text("Help -- Key functions", 10, 20);
   text("Space : Play/Pause", 10, 40);
   text("i : Toggle info", 10, 60);
   text("h : Toggle help", 10, 80);
   text("r : Restart", 10, 100);
   text("b : Toggle balls", 10, 120);
   text("f : Toggle fill", 10, 140);
   text("t : Switch torus type", 10, 160);
   text("d/D : Time step -/+", 10, 180);
   text("n/N : Noise level -/+", 10, 200);
   text("q/Q : Number of quads -/+", 10, 220);
   text("s/S : Overall scale -/+", 10, 240);
   text("c/C : Change computation factor -/+", 10, 260);
   text("p : Save screenshot", 10, 280);
   text("m : Press continuously to shoot a movie", 10, 300);

  }
  else if(t.showinfo){
    fill(0);
   text("Frame rate: " + int(frameRate), 10, 120);
   text("Time : " + timer, 10, 40);
   text("Time step : " + Delta/(N*N*N*N), 10, 60);
   text("Number of quads : " + int(N*N), 10, 80);
   text("Noise level : " + noise, 10, 100);
   switch(torus_type){
    case 0:
      text("Clifford torus", 10, 20);
      break;
    case 1:
      text("Chekanov torus", 10, 20);
      break;
    case 2:
      text("Complex torus type A", 10, 20);
      break;
   case 3:
      text("Complex torus type B", 10, 20);
      break;
   case 4:
      text("Tubular trefoil knot", 10, 20);
      break;
   case 5:
      text("Tubular (2,5)-knot", 10, 20);
      break;
   case 6:
      text("Tubular (3,5)-knot", 10, 20);
      break;
   case 7:
      text("Tubular (5,7)-knot", 10, 20);
      break;
   case 8:
      text("Spin torus on a trefoil knot", 10, 20);
      break;

   }
  text("Overall scale : " + general_scale, 10, 200);
  text("Minimal symplectic density: " + mu_min, 10, 160);
  text("Maximal symplectic density: " + mu_max, 10, 180);
  text("Flow terations per frame: " +  Boost_comp, 10, 140);
    }
    
 //check if a movie has to be recorded
 if(movie_flag)
   saveFrame("movie-########.tiff");
}

//Some toggle and incremental functions for parameters
void play_switch(){
       run = ! run;
  }

void showinfo_switch(){
    showinfo= ! showinfo; 
  }

void restart(){
   restart=true; 
  }
  
void fill_switch(){
    fill = ! fill; 
  }   
 
void balls_switch(){
   balls=!balls; 
  }

void help_switch(){
 help=!help; 
  }

void type_next(){
  torus_type++;
  if(torus_type==number_torus_types)
      torus_type=0;  
   restart();
}

void type_prev(){
  torus_type--;
  if(torus_type==-1)
      torus_type=number_torus_types-1;  
  restart();
}

void less_quads(){
 if (N >2) {
      N--;
       resetvertex();
       restart(); 
    } 
}

void more_quads(){
 if (N < N_MAX) {
      N++;
       resetvertex();
       restart(); 
    } 
}

void decrease_dt(){
    if (Delta_pow > -10) {
      Delta_pow--;
      Delta = pow(10.0,Delta_pow); 
    }
  
}
void increase_dt() {
    if (Delta_pow < 2) {
      Delta_pow++;
      Delta = pow(10.0,Delta_pow); 
    } 
}

void less_noise(){
     if (noise < .5)
      noise=0.0;
      else noise-=.5;
    restart(); 
}

void more_noise(){
     if (noise < 50.0)
      noise+=.5;
    restart(); 
}

void scale_down(){
     if (general_scale > 5) {
      general_scale-=4;
    }
}

void scale_up(){
     if (general_scale < 100) {
      general_scale+=4;
    }
}

void compute_down(){
    if (Boost_comp > 30) {
      Boost_comp-=10;
    }
    else if(Boost_comp>1)
      Boost_comp--;    
}

void compute_up(){
     if (Boost_comp >= 30) 
      Boost_comp+=10; 
      else Boost_comp++;
    }
void movie_rec(){
    movie_flag=true;
    }
void movie_stop(){
    movie_flag=false;
    }
    
    
}