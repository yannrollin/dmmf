More informations about this program and its mathematical motivations are available at the homepage  http://www.math.sciences.univ-nantes.fr/~rollin/index.php?page=flow

## Synopsis
The Discrete Moment Map Flow (DMMF) program provides an approximation of the flow defined for discrete surfaces in the four dimensional space computed and reprensented as a  realtime movie. The user can interact with various parameters of the flow and explore its behavior via keyboard and mouse inputs.

## Motivation

This program implements a particular evolution equation for quadrangular meshes in the four dimensional Euclidean space, associated to a quadrangulation of the torus. 
When running the program, the user may see the flow evolving in real time starting. As it is difficult to visualize a surface in dimension 4, we choose a radial projection on the 3-sphere followed by a stereographic projection onto the 3-dimensional Euclidean space. The faces of the quadrangular mesh come with a color gradient which depend on their symplectic density. In all experiments so far, we observe that the flow is converging toward a Lagrangian quadrangular mesh (with zero symplectic density).

This program provides many examples of discrete immersed surfaces in the four dimensional Euclidean space, at least from an experimental perspective. 


## Installation

The program is coded using Processing3. 
In order to execute the source code, you need Processing3P, which may be downloaded on this page: https://www.processing.org

The source code (i.e. every *.pde files) must be  placed in a directory named dmmf, located in your Processing sketchbook.

Once this is done, lauch Processing, go to the menu File->Sketchbook and open the dmmf sketch. Then run the program by pressing the "play" button of the processing IDLE.


## Contributors

Yann Rollin 
Contact: yann.rollin@univ-nantes.fr

## License

GNU GPL v3