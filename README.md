More informations about this program and its mathematical motivations are available at the homepage  http://www.math.sciences.univ-nantes.fr/~rollin/index.php?page=flow

## Synopsis
The Discrete Moment Map Flow (DMMF) program provides an approximation of the flow defined for discrete surfaces in the four dimensional space computed and reprensented as a  realtime movie. The user can interact with various parameters of the flow and explore its behavior via keyboard and mouse inputs.

## Motivation

The DMMF program implements a particular evolution equation for quadrangular meshes in the four dimensional Euclidean space, associated to a quadrangulation of the torus that came up in our paper Discrete Geometry and Isotropic Surfaces (F. Jauberteau, Y. Rollin, S. Tapie). 

When running the program, the user may see the flow evolving in real time. The program uses a collection of parametrized tori   of R^4 (like the Clifford torus, Chekanov torus, (p,q)-knotted tori, etc...) to define the inital mesh of the flow, at time 0. Since it is difficult to visualize a surface in dimension 4, we choose a radial projection on the 3-sphere followed by a stereographic projection onto the 3-dimensional Euclidean space. The faces of the quadrangular mesh come with a color gradient which depends on their symplectic density. In all experiments so far, we observe that the flow is converging toward a Lagrangian quadrangular mesh (with zero symplectic density).

Limits of the flow provide many examples of discrete immersed Lagrangian surfaces in the four dimensional Euclidean space, at least from an experimental perspective.


## Installation

The program is coded using the Processing language. In order to execute the source code, first install Processing3, which may be downloaded from the website: https://www.processing.org

**The source code (i.e. every \*.pde files) must be  placed in a directory named dmmf**, located in your Processing sketchbook.

Once this is done, launch the Processing IDLE, and go to the menu File->Sketchbook to open the dmmf sketch. Then run the program by pressing the "play" button of the processing IDLE.


## Contributors

Yann Rollin 

Contact: yann.rollin@univ-nantes.fr

## License

GNU GPL v3
