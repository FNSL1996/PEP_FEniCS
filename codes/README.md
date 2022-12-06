# Codes

This folder contains the codes used to benchmark our implementation, as well as the script for the reference model. Details on the poro-elasto-plastic code are detailed on the cryer_test.ipynb notebook. The main outline of the code is represented in Figure 1 through a diagram.

<figure>
<p align="center">
<img src="PEPOutline.png" style="background-color:white;width:50%">
</p>
<figcaption align = "center"><b> Figure 1: Pseudocode of poro-elasto-plastic implementation. FEniCS allows the user to select from a wide variety which finite elements, test and trial functions to use. In our implementation, we use a second and first-order continuous Lagrange element for the solid displacement and the fluid pressure, respectively. Meshes can be either created directly on FEniCS or imported from a mesh generator like Gmsh. System assembly is done automatically by FEniCS when inputting the residual form, the tangent operator, and the boundary conditions. von Mises’ correction mapping was implemented manually using FEniCS objects.</b></figcaption>
</figure>