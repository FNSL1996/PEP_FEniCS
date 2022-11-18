# Fluid flow migration and rock stress and deformation fields within a geothermal system due to a crustal fault: a numerical perspective

Felipe Sáez-Leiva (_[@FNSL1996](https://github.com/FNSL1996)_), Daniel E. Hurtado (_[@dehurtado](https://github.com/dehurtado)_), Muriel Gerbault (_[@murielgerbault](https://github.com/murielgerbault)_), Javiera Ruz-Ginouves(_[@JaviRG](https://github.com/JaviRG)_), Pablo Iturrieta (_[@pabloitu](https://github.com/pabloitu)_), José Cembrano.

Repository housing the codes developed for Sáez-Leiva et al., in review.

## Abstract

Geothermal systems are commonly spatially and temporally associated with volcanic complexes, which, in turn, are located nearby crustal fault systems. Such systems can alter fluid flow in its surroundings, potentially acting as barriers or conduits for fluids, depending mainly on its architecture and slip rate. However, the fundamental control a crustal fault exerts on fluid migration, rock deformation and stress state within a geothermal system is still poorly understood.  Numerical studies have the potential to help unraveling the nature of transient processes such as the interplay between a fault and a geothermal reservoir. Most of the previous efforts on volcanic and hydrothermal processes do not directly consider
fluids in their formulations, using often a purely mechanical approach. In this work, we present
a poro-elasto-plastic Finite Element Method (FEM) model to address the first-order, time-dependent control that a strike-slip crustal fault exerts on a nearby geothermal reservoir.  
For the model, we selected the Planch\'on-Peteroa geothermal system in the Southern Andes Volcanic Zone (SAVZ), for which the geometry and kinematics of the system is well-constrained from previous geological and geophysical studies.  We assess the emergence and diffusion
of fluid pressure domains because of fault slip, as well as the development of tensile/dilational and compressive/contractional domains in the surroundings of the fault. We also investigate
the spatial and temporal evolution of such domains resulting from changes in fault permeability and shear modulus, fluid viscosity, and rock rheology. Our
results show the appearance of negative and positive fluid pressure domains due to onset of tensile/dilational and compressive/contractional regions, respectively.
Mean stress and volumetric strain magnitude of the previously described domains range from around $\pm 10^{6}$ [MPa] and $\pm 10^{-4}$ [-]. Distinctive fluid pressure domains alter the trajectory of the reservoir fluids, causing increased migration to the eastern half of the fault, reaching a maximum fluid flux of 6 to 70 times greater than the stationary flux. Pressure-driven fluid diffusion over time causes fluid flow to return to the stationary state between weeks to months after fault slip. These results suggest that one first-order mechanism that controls fluid flow in this setting is a suction pump mechanism driven by fault slip, whose duration would heavily depend on fault permeability and fluid viscosity. The transient processes analyzed in this work highlight the importance of considering fluids directly in numerical modeling for volcano-tectonic studies. 

## Directories

- `benchmarks_data`: Benchmarks results in form of images and raw data
- `codes`: Codes used for poro-elasto-plastic simulation and validation
- `meshes`: Mesh files used for benchmarks and geothermal simulation

## Benchmarks

Two classical benchmarks were used to validate both the poroelastic and the elastoplastic behavior of the model

### Cryer's Problem

Formulated by Cryer (1963), this problem considers a poroelastic sphere of radius *r* subjected to a confining pressure *t*. Free flow condition is imposed on its surface. The traditional concept of this test is shown in Figure 1(a), while the setup used in this test is shown in Figure 1(b), which takes advantage of the symmetry of the original problem to improve computational performance. 

<figure>
<p align="center">
<img src="/benchmarks_data/cryer_setup.png" style="background-color:white;width:50%">
</p>
<figcaption align = "center"><b> Figure 1: (a) Traditional concept for Cryer's Problem. The analytical solution for normalized fluid pressure temporal evolution at the center of the sphere p is known. (b) Setup used in our benchmark. An octant of the sphere, with free-slip and no-flow conditions on the inner surface, is used in order to improve computational performance.</b></figcaption>
</figure>


Normalized fluid pressure results were plotted using [Paraview](https://www.paraview.org/), as shown in Figure 2, while temporal evolution of normalized fluid pressure was plotted against the analytical solution for two mesh sizes (Figure 3).

<figure>
<p align="center">
<img src="/benchmarks_data/cryer_evolution.gif" style="background-color:white;width:70%">
</p>
<figcaption align = "center"><b> Figure 2: Temporal evolution of normalized fluid pressure.</b></figcaption>
</figure>

<figure>
<p align="center">
<img src="/benchmarks_data/cryer_plot.png" style="background-color:white;width:80%">
</p>
<figcaption align = "center"><b> Figure 3: Temporal evolution of normalized fluid pressure at the sphere's center. Two mesh different meshes were used to test mesh dependency. The number of elements for the regular and the refined mesh was 2124 and 14843, respectively.</b></figcaption>
</figure>

### Thick-walled hollow cylinder loading test

In this benchmarks, an elasto-plastic hollow cylinder, confined in the axial direction, is subjected to an internal pressure *t*. As the pressure increases, plastic strain starts to propagate from the inner to the external surface. Analytical solutions, for every internal pressure, are known for:

- Plastic zone depth 
- Tangential stress 
- Radial stress

These solutions can be found in Nadai, 1950. The test setup is shown in Figure 4, and the results are shown in Figures 5, 6, 7, and 8. 

<figure>
<p align="center">
<img src="/benchmarks_data/cylinder_setup.png" style="background-color:white;width:50%">
</p>
<figcaption align = "center"><b> Figure 4: Setup of the hollow cylinder loading test in (a) 3D view, and (b) plan view. Only a quarter of the cylinder is modeled using the problem' symmetry, using the boundary conditions shown in (b). Note that dimensions are not to scale in (b).</b></figcaption>
</figure>

<figure>
<p align="center">
<img src="/benchmarks_data/cylinder_evolution.gif" style="background-color:white;width:60%">
</p>
<figcaption align = "center"><b> Figure 5: von Mises Stress evolution in cylinder wall with increasing internal pressure.</b></figcaption>
</figure>

<figure>
<p align="center">
<img src="/benchmarks_data/cylinder_plasticdepth.png" style="background-color:white;width:80%">
</p>
<figcaption align = "center"><b> Figure 6: Plastic depth evolution at increasing internal pressure. Plastification value of 1 indicates plastic strain development at that distance from the center of the cylinder. A value of 0 indicated elastic behavior.</b></figcaption>
</figure>

<figure>
<p align="center">
<img src="/benchmarks_data/cylinder_radial.png" style="background-color:white;width:80%">
</p>
<figcaption align = "center"><b> Figure 7: Radial stress across the cylinder wall at increasing internal pressure.</b></figcaption>
</figure>

<figure>
<p align="center">
<img src="/benchmarks_data/cylinder_tangential.png" style="background-color:white;width:80%">
</p>
<figcaption align = "center"><b> Figure 8: Tangential stress across the cylinder wall at increasing internal pressure.</b></figcaption>
</figure>

## Results

## Dependencies

- [FEniCS](https://fenicsproject.org/): 2019.1.0 (Legacy version)
- [numpy](https://numpy.org/)
