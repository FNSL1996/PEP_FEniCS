# Fluid flow migration and rock stress and deformation fields within a geothermal system due to a crustal fault: a numerical perspective

Felipe Sáez-Leiva (_[@FNSL1996](https://github.com/FNSL1996)_), Daniel E. Hurtado (_[@dehurtado](https://github.com/dehurtado)_), Muriel Gerbault (_[@murielgerbault](https://github.com/murielgerbault)_), Javiera Ruz-Ginouves (_[@JaviRG](https://github.com/JaviRG)_), Pablo Iturrieta (_[@pabloitu](https://github.com/pabloitu)_), José Cembrano.

Repository housing the codes developed for Sáez-Leiva et al., in review.

## Abstract

Geothermal systems are commonly genetically and spatially associated with volcanic complexes, which in turn, are located nearby crustal fault systems. Faults can alter fluid flow in their surroundings, potentially acting as barriers or conduits for fluids, depending on their architecture and slip-rate. However, this fundamental control on fluid migration is still poorly constrained. Most previous modeling efforts on volcanic and hydrothermal processes consider either only fluid flow in their formulations, or only a mechanical approach, and seldom a full, monolithic coupling between both. In this work, we present a poro-elasto-plastic Finite Element Method
(FEM) to address the first-order, time-dependent control that a strike-slip crustal fault exerts on a nearby geothermal reservoir. For the model setting, we selected the Planchón-Peteroa geothermal system in the Southern Andes Volcanic Zone (SAVZ), for which the geometry and kinematics of a potentially seismogenic fault and fluid reservoir is constrained from previous geological and geophysical studies. We assess the emergence and diffusion of fluid pressure domains due to fault slip, as well as the development of tensile/dilational and compressive/contractional domains in the fault’ surroundings. Mean stress and volumetric strain magnitudes in these domains range
between $\pm 1$ [MPa] and $\pm 10^{-4}$, respectively. Our results show the appearance of negative and positive fluid pressure domains in these dilational and contractional regions, respectively. We also investigate the spatial and temporal evolution of such domains resulting from changes in fault permeability and shear modulus, fluid viscosity, and rock rheology. These variations in fluid pressure alter the trajectory of the reservoir fluids, increasing migration to the eastern half of the fault, reaching a maximum fluid flux of 8 to 70 times the stationary flux. Pressure-driven fluid diffusion over time causes fluid flow to return to the stationary state between weeks to months after fault slip. These results suggest that the mechanism that exerts a first-order control is similar to a suction pump, whose duration heavily depends on fault permeability and fluid viscosity. We also show how a von Mises plasticity criterion locally enhances fluid flow. The transient process analyzed in this work highlights the importance of addressing the solid-fluid coupling in numerical models for volcano-tectonic studies.

## Directories

- `benchmarks_data`: Benchmarks results in form of images and raw data
- `codes`: Codes used for poro-elasto-plastic simulation and validation
- `meshes`: Mesh files used for benchmarks and geothermal simulation
- `flux_results`: Modules used in flux processing and notebook with plots used in article 

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

This repository includes the results of the 8 selected models of our article, plus 4 additional tests. Table 1 resumes the models' properties.

|Model Name|Rheology and Heterogeneity|Slip-rate [m/s]|$G_{f}$ [GPa]|$\kappa_{f}$ [m $^{2}$]|$\sigma_{y 0}$ [MPa]|$\mu$ [Pa $\cdot$ s]|
|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|
|Elastic_SR_1|**P-Elastic**, Homogeneous|1|6|1.9 $\cdot$10 $^{-13}$|*|10 $^{-3}$|
|SR_01|P-EPlas, Homogeneous|**0.1**|6|1.9 $\cdot$10 $^{-13}$|5|10 $^{-3}$|
|**SR_1**|P-EPlas, Homogeneous|1|6|1.9 $\cdot$10 $^{-13}$|5|10 $^{-3}$|
|SR_10|P-EPlas, Homogeneous|**10**|6|1.9 $\cdot$10 $^{-13}$|5|10 $^{-3}$|
|Yield_SR_1|P-EPlas, Homogeneous|1|6|1.9 $\cdot$10 $^{-13}$|**5**|10 $^{-3}$|
|Shear_SR_1|P-EPlas, **Compliant Fault**|1|6|1.9 $\cdot$10 $^{-13}$|5|10 $^{-3}$|
|Permeability_SR_1|P-EPlas, **Permeable Fault**|0.1|6|1.9 $\cdot$10 $^{-13}$|*|10 $^{-3}$|
|Viscosity_SR_1|P-EPlas, Homogeneous|1|6|1.9 $\cdot$10 $^{-13}$|4|**10 $^{-2}$**|

## Dependencies

- [FEniCS](https://fenicsproject.org/): 2019.1.0 (Legacy version)
- [numpy](https://numpy.org/)
