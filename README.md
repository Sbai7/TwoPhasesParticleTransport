# TwoPhasesParticleTransport

# What's New
* All code functions are now fully compatible with MATLAB (tested on MATLAB R2013a and later) 
* All functions/classes of this toolkit are documented 
* Live-script worked tutorials to get started are provided in the examples sub-folder 
* Dedicated web site which publishes the [worked tutorials online](https://sbai7.github.io/TwoPhasesParticleTransport/)
* Add a new function to drive the single-phase flow pressure solver 
* Fluid phase viscosity could be specified as a handle or anonymous function. This is useful to model coupled miscible flow and transport processes during EOR and aquifer remediation. This is also the case for coupled viscosity-dependent flow and heat transport in geothermal reservoirs when a cooled fluid is reinjected back into the host formation. 



![Alt text](pictures/spe36_results_1m.jpg?raw=true "")

This is a contributed code for numerical modelling of formation damage by two-phase particulate transport processes in heterogeneous porous media. 

It is written with standard [Matlab](http://www.mathworks.com/products/matlab/) high level scripting language. 

The aim of this code is to solve, numerically, the sequentially coupled set of 3D equations governing particlulate transport processes in two-phase flows by means of a standard Finite difference technique. This includes:

* The global pressure equation (for fractional two-phase flow formulation); 
* The saturation equation for slightly compressible flow; 
* A set of multispecies coupled particle convection-diffusion equations; 
* Particles kinetics ODE's;
* Cellwise porosity/permeability change (retroactive effects). 

*citing* 

Please kindly cite this code in your publications if it helps your research:


```
@article{sbai-2011,
  author  = {M. A. Sbai},
  author  = {M. Azaroual}
  title   = {{Numerical modelling of formation damage by two-phase particulate transport processes during CO2 injection in deep heterogeneous porous media}},
  journal = {Advances in Water Resources },
  volume  = {34},
  number  = {1},
  pages   = {62 - 82},
  year    = {2011},
  doi     = {http://dx.doi.org/10.1016/j.advwatres.2010.09.009},
  url     = {http://www.sciencedirect.com/science/article/pii/S0309170810001715},
}
```

# Copyleft 
This is an open source project, it is distributed under the [GPL v3](https://www.gnu.org/licenses/gpl-3.0.html). Anyone is free to use, learn, modify, develop or contribute to the code is welcome. Take a look at the contributing guidelines for starting to contribute to the project.

# Documentation 
Detailed documentation of the code is under preparation. Start by reading the accompagning paper to understand the theory behind it. Examples in this paper are included in the examples subfolder to promote *reproducible science*. 

# FAQ 

## Three-dimensional simulations 

A frequent question is related to 3D support since the code was tested only for one layer models as those presented in the accompagning paper. The code is indeed designed to run fully 3D models. This could be quite demanding computationally though. However, the code is optimally vectorized to achieve the best efficiency. Additional examples will be provided later for 3D geometries.  
