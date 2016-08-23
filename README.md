# TwoPhasesParticleTransport

![Alt text](pictures/spe36_results_1m.jpg?raw=true "")

This is a contributed code for numerical modelling of formation damage by two-phase particulate transport processes in heterogeneous porous media. 

It is written with standard [GNU Octave](https://www.gnu.org/software/octave/) (an open source [Matlab](http://www.mathworks.com/products/matlab/) equivalent) high level scripting language. 

The aim of this code is to solve, numerically, the sequentially coupled set of 3D equations governing particlulate transport processes in two-phase flows by means of a standard Finite Volume technique. This includes:

* The global pressure equation (fractional two-phase flow formulation); 
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
Detailed documentation of the code is under preparation. Start by reading the accompagning paper to understand the theory behind it. Examples in this paper are included in the code to promote *reproducible science*. They will be explained and posted to the project wiki pages as soon as possible. 






