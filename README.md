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

# FAQ 

## MATLAB support 

The most frequently asked question we get concerns [MATLAB](http://www.mathworks.com/products/matlab/) compatibility. The only incompatible portion of the code is related to the use of the DAE solver [DASPK](https://www.gnu.org/software/octave/doc/v4.0.1/Differential_002dAlgebraic-Equations.html) (from the particles kinetics routines) which is supported under OCTAVE but not in MATLAB. A simple fix consists on replacing the corresponding calls with a stiff ODE in MATLAB such as [ode15s](http://www.mathworks.com/help/matlab/ref/ode15s.html). We will test this and make a special branch for MATLAB support as soon as we could. 

## Three-dimensional simulations 

Another question is related to 3D support since the code was tested only for one layer models as those presented in the accompagning paper. The code is indeed designed to run fully 3D models. This could be quite demanding computationally though. However, the code is optimally vectorized to achieve the best efficiency. Additional speedup should be achieved by other means. 





