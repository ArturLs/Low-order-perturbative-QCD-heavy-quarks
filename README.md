# Low-order-perturbative-QCD-heavy-quarks


### Introduction

The low-order perturbative Quantum Chromodynamics (QCD) theory describes the interaction between heavy particles, such as charm and bottom quarks, and the strong interaction mediated by gluons. The differential cross section with respect to the squared transverse momentum, $\frac{d\sigma}{dP^2_t}$, is an important measure in this context, as it describes the probability of interaction between quark-antiquark pairs and the gluons of the surrounding quantum field. For charm and bottom quarks, the differential cross section can be calculated using perturbative techniques, taking into account contributions from different orders of strong coupling. These calculations are crucial for understanding and predicting the production processes of heavy particles in high-energy collisions, such as those performed in particle accelerators like the LHC.

-----------------------------------------------------------------------------
### Instructions

- execute the simulation:

``g++ -g -o `root-config --cflags` `lhapdf-config --cflags` `gsl-config --cflags` res.cpp `root-config --glibs` `lhapdf-config --libs` `gsl-config --libs` -I ./ -lstdc++ -lm -lMathMore -lLHAPDF -o res && ./res``

 - If you want to change the quark:

 If you want to change the quark, simply modify the first argument, `double params[pdim] = {4.5, pT};`, in line 364. For example: 
 
  `` Bottom: 4.5 GeV  Charm: 1.27 GeV ``


-----------------------------------------------------------------------------
### Simulation<br><br>






![Screenshot from 2023-05-18 23-34-54](https://github.com/ArturLs/Low-order-perturbative-QCD-heavy-quarks/blob/main/Bottom.png)
