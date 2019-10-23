milliQan Cosmic Muons Simulation
================================

Simulation of Cosmic muons entering the milliQan chamber from different radial distances. The chamber is approximated as a hemisphere of radius 2 meters under 58 meters of rock. The system is assumed to be
azimuthally symmetric. Flux is calculated based on the energy of the cosmic muon and the zenith angle at which it hits the surface of the Earth [\[1\]](#ref1). The program then uses a rock-throwing Monte Carlo
method to generate muons with the calcuated flux distribution.

System Requirements
-------------------

Python 2.7.11 and above (with scipy)

Configuration and Running
-------------------------

[csim.py](../master/csim.py) contains the simulation variables in the section third from the top. That can be modified to select the radial range at the surface where the muons will strike. You can also modify
the resolution at which the surface is sampled by changing the rstep variable. The total number of particles simulated is set by the npart variable. You can also make the program produce a graph of the Monte
Carlo simulation by uncommenting test(I,E_min) on line 98 in the annul function.

References
==========
<a name="ref1"></a>\[1\] [arXiv:1606.06907v3](https://arxiv.org/pdf/1606.06907.pdf)
