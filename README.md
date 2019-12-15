milliQan Cosmic Muons Generation
================================

This program generates Cosmic muons entering the milliQan chamber from different radial distances. The chamber is approximated as a hemisphere of radius 2 meters under 58 meters of rock (center of hemisphere is
60 meters under). The system is assumed to be azimuthally symmetric. Flux is calculated based on the energy of the cosmic muon and the zenith angle at which it hits the surface of the Earth \[1\]. The program 
then uses a rock-throwing Monte Carlo method (Von Neumann method) to generate muons with the calcuated flux distribution.

System Requirements
-------------------

Python 2.7.11 and above (with scipy)

Configuration and Running
-------------------------

[csim.py](../master/csim.py) contains the generation variables in the section third from the top. That can be modified to select the radial range at the surface where the muons will strike. You can also modify
the resolution at which the surface is sampled by changing the rstep variable. The minimum number of particles simulated is set by the npart variable. There is also a section for generating normalized histograms
of muon counts vs energy in each annulus. It creates and saves a numpy array file which can be opened by other programs to analyze. The file contains n-tuples of the form: (E, r, ϕ, θ’, ϕ’, f)

[analysis.py](../master/analysis.py) contains the code to construct plots comparing rates and surface rates at different annuli and different theta primes. It also constructs a histogram of the log E (at milliQan) of all the generated muons. It takes in the numpy array file created by [csim.py](../master/csim.py) for analysis.

[4vector.py](../master/4vector.py) contains the code that converts the n-tuples from [csim.py](../master/csim.py) into kinematic variables. It also plots the 3d scatter and vector plots of the muons when they hit the sphere. It outputs numpy array files for the milliQan kinematic variables as well as the surface kinematic variables. The output n-tuples are of the form: (x, y, z, E, px, py, pz), where the first three elements are the Cartesian coordinates of the muon (at the boundary of the milliQan chamber or at the surface), and the last four elements comprise the Energy-momentum 4 vector of the muons.

References
----------
\[1\] [arXiv:1606.06907v3](https://arxiv.org/pdf/1606.06907.pdf)
