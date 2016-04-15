DiffStress
==========

DiffStress is a Python-based diffraction analysis tool that 
analyze a *tensorial* stress pertains to a set of lattice spacings 
(interatomic spacings) obtained from a polycrystalline specimen 
using a diffractometer.


This python package allows one to analyze both experimental and 
model-predicted diffraction lattice strains in order to estimate 
a tensorial stress.

For the purpose of experimental d-spacing analysis, it accepts
the format of Proto data. Current software also enables to plot
various graphics for the purpose of data visualization by using
matplotlib and some in-house scripts written for matplotlib, which
is also available under current GitHub account (fork "mpl-lib").

In order to estimate the experimental/model-predicted stress,
current software requires the use of diffraction elastic constants written
in the format of 'sff', which is used in the 'PF' software
managed by Thomas Gnaeupel-Herold in NCNR, NIST. One may find the
template of the 'sff' file from the current software package.


This software was used for the following publications
-----------------------------------------------------
1. Uncertainty in flow stress measurements using X-ray diffraction for sheet metals subjected to large plastic deformations, Y Jeong, T Gnaeupel-Herold, M Iadicola, A Creuziger (Submitted to Journal of Applied Crystallography)
2. Multiaxial constitutive behavior of an interstitial-free steel: measurements through X-ray and digital image correlation, Y Jeong, M Iadicola, T Gnaeupel-Herold, A Creuziger (Accepted for publication in Acta Materialia)
3. Evaluation of biaxial flow stress based on Elasto-Viscoplastic Self-Consistent analysis of X- ray Diffraction Measurements, Y Jeong, T Gnaeupel-Herold, F Barlat, M Iadicola, A Creuziger, M-G Lee (2015) International Journal of Plasticity, 66, 103-118
