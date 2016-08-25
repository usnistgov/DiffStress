DiffStress
==========
DiffStress is a Python-based diffraction analysis tool that
analyzes the *tensorial* stress state present on a polycrystal sample.
The stress analysis utilizes a set of lattice-spacings (interatomic
spacings) obtained from a polycrystalline specimen using a diffractometer
(either X-ray or neutron)


This python package allows one to analyze both experimental diffraction
lattice strains (d-spacing) and model-emulated lattice strains (internal
elastic strain), e.g., those data resulting from EVPSC simulations.


For the purpose of experimental d-spacing analysis, it accepts the format
of Proto data. Current software can generate various graphics for the
purpose of data visualization by using matplotlib and some in-house
scripts written for matplotlib, which is also available under current
GitHub account (find "mpl-lib").


In order to estimate the experimental/model-predicted stress, current
software requires the use of diffraction elastic constants written
in the format of 'sff', which is used in the 'PF' software managed by
Thomas Gnaeupel-Herold in NCNR, NIST. One may find the template of the
'sff' file from the current software package.


Also, given set of diffraction data (say, from EVPSC code), one can
perturb the lattice strains in order to quantify the uncertainty in the
obtained stress by conducting Monte Carlo virtual experiments.


Application
===========

- Stress analysis from d-spacings

DiffStress can analyze the d-spacings collected in the form of Proto X-ray system.
One can obtain a fitting result in various orientation and nicely illustrate the 
goodness of fit. Below is an example of an interstitial-free steel, that shows the
undulation in d-spacing vs. sin2psi curve. With using experimental DECs obtained
at various orientation, one can obtain stress by least-square method.
The intensity can be used to 'estimate' the size (or height) of the diffraction peaks.


![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/protoEx.png)


- Monte Carlo experiments for uncertainty estimation

DiffStress comes with a feature that allows the Monte Carlo experiment
to quantify uncertainties present in flow-stress measurement technique
using in-situ diffraction experiments. Below image shows the resulting
internal elastic strain (that is corresponding to d-spacing vs. sin2psi
curve in a real in-situ diffraction experiment) and the fitting result
(like the conventional sin2psi method widely used in the community of
residual stress measurement).


![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/illu_1.png)



The stress measurement using the internal elastic strain (d-spacing)
should lead to a *self-consistent* result such that the weighted
average of stress (the macro stress) is equivalent to the stress
measured by this virtual diffraction experiment. Below figure shows that
the two flow stress curves (one with weighted average stress another
obtained by this virtual diffraction experiment) are equivalent.


![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/illu_1f.png)


In the real experiments, however, this self-consistency may be
challenged by various factors. For example, due to the counting
statistical nature and a finite period of exposure time in use when
determining the peak position (thus d-spacings), the d-spacing obtained
by a diffraction peak may contain a degree of uncertainty. Also, the peak
height (the total counts collected at a particular orientation) is
influenced by the presence of crystallographic (and its evolution w.r.t p
lastic deformation) thus further affecting the uncertainty in the measured
d-spacing. In DiffStress, one can mimic various types of uncertainties
existing in real experiments and can simulate the propagation of these
uncertainties to the final stress estimated using the virtual diffraction
experiments. The first procedure is to superimpose these uncertainties
to the internal-strain used in the virtual diffraction experiment. The
internal-strain is then *perturbed* by counting statistical error,
crystallographic texture, and incomplete measurements of diffraction
elastic constants, and the finite number of exposure to X-ray beam. By
using the *perturbed* internal strain, the statistical undertainty in the
d-spacing measured by diffraction peak can be mimicked. Finally, the
difference between the weighted-average stress and the one obtained by
the diffraction technique can quantify the propagated error to the stress
measured by the diffraction technique.


The below figure show an ensemble, in which the *perturbed* internal
strain is scattered and deviated from the 'fitting'.
![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/illu_2.png)


Based on these perturbed internal strain, the stress obtained by the
virtual experiment deviates from the 'weighted average' stress as below.
![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/illu_2f.png)


One can repeatedly conduct these virtual diffraction experiments to quantify
a statistically meaningful uncertainty, which is very difficult to obtain
in real experiments.


Contact information
===================

Youngung Jeong
youngung.jeong@gmail.com
younguj@clemson.edu

(2016 March -)
International Center for Automotive Research
Clemson University
(2014 March -)
Center for Automotive Lightweighting
National Institute of Standards and Technology


This software has been used for the following publications
==========================================================
------------------------------------------------------------------------
1. Uncertainty in flow stress measurements using X-ray diffraction for sheet metals subjected to large plastic deformations,
 **Y Jeong**, T Gnaeupel-Herold, M Iadicola, A Creuziger, Journal of Applied Crystallography (Accepted)
2. [Multiaxial constitutive behavior of an interstitial-free steel: measurements through X-ray and digital image correlation](http://dx.doi.org/10.1016/j.actamat.2016.04.013),
 **Y. Jeong**, T. Gnaeupel-Herold, M. Iadicola, A. Creuziger, Acta Materialia 112, 84-93 (2016)
3. [Evaluation of biaxial flow stress based on Elasto-Viscoplastic Self-Consistent analysis of X-ray Diffraction Measurements](http://dx.doi.org/10.1016/j.ijplas.2014.06.009),
 **Y Jeong**, T Gnaeupel-Herold, F Barlat, M Iadicola, A Creuziger, M-G Lee (2015) International Journal of Plasticity, 66, 103-118
