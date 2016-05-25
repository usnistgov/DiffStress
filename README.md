DiffStress
==========

DiffStress is a Python-based diffraction analysis tool that
analyzes the *tensorial* stress state present on a polycrystal sample.
The stress analysis utilizes a set of lattice-spacings (interatomic spacings)
obtained from a polycrystalline specimen using a diffractometer (either X-ray or neutron)

This python package allows one to analyze both experimental diffraction lattice strains and
model-emulated lattice strains, e.g., those data from EVPSC.

For the purpose of experimental d-spacing analysis, it accepts
the format of Proto data. Current software can generate
various graphics for the purpose of data visualization by using
matplotlib and some in-house scripts written for matplotlib, which
is also available under current GitHub account (find "mpl-lib").

In order to estimate the experimental/model-predicted stress,
current software requires the use of diffraction elastic constants written
in the format of 'sff', which is used in the 'PF' software
managed by Thomas Gnaeupel-Herold in NCNR, NIST. One may find the
template of the 'sff' file from the current software package.

Also, given set of diffraction data (say, from EVPSC code), one can perturb the lattice strains
in order to quantify the uncertainty in the obtained stress by conducting Monte Carlo virtual experiments.


Application
===========

DiffStress comes with a feature that allows Monte Carlo experiments to quantify uncertainties present in flow-stress measurement using
in-situ diffraction experiments to obtain flow stress curve. Below image shows the resulting
internal elastic strain (that is corresponding to d-spacing vs. sin2psi curve in a real experiment) and the fitting
to the elastic law (like the conventional sin2psi method widely used in the community of residual stress measurement).

![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/illu_1.png)


Of course, this method, within the model-predicted values of virtual diffraction experiments are *self-consistent*,
meaning that the weighted average of stress (the macro stress) is equivalent to the stress esimated by this virtual
diffraction experiments. Below figure shows that the two flow stress curves (one with weighted average stress another obtained
by this virtual diffraction experiment) are equivalent.


![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/illu_1f.png)


In the real experiments, however, this self-consistency may be challenged by various factors. For example, due to the counting 
statistical nature, your diffraction peak contains a certain degree of uncertainty when determining the correct d-spacing.
Also, your peak height is influenced by the presence of crystallographic (and its development w.r.t plastic deformation) thus further affecting the uncertainty in the peak location.
In DiffStress, one can mimic various types of uncertainties existing in real experiments and propagate these uncertainties to the
final stress estimation using the aforementioned virtual diffraction experiments. Ultimately, one could compare the two types of
stresses: the weighted-averaged stress and the stress from this virtual diffraction experiments using the 'perturbed' internal strains.


The below figure shows that the internal-spacing is scattered and deviated from the 'fitting' naturally.
![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/illu_2.png)


Based on these 'perturbed' internal strain, the stress by diffraction experiment deviates from the 'weighted average' stress as below.
![image of IF stee diffraction](https://github.com/usnistgov/DiffStress/blob/dev/images/illu_2f.png)


One can repeatedly conduct these virtual diffraction experiments to obtain a statistically meaning data - the Monte Carlo experiments.



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
----------------------------------------------------------
1. Uncertainty in flow stress measurements using X-ray diffraction for sheet metals subjected to large plastic deformations, Y Jeong, T Gnaeupel-Herold, M Iadicola, A Creuziger (Submitted to Journal of Applied Crystallography)
2. Multiaxial constitutive behavior of an interstitial-free steel: measurements through X-ray and digital image correlation, Y Jeong, M Iadicola, T Gnaeupel-Herold, A Creuziger (2016) Acta Materialia 112, 84-93
3. Evaluation of biaxial flow stress based on Elasto-Viscoplastic Self-Consistent analysis of X-ray Diffraction Measurements, Y Jeong, T Gnaeupel-Herold, F Barlat, M Iadicola, A Creuziger, M-G Lee (2015) International Journal of Plasticity, 66, 103-118
