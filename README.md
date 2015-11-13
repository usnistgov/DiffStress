DiffStress
==========

Diffraction analysis tool that estimates stress present in polycrystalline specimen.

This python package allows one to analyze both experimental and
model-predicted diffraction strains in order to estimate a tensorial
stress values estimated from a polycrystalline sample.

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
