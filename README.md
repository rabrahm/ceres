# ceres
A set of routines for processing data of echelle spectrographs

Author: Rafael Brahm (rbrahm@astro.puc.cl)

# About the code
The principal goal of the CERES routines is the developement of fully automated pipelines for the reduction, extraction and analysis of echelle data. CERES currently counts with dedicated pipelines for eleven different instruments, which are included in the repository: APO3.5m/ARCES, CAHA2.2m/CAFE, Euler1.2m/Coralie, DuPont2.5m/Echelle, MPG2.2m/FEROS, ESO1.0m/FIDEOS, ESO3.6m/HARPS, KECK10m/HIRES, MAGELLAN6.5m/MIKE, MALLEGAN6.5m/PFS, PUC0.5m/PUCHEROS. These processing recipes can be used as reference for developing pipelines for new instruments. A detailed decription of the structure and performance of these pipelines can be found in [Brahm et al 2016](http://adsabs.harvard.edu/abs/2016arXiv160705792B). The standard output of the pipelines is a cube fits with the optimally extracted, wavelength calibrated and instrumetal drift corrected for each observed spectrum. Additionally, CERES includes routines for the computation of precise radial velocities and bisector spans via the cross-correlation method, and an automated algorithm to obtain a relatively precise guess of the atmospheric parameters of the observed star.



