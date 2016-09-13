# ceres
A set of pipelines and routines for echelle spectrographs

Authors: Rafael Brahm, Andrés Jordán, Néstor Espinoza

# About the code
The principal goal of CERES is the developement of fully automated pipelines for the reduction, extraction and analysis of echelle spectrograph data. CERES currently counts with dedicated pipelines for eleven different instruments, which are included in the repository: APO3.5m/ARCES, CAHA2.2m/CAFE, Euler1.2m/Coralie, DuPont2.5m/Echelle, MPG2.2m/FEROS, ESO1.0m/FIDEOS, ESO3.6m/HARPS, KECK10m/HIRES, MAGELLAN6.5m/MIKE, MALLEGAN6.5m/PFS, PUC0.5m/PUCHEROS. These processing recipes can be used as a reference for developing pipelines for new instruments. A detailed decription of the structure and performance of these pipelines can be found in [Brahm et al 2016](http://arxiv.org/abs/1609.02279). The standard output of the pipelines is a fits cube with the optimally extracted, wavelength calibrated and instrumetal drift corrected spectrum for each of the science images. Additionally, CERES includes routines for the computation of precise radial velocities and bisector spans via the cross-correlation method, and an automated algorithm to obtain an estimate of the atmospheric parameters of the observed star.

# Usage
In order to run one of the pipelines, you have to be in the CERES directory where the pipeline for the instrument of interest lies and run the corresponding:


```
python INSTRUMENTpipe.py
```

code followed by the path to the raw data. For example, in the case of the FEROS pipeline you have to enter:

```
cd ceres/feros
python ferospipe.py /path/to/the/raw/data/to/be/processed
```

Additionally you can specify some options in order to modify some of the processing steps. The available options are:

    -avoid_plot     if activated, the code will not generate a pdf file with the plot of the computed CCF.
    -dirout         path to the directory where the pipeline products will be placed. The default path will be a
                    new directory with the same name that the input directory but followed by a '_red' suffix.
    -do_class       this option will enable the estimation of atmospheric parameters.
    -just_extract   if activated, the code will not compute the CCF and atmospheric parameters.
    -npools         number of CPU cores to be used by the code.
    -o2do           if you want to process just one particular science object you have to enter this option
                    followed by the name of the object.
    -reffile        name of the auxiliary file that is described below. The default is './reffile.txt', a file
                    located in the directory where the raw data is.
    -ofind          only valid for slit spectrographs. Name of the image that will be used to identify
                    and trace the echelle orders.

For example, if you want your output directory to bee foo/, you don't want to generate the CCF plots, you want to perform the spectral classification, you want to use 10 CPU cores, and you want to process only the data of the target called HATS-17, then you would do:

```
python ferospipe.py /path/to/the/raw/data/to/be/processed -dirout foo/ -avoid_plot -do_class -npools 10 -o2do HATS-17
```

The auxiliary file mentioned above (reffile) corresponds to a plain text file that can be used to give data for particular targets in order to improve the precision in the radial velocity estimation. The file should contain 8 comma-sparated colummns containing the following information:

    1-  name of the target as specified in the image header.
    2-  right ascension of the target (J2000) with format hh:mm:ss.
    3-  declination of the target (J2000) with format dd:mm:ss.
    4-  proper motion in RA [mas/yr].
    5-  proper motion in DEC [mas/yr].
    6-  integer (0 or 1). If 0, the code uses the cordinates given in the image header.
        If 1, the code uses the coordinates given in this file.
    7-  mask that will be used to compute the CCF. Allowed entries are G2, K5 and M2.
    8-  velocity width in km/s that is used to broaden the lines of the binary mask.
        It should be similar to the standard deviation of the Gaussian that is fitted to the CCF.

Here is an example of a reffile:

    HD157347,17:22:51.28809,-02:23:17.4297,49.39,-107.16,1,G2,4.0
    HD32147,05:00:48.99977,-05:45:13.2303,550.12,-1109.23,1,G2,4.0
    HD72673,08:32:51.49609,-31:30:03.0717,-1113.37,761.57,1,K5,4.0

If the pipeline doesn't find any reffile, it uses the coordinates found in the image header to compute the barycentric correction and uses the G2 mask to compute the CCF.

Additionally, there are other two auxiliary files that can be placed in the same directory of the raw images for improving the reductions. A file called 'bad_files.txt' can be used to list all the images of the raw directory that should be ignored by the pipeline. Each line of the bad_files file must have the name of the image without the complete path. Another file named 'moon_corr.txt' can be used to specify the images for which the CCF computation will include a double gaussian fit in order to correct the radial velocity for the effects of scattered moonlight contamination.

# Output
While the directory specified in the dirout directory will contain several intermediate reduction files, the final results will be placed in the directory 'proc/'. This directory can contain up to three types of files:

   1. .fits files with the extracted and wavelength calibrated spectra.
   1. .pdf files that show the CCF plots.
   1. a text file (results.txt) that contains a summary of the results of the reduction including the radial velocity measurements and the atmospheric parameters.

The .fits files are data cubes with dimensions (10,nords,npix), where nords are the total number of processed echelle orders, and npix in the number of pixels in the dispersion direction. The ten entries in the first dimension correspond to:

    0-  Wavelength
    1-  Extracted Flux
    2-  Measurement of the error in the extracted flux [1./sqrt(Var)]
    3-  Blaze corrected Flux
    4-  Measurement of the error in the blaze corrected flux
    5-  Continuum normalized flux
    6-  Measurement of the error in the continuum normalized flux
    7-  Estimated continuum
    8-  Signal-to-noise ratio
    9-  Continumm normalized flux multiplied by the derivative of the wavelength with respect to the pixels
    10- Corresponding error of the 9th entrance

The colummns in the results.txt file are the following:

    0-  Object name
    1-  MBJD
    2-  RV
    3-  error in RV
    4-  Bisector span
    5-  error in bisector span
    6-  instrument
    7-  pipeline
    8-  resolving power
    9-  Efective Temperture
    10- log(g)
    11- [Fe/H]
    12- v*sin(i)
    13- value of the continuum normalized CCF at it lowest point
    14- standard deviation of the gaussian fitted to the CCF
    15- Exposure time
    16- Signal to noise ratio at ~5150 \AA
    17- path to the CCF plot file

# Dependencies
CERES has been succesfully tested with python2.6 and python2.7 on OS X and Linux systems.

Some python packages need to be installed:

   1. python-numpy
   1. python-scipy
   1. python-matplotlib
   1. python-pyfits
   1. python-pycurl
   1. python-ephem
   1. python-rpy2 (depends on R)

All of the above can be installed using pip, e.g.:

```pip install numpy, scipy, matplotlib, pyfits, pycurl, ephem, rpy2```

Additionally:

   1. SWIG
   1. gcc, g++, gfortran
   1. gsl

are required. All of the above can be install using homebrew on OSX, e.g.:

```brew install swig, gcc, gsl```

# Installation
The python code does not need any previous installation. However, CERES uses some C, C++, and fortran routines that must be compiled in order to be called from python by the pipelines. Before installation it is necessary to check and modify some variables of the installation files of these codes.

First, it must be checked if the default f2py version is compatible with the default python version. If that is not the case, the file named "Proceso_f2py" in ceres/utils/CCF must be modified in order to point to the compatible f2py binary file.

CERES uses the [SSephem](http://www.cv.nrao.edu/~rfisher/Python/py_solar_system.html) package coupled to the [SOFA](http://www.iausofa.org/) C functions for computing the barycentric velocity corrections. A version of this package with minor modifications is included in the CERES repository. In order to install this package, it has to be checked that the ceres/utils/SSephem/Makefile points to the correct Python, SWIG, and numpy paths. Edit the lines:

```
PY_PREFIX := /usr/local
SWIG = /usr/local/bin/swig
```

SSephem requires also three separete files for running: the binary solar system ephemeris, a leap-second table, and a UT1 - UTC offset table. These three files are downloaded during the CERES installation. However, the latter two files should be periodically updated, which can be easily done by reinstalling CERES.

After checking the C, C++, and fortran installation files, CERES can be istalled with the following command:

```python install.py```

For uninstalling CERES you can enter:

```python clean.py```

