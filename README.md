# ceres
A set of routines for processing data of echelle spectrographs

Author: Rafael Brahm (rbrahm@astro.puc.cl)

# About the code
The principal goal of the CERES routines is the developement of fully automated pipelines for the reduction, extraction and analysis of echelle data. CERES currently counts with dedicated pipelines for eleven different instruments, which are included in the repository: APO3.5m/ARCES, CAHA2.2m/CAFE, Euler1.2m/Coralie, DuPont2.5m/Echelle, MPG2.2m/FEROS, ESO1.0m/FIDEOS, ESO3.6m/HARPS, KECK10m/HIRES, MAGELLAN6.5m/MIKE, MALLEGAN6.5m/PFS, PUC0.5m/PUCHEROS. These processing recipes can be used as reference for developing pipelines for new instruments. A detailed decription of the structure and performance of these pipelines can be found in [Brahm et al 2016](http://adsabs.harvard.edu/abs/2016arXiv160705792B). The standard output of the pipelines is a cube fits with the optimally extracted, wavelength calibrated and instrumetal drift corrected spectrum for each of the science images. Additionally, CERES includes routines for the computation of precise radial velocities and bisector spans via the cross-correlation method, and an automated algorithm to obtain a relatively precise guess of the atmospheric parameters of the observed star.

# Usage
In order to run one of the pipelines, you have to be in the CERES directory where the pipeline lives and run the correspondig python pipe code followed by the path to the raw data. For example in the case of the FEROS pipeline you have to enter:

    $ cd ceres/feros
    $ python ferospipe.py /path/to/the/raw/data/to/be/processed

Additionally you can add some options to the command in order to modify some of the processing steps. The allowed options are:

    -avoid_plot     if activated, the code will not generate the pdf file with the plot of the computed CCF.
    -dirout         path to the directory where the reductions will be placed. The default path will be a
                    new directory with the same name that the input directory but followed by a '_red' suffix.
    -do_class       this option will allow the code to perform the estimation of the atmospheric parameters.
    -just_extract   if activated, the code will not compute the CCF and atmospheric parameters.
    -npools         number of CPU cores to be used by the code.
    -o2do           if you want to process just one particular science object you have to enter this option
                    followed by the name of the object.
    -reffile        name of the auxiliary file that is described below. The default is a '/reffile.txt' file
                    located inside the directory with the raw data.
    -ofind          only in the case of slit spectrographs. Name of the image that will be used to identify
                    and trace the echelle orders.
    
For example, if you want your output directory to bee foo/, and you don't want to generate the CCF plots, and you want to perform the spectral classification, and you want to use 10 CPU cores, and you want to process only the data of the target called HATS-17, then you have to enter:

    $ python ferospipe.py /path/to/the/raw/data/to/be/processed -dirout foo/ -avoid_plot -do_class -npools 10 -o2do HATS-17
    
The auxiliary file mentioned above (reffile) corresponds to a plain text file that can be used to give specifications for particular targets that can improve the precision in the radial velocity computation. The file should contain 8 comma-sparated colummns. The colummns should contain the following information:

    1-  name of the target as specified in the image header.
    2-  right ascension of the target (J2000) with format hh:mm:ss.
    3-  declination of the target (J2000) with format dd:mm:ss.
    4-  proper motion in RA [mas/yr].
    5-  proper motion in DEC [mas/yr].
    6-  integer (0 or 1). If 0, the code uses the cordinates given in the image header.
        If 1, the code uses the coordinates given in this file.
    7-  mask that will be used to compute the CCF. Allowed entries are G2, K5 and M2.
    8-  velocity width in km/s that is used to broaden the lines of the binary mask.
        I should be similar to the standard deviation of the gaussian that is fitted to the CCF. 
        
Here is an example of a reffile:

    HD157347,17:22:51.28809,-02:23:17.4297,49.39,-107.16,1,G2,4.0
    HD32147,05:00:48.99977,-05:45:13.2303,550.12,-1109.23,1,G2,4.0
    HD72673,08:32:51.49609,-31:30:03.0717,-1113.37,761.57,1,K5,4.0
    
If the pipeline doesn't find any reffile, it uses the coordinates found in the image header to compute the barycentric correction and uses the G2 mask to compute the CCF.

Additionally, there are other two auxiliary files that can be placed in the same directory of the raw images for improving the reductions. A file called 'bad_files.txt' can be used to list all the images of the raw directory that should be ignored by the pipeline. Each line of the bad_files file must have the name of the image without the complete path. Another file named 'moon_corr.txt' can be used to specify the images for which the CCF computation will include a double gaussian fit in order to correct the radial velocity by scattered moonlight contamination.

# Outputs
While the directory specified in the dirout directory will contain several intermediate reduction files, the final results will be placed in the directory 'proc/'. This directory should contain three types of files:

1) .fits files with the extracted and wavelength calibrated spectra.

2) .pdf files that show the CCF plots.

3) a text file (results.txt) that contains a summary of the results of the reduction including the radial velocity measurements and the atmospheric parameters.

The .fits files are data cubes with dimensions (10,nords,npix), where nords are the total number of processed echelle orders, and npix in the number of pixels in the dispersion direction. The ten entrances in the first dimension correspond to:

    0-  Wavelength
    1-  Extracted Flux
    2-  Measurment of the error in the extracted flux 1./sqrt(Var)
    3-  Blaze corrected Flux
    4-  Measurment of the error in the blaze corrected flux
    5-  Continuum normalized flux
    6-  Measurment of the error in the continuum normalized flux
    7-  Estimated continuum
    8-  Signal to noise ratio
    9-  Continumm normalized flux multiplied by the derivative of the wavelength with respect to the pixels
    10- Corresponding error of the 9th entrance

The colummns in the results.txt file have the following meaning:

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
    
    



