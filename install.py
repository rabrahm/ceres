######################## Install Script v.1.0. #########################
#                                                                      #
# This script installs CERES in your computer                         #
#								       #
# UPDATES:							       #
#            							       #	
#          02-09-2016: Finished code. Working like a charm!            #
#                                                                      #  
########################################################################

p_name = "CERES"


import glob
import sys
import os
import shutil
import subprocess
sys.path.append("utils/SSEphem")
import update_ssephem

def CheckLibraries():
    print "     ----------------------------------------------------------"
    try:
      import numpy
    except ImportError:
      print "     ----------------------------------------------------------"
      print '     ERROR: '+p_name+' will not be installed in your system because' 
      print '             numpy is not installed in your system.'
      print '             To install it, go to: http://www.numpy.org/\n\n'
      sys.exit(1)
    print "     > Numpy is ok!"
    try:
      import scipy
    except ImportError:
      print "     ----------------------------------------------------------"
      print '     ERROR: '+p_name+' will not be installed in your system because'
      print '             scipy is not installed in your system.'
      print '             To install it, go to: http://www.scipy.org/\n\n'
      sys.exit(1)
    print "     > Scipy is ok!"
    try:
      import pyfits
    except ImportError:
      print "     ----------------------------------------------------------"
      print '     ERROR: '+p_name+' will not be installed in your system because'
      print '            pyfits is not installed in your system.'
      print '            To install it, go to: http://www.stsci.edu/institute/software_hardware/pyfits \n\n'
      sys.exit(1)
    print "     > Pyfits is ok!"
def getDirs(foldername):
    return os.walk(foldername).next()[1]

def spaced(input,space):
    fixed = False
    i = 0
    input = space+input
    while(not fixed):
        if(input[i:i+1] == '\n'):
           input = input[0:i+1]+space+input[i+1:]
           i = i + len(space)
        i = i + 1
        if(i == len(input)-1):
          fixed = True
    return input

def Build(directory):
    # We obtain al files and folders of the current directory...

    files_and_folders = glob.glob(directory+'/*')
    CFileFound = False
    SetupFileFound = False
    # ...and we check each folder or file:
    for cf in files_and_folders:
        # We search for files named Proceso_f2py, or for the setup.py - file.c
        # combination. If present, we build the process.
        pf2py = 'Proceso_f2py'
        stp = 'setup.py'
        if( cf[-len(pf2py):] == pf2py ):
           print "     > Fortran code found in directory "+directory+". Building..."
           cwd = os.getcwd()
           os.chdir(directory)
           subprocess.Popen('chmod u+wrx Proceso_f2py',shell = True).wait()
           subprocess.Popen('chmod g+wrx Proceso_f2py',shell = True).wait()
           subprocess.Popen('chmod o+wrx Proceso_f2py',shell = True).wait()
           p = subprocess.Popen('./Proceso_f2py',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
           p.wait()
           if(p.returncode != 0 and p.returncode != None):
             print "     ----------------------------------------------------------" 
             print "     > ERROR: "+p_name+" couldn't be installed."
             print "     > Problem building code in "+directory+". The error was:\n"
             out, err = p.communicate()
             print err
             print "     > If you can't solve the problem, please communicate"
             print "     > with the "+p_name+" team for help.\n\n"
             os.chdir(cwd)
             sys.exit()
           os.chdir(cwd)
           print "     >...done!"
        elif( cf[-len(stp):] == stp ):
           SetupFileFound = True
        elif( cf[-2:] == '.c' ):
           CFileFound = True
        if( SetupFileFound and CFileFound ):
           print "     > C code found in directory "+directory+". Building..."
           cwd = os.getcwd()
           os.chdir(directory)
           p = subprocess.Popen('{python} setup.py build'.format(python=sys.executable),stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
           p.wait()
           if(p.returncode != 0 and p.returncode != None):
             print "     ----------------------------------------------------------"
             print "     > ERROR: "+p_name+" couldn't be installed."
             print "     > Problem building code in "+directory+". The error was:\n"
             out, err = p.communicate()
             print spaced(err,"\t \t")
             print "     > If you can't solve the problem, please communicate"
             print "     > with the "+p_name+" team for help.\n \n"
             os.chdir(cwd)
             sys.exit()
           libfolder = getDirs('build/.')
           for name in libfolder:
               if(name[0:3]=='lib'):
                 filename = glob.glob('build/'+name+'/*')
                 shutil.copy2(filename[0],'.')
           shutil.rmtree('build')
           os.chdir(cwd)
           print '     >...done!'
        # If the current file or folder is a directory, we apply the same 
        # code to it:
        elif( os.path.isdir(cf) ):
           Build(cf)
print " \n\n                       "+p_name+" Installer v.1.0.  \n\n"
print " The "+p_name+" team is composed of: \n"
print "     - Rafael Brahm (PUC, rabrahm@astro.puc.cl)."
print "     - Andres Jordan (PUC, ajordan@astro.puc.cl)."
print "     - Nestor Espinoza (PUC, nespino@astro.puc.cl). \n"
print " DISCLAIMER: If you use this pipeline or part of it, please  "
print " acknowledge us and our current institutions. If you find any bugs,"
print " please contact us. \n"
print " 1.- Preparing to install. Checking if your system has the libraries"
print "     needed to compile the codes...\n"
CheckLibraries()
print "     ---------------------------------------------------------- \n"
print "     Done! All libraries checked. \n"
import numpy
print " 2.- Building processes...\n"
print "     ----------------------------------------------------------"
# First, we get all the directories in the current folder:
dirs = getDirs('utils/')
ndirs = []
for dire in dirs:
	ndirs.append('utils/'+dire)
dirs = ndirs
for directory in dirs:
    # To each directory, we apply the build function:
    if(directory == 'utils/SSEphem'):
       print "     > Installing SSEphem...\n \n"
       os.chdir('utils/SSEphem/SOFA')
       os.system('mkdir lib')
       os.system('mkdir include')
       p = subprocess.Popen('make',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
       p.wait()
       if(p.returncode != 0 and p.returncode != None):
          print '           Error in SSEphem installation! The error was:'
          out, err = p.communicate()
          print spaced(err,"\t \t")
          sys.exit()
       p = subprocess.Popen('make test',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
       p.wait()
       if(p.returncode != 0 and p.returncode != None):
          print '           Error in SSEphem installation! The error was:'
          out, err = p.communicate()
          print spaced(err,"\t \t")
          sys.exit()
       os.chdir('../')
       p = subprocess.Popen('make clean',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
       p.wait()
       if(p.returncode != 0 and p.returncode != None):
          print '           Error in SSEphem installation! The error was:'
          out, err = p.communicate()
          print spaced(err,"\t \t")
          sys.exit()
       p = subprocess.Popen('make',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
       p.wait()
       if(p.returncode != 0 and p.returncode != None):
          print '           Error in SSEphem installation! The error was:'
          out, err = p.communicate()
          print spaced(err,"\t \t")
          sys.exit()
       else:
          update_ssephem.SSEphemDownload()
          update_ssephem.LeapSecUpdate()
          update_ssephem.IersUpdate()

       os.chdir('../../')
    else:
       Build(directory)
print "     ---------------------------------------------------------- \n \n"
print "     Installation of "+p_name+" finished without problems! \n \n"
print "     Please read the README file in order to learn how to "
print "     use the routines. \n \n"

