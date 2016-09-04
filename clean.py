######################## Clean Script v.1.0. ###########################
#                                                                      #
# This script cleans all the subdirectories of the current directory   #
# os .so, .pyc, .something~ and build folders. USE WITH CAUTION;       #
# it was designed to be used in the ceres folder only.                #
#                                                                      #
#								       #
# UPDATES:							       #
#            							       #	
#          02-09-2016: Finished code. Working like a charm!            #
#                                                                      #  
########################################################################

import glob
import os
import shutil
import subprocess
def getDirs(foldername):
    return os.walk(foldername).next()[1]

def SeeknDestroy(directory):
    if(directory == 'utils/JPLEphemx'):
       os.chdir('utils/JPLEphemx')
       p = subprocess.Popen('make clean',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
       p.wait()
       os.chdir('SOFA')
       p = subprocess.Popen('make clean',stdout = subprocess.PIPE, stderr = subprocess.PIPE,shell = True)
       p.wait()
       os.chdir('lib')
       if(os.access('libsofa_c.a',os.F_OK)):
          os.remove('libsofa_c.a')
       os.chdir('../../')
    # We obtain al files and folders of the current directory...
    files_and_folders = glob.glob(directory+'/*')
    # ...and we check each folder or file:
    for cf in files_and_folders:
        # We search for .so, .pyc, .something~ files or build folders.
        # If we find them, we erase them:
        if( cf[-3:] == '.so' or cf[-4:] == '.pyc' or cf[-1:] == '~'):
           os.remove(cf)
        elif( cf[-5:] == 'build' ):
           shutil.rmtree(cf)
        # If the current file or folder is a directory, we apply the same 
        # SeeknDestroy code to it:
        elif( os.path.isdir(cf) ):
           SeeknDestroy(cf)
# First, we get all the directories in the current folder:

dirs2 = getDirs('utils')

for directory in dirs2:
    print directory
    # To each directory, we apply the SeeknDestroy function:
    SeeknDestroy('utils/'+directory)
