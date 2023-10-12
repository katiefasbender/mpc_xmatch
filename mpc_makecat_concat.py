#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# mpc_makecat_concat.py will create a finalized catalog of observations of minor planets
# from three Minor Planet Center (MPC) files: Numbered & Unnumbered Observations, and the
# Isolated Tracklet file.  The catalog is ordered by NSIDE=ring32 HEALPix.

# For the input directory, the script reads in all ring32 HEALPix fits files with MPC
# observations and combines them into one individual file for each ring32 pix.

# Imports
from astropy.table import Table,Column,vstack
import numpy as np
import subprocess
import sys
import os

# Main Code
if __name__=="__main__":

    # Inputs and Setup
    basedir = "/home/x25h971/catalogs/"
    rootdir = basedir+"mpc/ring32/"          # MPC mmts catalog, ordered by ring32 HEALPix
    dir = sys.argv[1]                        # ring32//1000 directory in rootdir
    unavails = basedir+"mpc/unavails.txt"    # file with a list of unreadable <ring32>_<mode>_<splitfile>.fits files

    # Get a list of files in the directory
    dir_flist = os.listdir(rootdir+dir)
    print("Concatenating "+str(len(dir_flist))+" files in ",rootdir+dir," directory")
    dir_r32s = [i.split("_")[0] for i in dir_flist]
    dir_flist = np.array(dir_flist) # filenames in directory
    dir_r32s = np.array(dir_r32s)   # ring32 HEALPix from filenames

    # For each ring32 pix, concatenate mmts from all files
    for r32 in np.unique(dir_r32s):
        r32_cat = Table() # initiate a table for r32 mmts
        fls = dir_r32s==r32
        r32_files = dir_flist[fls]
        # for each file, read in mmts and add to r32 table
        for fl in r32_files:
            if fl[0].isdigit():
                try:
                    fl_tab = Table.read(rootdir+dir+"/"+fl)
                    r32_cat = vstack([r32_cat,fl_tab])
                except: # if file is unreadable, save the filename
                    with open(unavails,'a') as fil:
                       fil.writelines(fl+"\n")
                       fil.close()
        # save r32 mmts to a nice new file
        r32_fname = rootdir+dir+"/"+str(r32)+"_mpc.fits.gz"
        #print(r32_fname)
        r32_cat.write(r32_fname,overwrite=True)
        print("ring32 "+str(r32)+" mmts written to "+r32_fname)
