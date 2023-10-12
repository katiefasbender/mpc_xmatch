#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# mpc_xmatch.py will read in the measurements from a healpix (NSIDE=RING32)
# included in the NSC CANFind tracklet mmt catalog, and in the Minor Planet
# Center's (MPC) 3 observation files:
#     (a) Numbered Observations file "NumObs.txt.gz"
#     (b) Unnumbered Observations file "UnnObs.txt.gz"
#     (c) Isolated Tracklet File "itf.txt.gz"
# and will crossmatch (xmatch) these measurements, returning a catalog with
# the NSC CANFind tracklet mmts and their matched MPC mmts.

# The MPC's 80-column line format is this:
#---------------------------
#   Columns     Format   Use
#---------------------------
#   For Minor Planets:   (assume this for all tracklets, initially )
#    1 -  5       A5     Packed minor planet number
#    6 - 12       A7     Packed provisional designation, or a temporary designation
#    13           A1     Discovery asterisk
#   For Comets:
#    1 -  4       I4     Periodic comet number
#    5            A1     Letter indicating type of orbit
#    6 - 12       A7     Provisional or temporary designation
#    13           X      Not used, must be blank
#   For Natural Satellites:
#    1            A1     Planet identifier [Only if numbered]
#    2 -  4       I3     Satellite number  [Only if numbered]
#    5            A1     "S"
#    6 - 12       A7     Provisional or temporary designation
#    13           X      Not used, must be blank
#---------------------------
#   Columns     Format   Use
#---------------------------
#   For all: (columns 14-80)
#    14            A1    Note 1
#    15            A1    Note 2
#    16 - 32             Date of observation
#    33 - 44             Observed RA (J2000.0)
#    45 - 56             Observed Dec (J2000.0)
#    57 - 65       9X    Must be blank
#    66 - 71    F5.2,A1  Observed magnitude and band
#    72 - 77       X     Must be blank
#    78 - 80       A3    Observatory code


#----------------------------------------------------------------------------
# Imports & Functions
#----------------------------------------------------------------------------
from astropy.coordinates import SkyCoord
from astropy.table import Table,Column,vstack
from astropy.time import Time
from astropy import units as u
import healpy as hp
import numpy as np
import os
import subprocess
import sys
import time
from mpc_makecat import *
#----------------------------------------------------------------------------
#def read_mpc80_line(line='     Hall2    C1999 06 05.03484 17 47 47.64 -25 27 24.3          16.6 R      706',mode="num"):

#----------------------------------------------------------------------------
# Main Code
#----------------------------------------------------------------------------
if __name__=="__main__":


    # --Set up inputs, out table, filenames, etc.--
    # ---------------------------------------------
    # output table (cat)
    #mpc_mmts=Table(names=("mjd","ra","dec","mag_augo","filter","obscode","name","line"),
    #              dtype=["float64","float64","float64","float64","U1","U3","U13","U90"])
    # inputs
    ring32 = int(sys.argv[1]) # ring32 HEALPix
    # repos & files
    basedir = "/home/x25h971/catalogs/" # for MSU tempest
    localdir = basedir+"files/"
    mpcdir = basedir+"mpc/ring32/"      # where the MPC observations catalog is stored
    cfdir = basedir+"canfind/dr2/ring32/"   # where CANFind NSC mmts are stored

    # --Read in MPC file --
    # ---------------------
    mpc_mmts = Table()
    mpc_file = mpcdir+str(ring32//1000)+"/"+str(ring32)+"_*.fits"
    mpc_files = subprocess.getoutput("ls -ltr "+mpc_file).split("\n")
    print("pulling MPC mmts from ",len(mpc_files)," ring32=",str(ring32)," MPC files")
    for mpcfl in mpc_files:
        mpcf = mpcfl.split(" ")[-1]
        mpc_mmts = vstack([mpc_mmts,Table.read(mpcf)])
    mpc_mjds = Time(mpc_mmts['mjd'],format='mjd',scale='utc')
    mpc_coords = SkyCoord(ra=mpc_mmts['ra'],dec=mpc_mmts['dec'],frame="icrs",unit="degree",obstime=mpc_mjds)

    # --Read in CANFind tracklets--
    #------------------------------
    # use mpc coords to get ring32 values
    cf_file = cfdir+str(ring32//1000)+"/"+str(ring32)+".fits.gz"
    print("pulling CANFind tracklet measurements from ",cf_file)
    cf_mmts = Table.read(cf_file)
    cf_mmts['matched_mpc_mmt'] = Column(length=len(cf_mmts),dtype="U92")
    cf_mjds = Time(cf_mmts['mjd'],format="mjd",scale="utc")
    cf_coords = SkyCoord(ra=cf_mmts['ra'],dec=cf_mmts['dec'],frame="icrs",unit="degree",obstime=cf_mjds)

    # --Cross-match MPC mmts with CANFind tracklet mmts--
    #----------------------------------------------------
    # ----space match----
    idx_mpc,idx_cf,d2d,d3d = cf_coords.search_around_sky(mpc_coords, 0.00015*u.deg) # 1.5e-4 degrees = 0.54 arcsec
    print(len(np.unique(mpc_mmts[idx_mpc]['name']))," tracklets matched from MPC's files")
    print(len(np.unique(cf_mmts[idx_cf]['tracklet_id']))," tracklets matched from NSC CF tracklet cat")
    # ----time match----
    good_matches = (abs(cf_mmts[idx_cf]['mjd']-mpc_mmts[idx_mpc]['mjd'])<0.0208) # 0.0208 days = 29.95 minutes
    print(len(np.unique(mpc_mmts[idx_mpc[good_matches]]['name']))," date-validated tracklet matches from MPC file")
    print(len(np.unique(cf_mmts[idx_cf[good_matches]]['tracklet_id']))," date-validated tracklet matches from NSC")
    # write out the CANFind mmts with any matched MPC mmts as new column
    cf_mmts['matched_mpc_mmt'][idx_cf[good_matches]] = mpc_mmts['line'][idx_mpc[good_matches]]
    cf_mmts.write(cf_file,overwrite=True)
#    mpc_mmts['matched_cf_measid'][idx_mpc[good_matches]] = cf_mmts['measid'][idx_cf[good_matches]]
#    mpc_mmts.write(mpc_file,overwrite=True)
    print("CANFind mmts & matches written to ",cf_file)
#    print("MPC mmts & matches written to ",mpc_file)
