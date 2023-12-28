#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# mpc_makecat.py will read in one of the x<letter><letter> files
# which each contain 1 million lines from one of three of the
# Minor Planet Center's (MPC) observation files:
#     (a) Numbered Observations file "NumObs.txt.gz"
#     (b) Unnumbered Observations file "UnnObs.txt.gz"
#     (c) Isolated Tracklet File "itf.txt.gz"
#     (d) Minor Planet Center Orbit database "MPCORB.DAT.gz" - https://minorplanetcenter.net/iau/MPCORB.html
# and will write them out to the correct HEALPix (NSIDE=RING32)
# file, with labeling as to which MPC file the mmts came from.

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
import gzip
import healpy as hp
import numpy as np
import os
import sys
import time
#----------------------------------------------------------------------------
def makedir(dir):
    '''makes a directory with name "dir" if it does not exist
    Arguments:
    ----------
    dir (str)
            directory name
    Returns:
    --------
        None; directory "dir" is created if not already there
    '''
    if not os.path.exists(dir):
        os.mkdir(dir)
#----------------------------------------------------------------------------
def read_mpc80_line(line='     Hall2    C1999 06 05.03484 17 47 47.64 -25 27 24.3          16.6 R      706',mode="num"):
    ''' read a line in MPC's 80-column format.  For refere
        *** indicates a value that will be included in the output astropy row.
    Arguments:
    ln (str)
        the line, in the Minor Planet Center's 80-column format.
    mode (str)
        default "num" for Numbered Observations
            or  "unn" for Unnumbered Observations
            or  "itf" for Isolated Tracklets
    Returns:
    row (list)
        for astropy table, following columns:
        mjd,ra,dec,obscode,mag,filter,obs_code,mmt/object_id,line
    '''
    try:
        if mode=="itf" and line[13]=="s": return("satellite!")
        elif mode=="num" and line[14]=="s": return("satellite!")
        if mode=="itf": # because why would the files all be in the same format??
            lns = line.split("         ")
            filler1 = "          "
            filler2 = "      "
            ln = lns[0][:56]+filler1+lns[1][-16:-9]+filler2+lns[1][-4:-1]
        else: ln = line
        mmt_id = ln[:13]      # Object/Measurement ID ***
        #----------------
        # Get string info
        #----------------
        # Time
        yyyy=ln[15:19]        # Year [yyyy]
        mm=ln[20:22]          # Month [mm]
        dd=ln[23:31]          # Day [dd]
        # Coordinates
        ra_hh=ln[32:34]       # RA hour [hh]
        ra_mm=ln[35:37]       # RA minute [mm]
        ra_ss=ln[38:44]       # RA second [ss.s]
        if ra_ss=="   0  ": ra_ss="00.000"
        dec_sign=ln[44]       # Dec sign ["" for +, or "-""]
        dec_dd=ln[45:47]      # Dec degrees [dd]
        dec_mm=ln[48:50]      # Dec minutes [mm]
        dec_ss=ln[51:56]      # Dec seconds [ss.s]
        if dec_ss=="     ": dec_ss="00.00"
        # And the rest
        mag=ln[65:70]         # Magnitude [mm.m]
        filt=ln[70]           # Filter ***
        obscode=ln[77:81]     # Observatory code [COD] ***
        #ra_err=ln[81:86]      # RA error []
        #dec_err=ln[87:92]     # Dec error []
        #----------------------------
        # Transform strings to values
        #----------------------------
        # Time
        ddd = float(dd) - int(float(dd))
        h = int(ddd*24)
        m = int(((ddd*24)-h)*60)
        s = ((((ddd*24)-h)*60)-m)*60
        mjd = Time({'year': int(yyyy), 'month': int(mm), 'day': int(float(dd)),'hour': h, 'minute': m, 'second': s},scale='utc')
        mjd = mjd.mjd             # MDJ ***
        # Coordinates
        c = SkyCoord(ra_hh+'h'+ra_mm+'m'+ra_ss+'s', dec_sign+dec_dd+'d'+dec_mm+'m'+dec_ss+'s', frame='icrs')
        ra= c.ra.value            # RA [degrees] ***
        dec = c.dec.value           # Dec [degrees] ***
        # Magnitude & Observatory code
        if mag=="     ": mag_val = 0.0
        else: mag_val=float(mag)  # Magnitude ***
        if mode=="itf": obscode = lns[1][-6:-3]
        else: obscode = ln[77:81] # Observatory code ***
        #-----------------
        # Add row to table
        #-----------------
        row=[mjd,ra,dec,mag_val,filt,obscode,mmt_id,line]
        return(row)
    except Exception as e:
        print(line)
        print(e)
        return([-9999,-9999,-9999,-9999,"","","",line])


#----------------------------------------------------------------------------
# Main Code
#----------------------------------------------------------------------------
if __name__=="__main__":


    # --Set up inputs, out table, filenames, etc.--
    # ---------------------------------------------
    # output table (cat)
    if mode=="orb":
        mpc_mmts=Table(names=
    else:
        mpc_mmts=Table(names=("mjd","ra","dec","mag_augo","filter","obscode","name","line"),
                      dtype=["float64","float64","float64","float64","U1","U3","U13","U90"])
    # inputs
    mode = sys.argv[1]     # num for NumOBs.txt.gz, unn for UnnObs.txt.gz, itf for itf.txt.gz, orb for MPCORB.DAT.gz
    filename = sys.argv[2] # some split filename x<X><X> where the X's are two lowercase letters
    # repos & files
    basedir = "/home/x25h971/catalogs/"
    localdir = basedir+"files/"
    listdir = basedir+"mpc/"+str(mode)+"/"
    outdir = basedir+"mpc/ring32/"
    print("Mode ",mode)

    # --Read in MPC file --
    # ---------------------
    # read in maxlines, starting from tner
    w = open(listdir+filename)
    print("Opening MPC split mode ",mode,", file ",listdir+filename)
    for line in w:
        #row = read_mpc80_line(line.decode(),mode=mode)
        if mode=="orb": row = read_mpc80_line(line))
        else: row = read_mpc80_line(line,mode=mode)
        if row!="satellite!":
            if row[0]!=-9999:
                mpc_mmts.add_row(row)
    # format the MPC mmts nicely, make it nice!
    mpc_good = mpc_mmts['mjd']!=-9999
    print(len(mpc_mmts[mpc_good])," measurements of ",len(np.unique(mpc_mmts[mpc_good]['name']))," tracklets from ",filename)
    mpc_mmts['matched_cf_measid'] = Column(length=len(mpc_mmts),dtype="U19")
    mpc_mmts['mpc_mode'] = Column(np.repeat(mode,len(mpc_mmts)),dtype="U3")
    mpc_mmts['ring32'] = Column(np.repeat(-9999,len(mpc_mmts)),dtype='int')
    mpc_mjds = Time(mpc_mmts[mpc_good]['mjd'], format='mjd', scale='utc')
    mpc_coords = SkyCoord(ra=mpc_mmts[mpc_good]['ra'],dec=mpc_mmts[mpc_good]['dec'],frame="icrs",unit="degree",obstime=mpc_mjds)
    ring32s = hp.ang2pix(32,mpc_coords.ra.value,mpc_coords.dec.value,lonlat=True,nest=False)
    ring32s_unique = np.unique(ring32s)
    mpc_mmts['ring32'][mpc_good] = ring32s


    # Open or create the relevant <ring32>.fits files, add mmts
    print("writing mmts to ring32 files...")
    for r32 in ring32s_unique:
        #print("gathering mmts in RING32 = ",r32)
        mpc_r32s = mpc_mmts[mpc_good][mpc_mmts[mpc_good]['ring32']==r32]
        outfile = outdir+str(r32//1000)+"/"+str(r32)+"_"+str(mode)+"_"+str(filename)+".fits"
        r32dat = Table()
        try:
            r32dat = vstack([r32dat,mpc_r32s])
            makedir(outdir+str(r32//1000)+"/")
            r32dat.write(outfile,overwrite=True)
            #print(len(mpc_r32s)," MPC mmts written to ",outfile)
        except Exception as e: print("r32 = ",r32,e)

