#!/usr/bin/env python

# for a list of names, get orbital elements from MPC's orbital elements file.


# Imports
from astropy.table import *
import numpy as np
import os
import subprocess
from orbit_funcs import * # copy this over from /home/x25h971/orbits/files/nsc_sso_orbits

# Functions




if __name__=="__main__":

    basedir = "/home/x25h971/catalogs/"
    mpcorbdir = basedir+"mpc/mpcorb/"
    xmatch_file = basedir+"canfind/dr2/ring32/cfdr2_xmatch_mpc_mmts.fits.gz"
    cf_mmts=Table.read(xmatch_file)
    mpc_xmatches = cf_mmts[cf_mmts['match_flag']==True]
    mpc_xmatches['ppn'] = [i[:7] for i in mpc_xmatches['matched_mpc_mmt'].filled("".zfill(92))]
    nums = Column(np.unique(np.array(mpc_xmatches['ppn'])[mpc_xmatches['mpc_mode']=="num"])) # numbered objects
    perm_nums = np.array(nums[np.array([len(i.strip())==5 for i in nums])])
    unpacked_perm_nums = np.array([int(unpack_desig(i)[0]) for i in perm_nums])
    prov_nums = np.array(nums[np.array([len(i.strip())==7 for i in nums])])
    unpackedprov_nums = np.array([unpack_desig(i)[0] for i in prov_nums])

    mpcorb_files = subprocess.getoutput("ls "+mpcorbdir+"*").split("\n")
    first_mmts = []
    last_mmts = []
    for fl in mpcorb_files:
        with open(fl) as f:
            first_mmts.append(unpack_desig(f.readline().strip('\n')[0:7].strip()))
        with open(fl, 'rb') as f:
            try:  # catch OSError in case of a one line file 
                f.seek(-2, os.SEEK_END)
                while f.read(1) != b'\n':
                    f.seek(-2, os.SEEK_CUR)
            except OSError:
                f.seek(0)
            last_mmts.append(unpack_desig(f.readline().decode()[0:7].strip()))
    lostats = np.array([i[1] for i in first_mmts])
    histats = np.array([i[1] for i in last_mmts])
    loranges = np.array([i[0] for i in first_mmts])
    hiranges = np.array([i[0] for i in last_mmts])

    # Search files for permanent numbered minor planets
    perm_num_ranges = lostats=="num_minor_planet"
    for lo,hi,ind in zip(loranges[perm_num_ranges],hiranges[perm_num_ranges],range(0,len(loranges[perm_num_ranges]))):
        print("reading file ",mpcorb_files[ind])
        lines = readlines(mpcorb_files[ind])
        if np.char.isdigit(hi): filenums = perm_nums[(unpacked_perm_nums<=int(hi)) & (unpacked_perm_nums>=int(lo))]
        else: filenums = perm_nums[unpacked_perm_nums>=int(lo)]
        print("numbered minor planets between ",lo," and ",hi," = ",len(filenums))
        lbool = sum([grep(lines,"^"+o,index=True) for o in filenums],[])
        lbool = [int(l) for l in lbool]
        lines = Column(lines)[lbool]
        print(len(lines)," = lines")
