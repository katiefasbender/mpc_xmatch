#!/usr/bin/env python


# for a list of names, get orbital elements from MPC's orbital elements file.

# Imports


# Functions

def read_mpcorb(line):

    ln = line
    hdr = ["Des'n","H","G","Epoch","M","Peri.","Node","Incl.","e","n","a"," ","Reference","#Obs","#Opp","rms","Pert1","Pert2","Computer"]
    inds_low = [0, 8,14,20,26,37,48,59,70,80, 92,105,107,117,123,137,142,146,150]
    inds_hi =  [7,13,19,25,35,46,57,68,79,91,103,107,116,122,126,141,145,149,160]

    row = []
    for i,j in zip(inds_low,inds_hi):
        row.append(ln[i:j])#.strip())
        #print(ln[i:j])#.strip())
    floats = [1,2,4,5,6,7,8,9,10,15]
    for f in floats:
        row[f] = float(row[f])
    ints = [13,14]
    for ii in ints:
        row[ii] = int(row[ii])
    return(row)


#if __name__=="__main__":



    #names =
