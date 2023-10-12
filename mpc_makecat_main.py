#!/usr/bin/env python

# AUTHOR:  Katie Fasbender
#          katiefasbender@montana.edu

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from argparse import ArgumentParser
from astropy.table import Table,Column
from astropy.io import fits
from dlnpyutils import utils as dln, coords
import logging
import numpy as np
import os
import socket
import subprocess
import sys
import time

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
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


def write_jscript(job_name,partition,cmd,outdir):
    '''writes a SLURM job script to "job_name.sh"
    Arguments:
    ----------
    job_name (str)
            name of job, job script file
    partition (str)
            node/partition the job will run on
    cmd (str)
            python command to run exposure
    outdir (str)
            base directory
    Returns:
    --------
    job_file (str)
            job filename the job script is written to
    '''
    job_file = outdir+job_name+".sh"
    # The following code writes lines to the "job_name.sh" file.
    # Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
    # All other lines are read by the shell
    with open(job_file,'w') as fh:
        fh.writelines("#!/bin/bash\n")
        if partition=="priority": fh.writelines("#SBATCH --account=priority-davidnidever\n")       #specify the account to use
        fh.writelines("#SBATCH --job-name="+job_name+"\n")               # job name
        fh.writelines("#SBATCH --output="+outdir+job_name+".out\n")      # output file (%j = jobid)
        fh.writelines("#SBATCH --error="+outdir+job_name+".err\n")       # error file
        fh.writelines("#SBATCH --partition="+partition+"\n")             # queue partition to run the job in
        fh.writelines("#SBATCH --mem=9000\n")                            # memory, set --mem with care!!!!!
        fh.writelines("#SBATCH --time=1:00:00\n")                        # Maximum job run time
        fh.writelines("module load Anaconda3\n")                         # load anaconda, needed for running python on tempest!
        fh.writelines("source activate $HOME/condaenv/\n")               # activate the conda environment
        fh.writelines(cmd+"\n")                                          # write python command to read in the mpc split file and shit
        fh.writelines("conda deactivate")
    return job_file

def sacct_cmd(job_name,outputs,c=False,m=False):
    '''parses the output of a sacct command, returning specified information
    Arguments:
    ----------
    job_name (str)
            you know what this is
    outputs (str list)
            a list of information to get with the sacct command
            see sacct manual page for options
    c (bool)
            if job is known to be completed, default=False
    m (bool)
            if multiple tasks in the job, default=False
    Returns:
    --------
    outputs (list)
            a list of the sacct outputs specified
    '''
    if len(outputs)>1: spref = "sacct -n -P --delimiter=',' --format "
    else: spref = "sacct -n -X --format "
    scommand = (''.join([spref]+[i+"," for i in outputs]))[:-1]+" --name "+job_name
    if c and not m: jindex=1
    elif c and m: jindex=2
    else: jindex=0
    job_info = (subprocess.getoutput(scommand).split("\n")[jindex]).split(",")
    jinfo = [i.strip() for i in job_info]
    print("sacct cmd = ",scommand)
    print("sacct output = ",jinfo)
    if len(outputs)>1 and len(jinfo)>1: return(jinfo)
    elif len(outputs)>1 and len(jinfo)<=1: return(["" for ou in outputs])
    else: return(jinfo[0])



#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # Setup
    #------
    # Initiate input arguments
    parser = ArgumentParser(description='Crossmatch MPC measurements with CANFind tracklet measurements')
    parser.add_argument('--hplist',type=str,nargs=1,help='')
    parser.add_argument('--njobs',type=int,nargs=1,default=1,help='Number of jobs to maintain across all partitions')
    parser.add_argument('--partitions',type=str,nargs=1,help='Comma-delimited list of slurm partitions')
    parser.add_argument('--uncomp', action='store_true', help='xmatch missed ring32 HEALPix')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()

    # Input
    hplist = args.hplist[0]
    hplist = Table.read(localdir+hplist)
    njobs = int(args.njobs[0])
    partitions = args.partitions[0].split(',')  # the slurm partitions to submit jobs to
    npar = len(partitions)                     # number of slurm partitions
    nchanpar = njobs//npar                     # number of channels per partition
    uncomp = args.uncomp

    # Establish necessary directories
    basedir = "/home/x25h971/catalogs/"
    localdir = basedir+"files/"
    mpcdir = basedir+"mpc/ring32"
    cfdir = basedir+"canfind/dr2/mmts/ring32"
    rundir = basedir+"mpc/runfiles/"
    outfiledir = basedir+"mpc/outfiles/"
    makedir(rundir)
    makedir(outfiledir)

    # Prepare job info structure
    #---------------------------
    # Make a table for the job information
    jstr = Table()
    jstr['splitfile'] = Column(splitfiles)
    jstr['mode'] = Column(np.repeat(str(mode),nsf))
    jstr['cmd'] = Column(dtype="U1000",length=nsf)
    jstr['submitted'] = Column(np.repeat(False,nsf))
    jstr['job_done'] = Column(np.repeat(False,nsf))
    jstr['jobname'] = Column(dtype="U100",length=nsf)
    jstr['jobid'] = Column(dtype="U20",length=nsf)
    jstr['jobstatus'] = Column(dtype="U20",length=nsf)
    jstr['cputime'] = Column(dtype="U20",length=nsf)
    jstr['maxrss'] = Column(dtype="U20",length=nsf)
    jstr['maxvmsize'] = Column(dtype="U20",length=nsf)
    # assign partitions
    chans = np.reshape([[i+"_"+str(pchan) for i in partitions] for pchan in range(0,nchanpar)],njobs)
    jstr['partition'] = [chans[(i-njobs*(i//njobs))] for i in range(0,nsf)]
    # save job info table
    runfile = rundir+'mpc_makecat_'+str(t0)+'.fits'
    Table(jstr).write(runfile,overwrite=True)

    # Start submitting jobs
    #----------------------
    print(len(jstr)," split MPC files to read from ",mode)
    jcount = 0 # increments when a job is submitted
    eflag = 0 # = 1 when all jobs have been submitted
    while eflag==0:
        for ch in chans:
            print("Starting checks on channel ",ch)
            # Get indices of last & next jobs to be submitted for this channel
            partition_ind = set(np.where(jstr['partition']==ch)[0])  # indices for this partition
            unsubmitted_ind = set(np.where(jstr['submitted']==0)[0])   # indices for unsubmitted exposures (to get tsub for this s$
            submitted_ind = set(np.where(jstr['submitted']==1)[0])     # indices for submitted exposures (to get lsub for last sub$
            sub = list(submitted_ind & partition_ind)
            unsub = list(unsubmitted_ind & partition_ind)
            if len(unsub)>0: submit_flag = 1 # 1 for submit a job to this channel, 0 for no!
            else: submit_flag = 0
            # --Check status of last job submitted, if applicable--
            if len(sub)>0: # if there has already been a submitted job,
                # get index & status of last job submitted
                print("Checking status of last job submitted")
                lsub = np.sort(sub)[-1]
                last_jname = jstr[lsub]['jobname']
                last_jstat,last_jid = sacct_cmd(last_jname,["state","jobid"],c=False,m=False)
                jstr['jobstatus'][lsub] = last_jstat
                jstr['jobid'][lsub] = last_jid
                print("status of job ",last_jname,last_jid," = ",last_jstat)
                if (last_jstat!="RUNNING" and last_jstat!="PENDING" and last_jstat!="REQUEUED"):
                    print("last job ",last_jname," is ",last_jstat,"!")
                    jstr['job_done'][lsub] = 1
                    if len(unsub)>0: submit_flag = 1
                    # --If last job was completed, get some info about it--
                    if last_jstat=="COMPLETED":
                        last_jinfo = sacct_cmd(last_jname,["cputimeraw","maxrss","maxvmsize"],c=True)
                        jstr['cputime'][lsub] = last_jinfo[0]
                        jstr['maxrss'][lsub] = last_jinfo[1]
                        jstr['maxvmsize'][lsub] = last_jinfo[2]
                # ---If last job is still running, skip to next channel--
                else:
                    print("job ",last_jname,last_jid," is ",last_jstat,", move on to next channel in 10 sec")
                    Table(jstr).write(runfile,overwrite=True)
                    submit_flag = 0
                    time.sleep(10)

            if submit_flag==1:
                # --Submit a new job--
                print("--Submitting new job--")
                jsub = np.sort(unsub)[0]
                # info of next job to submit
                sfname = jstr['splitfile'][jsub].split("/")[-1]
                cmd = "python "+localdir+"mpc_makecat.py "+str(mode)+" "+str(sfname)
                prt = jstr['partition'][jsub].split("_")[0]
                jstr['cmd'][jsub] = cmd
                jname = 'mpc_makecat_'+str(t0)+'_'+str(jcount)
                jstr['jobname'][jsub] = jname
                # write job script to file
                jfile=write_jscript(jname,prt,cmd,outfiledir)
                # submit job to slurm queue
                os.system("sbatch "+jfile)
                jstr['submitted'][jsub] = True
                jcount+=1
                # save job structure, sleep before checking/submitting again
                Table(jstr).write(runfile,overwrite=True)
                print("Job "+jname+" submitted to ",str(prt)," partition; sleeping for 10 sec")
                time.sleep(10)

            if len(jstr[jstr['job_done']==1])==nsf:
                Table(jstr).write(runfile,overwrite=True)
                eflag=1


