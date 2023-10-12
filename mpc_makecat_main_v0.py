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
        fh.writelines("#SBATCH --job-name="+job_name+"\n")        # job name
        fh.writelines("#SBATCH --output="+outdir+job_name+".out\n")      # output file (%j = jobid)
        fh.writelines("#SBATCH --error="+outdir+job_name+".err\n")       # error file
        fh.writelines("#SBATCH --partition="+partition+"\n")     # queue partition to run the job in
        fh.writelines("#SBATCH --ntasks=1\n")                    # for running in parallel
        fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
        fh.writelines("#SBATCH --cpus-per-task=10\n")
        fh.writelines("#SBATCH --mem=9000\n")                    # memory, set --mem with care!!!!! refer to hyalite quickstart guide
        fh.writelines("#SBATCH --time=6:00:00\n")               # Maximum job run time
        fh.writelines("module load Anaconda3\n")         # load anaconda, needed for running python on Hyalite!
        fh.writelines("source activate $HOME/condaenv/\n")
        fh.writelines(cmd+"\n")                                       # write python command to analyze exposure
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
    parser = ArgumentParser(description='Create a ring32-sorted catalog of MPC measurements')
    parser.add_argument('--mode',type=str,nargs=1,help='Which MPC file to read from (num for NumObs, unn for UnnObs, and itf)')
    args = parser.parse_args()

    # Start time, get hostname (should be tempest)
    t0 = time.time()

    # Input
    mode = args.mode[0]

    # Establish necessary directories
    basedir = "/home/x25h971/catalogs/"
    localdir = basedir+"files/"
    listdir = basedir+"mpc/"
    outdir = basedir+"mpc/ring32/"
    outfiledir = basedir+"mpc/outfiles/"
    checkfiledir = outfiledir+"checks/"
    for f in os.listdir(checkfiledir): os.remove(os.path.join(checkfiledir,f))
    makedir(outfiledir)
    killfile = listdir+"kill_"+str(mode)+".txt"
    if os.path.exists(killfile): os.remove(killfile)
    if mode=="num": listfile = listdir+"NumObs.txt.gz"
    elif mode=="unn": listfile = listdir+"UnnObs.txt.gz"
    elif mode=="itf": listfile = listdir+"itf.txt.gz"

    # Prepare job info structure
    #---------------------------
    jstr = Table(names=('tner','mode','done','cmd','submitted',
                        'jobname','jobid','jobstatus','outfile',
                        'cputime','maxrss','maxvmsize'),
                 dtype=('int','U3','bool','U1000','bool',
                        'U100','U20','U20','U1000',
                        'U20','U20','U20'))
    runfile = listdir+'runfiles/mpc_makecat_'+str(t0)+'.fits'
    Table(jstr).write(runfile,overwrite=True)

    # Start submitting jobs
    #----------------------
    maxlines = 100000 # in mpc_makecat.py, this is the number of lines to read from the MPC file
    jcount = 0 # increments when a job is submitted
    eflag = 0 # = 1 when killfile exists (all lines of MPC file have been read)
    while eflag==0:
        bindex = jcount*maxlines
        # --Check status of last job submitted, if applicable--
        if jcount!=0: # if there has already been a submitted job,
            # get index & status of last job submitted
            print("Checking status of last job submitted")
            lsub = -1
            last_jname = jstr[lsub]['jobname']
            last_jstat,last_jid = sacct_cmd(last_jname,["state","jobid"],c=False,m=False)
            jstr[lsub]['jobstatus'] = last_jstat
            jstr[lsub]['jobid'] = last_jid
            print("status of job ",last_jname,last_jid," = ",last_jstat)
            # ---If last job is still running: wait!--
            if (last_jstat=="RUNNING" or last_jstat=="PENDING" or last_jstat=="REQUEUED"):
                while (last_jstat=="RUNNING" or last_jstat=="PENDING" or last_jstat=="REQUEUED"):
                    last_jstat = sacct_cmd(last_jname,["state"])
                    print("job ",last_jname,last_jid," is ",last_jstat,", sleepin for 1 min")
                    time.sleep(60)
            print("last job ",last_jname," is ",last_jstat,"!")
            # --If last job was completed, get some info about it--
            if last_jstat=="COMPLETED":
                last_jinfo = sacct_cmd(last_jname,["cputimeraw","maxrss","maxvmsize"],c=True)
                jstr['cputime'][lsub] = last_jinfo[0]
                jstr['maxrss'][lsub] = last_jinfo[1]
                jstr['maxvmsize'][lsub] = last_jinfo[2]
                tner_names = subprocess.getoutput("ls -ltr "+checkfiledir+str(mode)+"_"+str(jcount*maxlines)+"_*.txt").split("\n")
                if ("No such" not in tner_names[0]):
                    if (tner_names[0].split("_")[-2]==tner_names[0].split("_")[-1].split(".")[0]): jstr['done'][lsub] = True
                    jstr['outfile'][lsub] = tner_names[0]
                    bindex = tner_names[0].split("_")[-3]
            # if the killfile exists, stop submitting jobs!
            if os.path.exists(killfile):
                Table(jstr).write(runfile,overwrite=True)
                print("all lines read from MPC file ",listfile,", exiting!")
                sys.exit(0)

        # --Submit a new job--
        print("--Submitting new job--")
        jstr.add_row([jcount*maxlines,mode,False,"",False,"","","","","","",""])
        jsub = -1
        # info of next job to submit
        cmd = "python "+localdir+"mpc_makecat.py "+str(mode)+" "+str(bindex-1)
        jstr['cmd'][jsub] = cmd
        jname = 'mpc_makecat_'+str(t0)+'_'+str(jcount)
        jstr['jobname'][jsub] = jname
        # write job script to file
        jfile=write_jscript(jname,"priority",cmd,outfiledir)
        # submit job to slurm queue
        os.system("sbatch "+jfile)
        jstr['submitted'][jsub] = True
        jcount+=1
        # save job structure, sleep before checking/submitting again
        Table(jstr).write(runfile,overwrite=True)
        print("Job "+jname+" submitted to priority partition; sleeping for 1 min")
        time.sleep(60)


