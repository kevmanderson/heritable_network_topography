#!/bin/python

import os
import sys
import glob
import numpy as np
import scipy.io as io
import pandas as pd
import nibabel as nb
import subprocess

# This script will:
#   1. Act as a wrapper for calling the "10_vertex_dice_matrix_slurm.R" script
#   2. will compute the subj-subj dice matrix for each vertex
#   3. submits chunks of vertices to be processed in parallel on a cluster


def submitSlurm(cmd, dependencies=None):
    if dependencies != None:
        # execute
        p = subprocess.Popen(['sbatch', '--dependency=afterany:' + ':'.join(dependencies), cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    else:
        p = subprocess.Popen(['sbatch', cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    # get slurm job id to set up job submission dependency
    job_id = str(out).split(' ')[-1].replace("\\n'", '')
    return(job_id)


def writeSlurm(slurm_file, partition, cmd, jobName, stime='6:00:00', nthreads=None, mem=None):
    '''
    Submit batch job to Yale Milgram cluster (SLURM)

    required arguments:
        - slurm_file        base filepath string for writing slurm command and output
                                (e.g. /gpfs/milgram/project/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/MRI/100024/slurm/100024_hello_world_cmd)
        - partition         short/long/scavenge
                                (e.g. short)
        - nthreads          up to 28 cpus on a node
                                (e.g. 4)
        - cmd           full string of the command to be written/run
                            (e.g. module load Apps/FREESURFER/5.3.0\\n print('Hello World'))
    '''
    slurm_name = slurm_file + '_slurm.txt'
    slurm_out  = slurm_file + '_slurmOut.txt'
    slurm_file = open(slurm_name, "w")
    slurm_file.write('#!/bin/bash\n')
    slurm_file.write('#SBATCH --partition=' +  partition + '\n')
    slurm_file.write('#SBATCH --output=' + slurm_out + '\n')
    slurm_file.write('#SBATCH --nodes=1\n')
    if mem != None:
        slurm_file.write('#SBATCH --mem=' + str(mem) + '\n')
    if nthreads != None:
        slurm_file.write('#SBATCH --ntasks=1 --cpus-per-task=' + str(nthreads) + '\n')
    slurm_file.write('#SBATCH --job-name=' + jobName + '\n')
    slurm_file.write('#SBATCH --time=' + stime + '\n')
    slurm_file.write(str(cmd))
    slurm_file.close()
    subprocess.call(['chmod', '0770', slurm_name])
    return(slurm_name)




base_dir    = '/gpfs/milgram/project/holmes/kma52/topo_herit'
script_path = os.path.join(base_dir, 'scripts/11_vertex_dice_matrix_slurm.R')
slurm_dir   = os.path.join(base_dir, 'slurm')
for chunk in np.arange(1,41):
    print(chunk)
    run_this = f'''source ~/oldRbashrc\n/gpfs/milgram/apps/hpc.rhel7/software/R/3.4.4-foss-2018a-X11-20180131/bin/Rscript {script_path} {chunk} 40 L'''
    print(run_this)
    script_file = os.path.join(slurm_dir, 'dice_mat_chunk{}'.format(chunk))

    cmd_path = writeSlurm(slurm_file=script_file, partition='long', nthreads=20, cmd=run_this, stime='48:00:00', jobName=str(chunk))
    submitSlurm(cmd_path, dependencies=None)




base_dir    = '/gpfs/milgram/project/holmes/kma52/topo_herit'
script_path = os.path.join(base_dir, 'scripts/11_vertex_dice_matrix_slurm.R')
slurm_dir   = os.path.join(base_dir, 'slurm')
#for chunk in np.arange(1,41):
for chunk in np.arange(21, 40):
    print(chunk)
    run_this = f'''source ~/oldRbashrc\n/gpfs/milgram/apps/hpc.rhel7/software/R/3.4.4-foss-2018a-X11-20180131/bin/Rscript {script_path} {chunk} 40 R'''
    print(run_this)
    script_file = os.path.join(slurm_dir, 'dice_mat_chunk{}'.format(chunk))

    cmd_path = writeSlurm(slurm_file=script_file, partition='long', nthreads=20, cmd=run_this, stime='48:00:00', jobName=str(chunk))
    submitSlurm(cmd_path, dependencies=None)










