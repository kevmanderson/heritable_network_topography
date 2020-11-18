#!/bin/python

import os
import subprocess
import numpy as np
import pandas as pd
import nibabel as nb
from scipy.io import loadmat


# submitSlurm/writeSlurm are general utilities for running jobs on our computre cluster
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





# set up directories
# -------
base_dir   = '/gpfs/milgram/project/holmes/kma52/topo_herit'
fs_symlink = os.path.join(base_dir, 'fs_symlink')
anat_dir   = '/gpfs/milgram/data/HCP/ANAT_PREPROCESS'
hcp_df     = pd.read_csv(os.path.join(base_dir, 'data/HCP/hcp_for_solar_allvars.csv'))


# Step 1: symlink freesurfer files
# -------
for sub in hcp_df.Subject:
    print(sub)
    sub      = str(sub)
    fs_dir   = os.path.join(anat_dir, sub, 'T1w', sub)
    sym_orig = os.path.join(fs_symlink, sub)
    sym_dest = os.path.join(fs_symlink, sub)

    os.symlink(fs_dir, sym_dest)



# Step 2: read matlab individualized data
# -------
mat_dat      = loadmat(os.path.join(base_dir, 'data/HCP/HCP_S1200_1029sub_17net_parcellation.mat'))
lh_label_df  = pd.DataFrame(mat_dat['lh_labels'])
rh_label_df  = pd.DataFrame(mat_dat['rh_labels'])
subject_list = [x[0][0] for x in mat_dat['subject_list']]


# path to workbench utiliites
wbc      = '/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command'
wb_short = '/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_shortcuts'
mesh_dir = os.path.join(base_dir, 'external/HCPpipelines/global/templates/standard_mesh_atlases/resample_fsaverage')

# color/net name file
label_table = '/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/kong_net17_labels.txt'


# Step 3: Resample fslr32k individualized parcellations to native freesurfer space
# -------
for idx,sub in enumerate(subject_list):
    print(idx)
    # create output dir to store resampled data
    fs_sub_dir = os.path.join(base_dir, 'fs_surfarea', sub)
    if not os.path.exists(fs_sub_dir):
        os.mkdir(fs_sub_dir)

    for hemi in ['lh','rh']:
        print(hemi)
        # indiv parcellation labels for cur sub
        if hemi == 'lh':
            HEMI='L'
            sub_labels = lh_label_df.iloc[:, idx].values
        elif hemi == 'rh':
            HEMI='R'
            sub_labels = rh_label_df.iloc[:, idx].values

        # Write fslr32k subject parcellatiions to a func.gii file
        gii_obj = nb.gifti.GiftiImage()
        gii_obj.add_gifti_data_array(nb.gifti.GiftiDataArray(sub_labels, datatype='NIFTI_TYPE_FLOAT32'))
        gii_file = os.path.join(fs_sub_dir, '{}_{}_indiv_parc_labels_32k.func.gii'.format(hemi, sub))
        nb.save(gii_obj, gii_file)

        # resample 32k group parcellation to native freesurfer space

        # freesurfer white/pial/sphere native space
        fs_white  = os.path.join(fs_symlink, sub, 'surf/{}.white'.format(hemi))
        fs_pial   = os.path.join(fs_symlink, sub, 'surf/{}.pial'.format(hemi))
        fs_sphere = os.path.join(fs_symlink, sub, 'surf/{}.sphere.reg'.format(hemi))
        fslr_32k_sphere = os.path.join(mesh_dir, 'fs_LR-deformed_to-fsaverage.{}.sphere.32k_fs_LR.surf.gii'.format(HEMI))

        # output and fslr32k group templates
        fs_mid      = os.path.join(fs_sub_dir, '{}.midthickness.surf.gii'.format(hemi))
        fs_32k_mid  = os.path.join(fs_sub_dir, '{}.{}.midthickness.32k_fs_LR.surf.gii'.format(sub, hemi))
        sphere_out  = os.path.join(fs_sub_dir, '{}.sphere.reg.surf.gii'.format(hemi))


        # spherical alignment of  native freesurfer to fslr32k
        resample_prep = f'''{wb_short} -freesurfer-resample-prep \\\n{fs_white} \\\n{fs_pial} \\\n{fs_sphere} \\\n{fslr_32k_sphere} \\\n{fs_mid} \\\n{fs_32k_mid} \\\n{sphere_out}'''
        cmd_long      = resample_prep

        # import parcellation label colors
        label_file  = os.path.join(fs_sub_dir, '{}_{}_indiv_parc_labels_32k.label.gii'.format(hemi, sub))
        wbc_label   = f'''{wbc} -metric-label-import {gii_file} {label_table} {label_file}'''
        cmd_long    = cmd_long + '\n\n' + wbc_label

        # resample individualized parcellation from fslr32k to native space
        fsavg_gii    = os.path.join(fs_sub_dir, '{}_{}_indiv_parc_labels_native_fs.label.gii'.format(hemi, sub))
        resample_cmd = f'''{wbc} -label-resample \\\n{label_file} \\\n{fslr_32k_sphere} \\\n{sphere_out} \\\nADAP_BARY_AREA \\\n{fsavg_gii} \\\n-area-surfs \\\n{fs_32k_mid} \\\n{fs_mid}'''
        cmd_long = cmd_long + '\n\n' + resample_cmd

        # convert label file to annot
        annot_cmd       = f'''export SUBJECTS_DIR=/gpfs/milgram/project/holmes/kma52/topo_herit/fs_symlink\n\n'''
        white_surf_gii  = os.path.join(fs_symlink, sub, 'surf/{}.white.surf.gii'.format(hemi))
        annot_out       = os.path.join(fs_sub_dir, '{}.fs_native_indiv_aparc.annot'.format(hemi))
        annot_cmd       = f'''{annot_cmd}mris_convert --annot {fsavg_gii} {white_surf_gii} {annot_out}\n\n'''
        cmd_long = cmd_long + '\n\n' + annot_cmd

        # run mris_anatomical stats to compute surface area in each individualized parcel
        stat_table  = os.path.join(fs_sub_dir, '{}.fs_native_indiv_aparc_stats.txt'.format(hemi))
        annot_stats = f'''mris_anatomical_stats -a {annot_out} -f {stat_table} -b {sub} {hemi}\n\n'''
        cmd_long = cmd_long + '\n\n' + annot_stats

        # run commands on cluster
        slurm_file = os.path.join(base_dir, 'slurm', '{}_{}_resample'.format(sub, hemi))
        cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', cmd=cmd_long, jobName='{}_{}_resample'.format(sub, hemi), nthreads=1)
        submitSlurm(cmd_path, dependencies=None)


# Step 4: Read freesurfer generative surface area estimates
# ---------------
surf_area_list = []
for idx, sub in enumerate(subject_list):
    print(idx)
    # create output dir to store resampled data
    fs_sub_dir = os.path.join(base_dir, 'fs_surfarea', sub)

    rh_stats          = pd.read_csv(os.path.join(fs_sub_dir, 'rh.fs_native_indiv_aparc_stats.txt'), comment='#', delim_whitespace=True, header=None)
    rh_stats.columns  = ['name', 'n_vert', 'surf_area', 'gray_vol', 'thick_avg', 'thick_std', 'mean_curv', 'gaus_curv', 'fold_in', 'curv_ind']
    rh_stats.name     = rh_stats.name.replace('???', 'None')
    rh_sa_row         = pd.DataFrame(rh_stats.surf_area).transpose()
    rh_sa_row.columns = list('rh_' + rh_stats.name)

    lh_stats          = pd.read_csv(os.path.join(fs_sub_dir, 'lh.fs_native_indiv_aparc_stats.txt'), comment='#', delim_whitespace=True, header=None)
    lh_stats.columns  = ['name', 'n_vert', 'surf_area', 'gray_vol', 'thick_avg', 'thick_std', 'mean_curv', 'gaus_curv', 'fold_in', 'curv_ind']
    lh_stats.name     = lh_stats.name.replace('???', 'None')
    lh_sa_row         = pd.DataFrame(lh_stats.surf_area).transpose()
    lh_sa_row.columns = list('lh_' + lh_stats.name)

    sa_row = pd.concat([lh_sa_row, rh_sa_row], 1)
    sa_row.insert(0, 'id', sub)
    surf_area_list.append(sa_row)


# save surface area dataframe
surf_area_df = pd.concat(surf_area_list)
surf_area_df.to_csv(os.path.join(base_dir, 'data/HCP/indiv_net_surfarea_native_freesurfer.csv'), index=None)







# Do the same as above, but with a different resampling technique (surf2surf)
# Just want to guard against potential mistakes/biases in surface alignment

wbc      = '/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command'
wb_short = '/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_shortcuts'
mesh_dir = os.path.join(base_dir, 'external/HCPpipelines/global/templates/standard_mesh_atlases/resample_fsaverage')

# color/net name file
label_table = '/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/kong_net17_labels.txt'

# Step 5: Resample native freesurfer surface area file to fslr32k space
# -------
for idx,sub in enumerate(subject_list):
    print(idx)
    # create output dir to store resampled data
    fs_sub_dir = os.path.join(base_dir, 'fs_surfarea', sub)
    if not os.path.exists(fs_sub_dir):
        os.mkdir(fs_sub_dir)

    for hemi in ['lh','rh']:
        print(hemi)
        # indiv parcellation labels for cur sub
        if hemi == 'lh':
            HEMI='L'
            sub_labels = lh_label_df.iloc[:, idx].values
            print(sub_labels)
        elif hemi == 'rh':
            HEMI='R'
            sub_labels = rh_label_df.iloc[:, idx].values

        # Write fslr32k subject parcellatiions to a func.gii file
        gii_obj = nb.gifti.GiftiImage()
        gii_obj.add_gifti_data_array(nb.gifti.GiftiDataArray(sub_labels, datatype='NIFTI_TYPE_FLOAT32'))
        gii_file = os.path.join(fs_sub_dir, '{}_{}_indiv_parc_labels_32k.func.gii'.format(hemi, sub))
        nb.save(gii_obj, gii_file)

        # import parcellation label colors
        label_file  = os.path.join(fs_sub_dir, '{}_{}_indiv_parc_labels_32k.label.gii'.format(hemi, sub))
        wbc_label   = f'''{wbc} -metric-label-import \\\n{gii_file} \\\n{label_table} \\\n{label_file}'''
        cmd_long    = wbc_label

        # convert native fs area file to fsaverage
        fs_area_out     = os.path.join(fs_sub_dir, '{}_{}_fsavg_area.nii.gz'.format(hemi, sub))
        native_to_fsavg = 'export SUBJECTS_DIR=/gpfs/milgram/project/holmes/kma52/topo_herit/fs_symlink\n'
        mris_preproc    = f'''mris_preproc --s {sub} --meas area --hemi {hemi} --target fsaverage --out {fs_area_out}'''
        cmd_long = cmd_long + '\n\n' + native_to_fsavg + mris_preproc

        # convert surface area to func.gii
        fs_area = os.path.join(fs_sub_dir, '{}_{}_fsavg_area.nii.gz'.format(hemi, sub))
        fs_area_gii  = os.path.join(fs_sub_dir, '{}_{}_fsavg_area.func.gii'.format(hemi, sub))
        fs_area_conv = f'''mri_convert {fs_area} {fs_area_gii}\n\n'''
        cmd_long = cmd_long + '\n\n' + fs_area_conv

        fsavg_sphere = os.path.join(mesh_dir, 'fsaverage_std_sphere.{}.164k_fsavg_{}.surf.gii'.format(HEMI, HEMI))
        fslr_sphere  = os.path.join(mesh_dir, 'fs_LR-deformed_to-fsaverage.{}.sphere.32k_fs_LR.surf.gii'.format(HEMI))
        cur_area     = os.path.join(mesh_dir, 'fsaverage.{}.midthickness_va_avg.164k_fsavg_{}.shape.gii'.format(HEMI, HEMI))
        new_area     = os.path.join(mesh_dir, 'fs_LR.{}.midthickness_va_avg.32k_fs_LR.shape.gii'.format(HEMI))
        fs_area_out  = os.path.join(fs_sub_dir, '{}_{}_fsavg_area.func.gii'.format(hemi, sub))
        fslr_area    = os.path.join(fs_sub_dir, '{}_{}_fslr32k_area.func.gii'.format(hemi, sub))

        print('FS sphere:\t\t{}'.format(fsavg_sphere))
        print('fslr32k sphere:\t\t{}'.format(fslr_sphere))
        print('FS area:\t\t{}'.format(cur_area))
        print('fslr area:\t\t{}'.format(new_area))
        print('FS area file:\t\t{}'.format(fs_area_out))
        print('FS new area file:\t{}'.format(fslr_area))

        # convert fsaverge to fslr32k
        resample_cmd = f'''{wbc} -metric-resample \\\n{fs_area_out} \\\n{fsavg_sphere} \\\n{fslr_sphere} \\\nADAP_BARY_AREA \\\n{fslr_area} \\\n-area-metrics \\\n{cur_area} \\\n{new_area}'''
        cmd_long = cmd_long + '\n\n' + resample_cmd
        print(cmd_long)

        # run commands on cluster
        slurm_file = os.path.join(base_dir, 'slurm', '{}_{}_v2_resample'.format(sub, hemi))
        cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', cmd=cmd_long, jobName='{}_{}_resample'.format(sub, hemi), nthreads=1)
        submitSlurm(cmd_path, dependencies=None)


label_df   = pd.read_table(label_table, header=None)
region_arr = np.array([x for x in label_df[0] if ' ' not in x])

# cifti concatenation and SA calc
# -------
sub_sa_list = []
for idx, sub in enumerate(subject_list):
    print(idx)
    # create output dir to store resampled data
    fs_sub_dir = os.path.join(base_dir, 'fs_surfarea', sub)

    hemi_sa_dict = pd.DataFrame()
    for hemi in ['lh', 'rh']:
        # indiv parcellation labels for cur sub
        if hemi == 'lh':
            HEMI = 'L'
            sub_labels = np.array(lh_label_df.iloc[:, idx].values)
        elif hemi == 'rh':
            HEMI = 'R'
            sub_labels = np.array(rh_label_df.iloc[:, idx].values)

        fslr_area = os.path.join(fs_sub_dir, '{}_{}_fslr32k_area.func.gii'.format(hemi, sub))
        fslr_obj  = nb.load(fslr_area)
        fslr_data = fslr_obj.agg_data()

        reg_arr = []
        for label in np.arange(1,18):
            roi_sum = np.sum(fslr_data[np.where(sub_labels == label)[0]])
            reg_arr.append(roi_sum)

        roi_mean_out = pd.DataFrame(reg_arr).transpose()
        roi_mean_out.columns = ['{}_{}_sa'.format(hemi, x) for x in region_arr]
        roi_mean_out.insert(0, 'id', sub)

        hemi_sa_dict = pd.concat([hemi_sa_dict, roi_mean_out], 1)

    sub_sa_list.append(hemi_sa_dict)

sub_sa_df = pd.concat(sub_sa_list)
sub_sa_df.to_csv(os.path.join(base_dir, 'data/HCP/method2_sa_fslr32k_freesurfer.csv'), index=None)



















