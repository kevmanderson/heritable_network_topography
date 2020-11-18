#!/bin/python

import os
import sys
import glob
import copy
import numpy as np
import scipy.io as sio
import pandas as pd
import nibabel as nb
import subprocess


# This script will:
#   1. Plot example individualized parcellation labels w/in a vertex for an MZ/UNR HCP pair



#############
# set up dirs
#############
out_dir  = '/gpfs/milgram/scratch/holmes/topo_herit/data/vert_dice'
roi_dir  = '/gpfs/milgram/project/holmes/kma52/topo_herit/data/HCP/roi_dir'
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'


# read individualized parcellations
mat_dat   = sio.loadmat(os.path.join(base_dir, 'data/HCP/HCP_S1200_1029sub_17net_parcellation.mat'))
lh_labels = mat_dat['lh_labels']
rh_labels = mat_dat['rh_labels']
subject_list = mat_dat['subject_list']
subject_list = [x[0][0] for x in subject_list]


# phenotype data
hcp_df = pd.read_csv(os.path.join(base_dir, 'data/HCP/hcp_wMRI_for_solar.csv'))
hcp_df_order = hcp_df.sort_values('id')

# select a vertex to plot as an example
vert_to_plot  = 13226 # 16302, 30721
lh_vert_label = copy.deepcopy(lh_labels)


# adjacent vertices
neighborhood_verts = os.path.join(roi_dir, 'adjacency_L_vertex{}_idxfrom0.txt'.format(vert_to_plot))
neigh_vert_df = pd.read_csv(neighborhood_verts)


# read dice matrix for example vertex
lh_vert_dice = pd.read_csv(os.path.join(out_dir, 'local_dice_L_vertex_{}_idxfrom0_radius10_matrix.csv'.format(vert_to_plot)), header=None)
lh_vert_dice = lh_vert_dice.values


# extract upper tri of dice matrix
all_dice     = lh_vert_dice[np.triu_indices(lh_vert_dice.shape[0], 1)]
all_dice_rev = np.array(sorted(all_dice, reverse=True))
all_dice     = np.array(sorted(all_dice, reverse=False))


# identify MZ twin pair with highest similarity
mz_fams = list(set(hcp_df_order.famid[hcp_df_order.Zyg_combined == 'MZ']))
mz_dice_df = list()
for fam in mz_fams:
    cur_fam = hcp_df_order.loc[hcp_df_order.famid == fam]
    cur_fam = cur_fam.loc[cur_fam.Zyg_combined == 'MZ']
    row_use = np.where(hcp_df_order.id == cur_fam.id.iloc[0])[0][0]
    col_use = np.where(hcp_df_order.id == cur_fam.id.iloc[1])[0][0]
    out_row = [str(cur_fam.id.iloc[0]), str(cur_fam.id.iloc[1]), lh_vert_dice[row_use, col_use]]
    mz_dice_df.append(out_row)
mz_dice_df = pd.DataFrame(mz_dice_df)
use_pair = mz_dice_df.loc[mz_dice_df[2] == np.max(mz_dice_df[2])]
mz_sub_a = use_pair[0]
mz_sub_b = use_pair[1]

high_row = np.where(hcp_df_order.id == int(mz_sub_a.values[0]))[0][0]
high_col  = np.where(hcp_df_order.id == int(mz_sub_b.values[0]))[0][0]


# just pick a high example subject pair
#match_pairs = np.where(lh_vert_dice == all_dice_rev[3])
#high_row = match_pairs[0][3]
#high_col = match_pairs[1][3]

hcp_df_order.iloc[high_row,0:15]
hcp_df_order.iloc[high_col,0:15]

# ..and low overlap pair
low_match_pairs = np.where(lh_vert_dice == all_dice[0])
low_row = low_match_pairs[0][0]
low_col = low_match_pairs[1][0]

plot_subjects = [high_row, high_col, low_row, low_col]
plot_ids = [hcp_df_order.id[x] for x in plot_subjects]



plot_vert = vert_to_plot
for idx in plot_subjects:
    #tri_idxs = np.triu_indices(n=lh_vert_dice.shape[0], k=1)

    plot_subj = hcp_df.id[idx]
    #plot_idx  = np.where(np.isin(hcp_df_order.id, str(plot_subj)))[0][0]
    bihemi_labels = np.concatenate([lh_labels[:,idx], rh_labels[:,idx]])

    idx_arr = np.arange(len(bihemi_labels))
    bihemi_labels[~np.isin(idx_arr, neigh_vert_df.vert.values)] = 0

    # read parcellation template
    dscalar_template = nb.load('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii')
    sub_labels = nb.Cifti2Image(dataobj=np.asarray([list(bihemi_labels)]), header=dscalar_template.header)
    out_path   = os.path.join(base_dir, 'figures/kma_{}_vert{}_sub_labels.dscalar.nii'.format(plot_subj, plot_vert))
    nb.save(sub_labels, out_path)
    out_path

    labelfile = '/gpfs/milgram/project/holmes/kma52/topo_herit/data/HCP/hcp_net17_lhrh_colortable.txt'
    dlabel_out = out_path.replace('dscalar.nii', 'dlabel.nii')
    add_colors = 'wb_command -cifti-label-import ' + out_path + ' ' + labelfile + ' ' + dlabel_out
    os.system(add_colors)


# make black ROI border
bihemi_labels[~np.isin(idx_arr, neigh_vert_df.vert.values)] = 0
bihemi_labels[np.isin(idx_arr, neigh_vert_df.vert.values)] = 1
sub_labels = nb.Cifti2Image(dataobj=np.asarray([list(bihemi_labels)]), header=dscalar_template.header)
out_path = os.path.join(base_dir, 'figures/vert{}_black_roi.dscalar.nii'.format(plot_vert))
nb.save(sub_labels, out_path)

labelfile = '/gpfs/milgram/project/holmes/kma52/topo_herit/data/HCP/black_hcp_net17_lhrh_colortable.txt'
dlabel_out = out_path.replace('dscalar.nii', 'dlabel.nii')
add_colors = 'wb_command -cifti-label-import ' + out_path + ' ' + labelfile + ' ' + dlabel_out
os.system(add_colors)

l_surf  = '/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'
l_label = dlabel_out.replace('.dlabel.nii', '_L.label.gii')
cifti_sep = 'wb_command -cifti-separate ' + dlabel_out + ' COLUMN -label CORTEX_LEFT ' + l_label
os.system(cifti_sep)
border_cmd = 'wb_command -label-to-border ' + l_surf + ' ' + l_label + ' ' + l_label.replace('.label.gii', '.border')
os.system(border_cmd)




# end
# end
# end










