#!/bin/python


import os
import glob
import shutil
import numpy as np
import nibabel as nb
import pandas as pd
import scipy.io as sio


def make_parcellation_ciis(sub, cii_out_dir, subject_list, lh_labels, rh_labels, net_names):
    print(sub)
    sub_idx = np.where(subject_list == sub)[0][0]

    # lh
    sub_lh = lh_labels[:,sub_idx]
    # rh
    sub_rh = rh_labels[:,sub_idx]
    sub_rh[sub_rh != 0] = sub_rh[sub_rh != 0] + 17
    bihemi_labels = np.concatenate([sub_lh, sub_rh])

    # read parcellation template
    dscalar_template = nb.load('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii')

    for idx in np.arange(1,18):
        sub_labels = nb.Cifti2Image(dataobj=np.asarray([bihemi_labels]), header=dscalar_template.header)
        tmp_labels = sub_labels.dataobj[0]
        tmp_labels[tmp_labels != idx] = 0
        write_labels = np.array([list(tmp_labels)])
        cur_name  = net_names[idx - 1]
        cii_write = nb.cifti2.Cifti2Image(dataobj=write_labels, header=sub_labels.header)
        out_path         = os.path.join(cii_out_dir, '{}_{}_kong_indiv_parcellation.dscalar.nii'.format(sub, cur_name))
        nb.save(cii_write, out_path)

        labelfile  = '/gpfs/milgram/project/holmes/kma52/topo_herit/data/HCP/hcp_net17_lhrh_colortable.txt'
        dlabel_out = out_path.replace('dscalar.nii','dlabel.nii')
        print(dlabel_out)
        add_colors = 'wb_command -cifti-label-import ' + out_path + ' ' + labelfile + ' ' + dlabel_out
        os.system(add_colors)

    return(dlabel_out)



def mz_get_most_similar_subs(mz_df, dice_mat):
    mz_dice_list = []
    for id in set(mz_df.mztwin):
        mz_subs   = mz_df.loc[mz_df.mztwin == id]
        keep_idxs = dice_mat.columns.isin(mz_subs.id.astype(str))
        mz_dice   = dice_mat.loc[keep_idxs, keep_idxs].iloc[0, 1]

        out_row = pd.Series({'dice' : np.array(mz_dice), 'sub_A':mz_subs.id.iloc[0], 'sub_B':mz_subs.id.iloc[1]})
        mz_dice_list.append(out_row)
    mz_df = pd.DataFrame(mz_dice_list)
    return(mz_df)



sub_a = '128329'
sub_b = '837964'
net_num = 6
def two_sub_overlap(sub_a, sub_b, net, net_num, subject_list, lh_labels, rh_labels):

    # identify subject indices for each sub in the pair
    sub_a_idx = np.where(subject_list == sub_a)[0][0]
    sub_b_idx = np.where(subject_list == sub_b)[0][0]

    # extract subject labels
    bihemi_a_labels = np.concatenate([lh_labels[:,sub_a_idx], rh_labels[:,sub_a_idx]])
    bihemi_b_labels = np.concatenate([lh_labels[:,sub_b_idx], rh_labels[:,sub_b_idx]])

    # set all labels not part of the current network to zero
    bihemi_a_labels[bihemi_a_labels != net_num] = 0
    bihemi_b_labels[bihemi_b_labels != net_num] = 0

    # combine subA/subB labels; 0=None; 1=SubA only; 2=SubB Only; 3=Both
    bihemi_a_labels[bihemi_a_labels == net_num] = 1
    bihemi_b_labels[bihemi_b_labels == net_num] = 2
    combined_labels = np.vstack((bihemi_a_labels, bihemi_b_labels)).sum(0)

    # load dscalar templace; plug in labels and write to disk
    cii_out_dir      = '/gpfs/milgram/project/holmes/kma52/topo_herit/figures'
    dscalar_template = nb.load('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti/Schaefer2018_200Parcels_17Networks_order.dlabel.nii')
    cii_write = nb.cifti2.Cifti2Image(dataobj=np.array([combined_labels]), header=dscalar_template.header)
    out_path  = os.path.join(cii_out_dir, '{}_{}_{}_kong_indiv_parcellation.dscalar.nii'.format(sub_a, sub_b, net))
    nb.save(cii_write, out_path)

    # plug in labels from pre-created file
    labelfile  = '/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/{}_colors.txt'.format(net)
    dlabel_out = out_path.replace('dscalar.nii', 'dlabel.nii')
    print(dlabel_out)
    add_colors = 'wb_command -cifti-label-import ' + out_path + ' ' + labelfile + ' ' + dlabel_out
    os.system(add_colors)


    # write file with only the overlap vertices (for border creation)
    combined_labels[combined_labels != 3] = 0
    cii_write = nb.cifti2.Cifti2Image(dataobj=np.array([combined_labels]), header=dscalar_template.header)
    out_path  = os.path.join(cii_out_dir, '{}_{}_{}_kong_indiv_parcellation_overlap_only.dscalar.nii'.format(sub_a, sub_b, net))
    nb.save(cii_write, out_path)

    # add black color
    labelfile  = '/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/{}_colors_border.txt'.format(net)
    dlabel_out = out_path.replace('.dscalar.nii', '_border.dlabel.nii')
    print(dlabel_out)
    add_colors = 'wb_command -cifti-label-import ' + out_path + ' ' + labelfile + ' ' + dlabel_out
    os.system(add_colors)

    l_surface  = '/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii'
    r_surface  = '/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii'
    make_border = f'''wb_command -cifti-label-to-border {dlabel_out} -border {l_surface} {dlabel_out.replace('.dlabel.nii', '_LH.border')} -border {r_surface} {dlabel_out.replace('.dlabel.nii', '_RH.border')}'''
    os.system(make_border)
    print(dlabel_out.replace('.dlabel.nii', '_LH.border'))



#####################################
# Set up required directories/paths
#####################################
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
vert_dir = os.path.join(base_dir, 'data/HCP/vert_similarity')
fig_dir  = os.path.join(base_dir, 'figures')
ref_dir  = os.path.join(base_dir, 'ref_files')



#####################################
# read individualized parcellations
#####################################
mat_dat   = sio.loadmat(os.path.join(base_dir, 'data/HCP/HCP_S1200_1029sub_17net_parcellation.mat'))
lh_labels = mat_dat['lh_labels']
rh_labels = mat_dat['rh_labels']
subject_list = mat_dat['subject_list']

labels    = sio.loadmat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = [labels['network_name'][0][i][0] for i in np.arange(0,17)]



##########################
# read HCP phenotype data
##########################
hcp_df   = pd.read_csv(os.path.join(base_dir, 'data/HCP/hcp_wMRI_for_solar.csv'))
mz_df    = hcp_df[~hcp_df.mztwin.isna()]
nonmz_df = hcp_df[hcp_df.mztwin.isna()]





##################################################
# Identify subjects with max/min Dice similarity
##################################################
net_names_idx = x + net_names

################
# Default B = 6
# Visual A = 8
################
mz_sub_list = []
dice_pairs  = []
for net in ['6', '9']:
    print(net)

    ################
    # Read HCP Dice coefficients for this network
    ################
    net_dice  = pd.read_csv(os.path.join(base_dir, 'data/HCP/hcp_bihemi_net_dice_net_{}_matrix.csv'.format(net)))
    keep_idxs = net_dice.columns.isin(mz_df.id.astype(str))
    mz_matrix = net_dice.loc[keep_idxs, keep_idxs]


    ################
    # find the most topographically similar twins
    ################
    mz_dice_df  = mz_get_most_similar_subs(mz_df=mz_df, dice_mat=mz_matrix)
    mz_dice_ordered = mz_dice_df.sort_values('dice', ascending=False)
    mz_sub_list = mz_sub_list + list(mz_dice_ordered.iloc[1, 1:3])
    high_pair   = mz_dice_ordered.head(1)
    dice_pairs.append((int(high_pair.sub_A), int(high_pair.sub_B)))
    print(mz_dice_ordered.head(1))


    ################
    # Subset to non-monozygotic twin dataframe to look for dissimilar twins
    ################
    nonmz_df = hcp_df[hcp_df.mztwin.isna()]
    nonmz_df = nonmz_df.reset_index()

    # non MZ matrix
    keep_idxs = net_dice.columns.isin(nonmz_df.id.astype(str))
    nonmz_matrix = net_dice.loc[keep_idxs, keep_idxs]
    nonmz_matrix = nonmz_matrix.reset_index()


    ################
    # flatten dice matrix
    ################
    df = nonmz_matrix.where(np.triu(np.ones(nonmz_matrix.shape)).astype(np.bool))
    df = df.stack().reset_index()
    df.columns = ['id', 'sub', 'dice']
    df = df.loc[df.dice != 1]
    df = df.loc[df.dice != 0]

    # fine minimum subject pair
    min_idx = df.loc[df.dice == df.dice.min()]
    sub_A = list(nonmz_df.id[min_idx.id])[0]
    sub_B = list(min_idx['sub'])[0]
    print(str(sub_A), str(sub_B))
    dice_pairs.append((str(sub_A), str(sub_B)))

    low_dice_subs = low_dice_subs + [str(sub_A), str(sub_B)]



############################################
# plot data for Default B high/low subjects
############################################
two_sub_overlap(sub_a=str(dice_pairs[0][0]), sub_b=str(dice_pairs[0][1]), net='DefaultB', net_num=6,
                subject_list=subject_list, lh_labels=lh_labels, rh_labels=rh_labels)

two_sub_overlap(sub_a=str(dice_pairs[1][0]), sub_b=str(dice_pairs[1][1]), net='DefaultB', net_num=6,
                subject_list=subject_list, lh_labels=lh_labels, rh_labels=rh_labels)



two_sub_overlap(sub_a=str(dice_pairs[2][0]), sub_b=str(dice_pairs[2][1]), net='VisualA', net_num=9,
                subject_list=subject_list, lh_labels=lh_labels, rh_labels=rh_labels)
two_sub_overlap(sub_a=str(dice_pairs[3][0]), sub_b=str(dice_pairs[3][1]), net='VisualA', net_num=9,
                subject_list=subject_list, lh_labels=lh_labels, rh_labels=rh_labels)


cii_out_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit/output/hcp/indiv_parcellations'
for sub in mz_sub_list + low_dice_subs:
    print(sub)
    cii_label = make_parcellation_ciis(sub=str(sub), cii_out_dir=cii_out_dir, subject_list=subject_list,
                                        lh_labels=lh_labels, rh_labels=rh_labels, net_names=net_names)


# non-monozygotic twin dataframe
nonmz_df = hcp_df[hcp_df.mztwin.isna()]
nonmz_df = nonmz_df.reset_index()

# non MZ matrix
keep_idxs = net_6.columns.isin(nonmz_df.id.astype(str))
nonmz_matrix = net_6.loc[keep_idxs, keep_idxs]
nonmz_matrix = nonmz_matrix.reset_index()

# flatten dice matrix
df = nonmz_matrix.where(np.triu(np.ones(nonmz_matrix.shape)).astype(np.bool))
df = df.stack().reset_index()
df.columns = ['id','sub','dice']
df = df.loc[df.dice != 1]
df = df.loc[df.dice != 0]

min_idx = df.loc[df.dice == df.dice.min()]
sub_A = list(nonmz_df.id[min_idx.id])[0]
sub_B = '281135'

for sub in [sub_A, sub_B]:
    cii_label = make_parcellation_ciis(sub=str(sub), cii_out_dir=cii_out_dir, subject_list=subject_list,
                                       lh_labels=lh_labels, rh_labels=rh_labels, net_names=net_names)









net_8     = pd.read_csv(os.path.join(base_dir, 'data/HCP/hcp_bihemi_net_dice_net_8_matrix.csv'))
keep_idxs = net_8.columns.isin(mz_df.id.astype(str))
mz_matrix = net_8.loc[keep_idxs, keep_idxs]


mz_net8_dice_df = mz_get_most_similar_subs(mz_df=mz_df, dice_mat=mz_matrix)
mz_net8_dice_ordered = mz_net8_dice_df.sort_values('dice', ascending=False)
mz_net8_dice_ordered.head()


cii_out_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit/output/hcp/indiv_parcellations'
cii_label = make_parcellation_ciis(sub=str(111312), cii_out_dir=cii_out_dir, subject_list=subject_list,
                                   lh_labels=lh_labels, rh_labels=rh_labels, net_names=net_names)
cii_label = make_parcellation_ciis(sub=str(144226), cii_out_dir=cii_out_dir, subject_list=subject_list,
                                   lh_labels=lh_labels, rh_labels=rh_labels, net_names=net_names)


nonmz_df = hcp_df[hcp_df.mztwin.isna()]
nonmz_df = nonmz_df.reset_index()
keep_idxs = net_6.columns.isin(nonmz_df.id.astype(str))
nonmz_matrix = net_6.loc[keep_idxs, keep_idxs]
nonmz_matrix = nonmz_matrix.reset_index()
df = nonmz_matrix.where(np.triu(np.ones(nonmz_matrix.shape)).astype(np.bool))
df = df.stack().reset_index()
df.columns = ['id','sub','dice']
df = df.loc[df.dice != 1]
df = df.loc[df.dice != 0]

min_idx = df.loc[df.dice == df.dice.min()]

sub_A = list(nonmz_df.id[min_idx.id])[0]
sub_B = list(min_idx['sub'])[0]


cii_label = make_parcellation_ciis(sub=str(sub_A), cii_out_dir=cii_out_dir, subject_list=subject_list,
                                   lh_labels=lh_labels, rh_labels=rh_labels, net_names=net_names)
cii_label = make_parcellation_ciis(sub=str(sub_B), cii_out_dir=cii_out_dir, subject_list=subject_list,
                                   lh_labels=lh_labels, rh_labels=rh_labels, net_names=net_names)










