#!/bin/python

import os
import glob
import numpy as np
import nibabel as nb
import pandas as pd
import scipy.stats
import scipy.io as sio
import multiprocessing as mp

# This script will:
#   1. Concatenate vertex-level h2-multi estimates created by "11_vertex_mat_h2_wrapper.m"
#   2. Create surface plots of vertex-wise heritability


# set up directories
# ------------
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
vert_dir = os.path.join(base_dir, 'data/HCP/vert_similarity')
fig_dir  = os.path.join(base_dir, 'figures')
vert_dir = '/gpfs/milgram/scratch/holmes/topo_herit/data/vert_dice'


# Read dice h2 estimates per vertex
# ------------
l_h2_files = glob.glob(os.path.join(vert_dir, 'L_*dice_network_topology_h2_rr.csv'))
l_df_list  = list()
for h2_in in l_h2_files:
    h2_df = pd.read_csv(h2_in)
    l_df_list.append(h2_df)
    print(h2_in)
    print(h2_df.shape)

# make into a dataframe
l_h2_df = pd.concat(l_df_list)
l_h2_df = l_h2_df.loc[l_h2_df.Var1.notna()]

# save the dataframe
write_df = os.path.join(vert_dir, 'all_verts_dice_L_network_topo_h2_rr.csv')
l_h2_df.to_csv(write_df, index=None)



#####
# RH
#####
r_h2_files = glob.glob(os.path.join(vert_dir, 'R_*dice_network_topology_h2_rr.csv'))
r_df_list = list()
for h2_in in r_h2_files:
    h2_df = pd.read_csv(h2_in)
    r_df_list.append(h2_df)
    print(h2_in)
    print(h2_df.shape)

# make into a data frame
r_h2_df = pd.concat(r_df_list)
r_h2_df = r_h2_df.loc[r_h2_df.Var1.notna()]

# save the dataframe
write_df = os.path.join(vert_dir, 'all_verts_dice_R_network_topo_h2_rr.csv')
r_h2_df.to_csv(write_df, index=None)



################################
# Plot vertex-level heritability
################################
dscalar_template = nb.load('/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/Schaefer2018_400Parcels_17Networks_order.dscalar.nii')
vert_num = len(dscalar_template.dataobj[0])/2


l_data_write = np.zeros(int(vert_num))
l_data_write[l_h2_df.Var5] = l_h2_df.Var1

r_data_write = np.zeros(int(vert_num))
r_data_write[r_h2_df.Var5] = r_h2_df.Var1

# combine LH/RH h2-multi values, save dscalar.nii
dat_arr    = np.array([list(np.concatenate([l_data_write, r_data_write]))])
sub_labels = nb.Cifti2Image(dataobj=dat_arr, header=dscalar_template.header)
out_path   = os.path.join(fig_dir, 'dice_vertex_heritability_vert_radius_rr.dscalar.nii')
nb.save(sub_labels, out_path)
out_path



mat_dat   = sio.loadmat(os.path.join(base_dir, 'data/HCP/HCP_S1200_1029sub_17net_parcellation.mat'))
net_names = [x[0] for x in labels['network_name'][0]]
net_dict  = {i+1:name for i,name in enumerate(net_names)}

lh_net_mode = pd.DataFrame(mat_dat['lh_labels']).mode(1)
rh_net_mode = pd.DataFrame(mat_dat['rh_labels']).mode(1)

r_h2_df = pd.read_csv(os.path.join(vert_dir, 'all_verts_dice_R_network_topo_h2_rr.csv'))
l_h2_df = pd.read_csv(os.path.join(vert_dir, 'all_verts_dice_L_network_topo_h2_rr.csv'))

l_h2_df['atlas_label']   = list(lh_net_mode[0][l_h2_df.Var5])
l_h2_df['atlas_parcels'] = list([net_dict[x] for x in l_h2_df.atlas_label])

r_h2_df['atlas_label']   = list(rh_net_mode[0][r_h2_df.Var5])
r_h2_df['atlas_parcels'] = list([net_dict[x] for x in r_h2_df.atlas_label])


# mean and of heteromodal/unimodal vertices
h2_df = pd.concat([l_h2_df, r_h2_df])
h2_df.Var1[h2_df.atlas_parcels.str.contains('Som|Vis|Aud')].mean()
h2_df.Var1[h2_df.atlas_parcels.str.contains('Som|Vis|Aud')].std()

h2_df.Var1[~h2_df.atlas_parcels.str.contains('Som|Vis|Aud')].mean()
h2_df.Var1[~h2_df.atlas_parcels.str.contains('Som|Vis|Aud')].std()

h2_df.to_csv(os.path.join(base_dir, 'data/HCP/h2_vert_df.csv'), index=None)


rsn_obj = nb.load('/gpfs/milgram/project/holmes/kma52/topo_herit/ref_files/null_lL_WG33/RSN-networks.32k_fs_LR.dlabel.nii')


# read network name information
labels    = sio.loadmat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = unlist(labels$network.name)


l_h2_df['atlas_label'] = group_atlas.dataobj[0][l_h2_df.Var5]
l_h2_df['atlas_parcels'] = [schaefer_dict[x] for x in l_h2_df.atlas_label]

r_h2_df['atlas_label'] = group_atlas.dataobj[0][r_h2_df.Var5 + 32492]
r_h2_df['atlas_parcels'] = [schaefer_dict[x] for x in r_h2_df.atlas_label]

uni_grep = ['SomMot', 'Vis','']

len(np.where(r_h2_df['atlas_label'] == 0)[0])


h2_lm_df = h2_df %>% group_by(Var4, atlas_parcels) %>% summarise(m_h2=mean(Var1))
h2_lm_df$type = ifelse(grepl('Vis|Som|Aud', h2_lm_df$atlas_parcels), 'uni', 'het')

# end
# end
# end
# end
