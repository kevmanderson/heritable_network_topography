library(tidyverse)
library(R.matlab)
library(foreach)
library(parallel)
library(aricode)


# This script will:
#   1. Read individualized parcellation data from Kong et al (2018)
#   2. In parallel, calculate dice coefficient between each subject.
#           Using overall network similarity, and separately for each network
#   3. Save multiple forms of the data
#       list: e.g. "hcp_dice_overlaps.Rdata"
#       long dataframe: e.g. "hcp_dice_overlaps_all_nets_df.csv"
#       matrix: e.g. "hcp_bihemi_net_1_matrix.csv"


dice_function = function(idx, uniq_combos, just_labels){

    # idx: int, required
    #       e.g. 100. will iterate over all possible pairwise comparisons
    # uniq_combos: matrix, required
    #       a matrix of subj-to-subj pairwise combinations
    # just_labels: matrix, required
    #       matrix of individualized parcellations; cols=59,411 non midline verices; rows=subjects

    print(idx)

    # subject A / subject B
    a = uniq_combos[1,idx]
    b = uniq_combos[2,idx]

    # network labels for subjects A/B
    x = just_labels[a,]
    y = just_labels[b,]

    # allocate matrix nrow=n_vertices; ncol=n_networks
    n_vert = length(x)
    arr    = 1:n_vert
    net_mat_a = matrix(0, n_vert, 17)
    net_mat_b = matrix(0, n_vert, 17)

    # convert array of network calls to matrix format
    net_mat_a[cbind(arr, x)] = 1
    net_mat_b[cbind(arr, y)] = 1

    # find network-wise overlap
    net_mult = net_mat_b*net_mat_a

    overlap_sums = colSums(net_mult)
    net_cts      = colSums(net_mat_a) + colSums(net_mat_b)
    net_dice     = (2*overlap_sums)/net_cts

    total_overlap = sum(rowSums(net_mult))/length(x)
    dice_df = data.frame(t(net_dice))
    colnames(dice_df) = paste0('dice_net_', 1:17)
    dice_df$subA = a
    dice_df$subB = b
    dice_df$overall_dice = total_overlap
    dice_df$idx  = idx

    return(dice_df)
}



# read individual parcellations from Kong et al 2018
# ------------------
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
mat_dat  = readMat(paste0(base_dir, '/data/HCP/HCP_S1200_1029sub_17net_parcellation.mat'))

# LH
lh_labels    = mat_dat$lh.labels
lh_labels_df = as.data.frame(t(lh_labels))
colnames(lh_labels_df) = paste0('lh_', colnames(lh_labels_df))
lh_labels_df$subject   = unlist(mat_dat$subject.list)

# RH
rh_labels    = mat_dat$rh.labels
rh_labels_df = as.data.frame(t(rh_labels))
colnames(rh_labels_df) = paste0('rh_', colnames(rh_labels_df))
rh_labels_df$subject   = unlist(mat_dat$subject.list)


# merge left and right hemisphere topology
bihemi_df = merge(x=lh_labels_df, y=rh_labels_df, by='subject')

# grab columns with network info
vert_cols = colnames(bihemi_df)[grep('lh|rh', colnames(bihemi_df))]
vert_sums = colSums(bihemi_df[vert_cols])

# identify non midline surface vertices
bihemi_use_df = bihemi_df[,which(vert_sums != 0)]




# calculate dice coefficients from network topology
# using all bihemispheric vertices
# ------------------
vert_use_cols = colnames(bihemi_use_df)[grepl('lh|rh', colnames(bihemi_use_df))]

# make dataframe with network labels (columns=vertices; rows=subjects)
just_labels     = as.matrix(bihemi_use_df[,vert_use_cols])
n_subs          = nrow(just_labels)
uniq_combos     = combn(1:n_subs, 2)
net_just_labels = just_labels




# BI-HEMISPHERIC
# compute the subject by subject dice in parallel
# ------------------
dice_list  = mclapply(1:ncol(uniq_combos), dice_function, uniq_combos=uniq_combos, just_labels=just_labels, mc.cores = 24)
save(x=dice_list, file=paste0(base_dir, '/data/HCP/hcp_dice_overlaps.Rdata'))

# load dice, convert to df, and save
load(file=paste0(base_dir, '/data/HCP/hcp_dice_overlaps.Rdata'), verbose=T)
#x = data.table::rbindlist(dice_list)
dice_df = do.call('rbind', dice_list)
write_csv(dice_df, paste0(base_dir, '/data/HCP/hcp_dice_overlaps_all_nets_df.csv'))



# read df version of the dice overlaps
# ------------------
dice_df = read.csv(paste0(base_dir, '/data/HCP/hcp_dice_overlaps_all_nets_df.csv'))

# convert from long to to matrix and write csv file for each network
# matrix format is required for later h2 estimation
net_measures = colnames(dice_df)[grep('dice', colnames(dice_df))]
for (dice_net in net_measures){
    write(dice_net,'')
    net_similarity = matrix(NA, 1029, 1029)
    net_similarity[cbind(dice_df$subA, dice_df$subB)] = dice_df[[dice_net]]
    net_similarity[cbind(dice_df$subB, dice_df$subA)] = dice_df[[dice_net]]

    diag(net_similarity) = 1
    net_similarity = as.data.frame(net_similarity)
    colnames(net_similarity) = bihemi_df$subject
    rownames(net_similarity) = bihemi_df$subject
    write_csv(net_similarity, paste0(base_dir, '/data/HCP/hcp_bihemi_net_',dice_net,'_matrix.csv'))
}



# Left Hemisphere
# calculate dice coefficients from network topology
# ------------------
lh_vert_use_cols = colnames(bihemi_use_df)[grepl('lh', colnames(bihemi_use_df))]
lh_just_labels   = as.matrix(bihemi_use_df[,lh_vert_use_cols])
n_subs           = nrow(lh_just_labels)
uniq_combos      = combn(1:n_subs, 2)


# calculate dice in parallel, save
lh_dice_list = mclapply(1:ncol(uniq_combos), dice_function, uniq_combos=uniq_combos, just_labels=lh_just_labels, mc.cores = 24)
save(x=lh_dice_list, file=paste0(base_dir, '/data/HCP/lh_hcp_dice_overlaps.Rdata'))

# reformat to dataframe
load(file=paste0(base_dir, '/data/HCP/lh_hcp_dice_overlaps.Rdata'), verbose=T)
lh_dice_df = do.call('rbind', lh_dice_list)
write_csv(lh_dice_df, paste0(base_dir, '/data/HCP/hcp_lh_dice_overlaps_all_nets_df.csv'))

# read overlap metrics
lh_dice_df   = read.csv(paste0(base_dir, '/data/HCP/hcp_lh_dice_overlaps_all_nets_df.csv'))
net_measures = colnames(lh_dice_df)[grep('dice', colnames(lh_dice_df))]
for (dice_net in net_measures){
    write(dice_net,'')
    net_similarity = matrix(NA, 1029, 1029)
    net_similarity[cbind(lh_dice_df$subA, lh_dice_df$subB)] = lh_dice_df[[dice_net]]
    net_similarity[cbind(lh_dice_df$subB, lh_dice_df$subA)] = lh_dice_df[[dice_net]]

    diag(net_similarity) = 1
    net_similarity = as.data.frame(net_similarity)
    colnames(net_similarity) = bihemi_df$subject
    rownames(net_similarity) = bihemi_df$subject
    write_csv(net_similarity, paste0(base_dir, '/data/HCP/hcp_lh_net_',dice_net,'_matrix.csv'))
}




# Right hemisphere
# calculate dice coefficients from network topology, using all bihemispheric vertices
# ------------------
rh_vert_use_cols = colnames(bihemi_use_df)[grepl('rh', colnames(bihemi_use_df))]
rh_just_labels   = as.matrix(bihemi_use_df[,rh_vert_use_cols])
n_subs      = nrow(rh_just_labels)
uniq_combos = combn(1:n_subs, 2)


rh_dice_list = mclapply(1:ncol(uniq_combos), dice_function, uniq_combos=uniq_combos, just_labels=rh_just_labels, mc.cores = 24)
save(x=rh_dice_list, file=paste0(base_dir, '/data/HCP/rh_hcp_dice_overlaps.Rdata'))
load(file=paste0(base_dir, '/data/HCP/rh_hcp_dice_overlaps.Rdata'), verbose=T)

rh_dice_df = do.call('rbind', rh_dice_list)
write_csv(rh_dice_df, paste0(base_dir, '/data/HCP/hcp_rh_dice_overlaps_all_nets_df.csv'))

# read overlap metrics
rh_dice_df = read.csv(paste0(base_dir, '/data/HCP/hcp_rh_dice_overlaps_all_nets_df.csv'))

net_measures = colnames(rh_dice_df)[grep('dice', colnames(rh_dice_df))]
for (dice_net in net_measures){
    write(dice_net,'')
    net_similarity = matrix(NA, 1029, 1029)
    net_similarity[cbind(rh_dice_df$subA, rh_dice_df$subB)] = rh_dice_df[[dice_net]]
    net_similarity[cbind(rh_dice_df$subB, rh_dice_df$subA)] = rh_dice_df[[dice_net]]

    diag(net_similarity) = 1
    net_similarity = as.data.frame(net_similarity)
    colnames(net_similarity) = bihemi_df$subject
    rownames(net_similarity) = bihemi_df$subject
    write_csv(net_similarity, paste0(base_dir, '/data/HCP/hcp_rh_net_',dice_net,'_matrix.csv'))
}

