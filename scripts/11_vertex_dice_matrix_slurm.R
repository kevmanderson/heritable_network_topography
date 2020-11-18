library(tidyverse)
library(foreach)
library(parallel)
library(aricode)
library(abind)
library(R.matlab)
library(data.table)
library(ragree)


dice_function = function(combo_idx, uniq_combos, neigh_labels){

    #if ((combo_idx %% 10000) == 0){
    #    print(combo_idx)
    #}

    # subject A / subject B
    a = uniq_combos[1,combo_idx]
    b = uniq_combos[2,combo_idx]

    # network labels for subjects A/B
    x = t(neigh_labels[a,])
    y = t(neigh_labels[b,])

    # allocate matrix nrow=n_vertices; ncol=n_networks
    n_vert    = length(x)
    arr       = 1:n_vert
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
    dice_df$idx  = combo_idx

    return(dice_df)
}




overlap_function = function(idx, just_labels, subj_arr, uniq_combos, hemi, radius=10, rerun_idx=FALSE){

    print(idx)
    vert_num = idx

    # directory to place dice coefficient matrices
    out_dir   = '/gpfs/milgram/scratch/holmes/topo_herit/data/vert_dice'

    # directory with vertex-wise ROI info
    roi_dir   = '/gpfs/milgram/project/holmes/kma52/topo_herit/data/HCP/roi_dir'

    if (hemi == 'L'){
        hemi_lower  = 'lh'
        hemi_labels = just_labels#[,grep('lh_', colnames(just_labels))]
        vertex = paste0('lh_V', as.character(vert_num))

    } else if (hemi == 'R'){
        hemi_lower  = 'rh'
        hemi_labels = just_labels#[,grep('rh_', colnames(just_labels))]
        vertex = paste0('rh_V', as.character(vert_num))

    }

    dice_write   = paste0(out_dir, '/local_dice_', hemi, '_vertex_', as.character(vert_num), '_idxfrom0_radius',as.character(radius),'.csv.gz')
    matrix_write = gsub('.csv.gz', '_matrix.csv', dice_write)
    if (file.exists(matrix_write)){
        return(NA)
    }

    neighborhood_verts = paste0(roi_dir, '/adjacency_', as.character(hemi), '_vertex', as.character(vert_num), '_idxfrom0.txt')
    if (file.exists(neighborhood_verts) == FALSE){
        return(NA)
    }
    neighborhood_df = read.csv(neighborhood_verts)
    neighborhood_df = neighborhood_df[neighborhood_df$distance <= radius,]

    # offset by 1 to match the column names
    neigh_cols   = paste0(hemi_lower, '_V', as.character(neighborhood_df$vert))
    neigh_cols   = neigh_cols[neigh_cols %in% colnames(hemi_labels)]
    neigh_labels = just_labels[neigh_cols]
    print(dim(neigh_labels))

    write('Calculating vertex dice...','')
    dice_list  = mclapply(1:ncol(uniq_combos), dice_function, uniq_combos=uniq_combos, neigh_labels=neigh_labels, mc.cores = 20)
    write('Dice Done!','')

    #dice_df    = do.call('rbind', dice_list)
    dice_dt    = rbindlist(dice_list)
    dice_dt$subA_id = subj_arr[dice_dt$subA]
    dice_dt$subB_id = subj_arr[dice_dt$subB]

    data.table::fwrite(dice_dt, file = dice_write)
    dice_df    = as.data.frame(dice_dt)


    # put into matrix format
    dice_mat = matrix(NA, max(dice_df$subB), max(dice_df$subB))
    dice_mat[cbind(dice_df$subA, dice_df$subB)] = dice_df$overall_dice
    dice_mat[cbind(dice_df$subB, dice_df$subA)] = dice_df$overall_dice
    diag(dice_mat) = 1

    # add subject info
    colnames(dice_mat) = subj_arr
    rownames(dice_mat) = subj_arr

    # order by subject number
    dice_mat_order = dice_mat[order(colnames(dice_mat)), order(colnames(dice_mat))]
    matrix_write   = gsub('.csv.gz', '_matrix.csv', dice_write)
    write.table(dice_mat_order, matrix_write, row.names=F, col.names=F, quote=F, sep=',')
}




#run_chunk   = 15


# process command line arguments
# --------------
args = commandArgs(trailingOnly=TRUE)
run_chunk = as.numeric(args[1])
num_chunk = as.numeric(args[2])
hemi      = as.character(args[3])



# read data
# --------------
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
mat_dat  = readMat(paste0(base_dir, '/data/HCP/HCP_S1200_1029sub_17net_parcellation.mat'))

# phenotype df
hcp_df = read_csv(paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))


# subj x subj combos for computation of dice
n_subs          = nrow(hcp_df)
uniq_combos     = combn(1:n_subs, 2)

# LH label information
lh_labels    = mat_dat$lh.labels
lh_labels_df = as.data.frame(t(lh_labels))
colnames(lh_labels_df) = paste0('lh_V', 0:(ncol(lh_labels_df)-1))
lh_labels_df$subject   = unlist(mat_dat$subject.list)


# RH label information
rh_labels    = mat_dat$rh.labels
rh_labels_df = as.data.frame(t(rh_labels))
colnames(rh_labels_df) =paste0('rh_V', 0:(ncol(rh_labels_df)-1))
rh_labels_df$subject   = unlist(mat_dat$subject.list)


# subset network labels to match those in HCP phenotype
lh_labels_df = lh_labels_df[lh_labels_df$subject %in% hcp_df$id,]
rh_labels_df = rh_labels_df[rh_labels_df$subject %in% hcp_df$id,]


# put dataframes in order
hcp_df       = hcp_df[order(hcp_df$id),]
lh_labels_df = lh_labels_df[order(lh_labels_df$subject),]
rh_labels_df = rh_labels_df[order(rh_labels_df$subject),]
which(lh_labels_df$subject != rh_labels_df$subject)


if (hemi == 'R'){
    # identify non midline surface vertices
    rh_vert_cols = colnames(rh_labels_df)[grep('rh', colnames(rh_labels_df))]
    rh_vert_sums = colSums(rh_labels_df[rh_vert_cols])

    # get rid of midlines
    rh_use_df = rh_labels_df[,which(rh_vert_sums != 0)]

    # compute the subject by subject dice in parallel
    chunk2 = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

    n_vert      = ncol(rh_use_df)
    vert_array  = as.numeric(gsub('rh_V', '', colnames(rh_use_df)))
    vert_chunks = chunk2(vert_array, num_chunk)
    run_idxs    = vert_chunks[[run_chunk]]
    lapply(run_idxs, overlap_function, hemi='R', subj_arr=rh_labels_df$subject, uniq_combos=uniq_combos, just_labels=rh_use_df)
}


if (hemi == 'L'){

    # identify non midline surface vertices
    lh_vert_cols = colnames(lh_labels_df)[grep('lh', colnames(lh_labels_df))]
    lh_vert_sums = colSums(lh_labels_df[lh_vert_cols])

    # get rid of midlines
    lh_use_df = lh_labels_df[,which(lh_vert_sums != 0)]

    # compute the subject by subject dice in parallel
    chunk2 = function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

    n_vert      = ncol(lh_use_df)
    vert_array  = as.numeric(gsub('lh_V', '', colnames(lh_use_df)))
    vert_chunks = chunk2(vert_array, num_chunk)
    run_idxs    = vert_chunks[[run_chunk]]
    lapply(run_idxs, overlap_function, hemi='L', subj_arr=lh_labels_df$subject, uniq_combos=uniq_combos, just_labels=lh_use_df)

}














