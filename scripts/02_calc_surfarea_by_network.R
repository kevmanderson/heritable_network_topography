library(tidyverse)
library(gifti)
library(R.matlab)
library(solarius)


# This script will:
#   1. Calculate surface area for each network and each individual
#   2. Reference each subject's vertex-level Freesurfer surface area surfaces
#   3. Save data as "surface_area_by_net.csv"


summarise_area_by_net = function(anat_dir, sub, lh_lab, rh_lab){


    avg_df = NULL
    L_subarea_file = read_gifti(paste0(anat_dir, sub, '/', sub,'.L.midthickness_areas.32k_fs_LR.func.gii'))
    for (idx in 1:17){
        # mean
        area_mean = mean(L_subarea_file$data$normal[which(lh_lab==idx)])
        avg_df[paste0('lh_', as.character(idx), '_mean_SA')] = area_mean
        # sum
        area_sum = sum(L_subarea_file$data$normal[which(lh_lab==idx)])
        avg_df[paste0('lh_', as.character(idx), '_sum_SA')] <- area_sum
    }

    R_subarea_file = read_gifti(paste0(anat_dir, sub, '/', sub,'.R.midthickness_areas.32k_fs_LR.func.gii'))
    for (idx in 1:17){
        # mean
        area_mean = mean(R_subarea_file$data$normal[which(rh_lab==idx)])
        avg_df[paste0('rh_', as.character(idx), '_mean_SA')] = area_mean
        # sum
        area_sum = sum(R_subarea_file$data$normal[which(rh_lab==idx)])
        avg_df[paste0('rh_', as.character(idx), '_sum_SA')] = area_sum
    }

    bihemi_sa  = rbind(L_subarea_file$data$normal, R_subarea_file$data$normal)
    bihemi_lab = c(lh_lab, rh_lab)
    for (idx in 1:17){
        # mean
        area_mean = mean(bihemi_sa[which(bihemi_lab==idx)])
        avg_df[paste0('bihemi_', as.character(idx), '_mean_SA')] = area_mean
        # mean
        area_sum = sum(bihemi_sa[which(bihemi_lab==idx)])
        avg_df[paste0('bihemi_', as.character(idx), '_sum_SA')] = area_sum
    }

    avg_df           = as.data.frame(t(avg_df))
    avg_df$Subject   = sub
    rownames(avg_df) = NULL
    return(avg_df)
}



# ---------------------------
# HCP phenotype/pedigree data
# ---------------------------
base_dir  = '/gpfs/milgram/project/holmes/kma52/topo_herit'
hcp_pheno = read_csv(paste0(base_dir, '/data/HCP/hcp_for_solar_allvars.csv'))
hcp_pheno$Subject = hcp_pheno$id


# -----------------------------
# read individual parcellations
# -----------------------------
mat_dat      = readMat(paste0(base_dir, '/data/HCP/HCP_S1200_1029sub_17net_parcellation.mat'))
lh_labels    = mat_dat$lh.labels
lh_labels_df = as.data.frame(t(lh_labels))
rh_labels    = mat_dat$rh.labels
rh_labels_df = as.data.frame(t(rh_labels))


# add subject numbers to the parcellation label dataframes
subj_list = as.numeric(unlist(mat_dat$subject.list))
lh_labels_df$Subject = as.numeric(subj_list)
rh_labels_df$Subject = as.numeric(subj_list)


sub_df = data.frame(Subject=subj_list, order_idx=1:length(subj_list))

# individual parcellations for each twin (no medial wall)
subj_indiparc   = subj_list[subj_list %in% hcp_pheno$Subject]
lh_lab_indiparc = lh_labels[,subj_list %in% hcp_pheno$Subject]
rh_lab_indiparc = rh_labels[,subj_list %in% hcp_pheno$Subject]


# ----------------------------------------------------
# merge phenotype with indiv part subject df (n=1,023)
# ----------------------------------------------------
hcp_df = merge(x=sub_df, y=hcp_pheno, by='Subject')
hcp_df = hcp_df[order(hcp_df$order_idx),]


# double check subject order
which(hcp_df$Subject != subj_indiparc)


# ---------------------------------------------
# get total surface area in each parcel/subject
# ---------------------------------------------
anat_dir       = paste0(base_dir, '/data/HCP/anat/')
subj_net_areas = list()
for (i in 1:length(hcp_df$Subject)){
  write(i,'')
  sub     = hcp_df$Subject[i]
  cur_row = summarise_area_by_net(anat_dir=anat_dir, sub=sub, lh_lab=lh_lab_indiparc[,i], rh_lab=rh_lab_indiparc[,i])
  subj_net_areas[[i]] = cur_row
}
subj_net_area_df         = as_tibble(do.call(rbind, subj_net_areas))
subj_net_area_df$Subject = as.character(subj_net_area_df$Subject)


out_file = paste0(base_dir, '/data/HCP/surface_area_by_net.csv')
write_csv(subj_net_area_df, out_file)

