library(tidyverse)
library(solarius)


# This script will:
#   1. Create the necessary csv file for multi-dimensional h2 analyses using the "h2_mat.m" functions
#       a. Dice phenotype similarity matrix (for each net/hemi)
#       b. Kinship matrix ("K.csv"), MZ=1, DZ/Sib=0.5, UNR=0.
#       c. Covariate dataframe ("covar.csv")
#       d. Family IDs ("F.csv"), necessary for bootstrapping error bars


# ---------
# read data
# ---------
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
hcp_df   = read_csv(paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))

dice_names = c(paste0('dice_net_',1:17), 'overall_dice')


# dice overlap, subjxsubj matrix
dice_name = dice_names[18]

for (hemi in c('lh','rh','bihemi')){
    for (dice_name in dice_names[1]){

        write(paste0(dice_name,'_',hemi),'')

        # sub by sub network similarity
        net_sim  = read_csv(paste0(base_dir, '/data/HCP/hcp_',hemi,'_net_',dice_name,'_matrix.csv'))
        net_sim  = net_sim[colnames(net_sim) %in% hcp_df$id, colnames(net_sim) %in% hcp_df$id]

        # make sure they are in the right order
        hcp_df_order = hcp_df[match(hcp_df$id, colnames(net_sim)),]
        which(hcp_df_order$id != colnames(net_sim))

        # create the kinship matrix for HCP
        kinship_big   = solarKinship2(hcp_df_order)

        # make sure kinship matrix is in right order, again
        kin_match_col = match(hcp_df_order$id, colnames(kinship_big))
        kin_match_row = match(hcp_df_order$id, rownames(kinship_big))
        kin_matrix    = kinship_big[kin_match_col, kin_match_row]
        kin_matrix    = as.data.frame(kin_matrix)


        # should be zero
        which(colnames(kin_matrix) != hcp_df_order$id)
        which(rownames(kin_matrix) != hcp_df_order$id)

        # make sure net similarity matrix is in right order, again
        net_match_col = match(hcp_df_order$id, colnames(net_sim))
        net_match_row = match(hcp_df_order$id, colnames(net_sim))
        net_matrix    = net_sim[net_match_col, net_match_row]
        net_matrix    = as.data.frame(net_matrix)

        # check order
        which(colnames(net_matrix) != hcp_df_order$id)
        which(colnames(net_matrix) != colnames(kin_matrix))
        which(colnames(net_matrix) != rownames(kin_matrix))


        # write data for input into Tian's heritability code
        out_dir    = paste0(base_dir, '/data/topology_heritability')

        # test with height
        #height_crossprod = tcrossprod(hcp_df_order$Height)
        #write_csv(x=as.tibble(height_crossprod), paste0(out_dir, '/height_P.csv'), col_names=F)

        write_csv(x=net_matrix, paste0(out_dir, '/',hemi,'_',dice_name,'_P.csv'), col_names=F)
    }
}

# write kinship matrix
write_csv(x=kin_matrix, paste0(out_dir, '/K.csv'))


write_csv(x=hcp_df_order['famid'], paste0(out_dir, '/F.csv'))


# format covariates
hcp_df_order$sex_bin        = ifelse(hcp_df_order$sex=='M',1,0)
hcp_df_order$Height[is.na(hcp_df_order$Height)] = mean(hcp_df_order$Height, na.rm=T)
hcp_df_order$Weight[is.na(hcp_df_order$Weight)] = mean(hcp_df_order$Weight, na.rm=T)
hcp_df_order$BMI[is.na(hcp_df_order$BMI)] = mean(hcp_df_order$BMI, na.rm=T)
hcp_df_order$age_2          = hcp_df_order$Age_in_Yrs^2
hcp_df_order$Age_in_Yrs_sex = hcp_df_order$Age_in_Yrs * hcp_df_order$sex_bin
hcp_df_order$age_2_sex      = hcp_df_order$age_2 * hcp_df_order$sex_bin
hcp_df_order$Ethnicity_bin = ifelse(grepl('Not', hcp_df_order$Ethnicity), 1, 0)


# write covariates
covars    = c('Age_in_Yrs', 'age2', 'age_2_sex', 'Height', 'Ethnicity_bin', 'Age_in_Yrs_sex', 'BMI', 'FS_IntraCranial_Vol_scale')
covar_mat = hcp_df_order[covars]
write_csv(x=covar_mat, paste0(out_dir, '/covar.csv'))


# height covariance
#height_pheno = hcp_df_order$Height %*% t(hcp_df_order$Height)
#write_csv(x=as.data.frame(height_cov), paste0(out_dir, '/Height_P.csv'), col_names=F)


# write covariates
#height_covars     = c('Age_in_Yrs', 'age2', 'age_2_sex', 'Ethnicity_bin', 'Age_in_Yrs_sex', 'FS_IntraCranial_Vol_scale')
#height_covars_mat = hcp_df_order[height_covars]
#write_csv(x=height_covars_mat, paste0(out_dir, '/Height_covar.csv'))










