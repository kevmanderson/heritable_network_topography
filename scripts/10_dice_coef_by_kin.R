library(tidyverse)
library(solarius)
library(R.matlab)
library(Cairo)


# This script will:
#   1. Plot Dice for an example vertex (Figure 4), split by MZ/DZ/SIB/UNR

####################
# set up directories
####################
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
out_dir  = '/gpfs/milgram/scratch/holmes/topo_herit/data/vert_dice'


# read dice of vert to plot
lh_vert_dice = read_csv(paste0(out_dir, '/local_dice_L_vertex_13226_idxfrom0_radius10_matrix.csv'), col_names=F)


# read individual parcellations
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
mat_dat  = readMat(paste0(base_dir, '/data/HCP/HCP_S1200_1029sub_17net_parcellation.mat'))


# read phenotype df
hcp_df = read_csv(paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))
hcp_df_order = hcp_df[order(hcp_df$id),]
colnames(lh_vert_dice) = hcp_df_order$id



###########################
# create HCP kinship matrix
###########################
kinship_big  = solarKinship2(hcp_df)

# make sure kinship matrix is in right order, again
kin_match_col = match(hcp_df_order$id, colnames(kinship_big))
kin_match_row = match(hcp_df_order$id, rownames(kinship_big))
kin_matrix    = kinship_big[kin_match_col, kin_match_row]
kin_matrix    = as.data.frame(kin_matrix)


# check counts
table(hcp_df_order$rel_group)
table(kin_matrix[upper.tri(kin_matrix)])


##############
# df of MZ IDs
##############
mz_twins = as.data.frame(which(kin_matrix == 1, arr.ind=T))
mz_twins = mz_twins[mz_twins$row != mz_twins$col,]
mz_twins$row_id = colnames(kin_matrix)[mz_twins$row]
mz_twins$col_id = colnames(kin_matrix)[mz_twins$col]
mz_twins$uniq_id = paste0(mz_twins$row_id, '_', mz_twins$col_id)
mz_twins$group = 'MZ'


##############
# df of DZ IDs
##############
sib_dz = as.data.frame(which(kin_matrix == 0.5, arr.ind=T))
sib_dz = sib_dz[sib_dz$row != sib_dz$col,]
sib_dz$row_id = colnames(kin_matrix)[sib_dz$row]
sib_dz$col_id = colnames(kin_matrix)[sib_dz$col]
sib_dz$uniq_id = paste0(sib_dz$row_id, '_', sib_dz$col_id)


##############
# df of SIB IDs
##############
sib_dz_df = merge(x=sib_dz, y=hcp_df_order[c('id','rel_group')], by.x='row_id', by.y='id')
sib_dz_df = merge(x=sib_dz_df, y=hcp_df_order[c('id','rel_group')], by.x='col_id', by.y='id')
sib_dz_df$group = ifelse(sib_dz_df$rel_group.x == 'DZ' & sib_dz_df$rel_group.y == 'DZ', 'DZ', 'SIB')
sib_dz2   = sib_dz_df
sib_dz2$uniq_id = paste0(sib_dz2$col_id, '_', sib_dz2$row_id)
sib_dz_df = rbind(sib_dz_df, sib_dz2)


##############
# df of UNR IDs
##############
unrelated = as.data.frame(which(kin_matrix == 0, arr.ind=T))
unrelated = unrelated[unrelated$row != unrelated$col,]
unrelated$row_id  = colnames(kin_matrix)[unrelated$row]
unrelated$col_id  = colnames(kin_matrix)[unrelated$col]
unrelated$uniq_id = paste0(unrelated$row_id, '_', unrelated$col_id)
unrelated$group = 'unrelated'


##########
# combine
##########
relatedness_df = rbind(mz_twins, sib_dz_df[colnames(mz_twins)], unrelated[colnames(mz_twins)])

# merge with dice coefficients
dice_in = as.data.frame(lh_vert_dice) #read.csv(paste0(base_dir, '/data/HCP/hcp_', hemi, '_net_',dice_net,'_matrix.csv'))

#upper_ids = which(upper.tri(dice_in))
dice_long_df = as.data.frame(which(upper.tri(dice_in), arr.ind=T))
dice_long_df$dice = dice_in[upper.tri(dice_in)]
dice_long_df$row_id = colnames(dice_in)[dice_long_df$row]
dice_long_df$col_id = colnames(dice_in)[dice_long_df$col]
dice_long_df$uniq_id = paste0(dice_long_df$row_id, '_', dice_long_df$col_id)

dice_related_df = merge(x=dice_long_df, y=relatedness_df, by='uniq_id')
dice_related_df$hemi = hemi
dice_related_df$net = dice_net
dice_related_df$hemi_net = paste0(hemi, '_', dice_net)
dice_related_df %>% group_by(group) %>% summarise(dice=mean(dice))

plot_df = dice_related_df
plot_df$group = factor(plot_df$group, levels=c('MZ','DZ','SIB','unrelated'))


##########
# plot
##########
plot_me = ggplot(data=plot_df, aes(x=group, y=dice, fill=group)) +
                    geom_boxplot(outlier.size=0.1) +
                    theme_classic() +
                    scale_y_continuous(limits=c(.2, 1), expand=c(0,0), breaks=seq(.2,1,.2)) +
                    scale_fill_manual(values=c('#464E57', '#79868F', '#B6BBC7', '#CED0D6')) +
                    scale_color_manual(values=c('#000000','#000000')) +
                    theme(legend.position='none', axis.text=element_text(color='black'),
                            axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                            axis.ticks.y=element_line(color='black'),
                            axis.title.x=element_blank())


CairoPDF('/gpfs/milgram/project/holmes/kma52/topo_herit/figures/vert13226_dice_by_group.pdf', width=2, height=2)
print(plot_me)
dev.off()

































