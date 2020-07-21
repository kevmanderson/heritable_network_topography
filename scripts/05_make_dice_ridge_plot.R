library(tidyverse)
library(solarius)
library(R.matlab)
library(ggridges)
library(Cairo)
library(corrplot)
library(ggridges)

# This script will:
#   1. Create ridge plot of inter-individual dice coefficients, goruped by network and hemisphere
#   2. Create a boxplot of dice coefficients, grouped by MZ/DZ/SIB/UNRELATED pairs


# ------------------------
# read HCP data/label info
# ------------------------
base_dir  = '/gpfs/milgram/project/holmes/kma52/topo_herit'
topo_dir  = paste0(base_dir, '/data/topology_heritability')
hcp_df    = read.csv(paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))

labels    = readMat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = unlist(labels$network.name)


# ------------------------------------------
# read dice matrices
# created by "01_hcp_dice_matrices.R" script
# ------------------------------------------
topo_dir = paste0(base_dir, '/data/topology_heritability')
dice_df  = NULL
for (net_num in 1:17){
    write(net_num,'')
    net_name = net_names[net_num]

    # left hemi dice
    lh_dice_path = paste0(topo_dir, '/lh_dice_net_',as.character(net_num),'_P.csv')
    lh_dice      = read.csv(lh_dice_path, header=F)
    lh_df        = data.frame(dice=lh_dice[upper.tri(lh_dice)], net_name=net_name, net_num=net_num, hemi='lh')

    # right hemi dice
    rh_dice_path = paste0(topo_dir, '/rh_dice_net_',as.character(net_num),'_P.csv')
    rh_dice      = read.csv(rh_dice_path, header=F)
    rh_df        = data.frame(dice=rh_dice[upper.tri(rh_dice)], net_name=net_name, net_num=net_num, hemi='rh')

    fill_row = data.frame(dice=NA, net_name=net_name, net_num=net_num, hemi='fill')

    # bihemi
    bihemi_df = rbind(lh_df, rh_df, fill_row)
    dice_df   = rbind(dice_df, bihemi_df)
}
dice_df     = dice_df %>% filter(hemi %in% c('lh','rh'))
dice_df$net = paste0(dice_df$hemi, '_', dice_df$net_num, '_relative_SA')
dice_avg_df = dice_df %>% group_by(net) %>% summarise(dice_mean=mean(dice))
dice_avg_df = dice_avg_df[dice_avg_df$dice_mean < 1,]


# plot unimodal before heteromodal networks
plot_order = c('VisualA','VisualB','VisualC','SomatomotorA','SomatomotorB','Auditory',
                'VentralAttentionA', 'VentralAttentionB', 'DorsalAttentionA', 'DorsalAttentionB',
                'ControlA','ControlB','ControlC',
                'DefaultA','DefaultB','DefaultC', 'TemporalParietal')

# adjust the plot_order variable to account for the "fill" rows
plot_full_order = NULL
for (net in plot_order){
    plot_full_order = c(plot_full_order, paste0(net, '_',  c('lh','rh','fill')))
}
dice_df$plot_nets = paste0(dice_df$net_name, '_', dice_df$hemi)
dice_df$plot_nets = factor(dice_df$plot_nets, levels=plot_full_order)

# plot variance - separate across hemispheres
ridge_plot = ggplot(data=dice_df, aes(x = dice, y = plot_nets)) +
                        geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3) +
                        scale_fill_gradientn(limits = c(0, 1),
                        colours = c("#C22A39",'#F7B494', "#F8F9FA", '#A2CFE5', "#297CBA"), name = 'dice') +
                        labs(title = paste0('HCP Dice')) +
                        theme_classic() +
                        scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
                        theme(axis.title.y=element_blank(),
                        axis.ticks = element_blank(),
                        text=element_text(size=10,  family="Arial", color='black'),
                        axis.text=element_text(size=10,  family="Arial", color='black'),
                        plot.title = element_text(hjust = 0.5),
                        legend.title=element_blank()) +
                        labs(title="HCP Dice Coefficients", x = 'Dice Coefficient')

# save plot
hist_file = paste0(base_dir, '/figures/hcp_dice_overlaps_ggridge_split_hemi.pdf')
CairoPDF(hist_file, width=4, height=4, family='Arial')
print(ridge_plot)
dev.off()
hist_file


# calculate average dice for each network/hemi
dice_avg_df = dice_df %>% group_by(net_name, hemi) %>% summarize(dice_mean=mean(dice), dice_var=sd(dice, na.rm=T))

# split networks into heteromodal/unimodal
dice_avg_df$id = paste0(dice_avg_df$net_name, '_', dice_avg_df$hemi)
dice_avg_df$modality = ifelse(grepl('Aud|Som|Vis', dice_avg_df$net_name), 'uni', 'het')

# linear model for het/uni difference in Dice
summary(lm(dice_mean ~ modality, data=dice_avg_df))
dice_avg_df %>% group_by(modality) %>% summarise(m_dice=mean(dice), sd_dice=sd(dice))




# ------------------------------------
# Create boxplot of MZ/DZ/SIB/UNR Dice
# ------------------------------------

# create the kinship matrix for HCP
hcp_df_order = hcp_df[order(hcp_df$id),]
kinship_big  = solarKinship2(hcp_df_order)

# make sure kinship matrix is in right order, again
kin_match_col = match(hcp_df_order$id, colnames(kinship_big))
kin_match_row = match(hcp_df_order$id, rownames(kinship_big))
kin_matrix    = kinship_big[kin_match_col, kin_match_row]
kin_matrix    = as.data.frame(kin_matrix)
#kin_matrix[lower.tri(kin_matrix)] = NA

table(hcp_df_order$rel_group)
table(kin_matrix[upper.tri(kin_matrix)])

# identify MZ pairs
mz_twins = as.data.frame(which(kin_matrix == 1, arr.ind=T))
mz_twins = mz_twins[mz_twins$row != mz_twins$col,]
mz_twins$row_id = colnames(kin_matrix)[mz_twins$row]
mz_twins$col_id = colnames(kin_matrix)[mz_twins$col]
mz_twins$uniq_id = paste0(mz_twins$row_id, '_', mz_twins$col_id)
mz_twins$group = 'MZ'

# identify DZ/sibling pairs
sib_dz = as.data.frame(which(kin_matrix == 0.5, arr.ind=T))
sib_dz = sib_dz[sib_dz$row != sib_dz$col,]
sib_dz$row_id = colnames(kin_matrix)[sib_dz$row]
sib_dz$col_id = colnames(kin_matrix)[sib_dz$col]
sib_dz$uniq_id = paste0(sib_dz$row_id, '_', sib_dz$col_id)

# merge MZ/DZ+sib DFs together
sib_dz_df = merge(x=sib_dz, y=hcp_df_order[c('id','rel_group')], by.x='row_id', by.y='id')
sib_dz_df = merge(x=sib_dz_df, y=hcp_df_order[c('id','rel_group')], by.x='col_id', by.y='id')
sib_dz_df$group = ifelse(sib_dz_df$rel_group.x == 'DZ' & sib_dz_df$rel_group.y == 'DZ', 'DZ', 'SIB')

# identify unrelated pairs
sib_dz2 = sib_dz_df
sib_dz2$uniq_id = paste0(sib_dz2$col_id, '_', sib_dz2$row_id)
sib_dz_df = rbind(sib_dz_df, sib_dz2)

unrelated = as.data.frame(which(kin_matrix == 0, arr.ind=T))
unrelated = unrelated[unrelated$row != unrelated$col,]
unrelated$row_id  = colnames(kin_matrix)[unrelated$row]
unrelated$col_id  = colnames(kin_matrix)[unrelated$col]
unrelated$uniq_id = paste0(unrelated$row_id, '_', unrelated$col_id)
unrelated$group = 'unrelated'

relatedness_df = rbind(mz_twins, sib_dz_df[colnames(mz_twins)], unrelated[colnames(mz_twins)])


# reformat to dataframe
dice_data_list = NULL
dice_nets = c(paste0('dice_net_',1:17), 'overall_dice')
dice_net  = dice_nets[1]


# read and merge dice coefficients for each network/hemi
for (dice_net in dice_nets){
    dice_bihemi_df = NULL
    write(dice_net,'')
    for (hemi in c('lh', 'rh')){
        dice_in = read.csv(paste0(base_dir, '/data/HCP/hcp_', hemi, '_net_',dice_net,'_matrix.csv'))
        colnames(dice_in) = gsub('X', '', colnames(dice_in))
        rownames(dice_in) = gsub('X', '', rownames(dice_in))

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

        dice_bihemi_df = rbind(dice_bihemi_df, dice_related_df)
    }
    dice_data_list[[dice_net]] = dice_bihemi_df
}


plot_df = dice_data_list[['overall_dice']]
plot_df$group = factor(plot_df$group, levels=c('MZ','DZ','SIB','unrelated'))

plot_df %>% group_by(group, hemi) %>% summarise(m_dice=mean(dice), sd_dice=sd(dice))


plot_me = ggplot(data=plot_df, aes(x=group, y=dice, fill=group, color=hemi)) +
                    geom_boxplot(outlier.size=0.1) +
                    theme_classic() +
                    scale_y_continuous(limits=c(.45,.9), expand=c(0,0), breaks=seq(0.45,.9,.15)) +
                    scale_fill_manual(values=c('#464E57', '#79868F', '#B6BBC7', '#CED0D6')) +
                    scale_color_manual(values=c('#000000','#000000')) +
                    theme(legend.position='none', axis.text=element_text(color='black'),
                            axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                            axis.ticks.y=element_line(color='black'),
                            axis.title.x=element_blank())

CairoPDF(paste0(base_dir, '/figures/overall_net_dice_by_group.pdf'), height=1.75, width=1.75)
print(plot_me)
dev.off()














