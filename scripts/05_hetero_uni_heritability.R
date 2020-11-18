library(tidyverse)
library(gifti)
library(R.matlab)
library(solarius)
library(ggpubr)
library(Cairo)


# This script will:
#   1. Boxplot/dot plot of SA heritability, split by hetero/unimodal networks


# network to color dataframe
hex_df = data.frame(VisualA='#781286',
                    VisualB='#ff0101',
                    VisualC='#7a8733',
                    SomatomotorA='#4682b2',
                    SomatomotorB='#2bcca2',
                    Auditory='#dcf8a5',
                    DorsalAttentionA='#4a9b3d',
                    DorsalAttentionB='#007610',
                    VentralAttentionA='#c43afb',
                    VentralAttentionB='#ff98d6',
                    ControlA='#e69423',
                    ControlB='#87324c',
                    ControlC='#778cb1',
                    DefaultA='#ffff02',
                    DefaultB='#cd3e50',
                    DefaultC='#000083',
                    TemporalParietal='#0929fa')

# read data
# -----------
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
hcp_df   = read_csv(paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))

# add vars
hcp_df$FS_IntraCranial_Vol_scale = as.numeric(scale(hcp_df$FS_IntraCranial_Vol))
hcp_df$age2 = hcp_df$Age_in_Yrs^2

# read network name information
labels    = readMat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = unlist(labels$network.name)


# read heritability of net surf area
# -----------
relSA_h2_df     = read_csv(paste0(base_dir, '/data/HCP/solar_surf_area_h2_by_net_rr.csv'))
relSA_h2_df$net = gsub('lh_|rh_|_sa', '', relSA_h2_df$trait)


# unimodal versus heteromodal
relSA_h2_df$modality = NA
relSA_h2_df$modality[which(grepl('Aud|Vis|Som', relSA_h2_df$net))] = 'unimodal'
relSA_h2_df$modality[which(!grepl('Aud|Vis|Som', relSA_h2_df$net))] = 'heteromodal'

table(relSA_h2_df$net, relSA_h2_df$modality)

relSA_h2_df %>% group_by(modality) %>% summarise(h2=mean(Var), h2_std=sd(Var))
summary(lm('Var~modality', relSA_h2_df))

# plot SA heritability, grouped by uni/hetero modal
boxplot = ggplot(data=relSA_h2_df, aes(x=modality, y=Var, color = modality)) +
                    geom_point() +
                    geom_boxplot(alpha=0.8) +
                    scale_color_manual(values=c('#82CFFD', '#00688B')) +
                    theme_classic() +
                    scale_y_continuous(limits=c(0, .6), breaks=seq(0, .6, .2), expand=c(0,0)) +
                    theme(axis.title.x=element_blank(),
                            text=element_text(size=10,  family="Arial", color='black'),
                            axis.text=element_text(size=10,  family="Arial", color='black'),
                            plot.title = element_text(hjust = 0.5),
                            legend.position='none') +
                    labs(title="Heritability", y = 'Heritability of Relative Surface Area')


hist_file = paste0(base_dir, '/figures/hcp_uni_hetero_modal_boxplot_rr.pdf')
CairoPDF(hist_file, width=1.75, height=1.75, family='Arial')
print(boxplot)
dev.off()
hist_file


# plot version 2
hex_df_t    = data.frame(color=t(hex_df), net=colnames(hex_df))
plot_df     = merge(x=relSA_h2_df, y=hex_df_t, all.x=T, by.x='net', by.y='net')
plot_df$net = factor(plot_df$net, levels=colnames(hex_df))
plot_df     = plot_df[order(plot_df$net),]

boxplot = ggplot(data=plot_df, aes(x=modality, y=Var, group = modality)) +
                    geom_point(aes(color=net, fill=net)) +
                    geom_boxplot(alpha=0.8) +
                    #scale_color_manual(values=c('#82CFFD', '#00688B')) +
                    theme_classic() +
                    scale_y_continuous(limits=c(0, .6), breaks=seq(0, .6, .2), expand=c(0,0)) +
                    theme(axis.title.x=element_blank(),
                            text=element_text(size=10,  family="Arial", color='black'),
                            axis.text=element_text(size=10,  family="Arial", color='black'),
                            plot.title = element_text(hjust = 0.5),
                            legend.position='none') +
                    labs(title="Heritability", y = 'Heritability of Relative Surface Area') +
                    scale_color_manual(values = unique(as.character(plot_df$color)))

hist_file = paste0(base_dir, '/figures/hcp_uni_hetero_modal_boxplot_v2_rr.pdf')
CairoPDF(hist_file, width=1.75, height=1.75, family='Arial')
print(boxplot)
dev.off()
hist_file















