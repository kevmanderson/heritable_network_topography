library(tidyverse)
library(gifti)
library(R.matlab)
library(solarius)
library(ggpubr)
library(Cairo)


# This script will:
#   1. Boxplot/dot plot of SA heritability, split by hetero/unimodal networks



# ---------
# read data
# ---------
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'
hcp_df   = read_csv(paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))

# add vars
hcp_df$FS_IntraCranial_Vol_scale = as.numeric(scale(hcp_df$FS_IntraCranial_Vol))
hcp_df$age2 = hcp_df$Age_in_Yrs^2

# read network name information
labels    = readMat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = unlist(labels$network.name)



# ----------------------------------
# read heritability of net surf area
# ----------------------------------
relSA_h2_df     = read_csv(paste0(base_dir, '/data/HCP/solar_relativeSA_h2_byNet.csv'))
relSA_h2_df$net = net_names[as.numeric(gsub('lh_|rh_|_relative_SA', '', relSA_h2_df$trait))]

# unimodal versus heteromodal
relSA_h2_df$modality = NA
relSA_h2_df$modality[which(grepl('Aud|Vis|Som', relSA_h2_df$net))] = 'unimodal'
relSA_h2_df$modality[which(!grepl('Aud|Vis|Som', relSA_h2_df$net))] = 'heteromodal'



relSA_h2_df %>% group_by(modality) %>% summarise(h2=mean(Var), h2_std=sd(Var))
summary(lm('Var~modality', relSA_h2_df))

# plot SA heritability, grouped by uni/hetero modal
boxplot = ggplot(data=relSA_h2_df, aes(x=modality, y=Var, color = modality)) +
                    geom_point() +
                    geom_boxplot(alpha=0.8) +
                    scale_color_manual(values=c('#82CFFD', '#00688B')) +
                    theme_classic() +
                    scale_y_continuous(limits=c(.2, .6), breaks=seq(.2, .6, .1), expand=c(0,0)) +
                    theme(axis.title.x=element_blank(),
                            text=element_text(size=10,  family="Arial", color='black'),
                            axis.text=element_text(size=10,  family="Arial", color='black'),
                            plot.title = element_text(hjust = 0.5),
                            legend.position='none') +
                    labs(title="Heritability", y = 'Heritability of Relative Surface Area')


hist_file = paste0(base_dir, '/figures/hcp_uni_hetero_modal_boxplot.pdf')
CairoPDF(hist_file, width=1.75, height=1.75, family='Arial')
print(boxplot)
dev.off()
hist_file





