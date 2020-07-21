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




covars = c('Age_in_Yrs', 'age2', 'age2*sex', 'Age_in_Yrs*sex', 'Height', 'BMI', 'FS_IntraCranial_Vol_scale')


# unimodal versus heteromodal
unimodal_nums = which(grepl('Aud|Vis|Som', net_names))
hetero_nums   = which(!grepl('Aud|Vis|Som', net_names))

# uni
hcp_df$lh_Unimodal_sum_SA = rowSums(hcp_df[paste0('lh_',unimodal_nums,'_sum_SA')])
hcp_df$rh_Unimodal_sum_SA = rowSums(hcp_df[paste0('rh_',unimodal_nums,'_sum_SA')])


# hetero
hcp_df$lh_Heteromodal_sum_SA = rowSums(hcp_df[paste0('lh_',hetero_nums,'_sum_SA')])
hcp_df$rh_Heteromodal_sum_SA = rowSums(hcp_df[paste0('rh_',hetero_nums,'_sum_SA')])

hcp_df$lh_Heteromodal_sum_SA_resid =  resid(lm('lh_Heteromodal_sum_SA~lh_total_SA', data=hcp_df))
hcp_df$rh_Heteromodal_sum_SA_resid =  resid(lm('rh_Heteromodal_sum_SA~lh_total_SA', data=hcp_df))

hcp_df$lh_Unimodal_sum_SA_resid =  resid(lm('lh_Unimodal_sum_SA~lh_total_SA', data=hcp_df))
hcp_df$rh_Unimodal_sum_SA_resid =  resid(lm('rh_Unimodal_sum_SA~lh_total_SA', data=hcp_df))



# calculate heritability of heteromodal versus unimodal cortex
covar_formula = paste0(covars, collapse=' + ')
df            = as.data.frame(hcp_df)
trait_arr     = c('lh_Unimodal_sum_SA', 'rh_Unimodal_sum_SA', 'lh_Heteromodal_sum_SA', 'rh_Heteromodal_sum_SA')
trait_arr     = c('lh_Unimodal_sum_SA_resid', 'rh_Unimodal_sum_SA_resid', 'lh_Heteromodal_sum_SA_resid', 'rh_Heteromodal_sum_SA_resid')
h2_hetUni_df  = NULL
for (trait in trait_arr){
    write(trait,'')
    df[[trait]] = scale(df[trait])
    formula   = as.formula(paste0(trait, ' ~ ', covar_formula))
    rhog      = solarPolygenic(formula, df)

    out_row       = rhog$vcf[1,]
    out_row$trait = trait
    h2_hetUni_df  = rbind(h2_hetUni_df, out_row)
}
h2_hetUni_df$hemi  = 'lh'
h2_hetUni_df$hemi[grep('rh', h2_hetUni_df$trait)] = 'rh'

h2_hetUni_df$split = 'uni'
h2_hetUni_df$split[grep('Het', h2_hetUni_df$trait)] = 'het'



# Mean SA - plot het/uni barplot
plot_me   = h2_hetUni_df[grep('mean', h2_hetUni_df$trait),]
plot_mean = ggplot(data=plot_me, aes(x=split, y=Var, fill=hemi)) +
                geom_bar(stat='identity', position=position_dodge(width=0.9)) +
                geom_errorbar(aes(ymin = Var-SE, ymax = Var+SE), width = 0.2, position=position_dodge(width=0.9)) +
                theme_classic() +
                scale_y_continuous(limits = c(0, .8), breaks=seq(0,.8,.2), expand = c(0,0)) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) +
                        theme(text = element_text(size = 10, colour = "black"),
                            axis.text.y = element_text(colour = "black"),
                            axis.text.x = element_text(colour = "black"),
                        plot.title = element_text(hjust = 0.5)) +
                scale_fill_manual(values=c('#D9E9F2', '#0F4D92')) + ggtitle('Mean SA')

fig_out = paste0(base_dir, '/figures/HetUni_h2_meanSA_barplot.pdf')
ggsave(plot=plot_mean, filename=fig_out, height=1.5, width=2)

