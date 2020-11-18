library(tidyverse)
library(gifti)
library(R.matlab)
library(solarius)
library(ggpubr)
library(Cairo)


# This script will:
#   1. Calculate relative surface area, by dividing indiv. network areas by total indiv. cortical surface area
#   2. SOLAR-based heritability of relative network size
#       save as "solar_relativeSA_h2_byNet.csv"
#   3. Plot LH/RH correlation of heritability estimates
#   4. Barplots of heritability, grouped by network and hemisphere


plot_hemi_h2_corr = function(h2_df, string, limits){
    #lh_traits = paste0('lh_', as.character(1:17), '_',  string)
    #rh_traits = paste0('rh_', as.character(1:17), '_',  string)

    lh_traits = h2_df$trait[grepl('lh', h2_df$trait)]
    rh_traits = h2_df$trait[grepl('rh', h2_df$trait)]

    lh = h2_df[h2_df$trait %in% lh_traits,]
    rh = h2_df[h2_df$trait %in% rh_traits,]
    hemi_cor = cor.test(lh$Var, rh$Var)

    hemi_cor_sp = cor.test(lh$Var, rh$Var, method='spearman')

    print(hemi_cor)
    print(hemi_cor_sp)

    colnames(lh) = paste0('lh_',colnames(lh))
    colnames(rh) = paste0('rh_',colnames(rh))
    bihemi = cbind(lh, rh)


    summary(bihemi$lh_Var)
    summary(bihemi$rh_Var)
    plot = ggplot(bihemi, aes(x=lh_Var, y=rh_Var)) +
                    geom_point(fill='lightskyblue', pch=21, size=3) +
                    geom_smooth(method='lm', se=FALSE, color='black', linetype='dashed', fullrange=T) +
                    theme_classic() +
                    scale_y_continuous(limits=c(limits[1],limits[2]), breaks=seq(limits[1],limits[2],limits[3]), expand=c(0,0)) +
                    scale_x_continuous(limits=c(limits[1],limits[2]), breaks=seq(limits[1],limits[2],limits[3]), expand=c(0,0)) +
                    ylab('Right Hemisphere Heritability (h2)') +
                    xlab('Left Hemisphere Heritability (h2)') +
                    theme(text = element_text(size = 10, colour = "black"),
                            axis.text.y = element_text(colour = "black"),
                            axis.text.x = element_text(colour = "black"),
                            plot.title = element_text(hjust = 0.5)) +
                    ggtitle(string)
    return(plot)
}


plot_h2_by_network = function(df, string, hex_df, limits){

    # add hemiphere information
    lh_labels   = as.character(df$trait[grep('lh_', df$trait)])
    rh_labels   = as.character(df$trait[grep('rh_', df$trait)])
    df$hemi     = ifelse(grepl('lh', df$trait), 'lh', 'rh')

    # merge heritability estimates with colors
    hex_df_t = data.frame(color=t(hex_df), net=colnames(hex_df))
    df       = merge(x=df, y=hex_df_t, all.x=T, by.x='net', by.y='net')

    # put networks in correct order
    df$net = factor(df$net, levels=colnames(hex_df))
    df     = df[order(df$net),]

    p = ggplot(df, aes(x=net, y=Var, alpha=hemi, fill=interaction(hemi,net))) +
                geom_bar(stat="identity", color="black", position=position_dodge(), width=.8) +
                geom_errorbar(aes(ymin=Var-SE, ymax=Var+SE), width=.4, position=position_dodge(.8)) +
                theme_classic() +
                scale_y_continuous(limits = c(limits[1], limits[2]), breaks=seq(limits[1], limits[2], limits[3]), expand = c(0,0)) +
                scale_fill_manual(values=as.character(df$color)) +
                theme(legend.position = "none") +
                scale_alpha_discrete(range  = c( 1, .5)) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) +
                        theme(text = element_text(size = 10, colour = "black"),
                            axis.text.y = element_text(colour = "black"), plot.title = element_text(hjust = 0.5),
                            axis.text.x = element_text(colour = "black")) + ggtitle(string)
    return(p)
}


calc_net_heritability = function(df, trait_list, covars, net_names){

    covar_formula = paste0(covars, collapse=' + ')

    #  trait = 'rh_SomatomotorA_sa'
    h2_df = NULL
    for (trait in trait_list){
        write(trait, '')
        hemi     = strsplit(trait, '_')[[1]][1]
        net_name = strsplit(trait, '_')[[1]][2]

        df[trait] = scale(df[trait])
        formula   = as.formula(paste0(trait, ' ~ ', covar_formula))
        rhog      = solarPolygenic(formula, df)

        out_row       = rhog$vcf[1,]
        print(out_row)
        out_row$net   = net_name
        out_row$trait = trait
        h2_df         = rbind(h2_df, out_row)
    }
    return(h2_df)
}


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


# read surface area network estimates and phenotypes
# ------------
base_dir    = '/gpfs/milgram/project/holmes/kma52/topo_herit'
hcp_twin_df = read_csv(paste0(base_dir, '/data/HCP/hcp_for_solar_allvars.csv'))

# read network name information
labels    = readMat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = unlist(labels$network.name)


# this isn't in the paper, but we double-checked that our method for calculating surface area is robust
# method #1 transformed each individualized parcellation from fs_LR32k to fsavg6 space, then calculated SA
# method #2 transformed SA data from favg6 into fs_LR32k space, then calculated network SA
# this code below shows that they both give basically the same data r-values>0.99

# surface area data, method #1
sa_df           = read_csv(paste0(base_dir, '/data/HCP/indiv_net_surfarea_native_freesurfer.csv'))
colnames(sa_df) = paste0('v1_', colnames(sa_df), '_sa')

# surface area data, method #2
sa_v2_df = read.csv(paste0(base_dir, '/data/HCP/method2_sa_fslr32k_freesurfer.csv'))
colnames(sa_v2_df) = paste0('v2_', colnames(sa_v2_df))

combo_sa_df = merge(x=sa_v2_df, y=sa_df, by.x='v2_id', by.y='v1_id_sa')
cor_arr = NULL
for (net in net_names){
    write(net,'')
    for (hemi in c('lh','rh')){
        cor_val = cor(combo_sa_df[[paste0('v1_', hemi, '_', net, '_sa')]], combo_sa_df[[paste0('v2_', hemi, '_', net, '_sa')]])
        cor_arr = c(cor_arr, cor_val)
    }
}
# high correspondence between each method
cor_arr


# use method #1, with SA values produced by mris_anatomical_stats
# ------------
sa_df = read_csv(paste0(base_dir, '/data/HCP/indiv_net_surfarea_native_freesurfer.csv'))
colnames(sa_df) = paste0(colnames(sa_df), '_sa')

# merge
hcp_df = merge(x=hcp_twin_df, y=sa_df, by.x='id', by.y='id_sa')


# format/create a few vars
hcp_df$FS_IntraCranial_Vol_scale = as.numeric(scale(hcp_df$FS_IntraCranial_Vol))
hcp_df$age2 = hcp_df$Age_in_Yrs^2
hcp_df$fs_aparc_total_sa_scale = as.numeric(scale(hcp_df$fs_aparc_total_sa))
hcp_df$Age_in_Yrs_scale = as.numeric(scale(hcp_df$Age_in_Yrs))
hcp_df$age2_scale = as.numeric(scale(hcp_df$age2))


# write a dataframe with all the info needed for SOLAR calculation
base_dir    = '/gpfs/milgram/project/holmes/kma52/topo_herit'

# only missing a single BMI value, so plug in median
hcp_df$BMI[is.na(hcp_df$BMI)] = median(hcp_df$BMI, na.rm=T)
hcp_df$Ethnicity_bin          = as.numeric(grepl('Not', hcp_df$Ethnicity))


# save
# save

write_csv(x=hcp_df, paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))

# save
# save


# covariates to include in heritability estimate
covars = c('Age_in_Yrs_scale',
            'sex',
            'age2_scale',
            'FS_IntraCranial_Vol_scale',
            'fs_aparc_total_sa_scale',
            'age2*sex',
            'Age_in_Yrs*sex')


hcp_df$sex_bin = ifelse(hcp_df$sex == 'M', 1, 0)
hcp_df$age_sex = hcp_df$Age_in_Yrs * hcp_df$sex_bin
hcp_df$age2_sex = hcp_df$age2 * hcp_df$sex_bin


# HERITABILITY -- Surface Area
# ------------
trait_list   = c(paste0('lh_', net_names, '_sa'), paste0('rh_', net_names, '_sa'))
relSA_h2_df  = calc_net_heritability(df=hcp_df, trait_list=trait_list, covars=covars, net_names=net_names)
out_file     = paste0(base_dir, '/data/HCP/solar_surf_area_h2_by_net_rr.csv')
write_csv(x=relSA_h2_df, out_file)
relSA_h2_df  = read.csv(paste0(base_dir, '/data/HCP/solar_surf_area_h2_by_net_rr.csv'))


hcp_df[paste0('lh_', net_names, '_sa')]


# plot heritability of surface area for each network
rel_h2_bynet_plot = plot_h2_by_network(df=relSA_h2_df, string='relative_SA', hex_df=hex_df, limits=c(0,.7,.1))
out_path = paste0(base_dir, '/figures/h2_relSA_h2_by_net_barplot_rr.pdf')
ggsave(plot=rel_h2_bynet_plot, out_path, width=4, height=2.5)
out_path


# correlate lh/rh h2
plot_relSA  = plot_hemi_h2_corr(h2_df=relSA_h2_df, string='relative_SA', limits=c(.2,.6,.1))
fig_out     = paste0(base_dir, '/figures/hemi_h2_relativeSA_corrplot.pdf')
CairoPDF(fig_out, height=4, width=4)
print(plot_relSA)
dev.off()
fig_out



# heritability of surfaces area, not controlling for eTIV
# just make sure the results are consistent
covars = c('Age_in_Yrs_scale',
            'sex',
            'age2_scale',
            'fs_aparc_total_sa_scale',
            'age2*sex',
            'Age_in_Yrs*sex')
relSA_h2_noicv_df = calc_net_heritability(df=hcp_df, trait_list=trait_list, covars=covars, net_names=net_names)
out_file     = paste0(base_dir, '/data/HCP/solar_surf_area_h2_by_net_noetiv_rr.csv')
write_csv(x=relSA_h2_noicv_df, out_file)

relSA_h2_noicv_df = read_csv(paste0(base_dir, '/data/HCP/solar_surf_area_h2_by_net_noetiv_rr.csv'))
#relSA_h2_df = read_csv(paste0(base_dir, '/data/HCP/solar_relativeSA_h2_byNet.csv'))

# doesn't matter if we control for ICV/eTIV or not
cor(relSA_h2_df$Var, relSA_h2_noicv_df$Var)



lh_df = relSA_h2_df[grep('lh_', relSA_h2_df$trait),]
rh_df = relSA_h2_df[grep('rh_', relSA_h2_df$trait),]

cor.test(lh_df$Var, rh_df$Var)
cor.test(lh_df$Var, rh_df$Var, method='spearman')





# same as above, but don't include fs total surface area
# just showing high level of convergence 
covars = c('Age_in_Yrs_scale',
            'sex',
            'age2_scale',
            'FS_IntraCranial_Vol_scale',
            'age2*sex',
            'Age_in_Yrs*sex')

lh_rel_sa = hcp_df[paste0('lh_', net_names, '_sa')]/hcp_df$lh_sa
rh_rel_sa = hcp_df[paste0('rh_', net_names, '_sa')]/hcp_df$rh_sa

colnames(lh_rel_sa) = paste0(colnames(lh_rel_sa), '_ratio')
colnames(rh_rel_sa) = paste0(colnames(rh_rel_sa), '_ratio')

rel_sa_df = cbind(lh_rel_sa, rh_rel_sa)
hcp_df = cbind(hcp_df, rel_sa_df)

trait_list = colnames(rel_sa_df)

ratio_df  = calc_net_heritability(df=hcp_df, trait_list=trait_list, covars=covars, net_names=net_names)


cor.test(ratio_df$Var, relSA_h2_df$Var, method='spearman')

