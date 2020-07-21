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
    lh_traits = paste0('lh_', as.character(1:17), '_',  string)
    rh_traits = paste0('rh_', as.character(1:17), '_',  string)

    lh = h2_df[h2_df$trait %in% lh_traits,]
    rh = h2_df[h2_df$trait %in% rh_traits,]
    hemi_cor = cor.test(lh$Var, rh$Var)

    hemi_cor_sp = cor.test(lh$Var, rh$Var, method='spearman')

    print(hemi_cor)

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
    df$hemi = 'lh'
    lh_labels  = paste0('lh_', as.character(1:17), '_', string)
    rh_labels  = paste0('rh_', as.character(1:17), '_', string)
    df = df[df$trait %in% c(lh_labels, rh_labels),]
    df$hemi[df$trait %in% rh_labels] = 'rh'

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

    covar_formula = paste0(covars[1:5], collapse=' + ')

    h2_df = NULL
    for (trait in trait_list){
        write(trait, '')
        hemi    = strsplit(trait, '_')[[1]][1]
        net_num = strsplit(trait, '_')[[1]][2]

        # naming information
        if ( net_num %in% as.character(1:17) ){
            net_num  = as.numeric(net_num)
            net_name = net_names[net_num]
        } else {
            net_num  = NA
            net_name = NA
        }

        df[trait] = scale(df[trait])
        formula   = as.formula(paste0(trait, ' ~ ', covar_formula))
        rhog      = solarPolygenic(formula, df)

        out_row       = rhog$vcf[1,]
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



# ----------------------------------------------
# read surface area network estimates/phenotypes
# ----------------------------------------------
base_dir    = '/gpfs/milgram/project/holmes/kma52/topo_herit'
out_file    = paste0(base_dir, '/data/HCP/surface_area_by_net.csv')
hcp_twin_df = read_csv(paste0(base_dir, '/data/HCP/hcp_for_solar_allvars.csv'))
sa_df       = read_csv(out_file)

# merge
hcp_df = merge(x=hcp_twin_df, y=sa_df, by.x='id', by.y='Subject')

# add vars
hcp_df$FS_IntraCranial_Vol_scale = as.numeric(scale(hcp_df$FS_IntraCranial_Vol))
hcp_df$age2 = hcp_df$Age_in_Yrs^2


# read network name information
labels    = readMat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = unlist(labels$network.name)


# ---------------------------------
# calc total surface area estimates
# ---------------------------------
hcp_df$total_SA    = rowSums(hcp_df[paste0('bihemi_', as.character(1:17), '_sum_SA')])
hcp_df$lh_total_SA = rowSums(hcp_df[paste0('lh_', as.character(1:17), '_sum_SA')])
hcp_df$rh_total_SA = rowSums(hcp_df[paste0('rh_', as.character(1:17), '_sum_SA')])

# get relative LH SA estimates
hcp_df[paste0('lh_', as.character(1:17), '_relative_SA')] = hcp_df[paste0('lh_', as.character(1:17), '_sum_SA')]/hcp_df$lh_total_SA
hcp_df[paste0('rh_', as.character(1:17), '_relative_SA')] = hcp_df[paste0('rh_', as.character(1:17), '_sum_SA')]/hcp_df$rh_total_SA

# write a dataframe with all the info needed for SOLAR calculation
base_dir    = '/gpfs/milgram/project/holmes/kma52/topo_herit'

# only missing a single BMI value, so just plug in median
hcp_df$BMI[is.na(hcp_df$BMI)] = median(hcp_df$BMI, na.rm=T)
hcp_df$Ethnicity_bin          = as.numeric(grepl('Not', hcp_df$Ethnicity))
write_csv(x=hcp_df, paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))


# covariates to include in heritability estimate
covars = c('Age_in_Yrs',
            'age2',
            'age2*sex',
            'Age_in_Yrs*sex',
            'Ethnicity_bin',
            'Height',
            'BMI',
            'FS_IntraCranial_Vol_scale')


# -------------------------------------
# HERITABILITY -- Relative Surface Area
# -------------------------------------
relative_sa_cols = c(paste0('lh_', as.character(1:17), '_relative_SA'), paste0('rh_', as.character(1:17), '_relative_SA'))
relSA_h2_df      = calc_net_heritability(df=hcp_df, trait_list=relative_sa_cols, covars=covars, net_names=net_names)
out_file    = paste0(base_dir, '/data/HCP/solar_relativeSA_h2_byNet.csv')
write_csv(x=relSA_h2_df, out_file)
relSA_h2_df = read_csv(paste0(base_dir, '/data/HCP/solar_relativeSA_h2_byNet.csv'))



# ---------------------------
# PLOT -- Relative SA h2 plot
# ---------------------------
plot_relSA  = plot_hemi_h2_corr(h2_df=relSA_h2_df, string='relative_SA', limits=c(.2,.6,.1))
fig_out     = paste0(base_dir, '/figures/hemi_h2_relativeSA_corrplot.pdf')
CairoPDF(fig_out, height=4, width=4)
print(plot_relSA)
dev.off()


# ----------------
# Relative SA plot
# ----------------
rel_h2_bynet_plot = plot_h2_by_network(df=relSA_h2_df, string='relative_SA', hex_df=hex_df, limits=c(0,.7,.1))

out_path = paste0(base_dir, '/figures/h2_relSA_h2_by_net_barplot.pdf')
ggsave(plot=rel_h2_bynet_plot, out_path, width=4, height=2.5)
out_path