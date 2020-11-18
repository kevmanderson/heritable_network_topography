library(tidyverse)
library(gifti)
library(R.matlab)
library(solarius)
library(ggpubr)
library(Cairo)
library(mgcv)
options(stringsAsFactors=F)


# read HCP pheno
# ------------
base_dir    = '/gpfs/milgram/project/holmes/kma52/topo_herit'
hcp_twin_df = read_csv(paste0(base_dir, '/data/HCP/hcp_for_solar_allvars.csv'))


# read surface area network estimates
# ------------
sa_df       = read_csv(paste0(base_dir, '/data/HCP/indiv_net_surfarea_native_freesurfer.csv'))
colnames(sa_df) = paste0(colnames(sa_df), '_sa')
sa_networks = colnames(sa_df)[grepl('lh_|rh_', colnames(sa_df))]
sa_networks = sa_networks[!grepl('None', sa_networks)]
sa_df$total_sa = rowSums(sa_df[grepl('lh_|rh_', colnames(sa_df))])


# merge
hcp_df = merge(x=hcp_twin_df, y=sa_df, by.x='id', by.y='id_sa')
hcp_df$sex_bin = ifelse(hcp_df$sex == 'M', 1, 0)



# as in Reardon/Seidlitz 2018, estimate general additive model for non-linear (allometric) scaling of network size
# --------------
hcp_df$sex = as.factor(hcp_df$sex)
gam_coefs_df = NULL
for (net in sa_networks){
    write(net, '')
    gam_model = gam(as.formula(paste0('log10(', net,') ~ s(Age_in_Yrs, by=sex) + log10(total_sa)')), data=hcp_df)

    hcp_df[paste0(net, '_allo_resid')] = gam_model$residuals
    print(print(gam_model$coefficients))

    gam_summary  = summary(gam_model)
    allo_beta    = as.numeric(gam_model$coefficients[2])
    gam_coefs_df = rbind(gam_coefs_df, data.frame(net=net, allo_beta=allo_beta))
}


# coefficients > 1 mean relative larger given head size (positive scaling)
gam_coefs_df         = gam_coefs_df[order(gam_coefs_df$allo_beta),]
gam_coefs_df$hemi    = ifelse(grepl('lh_', gam_coefs_df$net), 'lh', 'rh')
gam_coefs_df$network = unlist(lapply(gam_coefs_df$net, function(x) strsplit(x, '_')[[1]][[2]]))
head(gam_coefs_df)


# plot scaling coefficients
plot_order_df = gam_coefs_df %>% group_by(network) %>% summarise(m_allo_beta=mean(allo_beta))
plot_order    = plot_order_df$network[order(plot_order_df$m_allo_beta)]
gam_coefs_df$network = factor(gam_coefs_df$network, levels=plot_order)

gam_coefs_df = merge(x=gam_coefs_df, y=plot_order_df, by='network')

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

# color info
hex_df_t = data.frame(color=t(hex_df), net=colnames(hex_df))
hex_df_t = hex_df_t[match(hex_df_t$net, plot_order),]

# make sure the colors/networks are sorted proper
plot_df = merge(x=gam_coefs_df, y=hex_df_t, all.x=T, by.x='network', by.y='net')
plot_df$network = factor(plot_df$network, levels=plot_order)
plot_df = plot_df[order(plot_df$m_allo_beta),]


# plot allometric scaling coefficientfor each network
allo_plot = ggplot(data=plot_df, aes(x=network, y=allo_beta-1, alpha=hemi, fill=interaction(hemi,network))) +
                    geom_bar(stat="identity", color="black", position=position_dodge(), width=.8) +
                    coord_flip() +
                    theme_classic() +
                    geom_hline(yintercept=0) +
                    theme(axis.text.x = element_text(size=8, color='black'), axis.title.x=element_text(size=8, color='black'), axis.title.y=element_blank()) +
                        theme(text = element_text(size = 8, colour = "black"), legend.position='none',
                            axis.text.y = element_text(size=8, colour = "black"), plot.title = element_text(hjust = 0.5),
                            axis.text.x = element_text(size=8, colour = "black"), axis.ticks.x = element_line(color='black', size=.5)) +
                    scale_fill_manual(values=as.character(plot_df$color)) +
                    scale_y_continuous(limits=c(-.6,.6), expand=c(0,0), breaks=seq(-.6,.6,.3), labels=seq(-.6,.6,.3)+1) +
                    ylab('Allometric Scaling') +
                    scale_alpha_discrete(range = c(1, .5))

# save plot
out_path = paste0(base_dir, '/figures/indiv_network_allometric_scaling.pdf')
ggsave(plot=allo_plot, out_path, width=3, height=3)
out_path



# test whether allo coefs are different for uni/heteromodal networks
plot_df$type = ifelse(grepl('Vis|Som|Aud', plot_df$net), 'uni', 'het')
summary(lm('allo_beta ~ type', data=plot_df))
plot_df %>% group_by(type) %>% summarise(m_allo=mean(m_allo_beta), std_allo=sd(m_allo_beta))




# ------------
# ------------
# Heritability of allometric adjusted network size
# ------------
# ------------

calc_net_heritability = function(df, trait_list, covars){

    covar_formula = paste0(covars, collapse=' + ')

    h2_df = NULL
    for (trait in trait_list){
        write(trait, '')
        hemi     = strsplit(trait, '_')[[1]][1]
        net_name = strsplit(trait, '_')[[1]][2]

        df[trait] = scale(df[trait])
        formula   = as.formula(paste0(trait, ' ~ 1'))
        rhog      = solarPolygenic(formula, df)

        out_row       = rhog$vcf[1,]
        print(out_row)
        out_row$net   = net_name
        out_row$trait = trait
        h2_df         = rbind(h2_df, out_row)
    }
    return(h2_df)
}


# HERITABILITY -- Surface Area
# ------------
alloSA_h2_df  = calc_net_heritability(df=hcp_df, trait_list=paste0(sa_networks, '_allo_resid'), covars=covars)
out_file     = paste0(base_dir, '/data/HCP/solar_allometric_surf_area_h2_by_net.csv')
write_csv(x=alloSA_h2_df, out_file)

alloSA_h2_df  = read.csv(paste0(base_dir, '/data/HCP/solar_allometric_surf_area_h2_by_net.csv'))
relSA_h2_df   = read.csv(paste0(base_dir, '/data/HCP/solar_surf_area_h2_by_net_rr.csv'))


alloSA_h2_df$modality = ifelse(grepl('Som|Vis|Aud', alloSA_h2_df$net), 'uni', 'het')
alloSA_h2_df %>% group_by(modality) %>% summarise(m_h2=mean(Var))
summary(lm('Var ~ modality', data=alloSA_h2_df))


alloSA_h2_df$trait = gsub('_allo_resid', '', alloSA_h2_df$trait)
combo_df = merge(x=relSA_h2_df, y=alloSA_h2_df, by='trait')
cor.test(combo_df$Var.x, combo_df$Var.y, method='spearman')




# HERITABILITY -- Surface Area
# ------------

# plot SA heritability, grouped by uni/hetero modal
# plot version 2
hex_df_t    = data.frame(color=t(hex_df), net=colnames(hex_df))
plot_df     = merge(x=alloSA_h2_df, y=hex_df_t, all.x=T, by.x='net', by.y='net')
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

hist_file = paste0(base_dir, '/figures/hcp_uni_hetero_modal_boxplot_allometric_rr.pdf')
CairoPDF(hist_file, width=1.75, height=1.75, family='Arial')
print(boxplot)
dev.off()
hist_file




# plot correlation between SA h2, estimated with linear vs allometric adjustment
plot_df = data.frame(allo_h2=alloSA_h2_df$Var, h2=relSA_h2_df$Var)
hemi_cor_sp = cor.test(plot_df$allo_h2, plot_df$h2, method='spearman')

plot = ggplot(plot_df, aes(x=allo_h2, y=h2)) +
                geom_point(fill='lightgray', pch=21, size=3) +
                geom_smooth(method='lm', se=FALSE, color='black', linetype='dashed', fullrange=T) +
                theme_classic() +
                ylab('Right Hemisphere Heritability (h2)') +
                xlab('Left Hemisphere Heritability (h2)') +
                theme(text = element_text(size = 10, colour = "black"),
                        axis.text.y = element_text(colour = "black"),
                        axis.text.x = element_text(colour = "black"),
                        plot.title = element_text(hjust = 0.5)) +
                ggtitle('Allometric vs Linear Adjustment')
CairoPDF('/gpfs/milgram/project/holmes/kma52/topo_herit/figures/allo_vs_not.pdf', width=4, height=4)
print(plot)
dev.off()








