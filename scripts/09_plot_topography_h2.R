library(tidyverse)
library(R.matlab)


# This script will:
#   1. Read matlab-generated heritability estimates of network topography
#   2. Barplot of net heritbality (split by net/hemi)
#   3. Correlate LH/RH h2-multi estimates
#   4. Test for h2-multi differences between unimodal/hetero cortex


plot_h2_by_network = function(df, string, hex_df, limits){

    # merge heritability estimates with colors
    hex_df_t = data.frame(color=t(hex_df), net=colnames(hex_df))
    df       = merge(x=df, y=hex_df_t, all.x=T, by.x='net_names', by.y='net')

    # put networks in correct order
    df$net = factor(df$net_names, levels=colnames(hex_df))
    df     = df[order(df$net),]

    p = ggplot(df, aes(x=net, y=h2, alpha=hemi, fill=interaction(hemi,net))) +
                geom_bar(stat="identity", color="black", position=position_dodge(), width=.8) +
                geom_errorbar(aes(ymin=h2-jack_se, ymax=h2+jack_se), width=0.5, size=.25, color='black', position=position_dodge(width=0.8)) +
                #geom_errorbar(aes(ymin=Var-SE, ymax=Var+SE), width=.4, position=position_dodge(.8)) +
                theme_classic() +
                scale_y_continuous(limits = c(limits[1], limits[2]), breaks=seq(limits[1], limits[2], limits[3]), expand = c(0,0)) +
                scale_fill_manual(values=as.character(df$color)) +
                theme(legend.position = "none") +
                scale_alpha_discrete(range  = c( 1, .5)) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
                        theme(text = element_text(size = 10, colour = "black"),
                            axis.text.y = element_text(colour = "black"), plot.title = element_text(hjust = 0.5),
                            axis.text.x = element_text(colour = "black")) +
                ggtitle(string) +
                ylab('Heritability')
    return(p)
}


# color information for each network
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
                    TemporalParietal='#0929fa',
                    overall='#2B2D2F')



############
# read data
############
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit'

# read network name information
labels    = readMat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = unlist(labels$network.name)

# read h2-multi estimates
topology_h2    = read.csv(paste0(base_dir, '/data/topology_heritability/dice_network_topology_h2.csv'), header=T)
lh_topology_h2 = read.csv(paste0(base_dir, '/data/topology_heritability/lh_dice_network_topology_h2.csv'), header=T)
rh_topology_h2 = read.csv(paste0(base_dir, '/data/topology_heritability/rh_dice_network_topology_h2.csv'), header=T)

# combine LH/RH data
lh_topology_h2$net_names = c('overall', net_names[as.numeric(lh_topology_h2$network[2:18])])
rh_topology_h2$net_names = c('overall', net_names[as.numeric(rh_topology_h2$network[2:18])])
combined_df = rbind(lh_topology_h2, rh_topology_h2)



##################################
# h2-multi by heteromodal/unimodal
##################################
combined_netonly_df = combined_df[combined_df$net_names != 'overall',]
combined_netonly_df$modality = ifelse(grepl('Vis|Som|Aud',combined_netonly_df$net_names), 'unimodal', 'heteromodal')
combined_netonly_df %>% group_by(modality) %>% summarise(h2=mean(h2), stdev=sd(h2))

# test for differences between groups
lm(h2 ~ modality, data=combined_netonly_df)
summary(lm(h2 ~ modality, data=combined_netonly_df))



plot_obj = plot_h2_by_network(df=combined_df, string='Topology heritability', hex_df=hex_df, limits=c(0, .2, .05))
filename = paste0(base_dir, '/figures/h2_of_topology_heritability.pdf')
ggsave(filename=filename, plot=plot_obj, width=4, height=2.5)
filename


df_corr = merge(x=lh_topology_h2, y=rh_topology_h2, by='net_names')
df_corr = df_corr[!grepl('overall', df_corr$net_names),]
cor.test(df_corr$h2.x, df_corr$h2.y)
cor.test(df_corr$h2.x, df_corr$h2.y, method='spearman')





