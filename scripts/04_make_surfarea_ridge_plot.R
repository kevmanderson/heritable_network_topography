library(tidyverse)
library(solarius)
library(R.matlab)
library(ggridges)
library(Cairo)
library(corrplot)
library(ggridges)


# This script will:
#   1. Create ridge plot (Figure 1b) of relative surface area (SA) across HCP individuals
#   2. Calculate coefficient of variation for SA and plot by net/hemi (Figure 2c)
#   3. Standard devitation across individuals
#   4. Calculate hetero/uni network differences for CV and SD


# network colors for plotting
hex_df = data.frame(VisualA='#781286',
                    VisualB='#ff0101',
                    VisualC='#7a8733',
                    SomatomotorA='#4682b2',
                    SomatomotorB='#2bcca2',
                    Auditory='#dcf8a5',
                    VentralAttentionA='#c43afb',
                    VentralAttentionB='#ff98d6',
                    DorsalAttentionA='#4a9b3d',
                    DorsalAttentionB='#007610',
                    ControlA='#e69423',
                    ControlB='#87324c',
                    ControlC='#778cb1',
                    DefaultA='#ffff02',
                    DefaultB='#cd3e50',
                    DefaultC='#000083',
                    TemporalParietal='#0929fa',
                    overall='#2B2D2F')


# coefficient of variation
CV = function(arr){
    mean = mean(arr)
    sd   = sd(arr)
    cv   = (sd/mean)*100
    return(cv)
}


# ------------------------
# read HCP data/label info
# ------------------------
base_dir  = '/gpfs/milgram/project/holmes/kma52/topo_herit'
hcp_df    = read.csv(paste0(base_dir, '/data/HCP/hcp_wMRI_for_solar.csv'))
labels    = readMat('/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Kong2019_MSHBM/lib/group_priors/HCP_40/17network_labels.mat')
net_names = unlist(labels$network.name)


# set plot parameters
ct = 1
limit_arr = list(c(0,.16), c(0,5500), c(.5,2))


# -----------------------------
# bihemi 17 network relative SA
# -----------------------------
type    = '_relative_SA'
lh_cols = paste0('lh_', as.character(1:17), type)
rh_cols = paste0('rh_', as.character(1:17), type)
plot_me = hcp_df[c(lh_cols, rh_cols)]

# melt dataframe from wide to long
sa_melt_df = reshape2::melt(plot_me)
sa_melt_df$net_num = unlist(lapply(sa_melt_df$variable, function(x) strsplit(as.character(x), '_')[[1]][[2]]))
sa_melt_df$net_num = as.numeric(sa_melt_df$net_num)
sa_melt_df$net_name = net_names[sa_melt_df$net_num]
head(sa_melt_df)


# plot unimodal before heteromodal networks
plot_order = c('VisualA','VisualB','VisualC','SomatomotorA','SomatomotorB','Auditory',
                'VentralAttentionA', 'VentralAttentionB', 'DorsalAttentionA', 'DorsalAttentionB',
                'ControlA','ControlB','ControlC',
                'DefaultA','DefaultB','DefaultC', 'TemporalParietal')
sa_melt_df$net_name = factor(sa_melt_df$net_name, levels=rev(plot_order))
sa_melt_df$hemi     = ifelse(grepl('lh', sa_melt_df$variable), 'lh', 'rh')


# add empty rows to make space between ridges
tmp_row = sa_melt_df[1,]
for (net in 1:17){
    tmp_row$value = NA
    tmp_row$hemi  = 'fill'
    tmp_row$net_name = net_names[net]
    tmp_row$variable = paste0('fill_',net,'_relative_SA')
    sa_melt_df = rbind(sa_melt_df, tmp_row)
}
sa_melt_df$plot_order = paste0(sa_melt_df$net_name, '_', sa_melt_df$hemi)


# adjust the plot_order variable to account for the "fill" rows
plot_full_order = NULL
for (net in plot_order){
    plot_full_order = c(plot_full_order, paste0(net, '_',  c('lh','rh','fill')))
}
sa_melt_df$plot_order = factor(sa_melt_df$plot_order, levels=plot_full_order)


# plot variance - collapse across hemispheres
ridge_plot = ggplot(data=sa_melt_df, aes(x = value, y = net_name)) +
                        geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3) +
                        scale_fill_gradientn(limits = c(limit_arr[[ct]][1], limit_arr[[ct]][2]),
                        colours = c("#C22A39",'#F7B494', "#F8F9FA", '#A2CFE5', "#297CBA"), name = type) +
                        labs(title = paste0('HCP ', type)) +
                        theme_classic() +
                        scale_x_continuous(limits=c(limit_arr[[ct]][1], limit_arr[[ct]][2]), expand=c(0,0)) +
                        theme(axis.title.y=element_blank(),
                        axis.ticks = element_blank(),
                        text=element_text(size=10,  family="Arial", color='black'),
                        axis.text=element_text(size=10,  family="Arial", color='black'),
                        plot.title = element_text(hjust = 0.5),
                        legend.title=element_blank()) +
                        labs(title="HCP Relative Network Size", x = 'Proportion of Total Surface Area (%)')

# plot variance - separate across hemispheres
ridge_plot = ggplot(data=sa_melt_df, aes(x = value, y = plot_order)) +
                        geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3) +
                        scale_fill_gradientn(limits = c(limit_arr[[ct]][1], limit_arr[[ct]][2]),
                        colours = c("#C22A39",'#F7B494', "#F8F9FA", '#A2CFE5', "#297CBA"), name = type) +
                        labs(title = paste0('HCP ', type)) +
                        theme_classic() +
                        scale_x_continuous(limits=c(limit_arr[[ct]][1], limit_arr[[ct]][2]), expand=c(0,0)) +
                        theme(axis.title.y=element_blank(),
                        axis.ticks = element_blank(),
                        text=element_text(size=10,  family="Arial", color='black'),
                        axis.text=element_text(size=10,  family="Arial", color='black'),
                        plot.title = element_text(hjust = 0.5),
                        legend.title=element_blank()) +
                        labs(title="HCP Relative Network Size", x = 'Proportion of Total Surface Area (%)')


# ---------
# save plot
# ---------
hist_file = paste0(base_dir, '/figures/hcp_hemisplit_ggridge',type,'.pdf')
CairoPDF(hist_file, width=4, height=4, family='Arial')
print(ridge_plot)
dev.off()

ggsave(plot=ridge_plot, file=hist_file, width=5, height=5)




# --------------------------------------------------
# create a long df with relative SA of every subject
# --------------------------------------------------
type    = '_relative_SA'
lh_cols = paste0('lh_', as.character(1:17), type)
rh_cols = paste0('rh_', as.character(1:17), type)
plot_me = hcp_df[c('id', lh_cols, rh_cols)]

# wide to long
sa_melt_df2 = reshape2::melt(plot_me, id.vars='id')
sa_melt_df2 = merge(x=sa_melt_df2, y=hcp_df[c('id')], by='id')
sa_melt_df2$net_num = unlist(lapply(sa_melt_df2$variable, function(x) strsplit(as.character(x), '_')[[1]][[2]]))
sa_melt_df2$hemi = unlist(lapply(sa_melt_df2$variable, function(x) strsplit(as.character(x), '_')[[1]][[1]]))
sa_melt_df2$net_num = as.numeric(sa_melt_df2$net_num)
sa_melt_df2$net_name = net_names[sa_melt_df2$net_num]


# --------------------------------
# calc coef variation for lh/rh SA
# --------------------------------
lh_coef_var_net = mapply(CV, hcp_df[paste0('lh_', as.character(1:17), '_relative_SA')])
rh_coef_var_net = mapply(CV, hcp_df[paste0('rh_', as.character(1:17), '_relative_SA')])

# combine into dataframe
lh_cv_df = data.frame(coef_var=lh_coef_var_net, net_label=names(lh_coef_var_net), net_names=net_names, hemi='lh')
rh_cv_df = data.frame(coef_var=rh_coef_var_net, net_label=names(rh_coef_var_net), net_names=net_names, hemi='rh')
coef_var_df = rbind(lh_cv_df, rh_cv_df)
coef_var_df$net_names = factor(coef_var_df$net_names, levels=rev(plot_order))

head(coef_var_df)


# prepare for plotting
hex_df_t    = data.frame(color=t(hex_df), net=colnames(hex_df))
coef_var_df = merge(x=coef_var_df, y=hex_df_t, all.x=T, by.x='net_names', by.y='net')
coef_var_df$net_names = factor(coef_var_df$net_names, levels=colnames(hex_df))
coef_var_df     = coef_var_df[order(coef_var_df$net_names),]


# --------------------------------
# Plot CV by network (relative SA)
# --------------------------------
p = ggplot(coef_var_df, aes(x=net_names, y=coef_var, alpha=hemi, fill=interaction(hemi,net_names))) +
            geom_bar(stat="identity", color="black", position=position_dodge(), width=.8) +
            #geom_errorbar(aes(ymin=Var-SE, ymax=Var+SE), width=.4, position=position_dodge(.8)) +
            theme_classic() +
            scale_y_continuous(limits = c(0,40), breaks=seq(0,40,10), expand = c(0,0)) +
            scale_fill_manual(values=as.character(coef_var_df$color)) +
            theme(legend.position = "none") +
            scale_alpha_discrete(range  = c( 1, .5)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank()) +
                    theme(text = element_text(size = 10, colour = "black"),
                        axis.text.y = element_text(colour = "black"), plot.title = element_text(hjust = 0.5),
                        axis.text.x = element_text(colour = "black")) +
            coord_flip()


filename = paste0(base_dir, '/figures/variance_by_network.pdf')
ggsave(filename=filename, plot=p, width=2.5, height=3.5)
filename





# -----------------------------------
# get network-wise SD of relative SA
# -----------------------------------
plot_std_df = sa_melt_df2 %>%
                group_by(net_name, hemi) %>%
                summarise(net_size=mean(value), net_sd=sd(value), net_se=sd(value)/length(value))

plot_std_df = merge(x=plot_std_df, y=hex_df_t, all.x=T, by.x='net_name', by.y='net')
plot_std_df$net_name = factor(plot_std_df$net_name, levels=colnames(hex_df))

p = ggplot(plot_std_df, aes(x=net_name, y=net_sd, alpha=hemi, fill=interaction(hemi,net_name))) +
            geom_bar(stat="identity", color="black", position=position_dodge(), width=.8) +
            theme_classic() +
            scale_y_continuous(limits = c(0,.018), breaks=seq(0,.018,.016), expand = c(0,0)) +
            scale_fill_manual(values=as.character(coef_var_df$color)) +
            theme(legend.position = "none") +
            scale_alpha_discrete(range  = c( 1, .5)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.x=element_blank()) +
                    theme(text = element_text(size = 10, colour = "black"),
                        axis.text.y = element_text(colour = "black"), plot.title = element_text(hjust = 0.5),
                        axis.text.x = element_text(colour = "black")) +
            coord_flip()

filename = paste0(base_dir, '/figures/stdev_by_network.pdf')
ggsave(filename=filename, plot=p, width=2.5, height=3.5)
filename


# ----------------------------------
# test hetero/uni differences in SD
# ----------------------------------
plot_std_df$type = ifelse(grepl('Vis|Som|Aud', plot_std_df$net_name), 'uni', 'hetero')
summary(lm(formula = "net_sd ~ type", data=plot_std_df))

std_var_boxplot = ggplot(data=plot_std_df, aes(x=type, y=net_sd, color = type)) +
                    geom_point() +
                    geom_boxplot(alpha=0.8) +
                    scale_color_manual(values=c('#82CFFD', '#00688B')) +
                    theme_classic() +
                    scale_y_continuous(limits=c(0, 0.018), breaks=seq(0, 0.018, 0.006), expand=c(0,0)) +
                    theme(axis.title.x=element_blank(),
                            text=element_text(size=10,  family="Arial", color='black'),
                            axis.text=element_text(size=10, color='black'),
                            plot.title = element_text(hjust = 0.5),
                            legend.position='none') +
                    labs(y = 'stdev')


filename = paste0(base_dir, '/figures/net_sd_boxplot.pdf')
CairoPDF(filename, width=2, height=2)
print(std_var_boxplot)
dev.off()
filename




# ---------------------------------
# test hetero/uni differences in CV
# ---------------------------------
coef_var_df$type = ifelse(grepl('Vis|Som|Aud', coef_var_df$net_names), 'uni', 'hetero')
coef_var_df %>% group_by(type) %>% summarize(cv_mean=mean(coef_var), cv_sd = sd(coef_var))

# same, but for coefficient of variation
summary(lm(formula = "coef_var ~ type", data = coef_var_df))

# test again after removing outlier
summary(lm(formula = "coef_var ~ type", data = coef_var_df %>% filter(net_label != 'rh_11_relative_SA')))

coef_var_boxplot = ggplot(data=coef_var_df, aes(x=type, y=coef_var, color = type)) +
                    geom_point() +
                    geom_boxplot(alpha=0.8) +
                    scale_color_manual(values=c('#82CFFD', '#00688B')) +
                    theme_classic() +
                    scale_y_continuous(limits=c(10, 40), breaks=seq(10, 40, 10), expand=c(0,0)) +
                    theme(axis.title.x=element_blank(),
                            text=element_text(size=10,  family="Arial", color='black'),
                            axis.text=element_text(size=10, color='black'),
                            plot.title = element_text(hjust = 0.5),
                            legend.position='none') +
                    labs(y = 'Coef Var')


filename = paste0(base_dir, '/figures/coef_of_variance_boxplot.pdf')
CairoPDF(filename, width=2, height=2)
print(coef_var_boxplot)
dev.off()



