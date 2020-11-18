


% This script will:
%   1. Read pheno/covar/kinship matrices created by "07_net_topology_heritability.R"
%   2. Iterate over each network (bihemi/split by hemi) and calculate multi-dim heritability (using "h2_mat")


% set up paths/dirs
addpath('/gpfs/milgram/project/holmes/kma52/topo_herit/scripts')
base_dir = '/gpfs/milgram/project/holmes/kma52/topo_herit/data/topology_heritability'

% HCP kinship matrix
kin   = csvread(fullfile(base_dir, 'K.csv'), 1);

% HCP covariates
covar = csvread(fullfile(base_dir, 'covar.csv'), 1);

% Family IDs
fam = csvread(fullfile(base_dir, 'F.csv'), 1);

% calculate heritability of network topography
hemi_arr = {'lh','rh','bihemi'};
for (type = 1:3)
    hemi         = hemi_arr{type};
    hemi

    out_mat      = table();
    pheno        = csvread(fullfile(base_dir, [hemi '_overall_dice_P.csv']), 0);

    K=kin;
    P=pheno;
    X=covar;

    [h2, p_perm, jack_se] = h2_mat(pheno, kin, covar, 1000, fam);
    %[h2, p_perm] = h2_mat(pheno, kin, covar, 10);
    out_mat(1,:) = {h2, p_perm, jack_se, 'overall', hemi};
    out_mat.Properties.VariableNames = {'h2','pval','jack_se','network','hemi'}

    for (i = 1:17)
        disp(i)
        pheno = csvread(fullfile(base_dir, [hemi '_dice_net_', num2str(i),'_P.csv']), 0);
        [h2, p_perm, jack_se]   = h2_mat(pheno, kin, covar, 1000, fam);
        out_mat(i+1,:) = {h2, p_perm, jack_se, num2str(i), hemi};
    end

    out_path = fullfile(base_dir, [hemi '_dice_network_topology_h2.csv']);
    writetable(out_mat, out_path)
end

% end
% end
% end






% example using height phenotype

pheno = csvread(fullfile(base_dir, 'Height_P.csv'), 0);
covar = csvread(fullfile(base_dir, 'Height_covar.csv'), 1);

[h2, p_perm] = h2_mat(pheno, kin, covar, 10);


out_mat = {};
pheno = csvread(fullfile(base_dir, 'dice_net_12_P.csv'), 0);
[h2, p_perm] = h2_mat(pheno, kin, covar, 1000);
cur_row = [h2, p_perm];





