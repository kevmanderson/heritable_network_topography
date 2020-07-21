
% This script will:
%   1. Calculate multi-dim heritability at the vertex-level, using previously created dice sim matrices from "10_vertex_dice_matrix_slurm.R"
%   2. Given the number of vertices, it's best to chunk and run in parallel


% set up paths/dirs
addpath('/gpfs/milgram/project/holmes/kma52/topo_herit/scripts')
base_dir = '/gpfs/milgram/scratch/holmes/topo_herit/data/vert_dice';

% HCP kinship matrix
kin   = csvread(fullfile(base_dir, 'K.csv'), 1);

% HCP covariates
covar = csvread(fullfile(base_dir, 'C.csv'), 1);

start_list = {1,3001,6001,9001,12001,15001,18001,21001,24001,27001,30001};

hemi = 'L';
hemi = 'R';
for (i = 1:11);

    i = 1;
    out_mat = table();
    start   = start_list{i};
    stop    = start + 2999
    row     = 1;
    for (vert = start:stop);

        file_path = fullfile(base_dir, ['local_dice_' hemi '_vertex_' num2str(vert) '_idxfrom0_radius10_matrix.csv']);
        %file_path = fullfile(base_dir, ['local_dice_' hemi '_vertex_' num2str(vert) '_idxfrom0_radius5_matrix.csv']);
        disp(file_path)

        if (exist(file_path) == 2);
            pheno        = csvread(file_path, 0);
            if size(pheno,1) == size(pheno,2);
                [h2, p_perm] = h2_mat(pheno, kin, covar, 0);
                disp(h2)
                out_mat(row,:) = {h2, p_perm, 'overall', hemi, vert};
            end
        else
            out_mat(row,:) = {NaN, NaN, 'overall', hemi, vert};
        end
        row = row + 1;
    end
    out_path = fullfile(base_dir, [hemi '_' num2str(start) '_' num2str(stop) '_dice_network_topology_h2.csv']);
    %out_path = fullfile(base_dir, [hemi '_' num2str(start) '_' num2str(stop) '_dice_network_topology_h2_radius5.csv']);

    writetable(out_mat, out_path)
end






