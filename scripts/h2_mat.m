function [h2, p_perm, jack_se] = h2_mat(P, K, X, n_perm, F)
%%
%
% Heritability estimation from phenotypic and genetic similary mattrices
%
% Input:
% P: an Nsubj x Nsubj phenotypic similarity matrix
% K: an Nsubj x Nsubj genetic similarity matrix
% X: an Nsubj x Ncov matrix of covariates
% n_perm: number of permutations; set Nperm = 0 if permutation inference is not needed
%
% Output:
% h2: heritability estimate
% p_perm: nonparametric permutation p-value; if Nperm = 0, PermPval = NaN
%
%%

if nargin < 5
    disp('Not running block jack-knife')
    F = false;
end


K_orig=K;
X_orig=X;
P_orig=P;
F_orig=F;

[n_subj, n_cov] = size(X);
% -----
P0 = eye(n_subj) - X/(X'*X)*X';
[U,~,~] = svd(P0); U = U(:,1:n_subj-n_cov);
P = U'*P*U; K = U'*K*U;
n_subj = n_subj-n_cov;
%%
kappa = trace(K^2)/n_subj;
tau = trace(K)/n_subj;
vK = n_subj*(kappa-tau^2);

Qg = K-tau*eye(n_subj);
Qe = kappa*eye(n_subj)-tau*K;

tg = trace(Qg*P)/vK;
te = trace(Qe*P)/vK;
tp = tg+te;

h2 = max(min(tg/tp,1),0);
%%
if n_perm == 0
    p_perm = NaN;
    return
end

h2_perm = zeros(1,n_perm);
for s = 1:n_perm
    disp(['----- Permutation-', num2str(s), ' -----'])

    if s == 1
        K_perm = K;
    else
        subj_perm = randperm(n_subj);
        K_perm = K(subj_perm,subj_perm);
    end

    Qg_perm = K_perm-tau*eye(n_subj);
    Qe_perm = kappa*eye(n_subj)-tau*K_perm;

    tg_perm = trace(Qg_perm*P)/vK;
    te_perm = trace(Qe_perm*P)/vK;

    h2_perm(s) = max(min(tg_perm/(tg_perm+te_perm),1),0);
end

p_perm = sum(h2_perm>=h2)/n_perm;

jack_se = NaN;

if F ~= false;
    % block jack-knife
    fam_arr  = unique(F_orig);
    n_fam    = length(fam_arr)
    h2_jack  = zeros(1,n_fam);
    for f = 1:n_fam
        disp(['----- Block Jacknife-', num2str(f), ' -----'])

        K_jack = K_orig;
        P_jack = P_orig;
        X_jack = X_orig;

        use_idxs  = find(F_orig ~= f);
        subs_jack = length(use_idxs);

        % remove the family from the data
        K_jack = K_jack(use_idxs, use_idxs);
        X_jack = X_jack(use_idxs,:);
        P_jack = P_jack(use_idxs, use_idxs);

        [n_subj, n_cov] = size(X_jack);

        % -----
        P0      = eye(n_subj) - X_jack/(X_jack'*X_jack)*X_jack';
        [U,~,~] = svd(P0);
        U = U(:,1:n_subj-n_cov);
        P_jack = U'*P_jack*U;
        K_jack = U'*K_jack*U;
        n_subj = n_subj-n_cov;
        %%
        kappa = trace(K_jack^2)/n_subj;
        tau = trace(K_jack)/n_subj;
        vK = n_subj*(kappa-tau^2);

        Qg = K_jack-tau*eye(n_subj);
        Qe = kappa*eye(n_subj)-tau*K_jack;

        tg = trace(Qg*P_jack)/vK;
        te = trace(Qe*P_jack)/vK;
        tp = tg+te;

        h2_jack(f) = max(min(tg/tp,1),0);
        disp(h2_jack(f))
    end

    % variance of overall heritability, as in Ge 2017 (https://www.pnas.org/content/pnas/114/21/5521.full.pdf)
    A = (n_fam-1)/n_fam;
    B = sum((h2_jack - sum(h2_jack/n_fam)).^2)
    h2_var = A*B
    jack_se = sqrt(h2_var)
end

%%

