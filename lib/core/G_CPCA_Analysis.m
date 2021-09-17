function [G_CPCA] = G_CPCA_Analysis(Z, G, n, varmaxFlag)

switch nargin
    case 5
        G_test = G;
    case 4
        Z_test = Z;
        G_test = G;
    case 3
        varmaxFlag = false;
        Z_test = Z;
        G_test = G;
    case 2
        n = [];
        varmaxFlag = false;
        Z_test = Z;
        G_test = G;
end

% Code to run CPCA using G matrix

Z = (Z - nanmean(Z,1))./nanstd(Z,[],1); % Z-score the data 
Z(isnan(Z)) = 0;

if isempty(G)
    G_CPCA = runCPCA(Z,1,1,rank(Z),varmaxFlag); % Uncontrained PCA
    return;
end

if ~isempty(n)
    G_CPCA = runCPCA(Z,G,1,n,varmaxFlag); % PCA model
else
    G_CPCA = runCPCA(Z,G,1,min(rank(Z),rank(G)),varmaxFlag); % CPCA model
end
[G_CPCA.PCorr, G_CPCA.PCorr_pval] = corr(G, G_CPCA.scores_NU_GH); % Predictor correlations

% Save the Z and G matrices in the output struct

G_CPCA.G = G;
G_CPCA.Z = Z;
G_CPCA.n = n;

G_CPCA.varG = [];
G_CPCA.predictStruct = [];

end