function OutStruct = runCPCA(X, G, H, n, varmaxFlag)

% X - rows are observations (Points), columns are features (variables)
% X should be zscored along the columns
% G - Predictor variable matrix along 1st dim -rows
% H - Predictor variable matrix along 2nd dim - cols
% n - Select top n components

if nargin < 5
    varmaxFlag = false;
end

N = size(X,1);
if ~all(size(G) == 1) && all(size(H) == 1)
    C = pinv(G'*G)*G'*X; % Linear regression
    S = G*C; E = X - S;
elseif all(size(G) == 1) && ~all(size(H) == 1)
    B = X*H*pinv(H'*H); % Linear regression
    S = B*H'; E = X - S;
elseif all(size(G) == 1) && all(size(H) == 1)
    S = X; E = X - S;
else
    return;
end

temp = corrcoef(zscore(S));
temp(isnan(temp)) = 0;
[~, OutStruct.EigValGH] = eig(temp);
OutStruct.EigValGH = sort(diag(OutStruct.EigValGH),'descend');
[OutStruct.UGH, OutStruct.DGH, OutStruct.VGH] = svd(S,'econ'); % SVD for GC
[OutStruct.loadings_VD_sqrtN_GH, OutStruct.scores_NU_GH, OutStruct.loadings_VD_sqrtN_GH_inv, OutStruct.T] = findLoadingsScores(OutStruct.UGH, OutStruct.DGH, OutStruct.VGH, N, n, varmaxFlag);

if ~all(size(G) == 1) && all(size(H) == 1)
    temp = corrcoef(zscore(E));
    temp(isnan(temp)) = 0;
    [~, OutStruct.EigValE] = eig(temp);
    OutStruct.EigValE = sort(diag(OutStruct.EigValE),'descend');
    [OutStruct.UE, OutStruct.DE, OutStruct.VE] = svd(E,'econ'); % SVD for E
end

OutStruct.E = E;
OutStruct.S = S;

end


function [loadings_VD_sqrtN, scores_NU, loadings_VD_sqrtN_inv, T] = findLoadingsScores(U, D, V, N, n, varmaxFlag)

loadings_VD_sqrtN = V*diag(diag(D))/sqrt(N-1);
scores_NU = sqrt(N-1)*U;

loadings_VD_sqrtN_inv = pinv(loadings_VD_sqrtN)';

loadings_VD_sqrtN = loadings_VD_sqrtN(:,1:n);
loadings_VD_sqrtN_inv = loadings_VD_sqrtN_inv(:,1:n);
scores_NU = scores_NU(:,1:n);

T = NaN;

% Perform varimax rotation
if varmaxFlag && n>1
    [loadings_VD_sqrtN, T]= varimkn_beh(loadings_VD_sqrtN);
    loadings_VD_sqrtN_inv = loadings_VD_sqrtN_inv*T;
    scores_NU = scores_NU*T;
end

end