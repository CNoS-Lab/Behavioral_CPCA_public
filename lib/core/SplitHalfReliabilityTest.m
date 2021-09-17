function predictStruct = SplitHalfReliabilityTest(G_CPCA, nIter, load_cut_off, pred_cut_off)

Z = G_CPCA.Z;
G = G_CPCA.G;
n = G_CPCA.n;

predictStruct.nIter = nIter;
predictStruct.load_cut_off = load_cut_off;
predictStruct.pred_cut_off = pred_cut_off;
predictStruct.PCountp = zeros(size(G_CPCA.PCorr));
predictStruct.PCountn = zeros(size(G_CPCA.PCorr));
predictStruct.CompCount = zeros(1,size(G_CPCA.PCorr,2));
predictStruct.nos = 0;

for ii = 1:nIter
   cv_idx(:,ii) = crossvalind('KFold',size(Z,1),2); 
end

parfor ii = 1:nIter
    
    G_CPCA_CV_1 = G_CPCA_Analysis(Z(cv_idx(:,ii)~=1,:), G(cv_idx(:,ii)~=1,:), n, false);
    G_CPCA_CV_2 = G_CPCA_Analysis(Z(cv_idx(:,ii)~=2,:), G(cv_idx(:,ii)~=2,:), n, false);
    
    matchStruct(ii) = matchComp(G_CPCA, [G_CPCA_CV_1,G_CPCA_CV_2], load_cut_off, 'orthogonal');    
end

predictStruct = countPredictor(predictStruct, matchStruct);

predictStruct.CompReliabilityProportion = predictStruct.CompCount./predictStruct.nos;
for comp = 1:n
    if predictStruct.CompReliabilityProportion(comp) >= G_CPCA.params.comp_rel_prop_cutoff
        predictStruct.PredReliabilityProportion_p(:,comp) = (predictStruct.PCountp(:,comp))./predictStruct.CompCount(comp);
        predictStruct.PredReliabilityProportion_n(:,comp) = (predictStruct.PCountn(:,comp))./predictStruct.CompCount(comp);
    else
        predictStruct.PredReliabilityProportion_p(:,comp) = zeros(size(predictStruct.PCountp(:,comp)));
        predictStruct.PredReliabilityProportion_n(:,comp) = zeros(size(predictStruct.PCountn(:,comp)));
    end
end

end

function outStruct = matchComp(G_CPCA, G_CPCA_CV, load_cut_off, type)

a = G_CPCA_CV(1).loadings_VD_sqrtN_GH;
b = G_CPCA_CV(2).loadings_VD_sqrtN_GH;
x = G_CPCA.loadings_VD_sqrtN_GH;

%% Procrustes rotation

[A, Ta] = procrustes_todd(a, x, type);
[B, Tb] = procrustes_todd(b, x, type);

G_CPCA_CV(1).scores_NU_GH = G_CPCA_CV(1).scores_NU_GH*Ta;
G_CPCA_CV(1).PCorr = corr(G_CPCA_CV(1).scores_NU_GH, G_CPCA_CV(1).G)';

G_CPCA_CV(2).scores_NU_GH = G_CPCA_CV(2).scores_NU_GH*Tb;
G_CPCA_CV(2).PCorr = corr(G_CPCA_CV(2).scores_NU_GH, G_CPCA_CV(2).G)';

%%

r_mat = corr(A, B);
r_mat = diag(r_mat);
r = find(abs(r_mat) > load_cut_off);
idx_vals = [r, r];

mean_loadings = nan(size(x));

for jj = 1:length(r)
    mean_loadings(:,r(jj)) = (A(:,r(jj)) + B(:,r(jj)))/2;
end

r_mat2 = corr(x, mean_loadings);
r_mat2 = diag(r_mat2);
r2 = find(abs(r_mat2) > load_cut_off);

idx_vals2 = [r2, r2];

outStruct.idx_vals = idx_vals;
outStruct.idx_vals2 = idx_vals2;
outStruct.mean_loadings = mean_loadings;
outStruct.predictor_CV{1} = NaN(size(G_CPCA.PCorr));
outStruct.predictor_CV{2} = NaN(size(G_CPCA.PCorr));
outStruct.SampleCorrMat = corr(A, B);

for jj = 1:length(r2)
    outStruct.predictor_CV{1}(:,r2(jj)) = G_CPCA_CV(1).PCorr(:,r2(jj));    
    outStruct.predictor_CV{2}(:,r2(jj)) = G_CPCA_CV(2).PCorr(:,r2(jj));
end

end

function OutStruct = countPredictor(InStruct, matchStruct)

    OutStruct = InStruct;
    OutStruct.nos = length(matchStruct);
    cut_off = InStruct.pred_cut_off;
    
    pred_CV = cat(1,matchStruct.predictor_CV);
    
    pred_CV1 = cat(3,pred_CV{:,1});
    pred_CV2 = cat(3,pred_CV{:,2});
    
    
    idx = squeeze(~isnan(nanmean(cat(1,pred_CV1,pred_CV1),1)));
    OutStruct.CompCount = sum(idx,2)';
    
    pred_bool_p = pred_CV1 > cut_off & pred_CV2 > cut_off;
    OutStruct.PCountp = sum(pred_bool_p,3);
    pred_bool_n = pred_CV1 < -1*cut_off & pred_CV2 < -1*cut_off;
    OutStruct.PCountn = sum(pred_bool_n,3);
    
end

