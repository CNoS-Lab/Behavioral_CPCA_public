function Res = SelectCompLoadings(G_CPCA)

%comp_load = G_CPCA.loadings_VD_sqrtN_GH;
nIter = G_CPCA.predictStruct.nIter;
pred_cut_off = G_CPCA.predictStruct.pred_cut_off;
load_cut_off = G_CPCA.predictStruct.load_cut_off;
Z_dim = size(G_CPCA.Z,2);

pred_rel_p = G_CPCA.pred_rel_p;
pred_rel_n = G_CPCA.pred_rel_n;

Res = repmat(struct('pred_rel_p',cell(1,Z_dim),'pred_rel_n',cell(1,Z_dim),'pred_item_p_found',NaN(1,Z_dim),...
    'pred_item_n_found',NaN(1,Z_dim),'Sel_feat',cell(1,Z_dim),'n_feat',1:Z_dim,'least_idx',NaN(1,Z_dim)),[G_CPCA.n,1]);
Res = Res(:,1);

for comp = 1:G_CPCA.n
    
    rng('default');
    %comp_load_temp = comp_load(:,comp);
    
    Orig_pred_rel_map_p = pred_rel_p(:,comp);
    Orig_pred_rel_map_n = pred_rel_n(:,comp);
    
    if ~any(Orig_pred_rel_map_p + Orig_pred_rel_map_n)
        continue;
    end
    
    for ii = 1:Z_dim
        
        Rej_feat = ii;
        Sel_feat = setdiff(1:Z_dim, Rej_feat);
        
        predictStruct = SplitHalfReliabilityTest_SelFeat(G_CPCA, nIter, load_cut_off, pred_cut_off, Sel_feat);
        
        Res(comp).pred_rel_p{ii} = predictStruct.PredReliabilityProportion_p(:,comp);
        Res(comp).pred_rel_n{ii} = predictStruct.PredReliabilityProportion_n(:,comp);
        
        Res(comp).pred_item_p_found(ii) = sum(Orig_pred_rel_map_p & ...
            (Res(comp).pred_rel_p{ii} > G_CPCA.predrel_cutoff));
        Res(comp).pred_item_n_found(ii) = sum(Orig_pred_rel_map_n & ...
            (Res(comp).pred_rel_n{ii} > G_CPCA.predrel_cutoff));
        
        Res(comp).Sel_feat{ii} = Sel_feat;
        
    end
end

end