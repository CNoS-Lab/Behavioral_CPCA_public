function [G_CPCA] = BehCPCAgenerateReport(rpt_obj, temp, Z_index, G_index, Sub_index, params)

import mlreportgen.dom.*;
import mlreportgen.report.*;

%% Get Z matrix

Z = cell2mat(temp(Sub_index,Z_index));
Zlabel = temp(1,Z_index);

%% Get G Matrix

G = cell2mat(temp(Sub_index,G_index));
% G = [G ones(size(G,1),1)];
Glabel = temp(1,G_index);
% Glabel = [Glabel,{'Ones'}];

[Z, G] = removeNaN(Z, G);

Zdim = size(Z,2);
Gdim = size(G,2);
n_sub = size(Z,1);

if isempty(Z) || isempty(G) 
    P = Text('This analysis is not possible, because there is no data in either Z or G.');
    add(rpt_obj, P);
    return;
end

%% Perform initial CPCA

UnCon_PCA = G_CPCA_Analysis(Z,[]);
n_temp = min(rank(Z),rank(G));
G_CPCA_temp = G_CPCA_Analysis(Z,G,n_temp,params.varimaxFlag);

figure('Position',[0 0 1280 720]); % Scree plot
D = G_CPCA_temp.DGH;
D = diag(D).^2;
D = D./sum(D)*100;
k = length(D);
plot(1:k,D(1:k),'-O');
xlabel('No. of Components');
ylabel('Percentage of variance');
title('Scree Plot');
set(findall(gcf,'-property','FontSize'),'FontSize',15);

n = input('Enter the number of components: ');
close all;

params.load_cut_off = est_corr_thresholds(Zdim, params.n_bootstrap, params.p_val);
params.pred_cut_off = est_corr_thresholds(floor(n_sub/2), params.n_bootstrap, params.p_val);

%% Initial page for section

H = Heading4(Text('CPCA Analysis')); add(rpt_obj, H);
T = Text(sprintf('n = %d \n', n_sub)); add(rpt_obj, T);
T = Text(sprintf('%s \n', params.page_header)); add(rpt_obj, T);
T = Text(sprintf('Z-dim = %d \n', Zdim)); add(rpt_obj, T);
T = Text(sprintf('G-dim = %d \n', Gdim)); add(rpt_obj, T);
T = Text(sprintf('Estimated component loading correlation threshold = %.2f \n', params.load_cut_off)); add(rpt_obj, T);
T = Text(sprintf('Estimated predictor loading correlation threshold = %.2f \n\n', params.pred_cut_off)); add(rpt_obj, T);

br = PageBreak();
add(rpt_obj, br);

%% Redo the CPCA with n components

G_CPCA = G_CPCA_Analysis(Z,G,n,params.varimaxFlag);
G_CPCA.params = params;

%% Perform slpit-half reliability test

rng('default');
predictStruct = SplitHalfReliabilityTest(G_CPCA, params.nIter, params.load_cut_off, params.pred_cut_off);

%% Plot scree plot

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text('Scree plot')); add(rpt_obj, H);

fig = figure('Position',[0 0 1280 720],'visible','off'); % Scree plot
D = G_CPCA_temp.DGH;
D = diag(D).^2;
D = D./sum(D)*100;
k = length(D);
plot(1:k,D(1:k),'-O');
xlabel('No. of Components');
ylabel('Percentage of variance');
title('Scree Plot');
set(findall(gcf,'-property','FontSize'),'FontSize',15);

f(1) = Figure(fig); add(rpt_obj, f(1));

T = Text(sprintf('%d number of components were selected.',n));
add(rpt_obj, T);

br = PageBreak();
add(rpt_obj, br);

%% Build variance table

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text('Variance Table')); add(rpt_obj, H);

[varG, tbl] = buildVarianceTable(G_CPCA.DGH, G_CPCA.DE, UnCon_PCA.DGH, G_CPCA.n);

Tbl = MATLABTable(tbl);
Tbl.TableEntriesStyle = {Height('20pt')};

add(rpt_obj, Tbl);
br = PageBreak();
add(rpt_obj, br);

%% Variance significance test

H = Heading4(Text('Variance significance test')); add(rpt_obj, H);

[fig, G_CPCA.p] = Sig_Perm_Code_CP(G_CPCA.Z, G_CPCA.G, varG);
set(fig, 'Position',[0 0 1280 720]);     
set(findall(gcf,'-property','FontSize'),'FontSize',11);
f(end+1) = Figure(fig);
add(rpt_obj, f(end));

br = PageBreak();
add(rpt_obj, br);

%% Plot component loadings

for ii = 1:G_CPCA.n
    xStr{ii} = sprintf('C%d',ii);
end

H = Heading4(Text('Component Loadings'));
add(rpt_obj, H);

nz = 10;
temp = G_CPCA.loadings_VD_sqrtN_GH;
n_plots = ceil(size(temp,1)/nz);

for ii = 1:n_plots
    T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
    fig = figure('visible','off');
    if ii == n_plots
        h = plotLoadings(temp(nz*(ii-1)+1:end,:),xStr,Zlabel(nz*(ii-1)+1:end));
    else
        h = plotLoadings(temp(nz*(ii-1)+1:nz*ii,:),xStr,Zlabel(nz*(ii-1)+1:nz*ii));
    end
    h.Colormap = flipud(cbrewer('div','RdBu',64,'spline'));
    h.ColorLimits = [-1*max(abs(temp(:))) max(abs(temp(:)))];
    if n_plots == 1
        title(sprintf('Component Loadings'));
    else
        title(sprintf('Component Loadings plot no. %d of %d',ii,n_plots));
    end
    set(fig, 'Position',[0 0 700 900]);
    set(findall(gcf,'-property','FontSize'),'FontSize',12);
    f(end+1) = Figure(fig);
    add(rpt_obj, f(end));
    br = PageBreak();
    add(rpt_obj, br);
end

%% Plot predictor loadings

H = Heading4(Text('Predictor Loadings'));
add(rpt_obj, H);

nz = 10;
temp = G_CPCA.PCorr;
n_plots = ceil(size(temp,1)/nz);

for ii = 1:n_plots
    T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
    fig = figure('visible','off');
    if ii == n_plots
        h = plotLoadings(temp(nz*(ii-1)+1:end,:),xStr,Glabel(nz*(ii-1)+1:end));
    else
        h = plotLoadings(temp(nz*(ii-1)+1:nz*ii,:),xStr,Glabel(nz*(ii-1)+1:nz*ii));
    end
    h.Colormap = flipud(cbrewer('div','RdBu',64,'spline'));
    h.ColorLimits = [-1 1];
    if n_plots == 1
        title(sprintf('Predictor Loadings'));
    else
        title(sprintf('Predictor Loadings plot no. %d of %d',ii,n_plots));
    end
    set(fig, 'Position',[0 0 700 900]);
    set(findall(gcf,'-property','FontSize'),'FontSize',9);
    f(end+1) = Figure(fig);
    add(rpt_obj, f(end));
    br = PageBreak();
    add(rpt_obj, br);
end

%% Plot component reliability proportions

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text('Component reliability proportions'));
add(rpt_obj, H);

fig = figure('visible','off');
temp = predictStruct.CompReliabilityProportion;
h = plotLoadings(temp,xStr,'');
h.Colormap = cbrewer('seq','Reds',64,'spline');
h.ColorLimits = [0 1];
title('Component reliability proportions');
set(fig, 'Position',[0 0 1440 300]);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
f(end+1) = Figure(fig);
add(rpt_obj, f(end));

br = PageBreak();
add(rpt_obj, br);

comp_rel(1:n) = predictStruct.CompReliabilityProportion;

%% Predictor reliability threshold

T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
H = Heading4(Text('Predictor loading reliability proportion threshold'));
add(rpt_obj, H);

varimaxFlag = params.varimaxFlag;
nIter = params.nIter;
load_cut_off = params.load_cut_off;
pred_cut_off = params.pred_cut_off;

fprintf('Iter running = 0000');
for ii = 1:params.nIter
    G_temp = G(randperm(size(G,1)),:); Z_temp = Z;
    G_CPCA_val_temp = G_CPCA_Analysis(Z_temp,G_temp,n,varimaxFlag); % Perform CPCA 
    G_CPCA_val_temp.params = params;
    G_CPCA_val_predictStruct(ii) = SplitHalfReliabilityTest(G_CPCA_val_temp, nIter, load_cut_off, pred_cut_off);
    disp(ii);
end

temp = G_CPCA_val_predictStruct;
predRelp = cat(3,temp.PredReliabilityProportion_p);
predReln = cat(3,temp.PredReliabilityProportion_n);
predRel = abs(predRelp - predReln);

fig = figure('visible','off');
histogram(predRel);
hold on;
clear g;
g(1) = vline(quantile(predRel(:), 1-0.05)); g(1).LineWidth = 2; g(1).Color = 'r';
g(2) = vline(quantile(predRel(:), 1-0.01)); g(2).LineWidth = 2; g(2).Color = 'm';
g(3) = vline(quantile(predRel(:), 1-0.005)); g(3).LineWidth = 2; g(3).Color = 'g';
g(4) = vline(quantile(predRel(:), 1-0.001)); g(4).LineWidth = 2; g(4).Color = 'b';
legend(g,{'p = 0.05','p = 0.01','p = 0.005','p = 0.001'});
title('Predictor loading reliability proportion threshold');
set(fig, 'Position',[0 0 1280 720]);     
set(findall(gcf,'-property','FontSize'),'FontSize',15);

f(end+1) = Figure(fig);
add(rpt_obj, f(end));

G_CPCA.predictStruct = predictStruct;
G_CPCA.predrel_cutoff = quantile(predRel(:), 1-params.pred_rel_p_val);
G_CPCA.pred_rel_p = G_CPCA.predictStruct.PredReliabilityProportion_p > G_CPCA.predrel_cutoff;
G_CPCA.pred_rel_n = G_CPCA.predictStruct.PredReliabilityProportion_n > G_CPCA.predrel_cutoff;

G_CPCA.Pval_pred_rel_p = get_p_val(predRel(:), G_CPCA.predictStruct.PredReliabilityProportion_p);
G_CPCA.Pval_pred_rel_n = get_p_val(predRel(:), G_CPCA.predictStruct.PredReliabilityProportion_n);

T = Text(sprintf('Estimated Predictor loading reliability proportion threshold = %4f (at p <= %.4f) \n', G_CPCA.predrel_cutoff,params.pred_rel_p_val)); add(rpt_obj, T);

br = PageBreak();
add(rpt_obj, br);

%% Plot positive predictor reliability proportions

for ii = 1:length(xStr)
    xStr{ii} = strtrim(sprintf('%s \\newline %.2f',xStr{ii},...
        predictStruct.CompReliabilityProportion(ii)));
end

H = Heading4(Text(sprintf('Positive predictor loading reliability proportions (Thresh = %.4f)',G_CPCA.predrel_cutoff)));
add(rpt_obj, H);

nz = 10;
temp = predictStruct.PredReliabilityProportion_p;
n_plots = ceil(size(temp,1)/nz);

for ii = 1:n_plots
    T = Text(sprintf('Chapter: %s. Predictor loading reliability proportion threshold = %.4f\n', params.page_header,G_CPCA.predrel_cutoff)); add(rpt_obj, T);
    fig = figure('visible','off');
    if ii == n_plots
        h = plotLoadings2(temp(nz*(ii-1)+1:end,:),xStr,Glabel(nz*(ii-1)+1:end),G_CPCA.predrel_cutoff);
    else
        h = plotLoadings2(temp(nz*(ii-1)+1:nz*ii,:),xStr,Glabel(nz*(ii-1)+1:nz*ii),G_CPCA.predrel_cutoff);
    end
    colormap hot;
    h.Colormap = cbrewer('seq','Reds',64,'spline');
    h.CLim = [0 1]; %h.ColorLimits = [0 1];
    if n_plots == 1
        title(sprintf('Positive predictor loading reliability proportions'));
    else
        title(sprintf('Positive predictor loading reliability proportions plot no. %d of %d',ii,n_plots));
    end
    set(fig, 'Position',[0 0 700 900]);
    set(findall(gcf,'-property','FontSize'),'FontSize',9);
    f(end+1) = Figure(fig);
    add(rpt_obj, f(end));
    br = PageBreak();
    add(rpt_obj, br);
end

%% Plot negative predictor reliability proportions

H = Heading4(Text(sprintf('Negative predictor loading reliability proportions (Thresh = %.4f)',G_CPCA.predrel_cutoff)));
add(rpt_obj, H);

nz = 10;
temp = predictStruct.PredReliabilityProportion_n;
n_plots = ceil(size(temp,1)/nz);

for ii = 1:n_plots
    T = Text(sprintf('Chapter: %s. Predictor loading reliability proportion threshold = %.4f\n', params.page_header,G_CPCA.predrel_cutoff)); add(rpt_obj, T);
    fig = figure('visible','off');
    if ii == n_plots
        h = plotLoadings2(temp(nz*(ii-1)+1:end,:),xStr,Glabel(nz*(ii-1)+1:end),G_CPCA.predrel_cutoff);
    else
        h = plotLoadings2(temp(nz*(ii-1)+1:nz*ii,:),xStr,Glabel(nz*(ii-1)+1:nz*ii),G_CPCA.predrel_cutoff);
    end
    h.Colormap = cbrewer('seq','Reds',64,'spline');
    h.CLim = [0 1]; %h.ColorLimits = [0 1];
    if n_plots == 1
        title(sprintf('Negative predictor loading reliability proportions'));
    else
        title(sprintf('Negative predictor loading reliability proportions plot no. %d of %d',ii,n_plots));
    end
    set(fig, 'Position',[0 0 700 900]);
    set(findall(gcf,'-property','FontSize'),'FontSize',9);
    f(end+1) = Figure(fig);
    add(rpt_obj, f(end));
    br = PageBreak();
    add(rpt_obj, br);
end

%% Dominant Component Loadings

H = Heading4(Text(sprintf('Dominant component loadings'))); add(rpt_obj, H);
T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);

Res = SelectCompLoadings(G_CPCA);
G_CPCA.DomCompRes = Res;

n_rel_items_p = cat(1,Res.pred_item_p_found);
n_rel_items_n = cat(1,Res.pred_item_n_found);

pc_rel_items = (n_rel_items_p + n_rel_items_n)./(sum(G_CPCA.pred_rel_p + G_CPCA.pred_rel_n,1)');
[pc_rel_items, var_sort_idx] = sort(pc_rel_items,2);

n_feat = Res(1).n_feat;

fig = figure('visible','on');
plot(n_feat, pc_rel_items,'*-','LineWidth',2);
legend(sprintfc('C%d',1:n),'Location','bestoutside');
title('Dominant component loadings (Pos + Neg)');
ylabel('% reliable predictor items');
xlabel('Variable dominance rank #');
set(fig, 'Position',[0 0 1280 720]);     
set(findall(gcf,'-property','FontSize'),'FontSize',15);

fig = figure('visible','off');
plot(n_feat, pc_rel_items,'*-','LineWidth',2);
legend(sprintfc('C%d',1:n),'Location','bestoutside');
title('Dominant component loadings (Pos + Neg)');
ylabel('% reliable predictor items');
xlabel('Variable dominance rank #');
set(fig, 'Position',[0 0 1280 720]);     
set(findall(gcf,'-property','FontSize'),'FontSize',15);

f(end+1) = Figure(fig);
add(rpt_obj, f(end));
br = PageBreak();
add(rpt_obj, br);

%% Dominant loading plot for individual items

comp_rel = find(any(G_CPCA.pred_rel_p | G_CPCA.pred_rel_n, 1));
avg_rel_prop = NaN(length(n_feat),n);
sort_idx = NaN(length(n_feat),n);

for comp = comp_rel
    
    H = Heading4(Text(sprintf('Dominant component loadings - Comp %d',comp))); add(rpt_obj, H);
    T = Text(sprintf('Chapter: %s \n', params.page_header)); add(rpt_obj, T);
    
    if any(G_CPCA.pred_rel_p(:,comp) | G_CPCA.pred_rel_n(:,comp))
        
        idx = [G_CPCA.pred_rel_p(:,comp);G_CPCA.pred_rel_n(:,comp)];
        data_p = cat(2,Res(comp).pred_rel_p{:});
        data_p(isnan(data_p)) = 0;
        data_n = cat(2,Res(comp).pred_rel_n{:});
        data_n(isnan(data_n)) = 0;
        
        data = [data_p;-1*data_n];
        
        data_idx = data(idx,:);
        [avg_rel_prop(:,comp), sort_idx(:,comp)] = sort(nanmean(abs(data_idx), 1));
        xStr = sprintfc('%d',sort_idx(:,comp));
        
        Glabel_temp = [strcat(Glabel,'(+)'),strcat(Glabel,'(-)')];
        figure;
        h = plotLoadings(data(idx,sort_idx(:,comp)),xStr,Glabel_temp(idx));
        h.Colormap = flipud(cbrewer('div','RdBu',64,'spline'));
        h.ColorLimits = [-1 1];
        title(sprintf('Comp %d: Pred. load. rel. prop. (Pos + Neg)',comp));
        xlabel('Z Variable Removed #');
%         set(fig, 'Position',[0 0 1440 250]);
        set(findall(gcf,'-property','FontSize'),'FontSize',12);
        
        fig = figure('visible','off');
        h = plotLoadings(data(idx,sort_idx(:,comp)),xStr,Glabel_temp(idx));
        h.Colormap = flipud(cbrewer('div','RdBu',64,'spline'));
        h.ColorLimits = [-1 1];
        title(sprintf('Comp %d: Pred. load. rel. prop. (Pos + Neg)',comp));
        xlabel('Z Variable Removed #');
%         set(fig, 'Position',[0 0 1440 250]);
        set(findall(gcf,'-property','FontSize'),'FontSize',12);
        
        f(end+1) = Figure(fig);
        add(rpt_obj, f(end));
        
        T = Text([sprintf('Xlabel for comp %d = ',comp),sprintf('%d, ', sort_idx(1:end-1,comp)),sprintf('%d',sort_idx(end,comp))]);
        add(rpt_obj, T);
        
    end
    
    br = PageBreak();
    add(rpt_obj, br);
    
end

%% Dominant loading plot for individual items 2

n_feat = Res(1).n_feat;

fig = figure('visible','on');
plot(n_feat, avg_rel_prop,'*-','LineWidth',2);
legend(sprintfc('C%d',1:n),'Location','bestoutside');
title('Dominant component loadings (Pos + Neg)');
ylabel('Avg. Pred. load. rel. prop.');
xlabel('Variable dominance rank #');
set(fig, 'Position',[0 0 1280 720]);     
set(findall(gcf,'-property','FontSize'),'FontSize',15);

for comp = comp_rel
    fig = figure('visible','off');
    h = plot(n_feat, avg_rel_prop,'*-','LineWidth',2);
    set(h(setdiff(1:n,comp)),'Visible','off');
    legend(sprintfc('C%d',1:n),'Location','bestoutside');
    title(sprintf('Dominant component loadings (Pos + Neg) - Comp %d',comp));
    ylabel('Pred. load. rel. prop.');
    xlabel('Variable dominance rank #');
    set(fig, 'Position',[0 0 1280 720]);
    set(findall(gcf,'-property','FontSize'),'FontSize',15);
    
    f(end+1) = Figure(fig);
    add(rpt_obj, f(end));
    
    T = Text([sprintf('Xlabel for comp %d = ',comp),sprintf('%d, ', sort_idx(1:end-1,comp)),sprintf('%d',sort_idx(end,comp))]); 
    add(rpt_obj, T);
    
    br = PageBreak();
    add(rpt_obj, br);
end

%% Final summary

Add_more_summary_flag = true; 
summ_no = 1;

while Add_more_summary_flag
    
fprintf('Summary %d: \n',summ_no);

for c = 1:length(comp_rel)
    G_CPCA.Summary(summ_no,c).CompNo = comp_rel(c);
    G_CPCA.Summary(summ_no,c).CompRelProp = G_CPCA.predictStruct.CompReliabilityProportion(comp_rel(c));
    
    opt_comp_no = input(sprintf('Enter the number of dominant component loadings for comp %d: ',comp_rel(c))); 
    opt_comp_load_idx = sort_idx(1:opt_comp_no,comp_rel(c));
    
    for ii = 1:length(opt_comp_load_idx)
        G_CPCA.Summary(summ_no,c).CompLoadings(ii).Name = Zlabel{opt_comp_load_idx(ii)};
        G_CPCA.Summary(summ_no,c).CompLoadings(ii).Value = G_CPCA.loadings_VD_sqrtN_GH(opt_comp_load_idx(ii),comp_rel(c));
    end
    
    TempPredictStruct = SplitHalfReliabilityTest_SelFeat_no_regress(G_CPCA, nIter, load_cut_off, pred_cut_off, opt_comp_load_idx);
    G_CPCA.SelPredictStruct = TempPredictStruct;
    
    rel_pred_p_idx = find(TempPredictStruct.PredReliabilityProportion_p(:,comp_rel(c)) > G_CPCA.predrel_cutoff & G_CPCA.pred_rel_p(:,comp_rel(c)));
    
    for ii = 1:length(rel_pred_p_idx)
        G_CPCA.Summary(summ_no,c).PredLoadingsP(ii).Name = Glabel{rel_pred_p_idx(ii)};
        G_CPCA.Summary(summ_no,c).PredLoadingsP(ii).Value = G_CPCA.PCorr(rel_pred_p_idx(ii),comp_rel(c));
        G_CPCA.Summary(summ_no,c).PredLoadingsP(ii).Rel = G_CPCA.predictStruct.PredReliabilityProportion_p(rel_pred_p_idx(ii),comp_rel(c));
    end
    
    rel_pred_n_idx = find(TempPredictStruct.PredReliabilityProportion_n(:,comp_rel(c)) > G_CPCA.predrel_cutoff & G_CPCA.pred_rel_n(:,comp_rel(c)));
    
    for ii = 1:length(rel_pred_n_idx)
        G_CPCA.Summary(summ_no,c).PredLoadingsN(ii).Name = Glabel{rel_pred_n_idx(ii)};
        G_CPCA.Summary(summ_no,c).PredLoadingsN(ii).Value = G_CPCA.PCorr(rel_pred_n_idx(ii),comp_rel(c));
        G_CPCA.Summary(summ_no,c).PredLoadingsN(ii).Rel = G_CPCA.predictStruct.PredReliabilityProportion_n(rel_pred_n_idx(ii),comp_rel(c));
    end
    G_CPCA.Summary(summ_no,c).rel_item_pc = (length(rel_pred_p_idx)+length(rel_pred_n_idx))./(sum(G_CPCA.pred_rel_p(:,comp_rel(c)) + G_CPCA.pred_rel_n(:,comp_rel(c)),1));
    G_CPCA.Summary(summ_no,c).tot_rel_item = sum(G_CPCA.pred_rel_p(:,comp_rel(c)) + G_CPCA.pred_rel_n(:,comp_rel(c)),1);
    G_CPCA.Summary(summ_no,c).comp_load_rel_item = length(rel_pred_p_idx)+length(rel_pred_n_idx);
end

Add_more_summary_flag = input('Would you like to select different set of component loadings? 1 - Yes, 0 - No: ');
summ_no = summ_no + 1;

end

end