% Produces true if the percent of variance predictable from G is significant at level p, false otherwise
% Require G and Znorm (Z with standardized columns) to be already loaded in workspace

	% - randperm returns a permutation of G with rearranged rows, but does not rearrange the values in those rows
	% - p is significance level
	% - Num is the number of permutations we would like to compute
	% - count is the number of times perc_pred of the permutated G is greater than perc_pred of the original G, starts at 0

%Also outputs mean % of variance in Z predictable from G and standard deviation. The % of variance in Z predictable from G from all permutations is plotted in a histogram.
    
% =====================================================================================================================

% Inputs:

function [fig, significance] = Sig_Perm_Code_CP(Znorm, G, perc_pred_original)

p = 0.05;
Num = 10000;

% ======================================================================================================================

% Function:

%Setup:
perms=NaN(Num,1);
count = 0;
var_totalz=sum(var(Znorm));

%Permutations
for loop = 1:Num
	GP=G(randperm(size(G,1)),:); %Create list of permutations to rearrange rows of G

    % Add on 1s column
    g_const=[GP ones(size(G,1), 1)];
%     g_const=GP;
    
    %Regression to obtain C
    C = (pinv(g_const'* g_const))* g_const'* Znorm;
    GC = g_const * C;
    
    %SVD on GC
    [Ugc Dgc Vgc]=svd(GC, 'econ');
    
    %Variance Calculations
    var_totalgc=sum(diag(Dgc^2./(size(Znorm,1)-1))); %Variance in GC
    perc_pred_perm=(var_totalgc/var_totalz)*100;% percent of variance in Z predictable from G
	perms(loop)= perc_pred_perm; %Add the percent of variance in Z predictable from G to perms matrix

    %Increase count by 1 if the % variance in Z predictable by G (for this permutation of G) is greater than with the original G
    if perc_pred_perm > perc_pred_original
        count = 1 + count;
    end
%    disp(loop);
 end 

% Significant?
if (count / Num) > p

end
if (count / Num) <= p

end

%Outputs
significance = count / Num;
Means_perms = mean(perms);
STD_perms = std(perms);

fig = figure('visible','off');
histogram(perms)
xline(perc_pred_original,'r','LineWidth',2);
title(sprintf('Significance, p = %.4f',significance));

end


