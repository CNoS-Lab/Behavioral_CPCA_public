function [VarG, T] = buildVarianceTable(DG, DE, D, n)

V = (diag(D).^2); 
VG = (diag(DG).^2); 
VE = (diag(DE).^2);
VarG = sum(VG)/sum(V)*100;

Source = {'Total','Total (in %)','Between', '    Percent of Total','    Percent of between', 'Residue', '    Percent of total','    Percent of residue'};

VarNames = cell(1,n+2);
table_mat = cell(1,n+2);

Total_Variance = categorical(string(num2str(estimateRow(V, VG, VE)','%.3f')));
table_mat{1} = Total_Variance;
VarNames{1} = 'Total Variance';
for ii = 1:n
    table_mat{ii+1} = categorical(string(num2str(estimateRow(V, VG, VE, ii)','%.3f')));
    VarNames{ii+1} = num2str(ii);
end
table_mat{n+2} = categorical(string(num2str(estimateRow(V, VG, VE, 1:n)','%.3f')));
VarNames{n+2} = sprintf('Sum of %d comps',n);

temp = sprintfc('table_mat{%d},',1:length(table_mat));
table_str = cat(2,temp{:});

eval(sprintf('T = table(%s ''RowNames'', Source'', ''VariableNames'', VarNames);',table_str));

end

function [Row] = estimateRow(V, VG, VE, i)

if nargin > 3
    Row = [sum(V(i)), sum(V(i))/sum(V)*100, sum(VG(i)), sum(VG(i))/sum(V)*100, sum(VG(i))/sum(VG)*100, sum(VE(i)), sum(VE(i))/sum(V)*100, sum(VE(i))/sum(VE)*100];
else
    Row = [sum(V), sum(V)/sum(V)*100, sum(VG), sum(VG)/sum(V)*100, sum(VG)/sum(VG)*100, sum(VE), sum(VE)/sum(V)*100, sum(VE)/sum(VE)*100];
end

end