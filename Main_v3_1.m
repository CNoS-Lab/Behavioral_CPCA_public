%% This code performs 3 step CPCA on two sets of variables. 
%
% Make sure to
% 1. have all your data in one .csv file. 
% 2. encode all the missing data either using -99 or -88. 
% 3. Have or install the 'Bioinformatics' and 'the Statistics and Machine
% learning' toolbxes. 
%
%% Clearing workspace and adding necessary files to path

clear; clc; close all;
restoredefaultpath;
addpath(genpath('lib/'));
rng('default');

%% Stuff that you need to enter

data = readcell('Multisite_lab_data.csv'); % Insert your path/filename for your data here
Z_index_cell = {[14,17,20,32:40,31,47,48:51]}; % Enter the column numbers of the Z matrix (dependant variables)
G_index_cell = {[87:118,215:230]}; % Enter the column numbers of the G matrix (independant variables)
Sub_index_cell = {2:595}; % Enter the row numbers for subject groups on which CPCA analysis needs to be done as seen in .csv file

% Also change the Z, G, and Sub labels for each set that you enter above
% e.g. if your first set of Z variables are SSPI metrics, then change 'Z variable set 1' to 'SSPI'

Z_label = {'Cognitive Measures'};
G_label = {'Hallucinatory Experiences'};
Sub_label = {'Healthy Subjects'};

rpt_name = 'Multisite_dataset_CPCA_analysis_p_0.005_temp'; % Enter the filename of the report you would like to generate.
rpt_title = 'Multisite dataset: 3 stage CPCA Analysis v3.1'; % Enter the title of the report you would like to generate.

%% Set missing values to NaN

data_num_only = data(2:end,2:end);
data_num_only(cellfun(@(s)ischar(s),data_num_only)) = {NaN};
data_num_only(cellfun(@(s)ismissing(s),data_num_only)) = {NaN};
data_num_only(cellfun(@(s)isequal(s,-99),data_num_only)) = {NaN};
data_num_only(cellfun(@(s)isequal(s,-88),data_num_only)) = {NaN};

data(2:end,2:end) = data_num_only;

Glabel_temp = get_CAPS_LSHS_item_names();
data(1,87:118) = Glabel_temp{1};
data(1,215:230) = Glabel_temp{2};
data(1,1:end) = strrep(data(1,1:end),'_',' ');


%% Setting certain parameters for CPCA

params.comp_rel_prop_cutoff = 0.5;
params.nIter = 1000;
params.varimaxFlag = true;
params.p_val = 0.005;
params.pred_rel_p_val = 0.005;
params.n_bootstrap = 10000;

%% Report settings

import mlreportgen.dom.*;
import mlreportgen.report.*;

rpt = Report(rpt_name,'pdf'); 

tp = TitlePage; 
tp.Title = rpt_title; 
tp.Subtitle = sprintf(''); 
add(rpt,tp);

add(rpt,TableOfContents); 

%% Main interation 

for zi = 1:length(Z_index_cell)
    rpt_ch(zi) = Chapter; 
    rpt_ch(zi).Title = sprintf('Z matrix: %s',Z_label{zi}); 

    for gi = 1:length(G_index_cell)
        rpt_sec(zi,gi) = Section;
        rpt_sec(zi,gi).Title = sprintf('G matrix: %s',G_label{gi});

        for si = 1:length(Sub_index_cell)
            rpt_sub_sec(zi,gi,si) = Section;
            
            params.page_header = sprintf('Z matrix: %s, G matrix: %s, Subject set: %s',Z_label{zi},G_label{gi},Sub_label{si});
            G_CPCA(zi,gi,si) = BehCPCAgenerateReport(rpt_sub_sec(zi,gi,si), data, Z_index_cell{zi}, G_index_cell{gi}, ...
                Sub_index_cell{si}, params);
            
            rpt_sub_sec(zi,gi,si).Title = sprintf('Subject set: %s (n = %d, p = %0.2f)',...
                Sub_label{si},size(G_CPCA(zi,gi,si).Z,1),G_CPCA(zi,gi,si).p);
            
            add(rpt_sec(zi,gi), rpt_sub_sec(zi,gi,si));
            br = PageBreak();
            add(rpt_sec(zi,gi), br);
            
            pause(15);
            
        end
        add(rpt_ch(zi), rpt_sec(zi,gi));
    end
    add(rpt, rpt_ch(zi));
    close all;
end

save(sprintf('%s.mat',rpt_name),'G_CPCA');

%% add p-value table

for zi = 1:length(Z_index_cell)
    for gi = 1:length(G_index_cell)
        for si = 1:length(Sub_index_cell)
            table_mat{(zi-1)*length(Sub_index_cell)*length(G_index_cell)+(gi-1)*length(Sub_index_cell)+si,1} = categorical(string(Z_label{zi}));
            VarNames{1} = 'Z';
            table_mat{(zi-1)*length(Sub_index_cell)*length(G_index_cell)+(gi-1)*length(Sub_index_cell)+si,2} = categorical(string(G_label{gi}));
            VarNames{2} = 'G';
            table_mat{(zi-1)*length(Sub_index_cell)*length(G_index_cell)+(gi-1)*length(Sub_index_cell)+si,3} = categorical(string(Sub_label{si}(1)));
            VarNames{3} = 'Sub';
        end
    end
end

p = cat(1,G_CPCA(:).p);
temp = categorical(string(num2str(p(:),'%.2f')));
for jj = 1:length(temp)
    table_mat{jj,4} = temp(jj);
end
VarNames{4} = 'p-value';

temp = cat(1,G_CPCA.predictStruct);
comp_rel_prop_cell = {temp.CompReliabilityProportion};
max_comps = max(cellfun('length',comp_rel_prop_cell));

VarNames(5:4+max_comps) = sprintfc('C%d',1:max_comps);


for ii = 1:length(temp)
    for jj = 1:max_comps
        if jj <= length(comp_rel_prop_cell{ii})
            table_mat{ii,4+jj} = categorical(string(num2str(comp_rel_prop_cell{ii}(jj),'%.2f')));
        else
            table_mat{ii,4+jj} = categorical(string(num2str(NaN,'%.2f')));
        end
    end
end

temp = sprintfc('table_mat{%d},',1:length(table_mat));
table_str = cat(2,temp{:});

T = cell2table(table_mat);
T.Properties.VariableNames = VarNames;

H = Heading4(Text('Table with all significance values')); add(rpt, H);

Tbl = MATLABTable(T);
Tbl.TableEntriesStyle = {Height('20pt')};
add(rpt, Tbl);

%% Summary page

H = Heading1(Text('Summary: ')); add(rpt, H);

for zi = 1:length(Z_index_cell)
    for gi = 1:length(G_index_cell)
        for si = 1:length(Sub_index_cell)

            H = Heading1(Text(sprintf('Z: %s, G: %s, Sub: %s, p = %.3f, predictor loading reliability proportion threshold = %.4f',...
                Z_label{zi},G_label{gi},Sub_label{si},G_CPCA(zi,gi,si).p,G_CPCA(zi,gi,si).predrel_cutoff)));
            add(rpt, H);
            
            if isfield(G_CPCA(zi,gi,si), 'Summary')
                for summ_no = 1:size(G_CPCA(zi,gi,si).Summary,1)
                    
                    H = Heading2(Text(sprintf('Summary Number : %d',summ_no)));
                    add(rpt, H);
                    
                    for ii = 1:size(G_CPCA(zi,gi,si).Summary,2)
                        H = Heading3(Text(sprintf('Comp %d (Summary %d), Comp. Rel. Prop = %.3f',...
                            G_CPCA(zi,gi,si).Summary(summ_no,ii).CompNo, summ_no, G_CPCA(zi,gi,si).Summary(summ_no,ii).CompRelProp)));
                        add(rpt, H);
                        
                        CompLoadings = G_CPCA(zi,gi,si).Summary(summ_no,ii).CompLoadings;
                        H = Heading4(Text(sprintf('Dominant Component Loadings (r), %d selected', length(CompLoadings)))); add(rpt, H);
                        for jj = 1:length(CompLoadings)
                            T = Text(sprintf('%s (%.2f) \n', CompLoadings(jj).Name, CompLoadings(jj).Value)); add(rpt, T);
                        end
                        
                        H = Heading4(Text(sprintf('These %d component loadings above, predict %d out of the %d reliable predictor loadings',...
                            length(CompLoadings),G_CPCA(zi,gi,si).Summary(summ_no,ii).comp_load_rel_item, G_CPCA(zi,gi,si).Summary(summ_no,ii).tot_rel_item))); 
                        add(rpt, H);
                        
                        H = Heading4(Text(sprintf('Reliable positive predictor loadings (r, reliability proportion)'))); add(rpt, H);
                        if isfield(G_CPCA(zi,gi,si).Summary(summ_no,ii), 'PredLoadingsP')
                            PredLoadingsP = G_CPCA(zi,gi,si).Summary(summ_no,ii).PredLoadingsP;
                            for jj = 1:length(PredLoadingsP)
                                T = Text(sprintf('%s (%.2f, %.2f) \n', PredLoadingsP(jj).Name, ...
                                    PredLoadingsP(jj).Value, PredLoadingsP(jj).Rel));
                                add(rpt, T);
                            end
                        end
                        
                        H = Heading4(Text(sprintf('Reliable negative predictor loadings (r, reliability proportion)'))); add(rpt, H);
                        if isfield(G_CPCA(zi,gi,si).Summary(summ_no,ii), 'PredLoadingsN')
                            PredLoadingsN = G_CPCA(zi,gi,si).Summary(summ_no,ii).PredLoadingsN;
                            for jj = 1:length(PredLoadingsN)
                                T = Text(sprintf('%s (%.2f, %.2f) \n', PredLoadingsN(jj).Name, ...
                                    PredLoadingsN(jj).Value, PredLoadingsN(jj).Rel));
                                add(rpt, T);
                            end
                        end
                        br = PageBreak();
                        add(rpt, br);
                    end
                end
            else
                T = Text(sprintf('No reliable components!!')); add(rpt, T);
            end
            
        end
    end
end

close(rpt)
rptview(rpt)