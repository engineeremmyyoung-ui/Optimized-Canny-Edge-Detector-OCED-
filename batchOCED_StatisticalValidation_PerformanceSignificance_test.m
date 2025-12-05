%% batchOCED Statistical Performance Validation with Summary Sheet
clear; clc; close all;

%% Define metrics and modalities
metrics = {'Accuracy','MSE','MBE','PSNR','Entropy','Kappa_Index', ...
           'SSIM','Precision','F1_Score','Recall','Execution_Time'};
modalities = {'Xray','MRI','CT','Ultrasound'};
n_metrics = length(metrics);
n_modalities = length(modalities);

%% Load data from Excel interactively
[fileName,pathName] = uigetfile("C:\Users\User\Documents\OCED_Metrics.xlsx");
if isequal(fileName,0)
    error('No file selected. Exiting.');
end
fileName = fullfile(pathName,fileName);

% Determine number of samples per modality
sample_counts = zeros(1,n_modalities);
for mod = 1:n_modalities
    temp = readmatrix(fileName,'Sheet',modalities{mod});
    sample_counts(mod) = size(temp,1);
end
n_samples = max(sample_counts);

% Initialize data array
data = zeros(n_samples,n_metrics,n_modalities);
for mod = 1:n_modalities
    temp = readmatrix(fileName,'Sheet',modalities{mod});
    if size(temp,2) ~= n_metrics
        error('Number of columns in sheet %s does not match metrics list.', modalities{mod});
    end
    data(1:size(temp,1),:,mod) = temp;
end

%% Preallocate results
ANOVA_p = zeros(n_metrics,1);
KW_p = zeros(n_metrics,1);
ANOVA_sig = strings(n_metrics,1);
KW_sig = strings(n_metrics,1);

tTestPairs = strings(n_metrics, n_modalities*(n_modalities-1)/2);
tTestP = zeros(n_metrics, n_modalities*(n_modalities-1)/2);
tTestSig = strings(n_metrics, n_modalities*(n_modalities-1)/2);
CohenD = zeros(n_metrics, n_modalities*(n_modalities-1)/2);

TukeyPairs = strings(n_metrics, n_modalities*(n_modalities-1)/2);
TukeyP = zeros(n_metrics, n_modalities*(n_modalities-1)/2);
TukeySig = strings(n_metrics, n_modalities*(n_modalities-1)/2);

%% Statistical analysis
for m = 1:n_metrics
    metric_data = reshape(data(:,m,:),[],1);
    groups = repmat(1:n_modalities,n_samples,1);
    groups = groups(:);
    
    % ANOVA
    [p,~,stats] = anova1(metric_data, groups,'off');
    ANOVA_p(m) = p;
    ANOVA_sig(m) = sigLabel(p);
    
    % Tukey post hoc
    idxT = 1;
    if p < 0.05
        [c,~,~,gnames] = multcompare(stats,'Display','off','CType','tukey-kramer');
        for k = 1:size(c,1)
            i1 = str2double(gnames{c(k,1)});
            i2 = str2double(gnames{c(k,2)});
            TukeyPairs(m,idxT) = strcat(modalities{i1},' vs ',modalities{i2});
            TukeyP(m,idxT) = c(k,6);
            TukeySig(m,idxT) = sigLabel(c(k,6));
            idxT = idxT + 1;
        end
    end
    
    % Kruskal Wallis
    [p_kw,~,~] = kruskalwallis(metric_data, groups,'off');
    KW_p(m) = p_kw;
    KW_sig(m) = sigLabel(p_kw);
    
    % Pairwise t tests
    idx = 1;
    for i = 1:n_modalities-1
        for j = i+1:n_modalities
            x = data(:,m,i);
            y = data(:,m,j);
            [~,p_t] = ttest2(x,y);
            tTestP(m,idx) = p_t;
            tTestPairs(m,idx) = strcat(modalities{i},' vs ',modalities{j});
            tTestSig(m,idx) = sigLabel(p_t);
            
            pooled_std = sqrt(((numel(x)-1)*var(x) + (numel(y)-1)*var(y))/(numel(x)+numel(y)-2));
            CohenD(m,idx) = (mean(x)-mean(y))/pooled_std;
            idx = idx +1;
        end
    end
end

%% Compute means and standard deviations
Means = zeros(n_metrics,n_modalities);
Stds  = zeros(n_metrics,n_modalities);
for m = 1:n_metrics
    for mod = 1:n_modalities
        Means(m,mod) = mean(data(:,m,mod));
        Stds(m,mod)  = std(data(:,m,mod));
    end
end

%% Convert to tables
Means_Table = array2table(Means,'VariableNames',modalities,'RowNames',metrics);
Stds_Table  = array2table(Stds,'VariableNames',modalities,'RowNames',metrics);

%% Export to Excel
outFile = fullfile(pathName,'OCED_Statistical_Validation_Tukey.xlsx');

ANOVA_Table = table(metrics', ANOVA_p, ANOVA_sig, KW_p, KW_sig, ...
    'VariableNames', {'Metric','ANOVA_p','ANOVA_Significance','KW_p','KW_Significance'});
writetable(ANOVA_Table,outFile,'Sheet','ANOVA_KW');

Tukey_Table = table(TukeyPairs, TukeyP, TukeySig, ...
    'VariableNames', {'Pair','p_value','Significance'});
writetable(Tukey_Table,outFile,'Sheet','Tukey_PostHoc');

writetable(Means_Table,outFile,'Sheet','Means','WriteRowNames',true);
writetable(Stds_Table,outFile,'Sheet','StdDev','WriteRowNames',true);

tTest_Table = table(tTestPairs, tTestP, tTestSig, CohenD, ...
    'VariableNames', {'Pair','p_value','Significance','Cohens_d'});
writetable(tTest_Table,outFile,'Sheet','Pairwise_TTest');

%% Summary Sheet including Mean and StdDev
Summary = table('Size',[n_metrics, 2 + n_modalities*2], ...
                'VariableTypes',[{'string','string'}, repmat({'double'},1,n_modalities*2)], ...
                'VariableNames',[{'Metric','Significance'}, strcat(modalities,'_Mean'), strcat(modalities,'_Std')]);

for m = 1:n_metrics
    Summary.Metric(m) = metrics{m};
    Summary.Significance(m) = ANOVA_sig(m);
    
    for mod = 1:n_modalities
        Summary.(strcat(modalities{mod},'_Mean'))(m) = Means(m,mod);
        Summary.(strcat(modalities{mod},'_Std'))(m)  = Stds(m,mod);
    end
end

writetable(Summary,outFile,'Sheet','Summary','WriteRowNames',false);

disp(['All statistical results and Summary sheet exported to ', outFile]);

%% Helper function
function label = sigLabel(p)
    if p < 0.001
        label = '***';
    elseif p < 0.01
        label = '**';
    elseif p < 0.05
        label = '*';
    else
        label = '';
    end
end
