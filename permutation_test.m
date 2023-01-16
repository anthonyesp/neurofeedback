%%
% analyzing bipolar channels FC3-CP3, FCz-CPz, FC4-CP4 data from FlexEEG
% with cross-validation to evaluate statistically significant differences
% between real and randomized labels using a permutation test

close all;
clear;
clc
%%
addpath('C:\Users\HP\Desktop\Ricerca prof. Arpaia\FBCSP')

% ALGORITHM PARAMETERS
% Repeated Stratified K-Fold CV
k_folds = 5;                % number of folds
n_repts = 10;                % number of repetitions

% Initialize classes pairs
class1 = 1;                 % 1: left imagery
class2 = 2;                 % 2: right imagery

% Algorithm hyperparameters
mCSP = 2;                   % CSP components
kMIBIF = 5;                 % MIBIF components

% Training data selection
ch = 1:3;

% Time windows of the signal
tmin_extr = -1.00;
tmax_extr = 8.00;
tmin = 0.00;
tmax = 7.00;
tshift = 0.25;
twin = 2.00;
nwin = floor((tmax_extr-tmin_extr-twin)/tshift)+1;
norm_baseline = 1;

%% CHOOSE DATA TO ANALYZE IN CROSS-VALIDATION
[fileT,pathT] = uigetfile('*.mat', 'Select one or more files to analyze in CV', 'MultiSelect', 'on');
data_permutation = [];
accuracy = zeros(nwin,k_folds*n_repts);
info_data = split(pathT,'\');
subject = info_data(end-1);
block = info_data(end-2);
session = info_data(end-3);

% load data
if (ischar(fileT))
    nfil = 1;
else
    nfil = length(fileT);
end

imageryTot = [];
classTot = [];
for fls = 1:nfil
    try
        pathDATA = strcat(pathT,fileT{1,fls});
    catch
        pathDATA = strcat(pathT,fileT);
    end

    % Baseline extraction and normalization
    if norm_baseline == 1
        % Signal extraction
        [imagery, class_temp, fs, runs, ~ , trials] = extraction(pathDATA, [], ch, [], tmin_extr, tmax_extr);
        l_ch = length(ch);
        l_tr = length(trials);
        l_r = length(runs);
        imagery_baseline = zeros(l_ch,fs*(tmax_extr-tmin_extr),l_tr*l_r);
        for c_b = 1:l_tr*l_r
            imagery_baseline (:,:,c_b) = imagery(1+l_ch*(c_b-1):l_ch+l_ch*(c_b-1),:);
        end
        % Removing 100 ms before the cue
        EEG_LB_mean_rep=EEGbaseline(imagery_baseline,fs*(1.9+abs(tmin_extr)),fs*(2.0+abs(tmin_extr)));
        imagery_baseline = permute(EEG_LB_mean_rep,[1,3,2]);
        bas_norm = reshape(imagery_baseline,[size(imagery_baseline,1)*size(imagery_baseline,2),size(imagery_baseline,3)]);
        % Baseline removing
        imagery_temp = imagery-bas_norm;

        imageryTot = [imageryTot; imagery_temp];        %#ok
        classTot = [classTot; class_temp];              %#ok
    else
        [imagery_temp, class_temp, fs] = extraction(pathDATA, [], ch, [], tmin_extr, tmax_extr);

        imageryTot = [imageryTot; imagery_temp];        %#ok
        classTot = [classTot; class_temp];              %#ok
    end
end

accuracy_time_real = zeros(nwin,1);
std_time_real = zeros(nwin,1);
accuracy_time_perm = zeros(nwin,1);
std_time_perm = zeros(nwin,1);
for w = 1:nwin
    % Extraction
    wind_min = round(1+(tmin+tshift*(w-1))*fs);
    wind_max = round((tmin+twin+tshift*(w-1))*fs);
    imageryT = imageryTot(:,wind_min:wind_max);
    classT = classTot;

    % FILTER BANK
    if (exist('hd','var'))
        EEG = filterBank(imageryT,fs,hd);
    else
        [EEG, hd] = filterBank(imageryT,fs);
    end
    CLASS = classT;

    % Real labels
    [accuracy_time_real(w),std_time_real(w),accuracy(w,:)] = FBCSP_repeatCV_NBPW(EEG,CLASS,ch,mCSP,kMIBIF,class1,class2,k_folds,n_repts);
    % Permutated labels
    CLASS_perm = randsample(2,length(CLASS),'true');
    [accuracy_time_perm(w),std_time_perm(w)] = FBCSP_repeatCV_NBPW(EEG,CLASS_perm,ch,mCSP,kMIBIF,class1,class2,k_folds,n_repts);
end

figure
t_tot = 0:tshift:tmax;
boundedline(t_tot,accuracy_time_real,std_time_real,'b','transparency', 0.1,'alpha')
hold on
boundedline(t_tot,accuracy_time_perm,std_time_perm,'r','transparency', 0.1,'alpha')
xline(2,'-','Color','y');
xline(3,'-','Color','g');
xline(6,'-','Color','k');
axis([0 7 20 100])
xticks(0:1:7)
legend('real std','real meam','perm std','perm std','cue','MI','relax')
xlabel('time (s)')
ylabel('accuracy (%)')
title(strcat(subject,' block ',block,' session ',session))

%% WILCOXON SIGNED RANK TEST
% The Wilcoxon signed-rank test applies to the case of symmetric continuous distributions.
% Under these assumptions, the mean equals the median, and we can use this procedure
% to test the null hypothesis μ = μ0.
disp('WILCOXON SIGNED RANK TEST')
[p,h] = signrank(accuracy_time_real,accuracy_time_perm,'tail','right','method','approximate');
disp(['Mean accuracy real labels: ',num2str(mean(accuracy_time_real))])
disp(['Mean accuracy permutated labels: ',num2str(mean(accuracy_time_perm))])
disp(['p-value: ',num2str(p)])
disp(['Result of the hypothesis test (1= rejection of the null hypothesis): ',num2str(h)])

% Results
data_permutation.accuracy_time_real = accuracy_time_real;
data_permutation.std_time_real = std_time_real; 
data_permutation.accuracy_time_perm = accuracy_time_perm;
data_permutation.std_time_perm = std_time_perm;
data_permutation.p_value = p; 
data_permutation.noReject_0_Reject_1 = h; 