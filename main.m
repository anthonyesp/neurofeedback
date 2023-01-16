% analyzing bipolar channels FC3-CP3, FCz-CPz, FC4-CP4 data from FlexEEG
% with cross-validation

close all;
clear;
clc

% addpath of FBCSP algorithm

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
trials = 1:30;

%% CHOOSE DATA TO ANALYZE IN CROSS-VALIDATION
[fileT,pathT] = uigetfile('*.mat', 'Select one or more files to analyze in CV', 'MultiSelect', 'on');

% load data
if (ischar(fileT))
    nfil = 1;
else
    nfil = length(fileT);
end

% Time windows of the signal
tmin = 0.0;
tmax = 6.0;

tshift = 0.25;
twin = 2.00;

nwin = floor((tmax-tmin-twin)/tshift)+1;

accuracy_time = zeros(nwin,1);
std_time = zeros(nwin,1);

confusion_time = zeros(nwin,2,2);
conf_STD_time = zeros(nwin,2,2);
for w = 1:nwin
    imageryT = [];
    classT = [];
    for fls = 1:nfil
        try
            pathDATA = strcat(pathT,fileT{1,fls});
        catch
            pathDATA = strcat(pathT,fileT);      
        end
        
        [imagery_temp, class_temp, fs] = extraction(pathDATA, [], ch, [], tmin+tshift*(w-1), tmin+twin+tshift*(w-1));

        imageryT = [imageryT; imagery_temp];        %#ok
        classT = [classT; class_temp];              %#ok
    end

    % FILTER BANK
    if (exist('hd','var'))
        EEG = filterBank(imageryT,fs,hd);
    else
        [EEG, hd] = filterBank(imageryT,fs);
    end
    CLASS = classT;

    % reshape to have channels as first dimension
    nch = length(ch);
    EEG = reshape(EEG, [nch, size(EEG,1)/nch, size(EEG,2), size(EEG,3)]);
    CLASS = reshape(CLASS, [nch, size(CLASS,1)/nch]);

    % CROSS-VALIDATION
    % stratified k-fold partition with repeations
    cst = cvpartition(CLASS(1,:),'KFold', k_folds);
    cvPart.numTestSet = k_folds*n_repts;

    c = cst;
    for rep = 1:n_repts
        for ind = 1:k_folds
            cvPart.testInd{(rep-1)*k_folds+ind} = test(c,ind);
            cvPart.trainInd{(rep-1)*k_folds+ind} = training(c,ind);
        end
        c = repartition(c);
    end

    % cross-validation for each lambda
    accuracy_cross = zeros(cvPart.numTestSet,1);
    confusion_cross = zeros(2,2,cvPart.numTestSet);
    for k = 1:cvPart.numTestSet

        % Split dataset for train and test
        classEv = CLASS(:,cvPart.testInd{k});
        classTr = CLASS(:,cvPart.trainInd{k});

        eegEv = EEG(:,cvPart.testInd{k},:,:);
        eegTr = EEG(:,cvPart.trainInd{k},:,:);

        % Reshape
        classTr = reshape(classTr,[nch*size(classTr,2), 1]);
        eegTr = reshape(eegTr,[nch*size(eegTr,2), size(eegTr,3), size(eegTr,4)]);

        classEv = reshape(classEv,[nch*size(classEv,2), 1]);
        eegEv = reshape(eegEv,[nch*size(eegEv,2), size(eegEv,3), size(eegEv,4)]);

        % CSPcomposite + NBPW ALGORITHM
        % csp 2 tasks training
        [W1, W2] = CSPtrain(eegTr, classTr, ch, mCSP);
        [V, ~, ~, ~, Y] = CSPapply(eegTr, classTr, W1, W2,[],[]);

        % features selection
        mcsp = size(W1,2)/2;
        I = MIBIF(V,Y,1,mcsp,kMIBIF);

        % NBPW training
        f = V(:,I);
        cl = Y;

        % evaluation
        [Vev, ~, ~, ~, Y] = CSPapply(eegEv, classEv, W1, W2,[],[]);
        feval = Vev(:,I);
        nt = size(feval,1);
        label_eval = zeros(nt,1);
        for i = 1:nt
            pwx = NBPW(f, cl, feval(i,:), class1);
            if(pwx > 0.5)
                label_eval(i) = class1;
            else
                label_eval(i) = class2;
            end
        end
        accuracy_cross(k) = mean(label_eval == Y);
        confusion_cross(:,:,k) = confusionmat(Y,label_eval);
    end

    % save accuracy and standard deviation for each lambda
    accuracyCV = mean(accuracy_cross);
    stdCV = std(accuracy_cross);
    
    disp(['tmin = ',num2str(tmin+tshift*(w-1)),' tmax = ',num2str(tmin+twin+tshift*(w-1)),', CV result: ', num2str(accuracyCV), ' +/-' num2str(round(stdCV,2))])
    
    % save accuracy data
    accuracy_time(w) = accuracyCV;
    std_time(w) = stdCV;
    
    confusion_time(w,:,:) = mean(confusion_cross,3)/(length(Y)/2);
    conf_STD_time(w,:,:) = std(confusion_cross,[],3)/(length(Y)/2);
end

%% PLOT RESULTS
figure
subplot(1,2,1)
t = twin:tshift:tmax;
boundedline(t,100*accuracy_time,100*std_time,'*--')
xlabel('Time [s]')
ylabel('Accuracy [%]')
ylim([0 100])
grid

subplot(1,2,2)
plot(t,100*confusion_time(:,1,1),'*--')
hold on
plot(t,100*confusion_time(:,2,2),'r*--')
legend('lefh hand', 'right hand','Location','NorthWest')
xlabel('Time [s]')
ylabel('Accuracy [%]')
ylim([0 100])
grid