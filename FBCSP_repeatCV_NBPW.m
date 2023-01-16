function [accuracyCV,stdCV,accuracy] = FBCSP_repeatCV_NBPW(EEG,CLASS,ch,mCSP,kMIBIF,class1,class2,k_folds,n_repts)

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
accuracy = zeros(cvPart.numTestSet,1);
confusion = zeros(2,2,cvPart.numTestSet);
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
    accuracy(k) = mean(label_eval == Y)*100;
    confusion(:,:,k) = confusionmat(Y,label_eval);
end
accuracyCV = mean(accuracy);
stdCV = std(accuracy);
end