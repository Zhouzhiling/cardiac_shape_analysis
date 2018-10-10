%% 提取特征并输出
% classification 
% with lasso 
clear all;
rootpath = 'C:\Users\37908\Desktop\zzz\';
feature_file = 'all_features_50_50_50_24_Epi.mat';
F = load([rootpath, feature_file]);
feature = F.output;
dimention = length(size(feature));
numSpl = size(feature, dimention);
L = load([rootpath, 'label_3_24.txt']);
gnd = L;
accuracy = 0;

N=ndims(feature)-1;       %Order of the tensor sample
SZ=size(feature);	
numSpl=SZ(dimention);    % 23 train samples
MPCADADim=200;%200 most discriminative MPCA features for LDA
testQ=30;%Keep 97% variation in each mode
maxK=1;%One iteration only
fea3D = feature;
%    [tUs,odrIdx,TXmean,Wgt,LDAU] = MPCALDA(feature,gnd,MPCADADim,testQ,maxK);
%  [tUs, odrIdx, TXmean, Wgt]  = MPCA(TXX,gndTX,testQ,maxK);

[tUs, odrIdx, TXmean, Wgt] = MPCA(fea3D,-1,testQ,maxK);
fea3Dctr=fea3D-repmat(TXmean,[ones(1,N), numSpl]);%Centering
newfea = ttm(tensor(fea3Dctr),tUs,1:N);%MPCA projection
%Vectorization of the tensorial feature
newfeaDim = 1;
for t = 1:dimention-1
    newfeaDim=newfeaDim*size(newfea,t);
end
newfea=reshape(newfea.data,newfeaDim,numSpl)';%Note: Transposed
%P = ceil(length(odrIdx)*0.5);
%selfea=newfea(:,odrIdx(1:P));%Select the first "P" sorted features
selfea = newfea;

save ./classification/feature_24_case_24_feature.mat selfea

filename = 'feature_24_case_1000_feature.mat';
%% svm
clear
filename = './data/Data_5Neigh/feature_SUV.mat';
load(filename);
fea = feature;
caseNum = size(fea,1);
g = load('./label/Survival_Label.txt');
gnd = g';
% svm training and do leave-one-out test
acc_list = zeros(caseNum,1);

score_list = zeros(1, caseNum);
label_list = zeros(1,caseNum);
for i = 1:caseNum
    fea_known = fea; %fea_select;
    test = fea_known(i,:);
    gnd_test = gnd(i);
    
    fea_known(i,:) = [];    
    gnd_known = gnd;
    gnd_known(i) = [];
    svmStruct = fitcsvm(fea_known,gnd_known);
    score_list(i) = predict(svmStruct,test);
end
acc_list = (score_list==gnd);
accuracy = mean(acc_list)
sensitivity = sum(score_list==gnd & gnd==1) / sum(gnd==1)
specificity = sum(score_list==gnd & gnd==0) / sum(gnd==0)
%% lasso
clear
filename = 'feature_24_case_1000_feature.mat';
load(filename);
fea = selfea;
caseNum = size(fea,1);
gnd = load('label_2_24.txt');

%[b,fitinfo] = lasso(fea,gnd,'CV',5);
%lassoPlot(b,fitinfo,'PlotType','CV');

score_list = zeros(1,caseNum);
acc_list = zeros(1,caseNum);
for i = 1:caseNum
    i
    fea_train = fea; %fea_select;
    fea_test = fea_train(i,:);
    gnd_test = gnd(i);
    gnd_another = 1-gnd_test;
    fea_train(i,:) = [];    
    gnd_train = gnd;
    gnd_train(i) = [];
    
    [B, fitinfo] = lasso(fea_train,gnd_train,'Lambda',0.1,'CV',10);
    
    idxLambda1SE = fitinfo.Index1SE;
    coef = B(:,idxLambda1SE);
    coef0 = fitinfo.Intercept(idxLambda1SE);
    score_list(i) = fea_test*coef + coef0;
    if(abs(score_list(i)- gnd_test)< abs(score_list(i)- gnd_another))
        acc_list(i) = 1;
    end
end
lassoPlot(B,fitinfo,'PlotType','CV');
max_accuracy = mean(acc_list)
%% glm
clear
load(filename);
fea = selfea;
caseNum = size(fea,1);
gnd = load('label_2_24.txt');

%[b,fitinfo] = lasso(fea,gnd,'CV',5);
%lassoPlot(b,fitinfo,'PlotType','CV');

score_list = zeros(1,caseNum);
acc_list = zeros(1,caseNum);
for i = 1:caseNum
    i
    fea_train = fea; %fea_select;
    fea_test = fea_train(i,:);
    gnd_test = gnd(i);
    gnd_another = 1-gnd_test;
    fea_train(i,:) = [];    
    gnd_train = gnd;
    gnd_train(i) = [];
    
    b = glmfit(fea_train,gnd_train,'binomial');
    score_list(i) = glmval(b,fea_test,'probit');
    
    if(abs(score_list(i) - gnd_test)< abs(score_list(i) - gnd_another))
        acc_list(i) = 1;
    end
end
mean(acc_list)

%% random
clear
filename = './feature_24_case_200_feature.mat';
load(filename);
fea = selfea;
caseNum = size(fea,1);
g = load('./label_2_24.txt');
gnd = g';
score_list = zeros(1,caseNum);
acc_list = zeros(1,caseNum);
for i = 1:caseNum
    i
    fea_train = fea; %fea_select;
    fea_test = fea_train(i,:);
    gnd_test = gnd(i);
    fea_train(i,:) = [];    
    gnd_train = gnd;
    gnd_train(i) = [];
    
    Mdl = TreeBagger(200, fea_train, gnd_train, 'Method', 'classification', 'OOBPrediction','On');
    Yfit = predict(Mdl, fea_test);
    score_list(i) = str2num(Yfit{1});
    tmp = score_list(i) == gnd_test
    acc_list(i) = tmp;
end
accuracy = mean(acc_list)
sensitivity = sum(score_list==gnd & gnd==1) / sum(gnd==1)
specificity = sum(score_list==gnd & gnd==0) / sum(gnd==0)


    %% Boost
    filename = 'all_features_24_125000_Epi.mat';
    load(filename);
    fea = output(1:16,:);
    caseNum = size(fea,1);
    gnd = load('label_2_16.txt');
    score_list = zeros(1,caseNum);
    acc_list = zeros(2,caseNum);
    T = 50;
    for i = 1:caseNum
        i
        fea_train = fea; %fea_select;
        fea_test = fea_train(i,:);
        gnd_test = gnd(i);
        fea_train(i,:) = [];    
        gnd_train = gnd;
        gnd_train(i) = [];
        X = fea_train;
        Y = gnd_train;
        
        treeStump = templateTree('MaxNumSplits',1);
%        adaStump = fitcensemble(X,Y,'Method','AdaBoostM1','NumLearningCycles',T, ...
%            'Learners',treeStump);

        totalStump = fitcensemble(X,Y,'Method','TotalBoost','NumLearningCycles',T, ...
            'Learners',treeStump);

        lpStump = fitcensemble(X,Y,'Method','LPBoost','NumLearningCycles',T, ...
            'Learners',treeStump);
    
%        acc_list(1,i) = 1-loss(adaStump, fea_test, gnd_test, 'Lossfun', 'classiferror');
        acc_list(2,i) = 1-loss(totalStump, fea_test, gnd_test, 'Lossfun', 'classiferror');
        acc_list(3,i) = 1-loss(lpStump, fea_test, gnd_test, 'Lossfun', 'classiferror');
        
    end
    
    mean(acc_list,2)

    figure(1131);
        plot(resubLoss(adaStump,'Mode','Cumulative'));
        hold on
        plot(resubLoss(totalStump,'Mode','Cumulative'),'r');
        plot(resubLoss(lpStump,'Mode','Cumulative'),'g');
        xlabel('Number of stumps');
        ylabel('Training error');
        title('Training Error');
        legend('AdaBoost','TotalBoost','LPBoost','Location','NE');
        
        % cross validate the nesumbles
        cvlp = crossval(lpStump,'KFold',10);
        cvtotal = crossval(totalStump,'KFold',10);
        cvada = crossval(adaStump,'KFold',10);

        figure(1132);
        plot(kfoldLoss(cvada,'Mode','Cumulative'));
        hold on
        plot(kfoldLoss(cvtotal,'Mode','Cumulative'),'r');
        plot(kfoldLoss(cvlp,'Mode','Cumulative'),'g');
        hold off
        xlabel('Ensemble size');
        ylabel('Cross-validated error');
        legend('AdaBoost','TotalBoost','LPBoost','Location','NE');
        
%% LDA
        load(filename);
        fea = selfea;
        caseNum = size(fea,1);
        gnd = load('label_2_24.txt');
        score_list = zeros(1,caseNum);
        acc_list = zeros(1,caseNum);
        P = zeros(caseNum,2);
 for i = 1:caseNum
        i
        fea_train = fea; %fea_select;
        fea_test = fea_train(i,:);
        gnd_test = gnd(i);
        gnd_another = 1-gnd_test;
        fea_train(i,:) = [];    
        gnd_train = gnd;
        gnd_train(i) = [];
        
        W = LDA(fea_train,gnd_train);
        L = [ones(1,1) fea_test] * W';
        P(i,:) = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
        if(find(P(i,:)==0)-1 == gnd_test)
            acc_list(i) = 1;
        else
            acc_list(i) = 0;
        end
 end
        mean(acc_list)
   