clear
clc
close all
warning off;
LowRankRate = 0.05;
res = zeros(3 , 6);
path = './';
addpath(genpath(path));
addpath('Multi-Omic Data');
%Import the selected kernels
load GBM.mat
label = table2array(Response(:,2:end));
%Initial cluster centers
%GBM z=[19;50;79;167;179]
%BIC z=[34;42;80;88;94;96] 
%COAD z=[20;43;44;45;50]
%KRCCC z=[10;41;80;86]
%LSCC z=[37;39;80;88;90]
z=[19;50;79;167;179] ;
% 
%% initialization
%GBM最佳分类为5；BIC最佳分类为6；COAD最佳分类为6；KRCCC最佳分类为4；LSCC最佳分类为5
    CluNum =5;
    KerNum = size(KH,3);
    SampleNum = size(KH,1); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    KH = kcenter(KH);
    KH = knorm(KH);
    M = calculateM(KH);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qnorm = 2;
    %%%%%%%%---Average---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma0 = ones(KerNum,1)/KerNum;
    avgKer  = mycombFun(KH,gamma0);
    [KHL] = LocalKernelCalculation(KH  , 0.05);
    LargestIteration = 15;
    KernelWeight = ones(1 , KerNum) / KerNum;
    avgKer=sumKbeta(KHL,KernelWeight);
    BaseValue = norm(avgKer , 'fro')^2;
    RegularizationValue = 10^(-4) * BaseValue;
    Alpha_T = 0 * BaseValue;%默认为0
    Beta_T = 2.^[-8 , -2 , 2 , 6, 10];
    NNRate = [0.01,0.03,0.06,0.09 ,0.11,0.13,0.15,0.18,0.21];
    %   Beta_T = 2.^[-8:1:-2];
    %   NNRate = [0.01:0.01:0.1];
    j=1;
    %% Algorithm
    for IRate = 1 : length(NNRate)
    
        CRate = NNRate(IRate);
        [KHL , M1] = PreprocessFinal(KH , CRate);
        
        for IBeta = 1 : length(Beta_T)
            Beta_C = Beta_T(IBeta);
            InputKnum = max(2* CluNum , round(LowRankRate * SampleNum) );
            [Mu , Z , flag , TotalObj,Iter] = MKSSCERC(KHL, M1, RegularizationValue , Alpha_T , Beta_C , zeros(SampleNum) , LargestIteration , InputKnum);
            PON = TotalObj(1:end-1) - TotalObj(2 : end);
            DON = sum(PON < -10^(-10)) == 0;
            PI = Z > 0;
            Z = Z.*PI;
            %get clustering results
            [U] = baseline_spectral_onkernel( abs( (Z + Z') / 2) , CluNum);
            U_normalized=U./ repmat(sqrt(sum(U.^2, 2)), 1,CluNum);
            indx = litekmeans(U_normalized,CluNum, 'MaxIter',100, 'Start',z,'Replicates',30);
            group = num2str(indx);
            group1 = num2cell(group);
            %Survival analysis and predicts subtypes
            [p] = MatSurv(label(:,1),label(:,2),group1,'CensorLineLength',0,'LineWidth',1.2,'CensorLineWidth',1);
            P(j)=p;
            j=j+1;
            [minP,order]=min(P);
        end
    end
   
