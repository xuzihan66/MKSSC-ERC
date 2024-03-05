clear
clc
close all
warning off;


LowRankRate = 0.05;

res = zeros(3 , 6);
    
path = './';
addpath(genpath(path));
addpath('D:\ѧϰ\��˾���\SimpleMKKMcodes-master\Multi-Omic Data');

%dataName = 'flower17';%'flower17','flower102','proteinFold', 'CCV','UCI_DIGIT','caltech101_nTrain30_48'
%load(['flower17',dataName,'_Kmatrix'],'KH','Y');
load GBM.mat
Data{1,1} = ProgressData(Gene);
Data{1,2} = ProgressData(Methy);
Data{1,3} = ProgressData(Mirna);
label = table2array(Response(:,2:end));
KH(:,:,1)=kernel_matrix(Data{1,1},'RBF_kernel',2^(10));
KH(:,:,2)=kernel_matrix(Data{1,2},'RBF_kernel',2^(12));
KH(:,:,3)=kernel_matrix(Data{1,3},'RBF_kernel',2^(20));
z=[19;50;79;167;179];
% 
%% initialization
%GBM��ѷ���Ϊ3��BIC��ѷ���Ϊ5��Breast��ѷ���Ϊ3��COLON��ѷ���Ϊ3��Kidney��ѷ���Ϊ2��Lung��ѷ���Ϊ4;
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
    [H_normalized1] = mykernelkmeans(avgKer,CluNum);
    % %%%%%%%%%%---Single Best%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for p =1:KerNum
        [H_normalized2] = mykernelkmeans(KH(:,:,p),CluNum);
    end


    [KHL] = LocalKernelCalculation(KH  , 0.05);
    LargestIteration = 15;
    KernelWeight = ones(1 , KerNum) / KerNum;
    avgKer=sumKbeta(KHL,KernelWeight);
    BaseValue = norm(avgKer , 'fro')^2;
    RegularizationValue = 10^(-4) * BaseValue;
    Alpha_T = 0 * BaseValue;%Ĭ��Ϊ0
      Beta_T = 2.^[-8 , -2 , 2 , 6, 10];
      NNRate = [0.01,0.03,0.06,0.09 ,0.11,0.13,0.15,0.18,0.21];

 
%KRCCC������ʱ�Ĳ���
%      Beta_T = 2.^[-8:1:-2];
%      NNRate = [0.01:0.01:0.1];
       j=1;
%    Beta_T =2.^2;
%    NNRate =0.15;
%  j=1;
    %     LowRankRate = RankEstimation(avgKer , RegularizationValue , LowRankRate);
    
    GLMKC_ResZ = zeros(length(NNRate) , length(Beta_T) , 3);
    GLMKC_ResK = zeros(length(NNRate) , length(Beta_T) , 3);
    
    LMKC_ResZ = zeros(length(NNRate) , length(Beta_T) , 3);
    LMKC_ResK = zeros(length(NNRate) , length(Beta_T) , 3);
    
    %% Algorithm
    for IRate = 1 : length(NNRate)
        
        % Preprocessing (Calculating KHL, L , M1)
        CRate = NNRate(IRate);
        [KHL , M1] = LGMKCPreprocessFinal(KH , CRate);
        
        for IBeta = 1 : length(Beta_T)
            Beta_C = Beta_T(IBeta);
            
            InputKnum = max(2* CluNum , round(LowRankRate * SampleNum) );
            
            % Kernel weight and Z calculation
            [Mu , Z , flag , TotalObj,Iter] = LGMKCNew(KHL, M1, RegularizationValue , Alpha_T , Beta_C , zeros(SampleNum) , LargestIteration , InputKnum);
            
            PON = TotalObj(1:end-1) - TotalObj(2 : end);
            DON = sum(PON < -10^(-10)) == 0;
    
            PI = Z > 0;
            Z = Z.*PI;
            figure; % ����һ���µ�ͼ�δ���
            imagesc(Z); % ���ƾ��� A ����ͼ
             colorbar; % ��ʾ��ɫ��
              title('�����Ծ�����ͼ'); % ��ӱ���
             xlabel('����'); % X���ǩ
            ylabel('����'); % Y���ǩ

% ������ɫ�������ɫӳ��
               colormap('jet'); % ʹ����ɫ�������ɫӳ��

% �Ľ���ʾ
              axis square; % ʹX���Y�������ͬ�Ŀ̶ȳ���

% ���������������ÿ����Ԫ��
             grid on;
             ax = gca; % ��ȡ��ǰ������
             ax.GridColor = 'w'; % ������������ɫΪ��ɫ
             ax.GridLineStyle = '--'; % ������������ʽ
             ax.Layer = 'top'; % ȷ��������λ�ڶ���
            
%             [U] = baseline_spectral_onkernel( abs( (Z + Z') / 2) , CluNum);
%  
%             U_normalized=U./ repmat(sqrt(sum(U.^2, 2)), 1,CluNum);
%             indx = litekmeans(U_normalized,CluNum, 'MaxIter',100, 'Start',z,'Replicates',30);
%             group = num2str(indx);
%             group1 = num2cell(group);
%             [p] = MatSurv(label(:,1),label(:,2),group1,'CensorLineLength',0,'LineWidth',1.2,'CensorLineWidth',1);
%             P(j)=p;
%             j=j+1;
%             [minP,order]=min(P);
        end
    end
   