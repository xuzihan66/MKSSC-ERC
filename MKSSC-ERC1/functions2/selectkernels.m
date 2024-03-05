%%此代码用来选择核函数

clear
clc
close all
warning off;

path = './';
addpath(genpath(path));
addpath('D:\学习\多核聚类\SimpleMKKMcodes-master\Multi-Omic Data');
%dataName = 'flower17';%'flower17','flower102','proteinFold', 'CCV','UCI_DIGIT','caltech101_nTrain30_48'
%load(['flower17',dataName,'_Kmatrix'],'KH','Y');
load LSCC.mat
Data{1,1} = ProgressData(Gene);
Data{1,2} = ProgressData(Methy);
Data{1,3} = ProgressData(Mirna);
label = table2array(Response(:,2:end));
b=0.5;
z=[20;21;33;37;38] ;
%GBM最佳分类为3；Breast最佳分类为3；COLON最佳分类为3；Kidney最佳分类为2；Lung最佳分类为4；COAD最佳分类为4
numclass=5;
num=105;
for i=1:60
     %KH(:,:,i)=kernel_matrix(Data{1,1},'lin_kernel', 2^i);
     KH(:,:,i)=kernel_matrix(Data{1,2},'RBF_kernel', 2^(i));
  %  for j=1:3
  %      KT(:,:,j)=KH(:,:,i);
  %  end
  %  [H_normalized0,Sigma0,obj0] = mkkmeans_train(KT,numclass);
  %  U_normalized = H_normalized0 ./ repmat(sqrt(sum(H_normalized0.^2, 2)), 1,numclass);
    indx = litekmeans(KH(:,:,i),numclass, 'MaxIter',100,'Start',z, 'Replicates',10);
    group = num2str(indx);
    group = num2cell(group);
    [p] = MatSurv(label(:,1),label(:,2),group,'CensorLineLength',0);
    P(i)=p;
end

for i=1:60
    KH(:,:,i)=kernel_matrix(Data{1,3},'RBF_kernel', 2^(i));
    D=eigs(KH(:,:,i),numclass);
    r=sqrt(D(1))/sum(sqrt(D));
    div=-log(r)/log(numclass);
    acc=(1/P(i)-1/max(P))/(1/min(P)-1/max(P));
    F(i)=(b^2+1)*acc*div/(b^2*acc+div);
 
end

