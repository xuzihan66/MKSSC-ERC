%%This code is used to select base kernels
clear
clc
close all
warning off;
path = './';
addpath(genpath(path));
addpath('Multi-Omic Data');
load LSCC1.mat
Data{1,1} = ProgressData(Gene);
Data{1,2} = ProgressData(Methy);
Data{1,3} = ProgressData(Mirna);
label = table2array(Response(:,2:end));
b=0.5;
z=[20;21;33;37;38] ;
numclass=5;
num=105;
for i=1:60
     KH(:,:,i)=kernel_matrix(Data{1,2},'RBF_kernel', 2^(i));
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

