function [C,cl] = select_Initial_cluster_centers(J)
%C中存放聚类中心的值
%cl存放每一个点的类标
%输入数据，计算距离矩阵
[m, n, p] = size(J);
data = reshape(double(J), m*n, p);
dist = pdist2(data,data);

%计算截断误差dc
percent = 2.0;
N = size(dist,1);
position = round((N*(N-1)/2)*percent/100);
%提取dist的上三角阵，对角线均为0
tri_dist = triu(dist,1);
s_dist = sort(tri_dist(tri_dist~=0));
dc = s_dist(position);
%输出截断距离
fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);

%计算局部密度（利用高斯核）
rho = sum(exp(-(dist./dc).^2),2);
%求所有距离中的最大值
max_dist=max(max(dist));
%将 rho 按降序排列，ordrho存放原序，rho_sorted存放降序后rho
[rho_sorted,ordrho]=sort(rho,'descend');

% 处理 rho 值最大的数据点
delta(ordrho(1))=-1;
%最大的局部均值的点，nneigh为0
nneigh(ordrho(1))=0;
% 生成 delta 和 nneigh 数组
for k=2:N
   delta(ordrho(k))=max_dist;
   for p=1:k-1
     if(dist(ordrho(k),ordrho(p))<delta(ordrho(k)))
        delta(ordrho(k))=dist(ordrho(k),ordrho(p));
        nneigh(ordrho(k))=ordrho(p);    % 记录所有局部均值比ordrho(k)大的点中，距离最近的点的编号  
     end
   end
end
% 生成 rho 值最大数据点的 delta 值
delta(ordrho(1))=max(delta(:));

%画决策图，利用 rho 和 delta 画出一个所谓的“决策图”
disp('Generated file:DECISION GRAPH')
figure
subplot(1,2,1)
plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','r');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')

%计算gamma值
for i=1:N
  gamma(i)=rho(i)*delta(i);
end
[gamma_sorted,ordgamma]=sort(gamma,'descend');
subplot(1,2,2)
plot(gamma_sorted(:),'o','MarkerSize',5,'MarkerFaceColor','b','MarkerEdgeColor','r');
title ('Gamma Graph','FontSize',15.0)

%程序暂停，输入一个用于选取聚类个数的q，继续执行下面的命令
q = input('q=');

% 初始化 cluster 个数
NCLUST=0;
%cl(i)=j 表示第 i 号数据点归属于第 j 个 cluster
for b=1:N
  cl(b)=-1;
end 
% 统计数据点（即聚类中心）的个数
C = zeros(q-1,3);
for n=1:N
  if gamma(n) > gamma_sorted(q)
     NCLUST=NCLUST+1;
     cl(n)=NCLUST; % 第 n 号数据点是第 NCLUST类的聚类中心
     C(NCLUST,:) = data(n,:);
     icl(NCLUST)=n;%逆映射,第 NCLUST 个 cluster 的中心为第 n 号数据点
  end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);

%将其他数据点归类,借助ordrho和nneigh
disp('assignation');
for i=1:N
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end

%区分cluster_core和cluster_halo
for i=1:N
  halo(i)=cl(i);
end 
if (NCLUST>1)
  % 初始化数组 border_rho 为 0,每个 cluster 定义一个 border_rho 值
  for i=1:NCLUST
    border_rho(i)=0;
  end
  % 获取每一个 cluster 中平均密度的一个界 border_rho
  for i=1:N-1
    for j=i+1:N
      % 距离足够小但不属于同一个 cluster 的 i 和 j
      if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
        rho_aver=(rho(i)+rho(j))/2.; %% 取 i,j 两点的平均局部密度
        if (rho_aver>border_rho(cl(i)))
          border_rho(cl(i))=rho_aver;
        end
        if (rho_aver>border_rho(cl(j)))
          border_rho(cl(j))=rho_aver;
        end
      end
    end
  end
  % halo 值为 0 表示为 outlier(离群点)
  for i=1:N
    if (rho(i)<border_rho(cl(i)))
      halo(i)=0;
    end
  end 
end

%逐一处理每个 cluster
for i=1:NCLUST
  nc=0; %% 用于累计当前 cluster 中数据点的个数
  nh=0; %% 用于累计当前 cluster 中核心数据点的个数
  for j=1:N
    if (cl(j)==i)
      nc=nc+1;
    end
    if (halo(j)==i)
      nh=nh+1;       %离群点的halo为0
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh); 
end