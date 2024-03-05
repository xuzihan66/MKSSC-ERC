% 假设你有多组学数据存放在不同的矩阵中
% 具体数据矩阵命名为 data1, data2, data3, ...
data1=Data{1,1};
data2=Data{1,2};
data3=Data{1,3};
% 将多组学数据合并到一个大的矩阵中
combined_data = [data1, data2, data3];

% 对合并后的数据进行归一化（可根据具体情况选择合适的方法）
normalized_data = zscore(combined_data); 

% 计算降维后的主成分
[coeff, score, ~, ~, explained] = pca(normalized_data);

% 提取每个样本对应的癌症亚型标签（假设存放在 label 变量中）
label = cell2mat(cellfun(@str2double,group ,'UniformOutput',false));
labels = label(:);

% 根据癌症亚型标签将样本分类并分别绘制
unique_labels = unique(labels);
num_labels = length(unique_labels);
% % 
% %  figure;
% %  hold on;
% %  for i = 1:num_labels
% %      idx = labels == unique_labels(i);
% %      scatter3(score(idx, 1), score(idx, 2), score(idx, 3),'filled');
% % end
% % hold off;
% % 
% % 设置图例和标题
% %  legend('Subtype 1', 'Subtype 2', 'Subtype 3','Subtype 4','Subtype 5','Location','best');
% %  title('PCA of Multiple Omics Data for Cancer Subtype Prediction');

% 绘制PCA图
% figure;
% scatter3(score(:,1),score(:,2),score(:,3),[],labels,'filled');
% 
% 
% title('PCA of Multi-Omics Data for Cancer Subtype Prediction');
% 
% legend('Subtype 1', 'Subtype 2', 'Subtype 3','Subtype 4','Subtype 5','Location','best');
% 
figure;
scatter3(score(:,1), score(:,2), score(:,3), [], labels, 'filled');
title('PCA Scatter Plot of Subtype Prediction Results');

% 添加颜色条和图例
colormap(gca, jet(max(labels))); % 使用jet颜色映射
c = colorbar;
c.Label.String = 'Subtypes';









