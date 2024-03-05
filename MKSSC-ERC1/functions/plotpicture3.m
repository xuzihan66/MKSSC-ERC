% �������ж���ѧ���ݴ���ڲ�ͬ�ľ�����
% �������ݾ�������Ϊ data1, data2, data3, ...
data1=Data{1,1};
data2=Data{1,2};
data3=Data{1,3};
% ������ѧ���ݺϲ���һ����ľ�����
combined_data = [data1, data2, data3];

% �Ժϲ�������ݽ��й�һ�����ɸ��ݾ������ѡ����ʵķ�����
normalized_data = zscore(combined_data); 

% ���㽵ά������ɷ�
[coeff, score, ~, ~, explained] = pca(normalized_data);

% ��ȡÿ��������Ӧ�İ�֢���ͱ�ǩ���������� label �����У�
label = cell2mat(cellfun(@str2double,group ,'UniformOutput',false));
labels = label(:);

% ���ݰ�֢���ͱ�ǩ���������ಢ�ֱ����
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
% % ����ͼ���ͱ���
% %  legend('Subtype 1', 'Subtype 2', 'Subtype 3','Subtype 4','Subtype 5','Location','best');
% %  title('PCA of Multiple Omics Data for Cancer Subtype Prediction');

% ����PCAͼ
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

% �����ɫ����ͼ��
colormap(gca, jet(max(labels))); % ʹ��jet��ɫӳ��
c = colorbar;
c.Label.String = 'Subtypes';









