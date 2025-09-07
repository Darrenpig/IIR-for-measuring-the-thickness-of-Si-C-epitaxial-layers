function density_2D = density2D_KD(data,radius)
% 功能：利用KD树提取离散点2D密度特征
% 输入：data   - 原始数据(m*3)    
% 输出：planes - 拟合所得平面参数 
M = size(data,1);
density_2D = zeros(M,1);
idx = rangesearch(data(:,1:2),data(:,1:2),radius,'Distance','euclidean','NSMethod','kdtree');
for i = 1:M
    density_2D(i,1) = length(idx{i})/(pi*radius^2);
end
end
% Matlab点云处理与可视化第4期
% 基于KD树的离散点密度特征提取
% 公众号：阿昆的科研日常