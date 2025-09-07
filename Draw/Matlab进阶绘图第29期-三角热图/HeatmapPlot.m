% Matlab进阶绘图-三角热图
% 公众号：阿昆的科研日常

%% 数据准备
% 读取数据
load data.mat
% 数据矩阵
Z = data;
% 标签
xlb = {'Carb','Wt','Hp','Cyl','Disp','Qsec','Vs','Mpg','Drat','Gear'};
ylb = {'Carb','Wt','Hp','Cyl','Disp','Qsec','Vs','Mpg','Drat','Gear'};

%% 颜色定义
% TheColor函数获取方式：
% 公众号后台回复：TC
map = TheColor('sci',2061);
% map = flipud(map);

%% 图片尺寸设置（单位：厘米）
figureUnits = 'centimeters';
figureWidth = 15;
figureHeight = 15;

%% 窗口设置
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);

%% 三角热图绘制
% triheatmap(Z,map,xlb,ylb,'trid') % 下三角
triheatmap(Z,map,xlb,ylb,'triu') % 上三角

%% 细节优化
% 赋色
colormap(map)
% 字体和字号
set(gca, 'FontName', 'Arial', 'FontSize', 10)
% 背景颜色
set(gcf,'Color',[1 1 1])

%% 图片输出
exportgraphics(gca,'test.png','Resolution',300)