% Matlab进阶绘图-带伪彩图的特征渲染三维散点图
% 公众号：阿昆的科研日常

clear

%% 数据准备
% 读取数据
load data.mat
% 初始化绘图参数
sz = 10; % 散点尺寸
cellsize = 10; % 伪彩图格网尺寸
zpos = -10; % % 伪彩图z位置

%% 颜色定义
% TheColor函数获取方式：
% 公众号后台回复：TC
map = TheColor('sci',2073);
% map = flipud(map);

%% 图片尺寸设置（单位：厘米）
figureUnits = 'centimeters';
figureWidth = 15;
figureHeight = 13;

%% 窗口设置
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);

%% 带伪彩图的特征渲染三维散点图绘制
Scatter3withPcolor(data,10,cellsize,zpos)
% 标题、标签、视角
hTitle = title('Scatter3Feature with Pcolor');
hXLabel = xlabel('XAxis');
hYLabel = ylabel('YAxis');
hZLabel = zlabel('ZAxis');
view(-32,25)

%% 细节优化
% 赋色
colormap(map)
colorbar
% 坐标区调整
axis tight
set(gca, 'Box', 'on', ...                                                          % 边框
         'LineWidth', 1.5, 'GridLineStyle', '--',...                                   % 坐标轴线宽
         'layer','top',...
         'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on',...                            % 网格
         'TickDir', 'out', 'TickLength', [.01 .01], ...                             % 刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],'ZColor', [.1 .1 .1])          % 坐标轴颜色
% 字体和字号
set(gca, 'FontName', 'Arial', 'FontSize', 11)
set([hXLabel,hYLabel,hZLabel], 'FontName',  'Arial', 'FontSize', 11)
set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
% 背景颜色
set(gcf,'Color',[1 1 1])
set(gca,'Projection','Perspective');

%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = 'test';
print(figureHandle,[fileout,'.png'],'-r300','-dpng');