% Matlab进阶绘图-特征渲染的滑珠气泡图
% By：阿昆的科研日常

clear

%% 数据准备
% 读取数据
load data.mat
% 初始化绘图参数
X = data;
ylabels = p;
f1 = SZ;
f2 = SZ2;

%% 颜色定义
% TheColor函数获取方式：
% 公众号后台回复：TC
C = TheColor('sci',2068);

%% 图片尺寸设置（单位：厘米）
figureUnits = 'centimeters';
figureWidth = 18;
figureHeight = 14;

%% 窗口设置
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);

%% 特征渲染的滑珠气泡图绘制
[s,idx] = BubbleScatter(X,f1,f2);
bubblesize([5 20])
s.MarkerFaceAlpha = 0.8;
yl = ylabels(idx);
hTitle = title('BubbleScatter with feature');
hXLabel = xlabel('MeanDecreaseAccuracy');
hYLabel = ylabel('Product');

%% 细节优化
% 赋色
colormap(C)
% 坐标区调整
set(gca, 'Box', 'off', ...                                   % 边框
         'LineWidth', 1,...                                  % 线宽
         'XGrid', 'off', 'YGrid', 'on', ...                  % 网格
         'TickDir', 'out', 'TickLength', [.005 .005], ...    % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...       % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])        % 坐标轴颜色
% 坐标轴刻度调整
set(gca, 'YTick', 1:16,...
         'Ylim' , [0.5 15.5], ...
         'Xlim' , [-0.05 1.05], ...
         'XTick', 0:0.2:1,...
         'yticklabel',yl)
% 添加图例
blgd = bubblelegend('Style','vertical',...
    'BubbleSizeOrder','descending',...
    'box','off',...
    'Location','northeastoutside',...
    'NumBubbles',3,... ...
    'FontName', 'Helvetica',...
    'FontSize', 9);
colorbar('Position',[0.79,0.11,0.03,0.4])
% 字体和字号
set(gca, 'FontName', 'Arial', 'FontSize', 9)
set([hXLabel, hYLabel], 'FontSize', 11, 'FontName', 'Arial')
set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
% 背景颜色
set(gcf,'Color',[1 1 1])
% 添加上、右框线
xc = get(gca,'XColor');
yc = get(gca,'YColor');
unit = get(gca,'units');
ax = axes( 'Units', unit,...
           'Position',get(gca,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor',xc,...
           'YColor',yc);
set(ax, 'linewidth',1,...
        'XTick', [],...
        'YTick', []);

%% 图片输出
figW = figureWidth;
figH = figureHeight;
set(figureHandle,'PaperUnits',figureUnits);
set(figureHandle,'PaperPosition',[0 0 figW figH]);
fileout = 'test';
print(figureHandle,[fileout,'.png'],'-r300','-dpng');