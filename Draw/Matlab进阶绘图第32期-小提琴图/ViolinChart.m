% Matlab进阶绘图-小提琴图
% 公众号：阿昆的科研日常

%% 数据准备
% 读取数据
load data.mat
% 初始化绘图参数
data = Y;

%% 颜色定义
% TheColor函数获取方式：
% 公众号后台回复：TC
map = TheColor('sci',1006);
C = map(1:5,1:3);

%% 图片尺寸设置（单位：厘米）
figureUnits = 'centimeters';
figureWidth = 15;
figureHeight = 13;

%% 窗口设置
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);

%% 小提琴图绘制
violin(data,'x',1:size(data,2),...      % x坐标刻度
            'facecolor',C, ...          % 面颜色
            'edgecolor','k',...         % 轮廓颜色
            'facealpha',0.8,...         % 透明度
            'bw',0.2,...                % 核函数带宽Kernel bandwidth
            'mc',[],...                 % 平均值线颜色
            'medc','k');                % 中位数线颜色
hTitle = title('Violin plot');
hXLabel = xlabel('Samples');
hYLabel = ylabel('\Delta [yesno^{-2}]');

%% 细节优化
% 坐标轴美化
set(gca, 'Box', 'off', ...                                        % 边框
         'LineWidth',1,...                                        % 线宽
         'XGrid', 'off', 'YGrid', 'off', ...                      % 网格
         'TickDir', 'out', 'TickLength', [.005 .005], ...         % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...            % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...          % 坐标轴颜色
         'XTick', 1:5,...                                         % 坐标区刻度、范围
         'Xticklabels',{'Sample1','Sample2','Sample3','Sample4','Sample5'},...
         'XLim', [0 6],...
         'YLim',[-3.5 6])
set(gca, 'FontName', 'Arial', 'FontSize', 10)
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
print(figureHandle,[fileout,'.png'],'-r600','-dpng');