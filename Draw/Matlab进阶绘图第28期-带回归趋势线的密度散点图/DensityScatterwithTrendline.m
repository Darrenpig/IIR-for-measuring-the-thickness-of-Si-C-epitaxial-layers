% Matlab进阶绘图-带回归趋势线的密度散点图
% 公众号：阿昆的科研日常

%% 数据准备
% 读取数据
load data.mat
% 初始化绘图参数
data = [x,y];
% 密度计算
radius = 1.5; % 定义半径
density_2D = density2D_KD(data(:,1:2),radius); % 2D平面密度
% 回归趋势线生成
xq = min(data(:,1)):0.1:max(data(:,1));
p = polyfit(data(:,1),data(:,2),1);
x1 = linspace(min(data(:,1)),max(data(:,1)),100);
y1 = polyval(p,x1);

%% 颜色定义
% TheColor函数获取方式：
% 公众号后台回复：TC
map = TheColor('sci',2064);
map = flipud(map);

%% 图片尺寸设置（单位：厘米）
figureUnits = 'centimeters';
figureWidth = 12;
figureHeight = 9;

%% 窗口设置
figureHandle = figure;
set(gcf, 'Units', figureUnits, 'Position', [0 0 figureWidth figureHeight]);
hold on

%% 带拟合线的密度散点图绘制
scatter(data(:,1), data(:,2), 5, density_2D, 'filled')
p1 = plot(x1,y1,'LineStyle',':','LineWidth',2,'Color','k');
hTitle = title('Satellite-derived bathymetry');
hXLabel = xlabel('ICESat-2 bathymetric points in depth (m)');
hYLabel = ylabel('Estimated depth (m)');

%% 细节优化
% 赋色
colormap(map)
colorbar
% 坐标轴美化
set(gca, 'Box', 'off', ...                                        % 边框
         'LineWidth',1,...                                        % 线宽
         'XGrid', 'on', 'YGrid', 'on', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.005 .005], ...         % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...            % 小刻度
         'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...          % 坐标轴颜色
         'XTick', 0:40:160,...                                    % 坐标区刻度、范围
         'XLim', [0 160],...
         'YTick', 0:40:160,...
         'YLim', [0 160])
hLegend = legend(p1, ...
                 'Regression line',...
                 'Location', 'northwest'); 
set(gca, 'FontName', 'Arial', 'FontSize', 10)
set([hLegend, hXLabel, hYLabel], 'FontSize', 11, 'FontName', 'Arial')
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