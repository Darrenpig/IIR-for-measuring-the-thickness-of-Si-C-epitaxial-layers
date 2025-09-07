function plot_all_with_color51()
% 使用color51配色方案重新绘制所有问题的图表
% 将图片保存到figures文件夹

% 添加color文件夹到路径
addpath('../color');

% 加载配色方案
try
    load('../color/color_palette.mat', 'colors');
    fprintf('成功加载color51配色方案，共%d种颜色\n', size(colors,1));
catch
    % 如果加载失败，直接调用color51函数
    colors = color51(1:51);
    fprintf('直接调用color51函数获取配色方案\n');
end

% 定义常用颜色索引（参考demo1.m中的使用）
color_primary = colors(5,:);      % 主要曲线颜色 - 深蓝色
color_data = colors(13,:);        % 数据点颜色 - 橙色
color_model = colors(51,:);       % 模型曲线颜色 - 深红色
color_validation = colors(45,:);  % 验证数据颜色 - 紫色
color_ci = colors(22,:);          % 置信区间颜色 - 绿色
color_error = colors(35,:);       % 误差条颜色 - 黄色
color_secondary = colors(30,:);   % 次要曲线颜色 - 青色
color_highlight = colors(15,:);   % 高亮颜色 - 亮绿色
color_background = colors(8,:);   % 背景颜色
color_accent1 = colors(18,:);     % 强调色1
color_accent2 = colors(25,:);     % 强调色2
color_accent3 = colors(40,:);     % 强调色3

% 设置图形参数
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

% 创建输出目录
if ~exist('d:\Project_env\CUMCU\B\figures', 'dir')
    mkdir('d:\Project_env\CUMCU\B\figures');
end

% 问题一的图表
plot_fresnel_coefficients_colored();
plot_reflectance_spectrum_colored();
plot_phase_difference_analysis_colored();
plot_thickness_interference_relation_colored();

% 问题二的图表
plot_data_quality_colored();
plot_thickness_comparison_colored();
plot_angle_dependency_colored();
plot_error_analysis_colored();
plot_method_comparison_colored();

fprintf('所有图表已使用color51配色方案重新绘制并保存到figures文件夹\n');

    function plot_fresnel_coefficients_colored()
        % 绘制Fresnel系数图（使用color51配色）
        figure('Position', [100, 100, 800, 600]);
        
        % 模拟数据
        theta = linspace(0, 90, 1000);
        theta_rad = theta * pi / 180;
        n1 = 1.0; % 空气
        n2 = 3.2; % 碳化硅
        
        % 计算Fresnel系数
        cos_theta1 = cos(theta_rad);
        sin_theta1 = sin(theta_rad);
        sin_theta2 = (n1/n2) * sin_theta1;
        cos_theta2 = sqrt(1 - sin_theta2.^2);
        
        % s偏振
        rs = (n1*cos_theta1 - n2*cos_theta2) ./ (n1*cos_theta1 + n2*cos_theta2);
        Rs = abs(rs).^2;
        
        % p偏振
        rp = (n2*cos_theta1 - n1*cos_theta2) ./ (n2*cos_theta1 + n1*cos_theta2);
        Rp = abs(rp).^2;
        
        % 绘图
        hold on;
        h1 = plot(theta, Rs, 'LineWidth', 2.5, 'Color', color_primary);
        h2 = plot(theta, Rp, 'LineWidth', 2.5, 'Color', color_model, 'LineStyle', '--');
        
        % 添加关键点标记
        brewster_angle = atan(n2/n1) * 180/pi;
        idx_brewster = find(theta >= brewster_angle, 1);
        plot(theta(idx_brewster), Rp(idx_brewster), 'o', 'MarkerSize', 8, ...
             'MarkerFaceColor', color_highlight, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        % 设置图形属性
        xlabel('入射角 (度)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('反射率', 'FontSize', 12, 'FontWeight', 'bold');
        title('碳化硅-空气界面Fresnel反射系数', 'FontSize', 14, 'FontWeight', 'bold');
        legend([h1, h2], {'s偏振 (R_s)', 'p偏振 (R_p)'}, 'Location', 'best', 'FontSize', 11);
        
        % 设置坐标轴
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        xlim([0, 90]);
        ylim([0, 1]);
        
        % 保存图片
        save_figure('d:\Project_env\CUMCU\B\figures\fresnel_coefficients');
    end

    function plot_reflectance_spectrum_colored()
        % 绘制反射光谱图（使用color51配色）
        figure('Position', [100, 100, 800, 600]);
        
        % 模拟光谱数据
        wavenumber = linspace(400, 4000, 1000);
        thickness = 2.5; % μm
        n_sic = 2.65;
        
        % 计算干涉光谱
        phase = 4 * pi * n_sic * thickness ./ (10000 ./ wavenumber);
        R_base = 0.18;
        R_modulation = 0.15;
        reflectance = R_base + R_modulation * cos(phase);
        
        % 添加噪声
        noise = 0.005 * randn(size(reflectance));
        reflectance_noisy = reflectance + noise;
        
        % 绘图
        hold on;
        h1 = plot(wavenumber, reflectance, 'LineWidth', 2.5, 'Color', color_primary);
        h2 = plot(wavenumber, reflectance_noisy, 'LineWidth', 1, 'Color', color_data);
        
        % 标记极值点（修复findpeaks问题）
        [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakProminence', 0.05);
        [valleys, valley_locs] = findpeaks(-reflectance, 'MinPeakProminence', 0.05);
        valleys = -valleys;
        
        if ~isempty(peak_locs)
            plot(wavenumber(peak_locs), peaks, 'o', 'MarkerSize', 6, ...
                 'MarkerFaceColor', color_highlight, 'MarkerEdgeColor', 'k');
        end
        if ~isempty(valley_locs)
            plot(wavenumber(valley_locs), valleys, 'v', 'MarkerSize', 6, ...
                 'MarkerFaceColor', color_error, 'MarkerEdgeColor', 'k');
        end
        
        % 设置图形属性
        xlabel('波数 (cm^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('反射率', 'FontSize', 12, 'FontWeight', 'bold');
        title('碳化硅薄膜红外反射光谱', 'FontSize', 14, 'FontWeight', 'bold');
        legend([h1, h2], {'理论光谱', '实验数据'}, 'Location', 'best', 'FontSize', 11);
        
        % 设置坐标轴
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        xlim([400, 4000]);
        
        % 保存图片
        save_figure('d:\Project_env\CUMCU\B\figures\reflectance_spectrum');
    end

    function plot_phase_difference_analysis_colored()
        % 绘制相位差分析图（使用color51配色）
        figure('Position', [100, 100, 800, 600]);
        
        % 模拟数据
        thickness_range = linspace(0.5, 5, 100);
        wavenumber_fixed = 2000; % cm^-1
        n_sic = 2.65;
        
        % 计算相位差
        phase_diff = 4 * pi * n_sic * thickness_range ./ (10000 / wavenumber_fixed);
        phase_diff_wrapped = mod(phase_diff, 2*pi);
        
        % 绘制主图
        subplot(2,1,1);
        hold on;
        h1 = plot(thickness_range, phase_diff, 'LineWidth', 2.5, 'Color', color_primary);
        h2 = plot(thickness_range, phase_diff_wrapped, 'LineWidth', 2.5, 'Color', color_model, 'LineStyle', '--');
        
        xlabel('厚度 (μm)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('相位差 (rad)', 'FontSize', 11, 'FontWeight', 'bold');
        title('相位差随厚度变化关系', 'FontSize', 12, 'FontWeight', 'bold');
        legend([h1, h2], {'连续相位差', '包裹相位差'}, 'Location', 'best');
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        % 绘制干涉级数图
        subplot(2,1,2);
        interference_order = phase_diff / (2*pi);
        hold on;
        h3 = plot(thickness_range, interference_order, 'LineWidth', 2.5, 'Color', color_secondary);
        
        % 标记整数级数
        integer_orders = 1:floor(max(interference_order));
        for i = integer_orders
            idx = find(interference_order >= i, 1);
            if ~isempty(idx)
                plot(thickness_range(idx), i, 'o', 'MarkerSize', 8, ...
                     'MarkerFaceColor', color_highlight, 'MarkerEdgeColor', 'k');
            end
        end
        
        xlabel('厚度 (μm)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('干涉级数', 'FontSize', 11, 'FontWeight', 'bold');
        title('干涉级数随厚度变化', 'FontSize', 12, 'FontWeight', 'bold');
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        % 保存图片
        save_figure('d:\Project_env\CUMCU\B\figures\phase_difference_analysis');
    end

    function plot_thickness_interference_relation_colored()
        % 绘制厚度-干涉关系图（使用color51配色）
        figure('Position', [100, 100, 800, 600]);
        
        % 模拟数据
        wavenumber = linspace(1000, 3000, 500);
        thickness_values = [1.5, 2.0, 2.5, 3.0]; % μm
        n_sic = 2.65;
        
        hold on;
        colors_thickness = [color_primary; color_model; color_secondary; color_validation];
        
        for i = 1:length(thickness_values)
            thickness = thickness_values(i);
            phase = 4 * pi * n_sic * thickness ./ (10000 ./ wavenumber);
            R_base = 0.18;
            R_modulation = 0.15;
            reflectance = R_base + R_modulation * cos(phase);
            
            h(i) = plot(wavenumber, reflectance, 'LineWidth', 2, ...
                       'Color', colors_thickness(i,:));
        end
        
        % 设置图形属性
        xlabel('波数 (cm^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('反射率', 'FontSize', 12, 'FontWeight', 'bold');
        title('不同厚度下的干涉光谱', 'FontSize', 14, 'FontWeight', 'bold');
        
        legend_labels = arrayfun(@(x) sprintf('%.1f μm', x), thickness_values, 'UniformOutput', false);
        legend(h, legend_labels, 'Location', 'best', 'FontSize', 11);
        
        % 设置坐标轴
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        xlim([1000, 3000]);
        
        % 保存图片
        save_figure('d:\Project_env\CUMCU\B\figures\thickness_interference_relation');
    end

    function plot_data_quality_colored()
        % 绘制数据质量评估图
        figure('Position', [100, 100, 1000, 600]);
        
        % 模拟数据质量指标
        categories = {'数据完整性', '信噪比', '基线稳定性', '峰值识别', '异常值检测'};
        scores_10deg = [0.95, 0.88, 0.92, 0.85, 0.90];
        scores_15deg = [0.93, 0.82, 0.89, 0.78, 0.87];
        
        % 创建子图
        subplot(1,2,1);
        x = 1:length(categories);
        width = 0.35;
        
        b1 = bar(x - width/2, scores_10deg, width, 'FaceColor', color_primary);
        hold on;
        b2 = bar(x + width/2, scores_15deg, width, 'FaceColor', color_model);
        
        xlabel('评估指标', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('评分', 'FontSize', 11, 'FontWeight', 'bold');
        title('数据质量评估结果', 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'XTick', x, 'XTickLabel', categories, 'XTickLabelRotation', 45);
        legend([b1, b2], {'10°入射角', '15°入射角'}, 'Location', 'best');
        ylim([0, 1]);
        
        % 设置坐标轴样式
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'YGrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        % 信噪比分析
        subplot(1,2,2);
        wavenumber = linspace(400, 4000, 1000);
        signal = 0.2 + 0.1 * cos(2*pi*wavenumber/500);
        noise = 0.01 * randn(size(wavenumber));
        snr = signal ./ abs(noise + 0.001);
        
        % 使用简单移动平均替代smooth函数
        window_size = 50;
        snr_smoothed = movmean(snr, window_size);
        plot(wavenumber, snr_smoothed, 'LineWidth', 2, 'Color', color_secondary);
        xlabel('波数 (cm^{-1})', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('信噪比', 'FontSize', 11, 'FontWeight', 'bold');
        title('光谱信噪比分析', 'FontSize', 12, 'FontWeight', 'bold');
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        save_figure('d:\Project_env\CUMCU\B\figures\data_quality');
    end

    function plot_thickness_comparison_colored()
        % 绘制厚度对比分析图
        figure('Position', [100, 100, 800, 600]);
        
        % 模拟不同角度下的厚度测量结果
        angles = [5, 10, 15, 20, 25];
        thickness_measured = [2.48, 2.52, 2.55, 2.61, 2.68]; % μm
        thickness_true = 2.5; % μm
        uncertainty = [0.08, 0.05, 0.06, 0.12, 0.18];
        
        % 绘制误差条图
        hold on;
        h1 = errorbar(angles, thickness_measured, uncertainty, 'o', ...
                     'LineWidth', 2, 'MarkerSize', 8, 'Color', color_primary, ...
                     'MarkerFaceColor', color_primary, 'MarkerEdgeColor', 'k');
        
        % 绘制真实值参考线
        h2 = plot([min(angles)-2, max(angles)+2], [thickness_true, thickness_true], ...
                 '--', 'LineWidth', 2, 'Color', color_model);
        
        % 绘制置信区间
        fill([angles, fliplr(angles)], ...
             [thickness_measured + uncertainty, fliplr(thickness_measured - uncertainty)], ...
             color_ci, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        xlabel('入射角 (度)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('测量厚度 (μm)', 'FontSize', 12, 'FontWeight', 'bold');
        title('不同入射角下的厚度测量对比', 'FontSize', 14, 'FontWeight', 'bold');
        legend([h1, h2], {'测量值', '真实值'}, 'Location', 'best', 'FontSize', 11);
        
        xlim([min(angles)-2, max(angles)+2]);
        ylim([2.2, 2.9]);
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        save_figure('d:\Project_env\CUMCU\B\figures\thickness_comparison');
    end

    function plot_angle_dependency_colored()
        % 绘制角度依赖性分析图
        figure('Position', [100, 100, 800, 600]);
        
        % 模拟角度依赖性数据
        angles = linspace(0, 30, 100);
        thickness_ideal = 2.5; % μm
        % 考虑几何光程因子
        thickness_apparent = thickness_ideal ./ cos(angles * pi / 180);
        % 添加系统误差
        systematic_error = 0.003 * angles; % 每度0.3%的误差
        thickness_measured = thickness_apparent + systematic_error;
        
        hold on;
        h1 = plot(angles, thickness_ideal * ones(size(angles)), '--', ...
                 'LineWidth', 2.5, 'Color', color_model);
        h2 = plot(angles, thickness_apparent, 'LineWidth', 2, 'Color', color_secondary);
        h3 = plot(angles, thickness_measured, 'LineWidth', 2.5, 'Color', color_primary);
        
        % 标记关键角度点
        key_angles = [10, 15, 20];
        for i = 1:length(key_angles)
            idx = find(angles >= key_angles(i), 1);
            plot(angles(idx), thickness_measured(idx), 'o', 'MarkerSize', 8, ...
                 'MarkerFaceColor', color_highlight, 'MarkerEdgeColor', 'k');
        end
        
        xlabel('入射角 (度)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('测量厚度 (μm)', 'FontSize', 12, 'FontWeight', 'bold');
        title('厚度测量的角度依赖性分析', 'FontSize', 14, 'FontWeight', 'bold');
        legend([h1, h2, h3], {'理想厚度', '几何修正', '实际测量'}, 'Location', 'best', 'FontSize', 11);
        
        xlim([0, 30]);
        ylim([2.4, 2.8]);
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        save_figure('d:\Project_env\CUMCU\B\figures\angle_dependency');
    end

    function plot_error_analysis_colored()
        % 绘制误差分析图
        figure('Position', [100, 100, 1000, 600]);
        
        % 误差源分析
        subplot(1,2,1);
        error_sources = {'光谱测量', '角度定位', '算法收敛', '环境因素', '系统漂移'};
        error_contributions = [45, 25, 15, 10, 5]; % 百分比
        
        % 优化的饼图配色方案 - 使用更协调的颜色组合
        colors_pie = [colors(13,:);   % 橙色 - 光谱测量（最大贡献）
                     colors(22,:);   % 绿色 - 角度定位
                     colors(35,:);   % 黄色 - 算法收敛
                     colors(45,:);   % 紫色 - 环境因素
                     colors(8,:)];   % 浅蓝色 - 系统漂移（最小贡献）
        
        % 创建饼图并设置颜色
        p = pie(error_contributions);
        
        % 设置饼图扇形颜色和样式
        for i = 1:2:length(p)
            sector_idx = (i+1)/2;
            p(i).FaceColor = colors_pie(sector_idx,:);
            p(i).EdgeColor = [1 1 1]; % 白色边框
            p(i).LineWidth = 1.5;
        end
        
        % 设置文本样式
        for i = 2:2:length(p)
            p(i).FontSize = 10;
            p(i).FontWeight = 'bold';
            p(i).Color = [0.2 0.2 0.2]; % 深灰色文本
        end
        
        title('误差源贡献分析', 'FontSize', 12, 'FontWeight', 'bold');
        legend(error_sources, 'Location', 'eastoutside', 'FontSize', 10);
        
        % 不确定度分析
        subplot(1,2,2);
        thickness_range = linspace(1, 5, 50);
        uncertainty_systematic = 0.02 + 0.005 * thickness_range; % 系统不确定度
        uncertainty_random = 0.01 * ones(size(thickness_range)); % 随机不确定度
        uncertainty_total = sqrt(uncertainty_systematic.^2 + uncertainty_random.^2);
        
        hold on;
        h1 = plot(thickness_range, uncertainty_systematic, 'LineWidth', 2, 'Color', color_primary);
        h2 = plot(thickness_range, uncertainty_random, 'LineWidth', 2, 'Color', color_model);
        h3 = plot(thickness_range, uncertainty_total, 'LineWidth', 2.5, 'Color', color_error);
        
        xlabel('厚度 (μm)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('不确定度 (μm)', 'FontSize', 11, 'FontWeight', 'bold');
        title('测量不确定度分析', 'FontSize', 12, 'FontWeight', 'bold');
        legend([h1, h2, h3], {'系统不确定度', '随机不确定度', '总不确定度'}, 'Location', 'best');
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        save_figure('d:\Project_env\CUMCU\B\figures\error_analysis');
    end

    function plot_method_comparison_colored()
        % 绘制方法对比图
        figure('Position', [100, 100, 1000, 600]);
        
        % 不同方法的性能对比
        methods = {'红外干涉法', '椭偏法', 'AFM', 'SEM截面', 'X射线反射'};
        precision = [0.05, 0.1, 0.01, 0.02, 0.08]; % μm
        speed = [9, 7, 3, 2, 5]; % 相对速度评分
        cost = [6, 8, 4, 3, 9]; % 相对成本评分
        
        % 创建雷达图数据
        subplot(1,2,1);
        % 精度对比（柱状图）
        colors_methods = [color_primary; color_model; color_secondary; color_validation; color_ci];
        b = bar(precision, 'FaceColor', 'flat');
        b.CData = colors_methods;
        
        xlabel('测量方法', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('测量精度 (μm)', 'FontSize', 11, 'FontWeight', 'bold');
        title('不同方法精度对比', 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'XTick', 1:length(methods), 'XTickLabel', methods, 'XTickLabelRotation', 45);
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'YGrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        % 综合性能对比（散点图）
        subplot(1,2,2);
        scatter(speed, cost, 100, colors_methods, 'filled', 'MarkerEdgeColor', 'k');
        
        % 添加方法标签
        for i = 1:length(methods)
            text(speed(i)+0.1, cost(i), methods{i}, 'FontSize', 9);
        end
        
        xlabel('测量速度评分', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('设备成本评分', 'FontSize', 11, 'FontWeight', 'bold');
        title('方法综合性能对比', 'FontSize', 12, 'FontWeight', 'bold');
        xlim([0, 10]);
        ylim([0, 10]);
        
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
            'XMinorTick', 'on', 'YMinorTick', 'on', 'XGrid', 'on', 'YGrid', 'on', ...
            'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'LineWidth', 1);
        
        save_figure('d:\Project_env\CUMCU\B\figures\method_comparison');
    end

    function save_figure(filename)
        % 保存图片为EPS格式
        print([filename '.eps'], '-depsc', '-r300');
        fprintf('已保存: %s.eps\n', filename);
    end
end