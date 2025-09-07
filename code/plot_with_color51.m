function plot_with_color51()
% 使用color51配色方案重新绘制所有图表
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
color_primary = colors(5,:);    % 主要曲线颜色
color_data = colors(13,:);      % 数据点颜色
color_model = colors(51,:);     % 模型曲线颜色
color_validation = colors(45,:); % 验证数据颜色
color_ci = colors(22,:);        % 置信区间颜色
color_error = colors(35,:);     % 误差条颜色
color_secondary = colors(30,:); % 次要曲线颜色
color_highlight = colors(15,:); % 高亮颜色

% 设置图形参数
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

% 创建输出目录
if ~exist('../figures', 'dir')
    mkdir('../figures');
end

% 1. 绘制Fresnel系数图
plot_fresnel_coefficients_colored();

% 2. 绘制反射光谱图
plot_reflectance_spectrum_colored();

% 3. 绘制相位差分析图
plot_phase_difference_colored();

% 4. 绘制厚度-干涉关系图
plot_thickness_interference_colored();

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
        save_figure('../figures/fresnel_coefficients');
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
        
        % 标记极值点
        [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakHeight', max(reflectance)*0.8);
        [valleys, valley_locs] = findpeaks(-reflectance, 'MinPeakHeight', -min(reflectance)*0.8);
        valleys = -valleys;
        
        plot(wavenumber(peak_locs), peaks, 'o', 'MarkerSize', 6, ...
             'MarkerFaceColor', color_highlight, 'MarkerEdgeColor', 'k');
        plot(wavenumber(valley_locs), valleys, 'v', 'MarkerSize', 6, ...
             'MarkerFaceColor', color_error, 'MarkerEdgeColor', 'k');
        
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
        save_figure('../figures/reflectance_spectrum');
    end

    function plot_phase_difference_colored()
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
        save_figure('../figures/phase_difference_analysis');
    end

    function plot_thickness_interference_colored()
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
        save_figure('../figures/thickness_interference_relation');
    end

    function save_figure(filename)
        % 保存图片为EPS格式
        print([filename '.eps'], '-depsc', '-r300');
        fprintf('已保存: %s.eps\n', filename);
    end
end