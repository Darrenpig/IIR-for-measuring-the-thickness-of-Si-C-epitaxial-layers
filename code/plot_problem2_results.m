function plot_problem2_results()
% PLOT_PROBLEM2_RESULTS - 绘制问题二的分析结果图表
%
% 该函数读取问题二的计算结果并生成多种分析图表
% 包括厚度对比图、角度依赖性分析、数据质量评估等
%
% 输出格式: PDF, EPS, PNG
% 样式: Nature期刊风格，中文标注
%
% 作者: CUMCU数学建模团队
% 日期: 2024

    clc; close all;
    
    % 添加color目录到路径
    addpath('../color');
    
    fprintf('\n=== 开始绘制问题二结果图表 ===\n');
    
    % 设置图表样式
    setup_plot_style();
    
    % 创建输出目录
    output_dir = '../figures';
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    try
        % 加载计算结果
        fprintf('1. 加载计算结果...\n');
        load('results/problem2_results.mat', 'problem2_results');
        
        % 检查数据结构
        fprintf('数据结构检查完成\n');
        fprintf('角度数量: %d\n', length(problem2_results.angles));
        fprintf('厚度范围: %.1f - %.1f μm\n', min(problem2_results.thickness_values), max(problem2_results.thickness_values));
        
        % 1. 厚度计算结果对比图
        fprintf('2. 绘制厚度计算结果对比图...\n');
        plot_thickness_comparison(problem2_results, output_dir);
        
        % 2. 角度依赖性分析图
        fprintf('3. 绘制角度依赖性分析图...\n');
        plot_angle_dependency(problem2_results, output_dir);
        
        % 3. 计算方法对比图
        fprintf('4. 绘制计算方法对比图...\n');
        plot_method_comparison(problem2_results, output_dir);
        
        % 4. 数据质量评估图
        fprintf('5. 绘制数据质量评估图...\n');
        plot_data_quality(problem2_results, output_dir);
        
        % 5. 误差分析图
        fprintf('6. 绘制误差分析图...\n');
        plot_error_analysis(problem2_results, output_dir);
        
        fprintf('\n=== 所有图表绘制完成 ===\n');
        fprintf('图表保存位置: %s\n', output_dir);
        
    catch ME
        fprintf('绘图过程中出现错误: %s\n', ME.message);
        fprintf('错误位置: %s (第%d行)\n', ME.stack(1).name, ME.stack(1).line);
    end
    
end

function setup_plot_style()
% 设置Nature期刊风格的图表样式
    
    % 设置默认字体
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultAxesFontSize', 12);
    set(0, 'DefaultTextFontName', 'Arial');
    set(0, 'DefaultTextFontSize', 12);
    
    % 设置图表尺寸
    set(0, 'DefaultFigurePosition', [100, 100, 800, 600]);
    set(0, 'DefaultFigurePaperPositionMode', 'auto');
    
    % 使用color51配色方案
    try
        colors = [color51(5); color51(13); color51(51); color51(45); color51(22); color51(30)];
        set(0, 'DefaultAxesColorOrder', colors);
    catch
        % 如果color51不可用，使用默认配色
        colors = [0.2, 0.4, 0.8;    % 蓝色
                  0.8, 0.2, 0.2;    % 红色
                  0.2, 0.8, 0.2;    % 绿色
                  0.8, 0.6, 0.2;    % 橙色
                  0.6, 0.2, 0.8;    % 紫色
                  0.2, 0.8, 0.8];   % 青色
        set(0, 'DefaultAxesColorOrder', colors);
    end
    
end

function plot_thickness_comparison(results, output_dir)
% 绘制厚度计算结果对比图
    
    figure('Name', '厚度计算结果对比');
    
    angles = results.angles;
    thickness_values = results.thickness_values;
    
    % 主图：厚度随角度变化
    subplot(2, 2, 1);
    plot(angles, thickness_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('入射角 (°)');
    ylabel('厚度 (μm)');
    title('厚度随入射角变化');
    grid on;
    
    % 添加平均值线
    hold on;
    yline(results.mean_thickness, '--r', 'LineWidth', 1.5, ...
          'Label', sprintf('平均值: %.0f μm', results.mean_thickness));
    
    % 子图：厚度分布直方图
    subplot(2, 2, 2);
    histogram(thickness_values, 'NumBins', 5, ...
              'FaceColor', [0.2, 0.4, 0.8], 'FaceAlpha', 0.7);
    xlabel('厚度 (μm)');
    ylabel('频次');
    title('厚度分布直方图');
    
    % 子图：厚度误差棒图
    subplot(2, 2, 3);
    errorbar(angles, thickness_values, results.std_thickness * ones(size(angles)), ...
             'o-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('入射角 (°)');
    ylabel('厚度 (μm)');
    title('厚度测量不确定度');
    grid on;
    
    % 子图：统计信息
    subplot(2, 2, 4);
    stats_data = [results.mean_thickness; results.std_thickness; ...
                  results.thickness_range(2) - results.thickness_range(1)];
    stats_labels = {'平均值', '标准差', '范围'};
    bar(stats_data, 'FaceColor', [0.8, 0.2, 0.2], 'FaceAlpha', 0.7);
    set(gca, 'XTickLabel', stats_labels);
    ylabel('厚度 (μm)');
    title('统计信息汇总');
    
    % 调整布局
    sgtitle('问题二：SiC外延层厚度计算结果对比', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 保存图表
    save_figure(gcf, fullfile(output_dir, 'thickness_comparison'), {'pdf', 'eps', 'png'});
    
end

function plot_angle_dependency(results, output_dir)
% 绘制角度依赖性分析图
    
    figure('Name', '角度依赖性分析');
    
    angles = results.angles;
    thickness_values = results.thickness_values;
    
    % 主图：厚度-角度关系
    subplot(2, 1, 1);
    plot(angles, thickness_values, 'o-', 'LineWidth', 3, 'MarkerSize', 10, ...
         'Color', [0.2, 0.4, 0.8], 'MarkerFaceColor', [0.2, 0.4, 0.8]);
    
    % 拟合趋势线
    p = polyfit(angles, thickness_values, 1);
    fit_line = polyval(p, angles);
    hold on;
    plot(angles, fit_line, '--r', 'LineWidth', 2, ...
         'DisplayName', sprintf('拟合线: y = %.2fx + %.2f', p(1), p(2)));
    
    xlabel('入射角 (°)', 'FontSize', 14);
    ylabel('厚度 (μm)', 'FontSize', 14);
    title('厚度随入射角的变化趋势', 'FontSize', 16);
    legend('Location', 'best');
    grid on;
    
    % 计算相关系数
    R = corrcoef(angles, thickness_values);
    correlation = R(1, 2);
    
    % 添加文本注释
    text(min(angles) + 0.1 * range(angles), max(thickness_values) - 0.1 * range(thickness_values), ...
         sprintf('相关系数: %.3f', correlation), 'FontSize', 12, ...
         'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % 子图：角度对厚度的影响程度
    subplot(2, 1, 2);
    angle_effect = abs(thickness_values - mean(thickness_values));
    bar(angles, angle_effect, 'FaceColor', [0.8, 0.6, 0.2], 'FaceAlpha', 0.7);
    xlabel('入射角 (°)', 'FontSize', 14);
    ylabel('厚度偏差 (μm)', 'FontSize', 14);
    title('各角度下的厚度偏差', 'FontSize', 16);
    grid on;
    
    % 调整布局
    sgtitle('问题二：入射角对厚度测量的影响分析', 'FontSize', 18, 'FontWeight', 'bold');
    
    % 保存图表
    save_figure(gcf, fullfile(output_dir, 'angle_dependency'), {'pdf', 'eps', 'png'});
    
end

function plot_method_comparison(results, output_dir)
% 绘制计算方法对比图
    
    figure('Name', '计算方法对比');
    
    % 提取各方法的结果
    angles = results.angles;
    n_angles = length(angles);
    
    % 初始化数据矩阵
    extrema_results = [];
    phase_results = [];
    fft_results = [];
    
    for i = 1:n_angles
        angle_key = sprintf('angle_%d', angles(i));
        if isfield(results, angle_key)
            angle_data = results.(angle_key);
            
            if isfield(angle_data, 'extrema_method') && angle_data.extrema_method.success
                extrema_results = [extrema_results, angle_data.extrema_method.thickness];
            else
                extrema_results = [extrema_results, NaN];
            end
            
            if isfield(angle_data, 'phase_fitting') && angle_data.phase_fitting.success
                phase_results = [phase_results, angle_data.phase_fitting.thickness];
            else
                phase_results = [phase_results, NaN];
            end
            
            if isfield(angle_data, 'fft_method') && angle_data.fft_method.success
                fft_results = [fft_results, angle_data.fft_method.thickness];
            else
                fft_results = [fft_results, NaN];
            end
        end
    end
    
    % 主图：各方法结果对比
    subplot(2, 2, 1);
    hold on;
    if ~all(isnan(extrema_results))
        plot(angles, extrema_results, 'o-', 'LineWidth', 2, 'DisplayName', '极值法');
    end
    if ~all(isnan(phase_results))
        plot(angles, phase_results, 's-', 'LineWidth', 2, 'DisplayName', '相位拟合法');
    end
    if ~all(isnan(fft_results))
        plot(angles, fft_results, '^-', 'LineWidth', 2, 'DisplayName', 'FFT法');
    end
    
    xlabel('入射角 (°)');
    ylabel('厚度 (μm)');
    title('各计算方法结果对比');
    legend('Location', 'best');
    grid on;
    
    % 子图：方法成功率
    subplot(2, 2, 2);
    success_rates = [];
    method_names = {};
    
    if ~all(isnan(extrema_results))
        success_rates = [success_rates, sum(~isnan(extrema_results)) / n_angles * 100];
        method_names = [method_names, '极值法'];
    end
    if ~all(isnan(phase_results))
        success_rates = [success_rates, sum(~isnan(phase_results)) / n_angles * 100];
        method_names = [method_names, '相位拟合法'];
    end
    if ~all(isnan(fft_results))
        success_rates = [success_rates, sum(~isnan(fft_results)) / n_angles * 100];
        method_names = [method_names, 'FFT法'];
    end
    
    bar(success_rates, 'FaceColor', [0.2, 0.8, 0.2], 'FaceAlpha', 0.7);
    set(gca, 'XTickLabel', method_names);
    ylabel('成功率 (%)');
    title('各方法计算成功率');
    ylim([0, 105]);
    
    % 子图：方法精度对比（标准差）
    subplot(2, 2, 3);
    method_stds = [];
    if ~all(isnan(extrema_results))
        method_stds = [method_stds, nanstd(extrema_results)];
    end
    if ~all(isnan(phase_results))
        method_stds = [method_stds, nanstd(phase_results)];
    end
    if ~all(isnan(fft_results))
        method_stds = [method_stds, nanstd(fft_results)];
    end
    
    bar(method_stds, 'FaceColor', [0.8, 0.2, 0.8], 'FaceAlpha', 0.7);
    set(gca, 'XTickLabel', method_names);
    ylabel('标准差 (μm)');
    title('各方法精度对比');
    
    % 子图：综合评分
    subplot(2, 2, 4);
    % 简单的综合评分：成功率权重0.6，精度权重0.4
    if ~isempty(success_rates) && ~isempty(method_stds)
        normalized_success = success_rates / max(success_rates);
        normalized_precision = 1 - (method_stds / max(method_stds));
        composite_scores = 0.6 * normalized_success + 0.4 * normalized_precision;
        
        bar(composite_scores * 100, 'FaceColor', [0.6, 0.2, 0.8], 'FaceAlpha', 0.7);
        set(gca, 'XTickLabel', method_names);
        ylabel('综合评分');
        title('方法综合评价');
        ylim([0, 105]);
    end
    
    % 调整布局
    sgtitle('问题二：厚度计算方法对比分析', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 保存图表
    save_figure(gcf, fullfile(output_dir, 'method_comparison'), {'pdf', 'eps', 'png'});
    
end

function plot_data_quality(results, output_dir)
% 绘制数据质量评估图
    
    figure('Name', '数据质量评估');
    
    % 获取输入数据
    wavenumber = results.input_data.wavenumber;
    reflectance = results.input_data.reflectance;
    
    % 子图1：原始数据展示
    subplot(2, 3, 1);
    plot(wavenumber, reflectance, 'b-', 'LineWidth', 1);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率 (%)');
    title('原始光谱数据');
    grid on;
    
    % 子图2：数据分布直方图
    subplot(2, 3, 2);
    histogram(reflectance, 30, 'FaceColor', [0.2, 0.4, 0.8], 'FaceAlpha', 0.7);
    xlabel('反射率 (%)');
    ylabel('频次');
    title('反射率分布');
    
    % 子图3：数据完整性
    subplot(2, 3, 3);
    data_completeness = sum(~isnan(reflectance)) / length(reflectance) * 100;
    data_range_coverage = (max(wavenumber) - min(wavenumber)) / 4000 * 100;  % 假设理想范围4000 cm^-1
    
    quality_metrics = [data_completeness, data_range_coverage];
    quality_labels = {'数据完整性', '波数覆盖率'};
    
    bar(quality_metrics, 'FaceColor', [0.2, 0.8, 0.2], 'FaceAlpha', 0.7);
    set(gca, 'XTickLabel', quality_labels);
    ylabel('百分比 (%)');
    title('数据质量指标');
    ylim([0, 105]);
    
    % 子图4：信噪比估计
    subplot(2, 3, 4);
    % 简单的信噪比估计：信号均值/噪声标准差
    signal_mean = mean(reflectance);
    noise_std = std(diff(reflectance));  % 使用差分估计噪声
    snr_estimate = signal_mean / noise_std;
    
    bar(snr_estimate, 'FaceColor', [0.8, 0.6, 0.2], 'FaceAlpha', 0.7);
    ylabel('信噪比');
    title('估计信噪比');
    set(gca, 'XTickLabel', {'SNR'});
    
    % 子图5：数据平滑度
    subplot(2, 3, 5);
    smoothness = 1 / (std(diff(reflectance, 2)) + eps);  % 二阶差分的倒数作为平滑度指标
    
    bar(smoothness, 'FaceColor', [0.8, 0.2, 0.8], 'FaceAlpha', 0.7);
    ylabel('平滑度指标');
    title('数据平滑度');
    set(gca, 'XTickLabel', {'平滑度'});
    
    % 子图6：异常值检测
    subplot(2, 3, 6);
    % 使用3σ准则检测异常值
    mean_val = mean(reflectance);
    std_val = std(reflectance);
    outliers = abs(reflectance - mean_val) > 3 * std_val;
    outlier_rate = sum(outliers) / length(reflectance) * 100;
    
    pie([100 - outlier_rate, outlier_rate], {'正常数据', '异常值'});
    title(sprintf('异常值检测\n异常率: %.2f%%', outlier_rate));
    
    % 调整布局
    sgtitle('问题二：数据质量评估报告', 'FontSize', 16, 'FontWeight', 'bold');
    
    % 保存图表
    save_figure(gcf, fullfile(output_dir, 'data_quality'), {'pdf', 'eps', 'png'});
    
end

function plot_error_analysis(results, output_dir)
% 绘制误差分析图
    
    figure('Name', '误差分析');
    
    angles = results.angles;
    thickness_values = results.thickness_values;
    mean_thickness = results.mean_thickness;
    std_thickness = results.std_thickness;
    
    % 子图1：误差分布直方图
    subplot(2, 3, 1);
    errors = thickness_values - mean_thickness;
    try
        h1 = histogram(errors, 'BinWidth', std_thickness/3, ...
                      'FaceColor', color51(5), 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    catch
        h1 = histogram(errors, 'BinWidth', std_thickness/3, ...
                      'FaceColor', [0.2, 0.4, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    end
    xlabel('误差 (μm)', 'FontSize', 12);
    ylabel('频次', 'FontSize', 12);
    title('误差分布直方图', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 添加正态分布拟合
    hold on;
    x_fit = linspace(min(errors), max(errors), 100);
    y_fit = normpdf(x_fit, 0, std_thickness) * length(errors) * (max(errors) - min(errors)) / 10;
    try
        plot(x_fit, y_fit, 'Color', color51(13), 'LineWidth', 3, 'DisplayName', '正态拟合');
    catch
        plot(x_fit, y_fit, 'r-', 'LineWidth', 3, 'DisplayName', '正态拟合');
    end
    legend('Location', 'best', 'FontSize', 10);
    grid on; grid minor;
    
    % 子图2：置信区间
    subplot(2, 3, 2);
    confidence_levels = [68, 95, 99.7];  % 1σ, 2σ, 3σ
    confidence_intervals = [1, 2, 3] * std_thickness;
    
    try
        colors_bar = [color51(22); color51(30); color51(45)];
        b = bar(confidence_intervals, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        for i = 1:length(confidence_intervals)
            b.CData(i,:) = colors_bar(i,:);
        end
    catch
        bar(confidence_intervals, 'FaceColor', [0.2, 0.8, 0.2], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    end
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.1f%%', x), confidence_levels, 'UniformOutput', false));
    xlabel('置信水平', 'FontSize', 12);
    ylabel('置信区间 (μm)', 'FontSize', 12);
    title('测量置信区间', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; grid minor;
    
    % 子图3：相对误差
    subplot(2, 3, 3);
    relative_errors = abs(errors) / mean_thickness * 100;
    try
        plot(angles, relative_errors, 'o-', 'Color', color51(51), 'LineWidth', 2.5, ...
             'MarkerSize', 8, 'MarkerFaceColor', color51(51), 'MarkerEdgeColor', 'white', 'LineWidth', 1.5);
    catch
        plot(angles, relative_errors, 'o-', 'Color', [0.2, 0.8, 0.2], 'LineWidth', 2.5, ...
             'MarkerSize', 8, 'MarkerFaceColor', [0.2, 0.8, 0.2], 'MarkerEdgeColor', 'white');
    end
    xlabel('入射角 (°)', 'FontSize', 12);
    ylabel('相对误差 (%)', 'FontSize', 12);
    title('相对误差随角度变化', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; grid minor;
    
    % 子图4：累积误差分布
    subplot(2, 3, 4);
    sorted_errors = sort(abs(errors));
    cumulative_prob = (1:length(sorted_errors)) / length(sorted_errors) * 100;
    try
        plot(sorted_errors, cumulative_prob, 'Color', color51(13), 'LineWidth', 3);
        hold on;
        % 添加关键百分位点标记
        percentiles = [50, 90, 95, 99];
        for p = percentiles
            idx = round(p/100 * length(sorted_errors));
            if idx > 0 && idx <= length(sorted_errors)
                plot(sorted_errors(idx), p, 'o', 'Color', color51(45), ...
                     'MarkerSize', 8, 'MarkerFaceColor', color51(45), 'MarkerEdgeColor', 'white');
            end
        end
    catch
        plot(sorted_errors, cumulative_prob, 'b-', 'LineWidth', 3);
    end
    xlabel('绝对误差 (μm)', 'FontSize', 12);
    ylabel('累积概率 (%)', 'FontSize', 12);
    title('误差累积分布函数', 'FontSize', 14, 'FontWeight', 'bold');
    grid on; grid minor;
    
    % 子图5：不确定度分析
    subplot(2, 3, 5);
    % 不同来源的不确定度（示例）
    uncertainty_sources = {'测量重复性', '系统误差', '环境因素', '算法误差'};
    uncertainty_values = [std_thickness * 0.4, std_thickness * 0.3, ...
                         std_thickness * 0.2, std_thickness * 0.1];
    
    try
        pie_colors = [color51(5); color51(22); color51(30); color51(51)];
        p = pie(uncertainty_values, uncertainty_sources);
        % 设置饼图颜色
        for i = 1:2:length(p)
            p(i).FaceColor = pie_colors((i+1)/2,:);
            p(i).EdgeColor = 'white';
            p(i).LineWidth = 1.5;
        end
        % 设置文本样式
        for i = 2:2:length(p)
            p(i).FontSize = 10;
            p(i).FontWeight = 'bold';
        end
    catch
        pie(uncertainty_values, uncertainty_sources);
    end
    title('不确定度来源分析', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 子图6：误差统计摘要
    subplot(2, 3, 6);
    error_stats = [mean(abs(errors)), std(errors), max(abs(errors)), min(abs(errors))];
    error_labels = {'平均绝对误差', '误差标准差', '最大误差', '最小误差'};
    
    try
        % 使用不同颜色区分不同统计量
        bar_colors = [color51(5); color51(22); color51(45); color51(30)];
        b = bar(error_stats, 'FaceAlpha', 0.8, 'EdgeColor', 'none');
        b.CData = bar_colors;
        
        % 添加数值标签
        for i = 1:length(error_stats)
            text(i, error_stats(i) + max(error_stats)*0.02, ...
                 sprintf('%.3f', error_stats(i)), ...
                 'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        end
    catch
        bar(error_stats, 'FaceColor', [0.6, 0.2, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    end
    
    set(gca, 'XTickLabel', error_labels);
    ylabel('误差 (μm)', 'FontSize', 12);
    title('误差统计摘要', 'FontSize', 14, 'FontWeight', 'bold');
    xtickangle(45);
    grid on; grid minor;
    
    % 调整布局
    sgtitle('问题二：测量误差与不确定度分析', 'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.2, 0.2, 0.2]);
    
    % 保存图表
    save_figure(gcf, fullfile(output_dir, 'error_analysis'), {'pdf', 'eps', 'png'});
    
end



function save_figure(fig, filename, formats)
% 保存图表为多种格式
    
    for i = 1:length(formats)
        format = formats{i};
        full_filename = sprintf('%s.%s', filename, format);
        
        switch lower(format)
            case 'pdf'
                print(fig, full_filename, '-dpdf', '-r300');
            case 'eps'
                print(fig, full_filename, '-depsc', '-r300');
            case 'png'
                print(fig, full_filename, '-dpng', '-r300');
            otherwise
                warning('不支持的格式: %s', format);
        end
    end
    
    fprintf('图表已保存: %s\n', filename);
    
end