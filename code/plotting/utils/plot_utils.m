function varargout = plot_utils(varargin)
% PLOT_UTILS - 绘图工具函数集合
%
% 提供各种光谱数据和分析结果的可视化功能
% 包括光谱绘图、干涉条纹分析、厚度计算结果展示等
%
% 主要函数:
%   plot_spectrum - 绘制光谱数据
%   plot_interference_analysis - 绘制干涉分析结果
%   plot_thickness_results - 绘制厚度计算结果
%   plot_multi_beam_analysis - 绘制多光束干涉分析
%   plot_comparison - 绘制对比图
%


    if nargin == 0
        fprintf('绘图工具函数库已加载\n');
        fprintf('可用函数: plot_spectrum, plot_interference_analysis, plot_thickness_results\n');
        return;
    end
    
    % 根据输入参数调用相应的绘图函数
    func_name = varargin{1};
    switch func_name
        case 'spectrum'
            varargout{1} = plot_spectrum_internal(varargin{2:end});
        case 'interference'
            varargout{1} = plot_interference_analysis_internal(varargin{2:end});
        case 'thickness'
            varargout{1} = plot_thickness_results_internal(varargin{2:end});
        case 'multi_beam'
            varargout{1} = plot_multi_beam_analysis_internal(varargin{2:end});
        case 'comparison'
            varargout{1} = plot_comparison_internal(varargin{2:end});
        otherwise
            error('Unknown plot function: %s', func_name);
    end
end

%% 独立的光谱绘图函数
function fig_handle = plot_spectrum(spectrum_data, options)
% PLOT_SPECTRUM - 绘制光谱数据
%
% 输入参数:
%   spectrum_data - 光谱数据结构体或cell数组（多个数据集）
%   options - 绘图选项结构体
%
% 输出参数:
%   fig_handle - 图形句柄

    if nargin < 2
        options = struct();
    end
    
    % 默认选项
    default_options = struct(...
        'title', '光谱数据', ...
        'xlabel', '波数 (cm^{-1})', ...
        'ylabel', '反射率 (%)', ...
        'grid', true, ...
        'legend', true, ...
        'line_width', 1.5, ...
        'marker_size', 4, ...
        'colors', [], ...
        'save_path', 'd:\Project_env\CUMCU\B\figures\', ...
        'figure_size', [800, 600], ...
        'font_size', 12, ...
        'show_peaks', false, ...
        'show_statistics', false
    );
    
    % 合并选项
    options = merge_options(default_options, options);
    
    % 创建图形
    fig_handle = figure('Position', [100, 100, options.figure_size]);
    set(fig_handle, 'Color', 'white');
    
    % 处理单个或多个数据集
    if ~iscell(spectrum_data)
        spectrum_data = {spectrum_data};
    end
    
    num_datasets = length(spectrum_data);
    
    % 设置颜色
    if isempty(options.colors)
        colors = lines(num_datasets);
    else
        colors = options.colors;
    end
    
    hold on;
    legend_entries = {};
    
    for i = 1:num_datasets
        data = spectrum_data{i};
        
        % 绘制主曲线
        plot(data.wavenumber, data.reflectance, ...
            'Color', colors(i, :), ...
            'LineWidth', options.line_width, ...
            'DisplayName', sprintf('%s (%.1f°)', data.material, data.angle));
        
        legend_entries{end+1} = sprintf('%s (%.1f°)', data.material, data.angle);
        
        % 显示峰值
        if options.show_peaks
            [peaks, peak_locs] = findpeaks(data.reflectance, data.wavenumber, ...
                'MinPeakHeight', mean(data.reflectance) + std(data.reflectance), ...
                'MinPeakDistance', (max(data.wavenumber) - min(data.wavenumber)) / 20);
            
            if ~isempty(peaks)
                plot(peak_locs, peaks, 'o', ...
                    'Color', colors(i, :), ...
                    'MarkerSize', options.marker_size, ...
                    'MarkerFaceColor', colors(i, :), ...
                    'HandleVisibility', 'off');
                
                % 标注峰值
                for j = 1:length(peaks)
                    text(peak_locs(j), peaks(j) + max(data.reflectance) * 0.02, ...
                        sprintf('%.1f', peak_locs(j)), ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', options.font_size - 2, ...
                        'Color', colors(i, :));
                end
            end
        end
    end
    
    hold off;
    
    % 设置坐标轴
    xlabel(options.xlabel, 'FontSize', options.font_size);
    ylabel(options.ylabel, 'FontSize', options.font_size);
    title(options.title, 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    
    if options.grid
        grid on;
        grid minor;
    end
    
    if options.legend && num_datasets > 1
        legend(legend_entries, 'Location', 'best', 'FontSize', options.font_size - 1);
    end
    
    % 设置坐标轴范围
    xlim([min(cellfun(@(x) min(x.wavenumber), spectrum_data)), ...
          max(cellfun(@(x) max(x.wavenumber), spectrum_data))]);
    
    % 显示统计信息
    if options.show_statistics
        add_statistics_text(spectrum_data, options.font_size);
    end
    
    % 保存图形
    if ~isempty(options.save_path)
        saveas(fig_handle, options.save_path);
        fprintf('图形已保存至: %s\n', options.save_path);
    end
end

%% 内部光谱绘图函数
function fig_handle = plot_spectrum_internal(spectrum_data, options)
    fig_handle = plot_spectrum(spectrum_data, options);
end

%% 干涉分析结果绘图
function fig_handle = plot_interference_analysis(analysis_results, options)
% PLOT_INTERFERENCE_ANALYSIS - 绘制干涉分析结果
%
% 输入参数:
%   analysis_results - 干涉分析结果结构体
%   options - 绘图选项

    if nargin < 2
        options = struct();
    end
    
    default_options = struct(...
        'subplot_layout', [2, 2], ...
        'figure_size', [1200, 800], ...
        'font_size', 10, ...
        'save_path', ''
    );
    
    options = merge_options(default_options, options);
    
    fig_handle = figure('Position', [50, 50, options.figure_size]);
    set(fig_handle, 'Color', 'white');
    
    % 子图1: 原始光谱数据
    subplot(options.subplot_layout(1), options.subplot_layout(2), 1);
    plot(analysis_results.wavenumber, analysis_results.reflectance, 'b-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
    ylabel('反射率 (%)', 'FontSize', options.font_size);
    title('原始光谱数据', 'FontSize', options.font_size + 1);
    grid on;
    
    % 子图2: 理论拟合对比
    if isfield(analysis_results, 'theoretical_reflectance')
        subplot(options.subplot_layout(1), options.subplot_layout(2), 2);
        plot(analysis_results.wavenumber, analysis_results.reflectance, 'b-', 'LineWidth', 1.5, 'DisplayName', '实测');
        hold on;
        plot(analysis_results.wavenumber, analysis_results.theoretical_reflectance, 'r--', 'LineWidth', 1.5, 'DisplayName', '理论');
        hold off;
        xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
        ylabel('反射率 (%)', 'FontSize', options.font_size);
        title('理论拟合对比', 'FontSize', options.font_size + 1);
        legend('Location', 'best');
        grid on;
    end
    
    % 子图3: 相位分析
    if isfield(analysis_results, 'phase')
        subplot(options.subplot_layout(1), options.subplot_layout(2), 3);
        plot(analysis_results.wavenumber, analysis_results.phase, 'g-', 'LineWidth', 1.5);
        xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
        ylabel('相位 (rad)', 'FontSize', options.font_size);
        title('相位分析', 'FontSize', options.font_size + 1);
        grid on;
    end
    
    % 子图4: 残差分析
    if isfield(analysis_results, 'residual')
        subplot(options.subplot_layout(1), options.subplot_layout(2), 4);
        plot(analysis_results.wavenumber, analysis_results.residual, 'm-', 'LineWidth', 1.5);
        xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
        ylabel('残差', 'FontSize', options.font_size);
        title('拟合残差', 'FontSize', options.font_size + 1);
        grid on;
        
        % 添加零线
        hold on;
        plot([min(analysis_results.wavenumber), max(analysis_results.wavenumber)], [0, 0], 'k--', 'LineWidth', 0.5);
        hold off;
    end
    
    % 总标题
    sgtitle('干涉分析结果', 'FontSize', options.font_size + 3, 'FontWeight', 'bold');
    
    % 保存图形
    if ~isempty(options.save_path)
        saveas(fig_handle, options.save_path);
        fprintf('干涉分析图已保存至: %s\n', options.save_path);
    end
end

%% 内部干涉分析绘图函数
function fig_handle = plot_interference_analysis_internal(analysis_results, options)
    fig_handle = plot_interference_analysis(analysis_results, options);
end

%% 厚度计算结果绘图
function fig_handle = plot_thickness_results(thickness_results, options)
% PLOT_THICKNESS_RESULTS - 绘制厚度计算结果
%
% 输入参数:
%   thickness_results - 厚度计算结果结构体或数组
%   options - 绘图选项

    if nargin < 2
        options = struct();
    end
    
    default_options = struct(...
        'figure_size', [1000, 600], ...
        'font_size', 12, ...
        'show_error_bars', true, ...
        'show_statistics', true, ...
        'save_path', ''
    );
    
    options = merge_options(default_options, options);
    
    fig_handle = figure('Position', [100, 100, options.figure_size]);
    set(fig_handle, 'Color', 'white');
    
    % 处理单个或多个结果
    if ~iscell(thickness_results)
        thickness_results = {thickness_results};
    end
    
    num_results = length(thickness_results);
    
    % 提取厚度值和标签
    thickness_values = zeros(num_results, 1);
    labels = cell(num_results, 1);
    error_values = zeros(num_results, 1);
    
    for i = 1:num_results
        result = thickness_results{i};
        thickness_values(i) = result.thickness;
        
        if isfield(result, 'method')
            labels{i} = sprintf('%s\n%.1f°', result.method, result.angle);
        else
            labels{i} = sprintf('结果%d\n%.1f°', i, result.angle);
        end
        
        if isfield(result, 'uncertainty')
            error_values(i) = result.uncertainty;
        elseif isfield(result, 'std_error')
            error_values(i) = result.std_error;
        else
            error_values(i) = 0;
        end
    end
    
    % 绘制柱状图
    if options.show_error_bars && any(error_values > 0)
        bar_handle = bar(1:num_results, thickness_values, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'black');
        hold on;
        errorbar(1:num_results, thickness_values, error_values, 'k.', 'LineWidth', 1.5, 'MarkerSize', 10);
        hold off;
    else
        bar_handle = bar(1:num_results, thickness_values, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'black');
    end
    
    % 设置坐标轴
    set(gca, 'XTickLabel', labels, 'FontSize', options.font_size - 1);
    xlabel('测量条件/方法', 'FontSize', options.font_size);
    ylabel('厚度 (μm)', 'FontSize', options.font_size);
    title('外延层厚度计算结果', 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    
    % 添加数值标签
    for i = 1:num_results
        text(i, thickness_values(i) + max(thickness_values) * 0.02, ...
            sprintf('%.3f μm', thickness_values(i)), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', options.font_size - 1, ...
            'FontWeight', 'bold');
    end
    
    grid on;
    grid minor;
    
    % 显示统计信息
    if options.show_statistics && num_results > 1
        mean_thickness = mean(thickness_values);
        std_thickness = std(thickness_values);
        
        text_str = sprintf('平均厚度: %.3f μm\n标准偏差: %.3f μm\n变异系数: %.2f%%', ...
            mean_thickness, std_thickness, std_thickness/mean_thickness*100);
        
        text(0.02, 0.98, text_str, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', options.font_size - 1, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    
    % 保存图形
    if ~isempty(options.save_path)
        saveas(fig_handle, options.save_path);
        fprintf('厚度结果图已保存至: %s\n', options.save_path);
    end
end

%% 内部厚度结果绘图函数
function fig_handle = plot_thickness_results_internal(thickness_results, options)
    fig_handle = plot_thickness_results(thickness_results, options);
end

%% 多光束干涉分析绘图
function fig_handle = plot_multi_beam_analysis(multi_beam_results, options)
% PLOT_MULTI_BEAM_ANALYSIS - 绘制多光束干涉分析结果
%
% 输入参数:
%   multi_beam_results - 多光束干涉分析结果
%   options - 绘图选项

    if nargin < 2
        options = struct();
    end
    
    default_options = struct(...
        'figure_size', [1200, 900], ...
        'font_size', 10, ...
        'save_path', ''
    );
    
    options = merge_options(default_options, options);
    
    fig_handle = figure('Position', [50, 50, options.figure_size]);
    set(fig_handle, 'Color', 'white');
    
    % 子图1: 双光束 vs 多光束对比
    subplot(2, 3, 1);
    plot(multi_beam_results.wavenumber, multi_beam_results.measured_reflectance, 'b-', 'LineWidth', 1.5, 'DisplayName', '实测');
    hold on;
    if isfield(multi_beam_results, 'two_beam_theory')
        plot(multi_beam_results.wavenumber, multi_beam_results.two_beam_theory, 'r--', 'LineWidth', 1.5, 'DisplayName', '双光束理论');
    end
    if isfield(multi_beam_results, 'multi_beam_theory')
        plot(multi_beam_results.wavenumber, multi_beam_results.multi_beam_theory, 'g--', 'LineWidth', 1.5, 'DisplayName', '多光束理论');
    end
    hold off;
    xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
    ylabel('反射率', 'FontSize', options.font_size);
    title('理论模型对比', 'FontSize', options.font_size + 1);
    legend('Location', 'best');
    grid on;
    
    % 子图2: 精细度分析
    if isfield(multi_beam_results, 'finesse')
        subplot(2, 3, 2);
        plot(multi_beam_results.wavenumber, multi_beam_results.finesse, 'm-', 'LineWidth', 1.5);
        xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
        ylabel('精细度', 'FontSize', options.font_size);
        title('Fabry-Perot精细度', 'FontSize', options.font_size + 1);
        grid on;
    end
    
    % 子图3: 反射率分布
    if isfield(multi_beam_results, 'reflection_orders')
        subplot(2, 3, 3);
        orders = multi_beam_results.reflection_orders;
        contributions = multi_beam_results.order_contributions;
        
        bar(orders, contributions, 'FaceColor', [0.7, 0.3, 0.3]);
        xlabel('反射次数', 'FontSize', options.font_size);
        ylabel('贡献度', 'FontSize', options.font_size);
        title('各阶反射贡献', 'FontSize', options.font_size + 1);
        grid on;
    end
    
    % 子图4: 相位关系
    if isfield(multi_beam_results, 'phase_difference')
        subplot(2, 3, 4);
        plot(multi_beam_results.wavenumber, multi_beam_results.phase_difference, 'c-', 'LineWidth', 1.5);
        xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
        ylabel('相位差 (rad)', 'FontSize', options.font_size);
        title('相位差分析', 'FontSize', options.font_size + 1);
        grid on;
    end
    
    % 子图5: 误差分析
    if isfield(multi_beam_results, 'two_beam_error') && isfield(multi_beam_results, 'multi_beam_error')
        subplot(2, 3, 5);
        plot(multi_beam_results.wavenumber, multi_beam_results.two_beam_error, 'r-', 'LineWidth', 1.5, 'DisplayName', '双光束误差');
        hold on;
        plot(multi_beam_results.wavenumber, multi_beam_results.multi_beam_error, 'g-', 'LineWidth', 1.5, 'DisplayName', '多光束误差');
        hold off;
        xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
        ylabel('拟合误差', 'FontSize', options.font_size);
        title('模型误差对比', 'FontSize', options.font_size + 1);
        legend('Location', 'best');
        grid on;
    end
    
    % 子图6: 厚度修正效果
    if isfield(multi_beam_results, 'thickness_correction')
        subplot(2, 3, 6);
        correction = multi_beam_results.thickness_correction;
        
        bar([1, 2], [correction.original_thickness, correction.corrected_thickness], ...
            'FaceColor', [0.4, 0.7, 0.4]);
        set(gca, 'XTickLabel', {'修正前', '修正后'});
        ylabel('厚度 (μm)', 'FontSize', options.font_size);
        title('厚度修正效果', 'FontSize', options.font_size + 1);
        
        % 添加数值标签
        text(1, correction.original_thickness + 0.1, sprintf('%.3f', correction.original_thickness), ...
            'HorizontalAlignment', 'center', 'FontSize', options.font_size);
        text(2, correction.corrected_thickness + 0.1, sprintf('%.3f', correction.corrected_thickness), ...
            'HorizontalAlignment', 'center', 'FontSize', options.font_size);
        
        grid on;
    end
    
    % 总标题
    sgtitle('多光束干涉分析结果', 'FontSize', options.font_size + 3, 'FontWeight', 'bold');
    
    % 保存图形
    if ~isempty(options.save_path)
        saveas(fig_handle, options.save_path);
        fprintf('多光束分析图已保存至: %s\n', options.save_path);
    end
end

%% 内部多光束分析绘图函数
function fig_handle = plot_multi_beam_analysis_internal(multi_beam_results, options)
    fig_handle = plot_multi_beam_analysis(multi_beam_results, options);
end

%% 对比图绘制
function fig_handle = plot_comparison(data_sets, labels, options)
% PLOT_COMPARISON - 绘制多个数据集的对比图
%
% 输入参数:
%   data_sets - 数据集cell数组
%   labels - 标签cell数组
%   options - 绘图选项

    if nargin < 3
        options = struct();
    end
    
    default_options = struct(...
        'figure_size', [1000, 700], ...
        'font_size', 12, ...
        'colors', [], ...
        'line_styles', {'-', '--', '-.', ':'}, ...
        'save_path', '', ...
        'normalize', false
    );
    
    options = merge_options(default_options, options);
    
    fig_handle = figure('Position', [100, 100, options.figure_size]);
    set(fig_handle, 'Color', 'white');
    
    num_sets = length(data_sets);
    
    if isempty(options.colors)
        colors = lines(num_sets);
    else
        colors = options.colors;
    end
    
    hold on;
    
    for i = 1:num_sets
        data = data_sets{i};
        
        % 数据归一化
        if options.normalize
            y_data = (data.reflectance - min(data.reflectance)) / ...
                     (max(data.reflectance) - min(data.reflectance));
        else
            y_data = data.reflectance;
        end
        
        line_style = options.line_styles{mod(i-1, length(options.line_styles)) + 1};
        
        plot(data.wavenumber, y_data, line_style, ...
            'Color', colors(i, :), ...
            'LineWidth', 1.5, ...
            'DisplayName', labels{i});
    end
    
    hold off;
    
    xlabel('波数 (cm^{-1})', 'FontSize', options.font_size);
    if options.normalize
        ylabel('归一化反射率', 'FontSize', options.font_size);
    else
        ylabel('反射率 (%)', 'FontSize', options.font_size);
    end
    
    title('光谱数据对比', 'FontSize', options.font_size + 2, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', options.font_size - 1);
    grid on;
    grid minor;
    
    % 保存图形
    if ~isempty(options.save_path)
        saveas(fig_handle, options.save_path);
        fprintf('对比图已保存至: %s\n', options.save_path);
    end
end

%% 内部对比图绘制函数
function fig_handle = plot_comparison_internal(data_sets, labels, options)
    fig_handle = plot_comparison(data_sets, labels, options);
end

%% 辅助函数：合并选项
function merged_options = merge_options(default_options, user_options)
    merged_options = default_options;
    
    if ~isempty(user_options)
        field_names = fieldnames(user_options);
        for i = 1:length(field_names)
            merged_options.(field_names{i}) = user_options.(field_names{i});
        end
    end
end

%% 辅助函数：添加统计信息文本
function add_statistics_text(spectrum_data, font_size)
    if length(spectrum_data) == 1
        data = spectrum_data{1};
        stats_text = sprintf('数据点数: %d\n波数范围: %.1f - %.1f cm^{-1}\n反射率: %.2f ± %.2f%%', ...
            data.num_points, min(data.wavenumber), max(data.wavenumber), ...
            mean(data.reflectance), std(data.reflectance));
        
        text(0.02, 0.98, stats_text, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', font_size - 2, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
end