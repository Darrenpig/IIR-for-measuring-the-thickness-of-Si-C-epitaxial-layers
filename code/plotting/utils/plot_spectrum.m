function fig_handle = plot_spectrum(spectrum_data, options)
% PLOT_SPECTRUM - 绘制光谱数据
%
% 输入参数:
%   spectrum_data - 光谱数据结构体或cell数组（多个数据集）
%   options - 绘图选项结构体
%
% 输出参数:
%   fig_handle - 图形句柄
%


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
        'save_path', '', ...
        'figure_size', [800, 600], ...
        'font_size', 12, ...
        'show_peaks', false, ...
        'show_statistics', false ...
    );
    
    % 合并选项
    options = merge_options(default_options, options);
    
    % 创建图形
    fig_handle = figure('Position', [100, 100, options.figure_size(1), options.figure_size(2)]);
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
            'Color', colors(i,:), ...
            'LineWidth', options.line_width, ...
            'DisplayName', sprintf('数据集 %d', i));
        
        legend_entries{end+1} = sprintf('数据集 %d', i);
        
        % 显示峰值（如果需要）
        if options.show_peaks
            [peaks, locs] = findpeaks(data.reflectance, 'MinPeakHeight', max(data.reflectance)*0.1);
            if ~isempty(peaks)
                plot(data.wavenumber(locs), peaks, 'o', ...
                    'Color', colors(i,:), ...
                    'MarkerSize', options.marker_size, ...
                    'MarkerFaceColor', colors(i,:));
            end
        end
    end
    
    % 设置图形属性
    xlabel(options.xlabel, 'FontSize', options.font_size);
    ylabel(options.ylabel, 'FontSize', options.font_size);
    title(options.title, 'FontSize', options.font_size + 2);
    
    if options.grid
        grid on;
    end
    
    if options.legend && num_datasets > 1
        legend(legend_entries, 'Location', 'best');
    end
    
    % 显示统计信息（如果需要）
    if options.show_statistics
        add_statistics_text(spectrum_data{1}, gca);
    end
    
    % 保存图形（如果指定路径）
    if ~isempty(options.save_path)
        % 确保输出目录存在
        output_dir = 'd:\Project_env\CUMCU\B\figures';
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        % 如果save_path是相对路径，则使用新的输出目录
        if ~isabs(options.save_path)
            full_save_path = fullfile(output_dir, options.save_path);
        else
            full_save_path = options.save_path;
        end
        
        saveas(fig_handle, full_save_path);
    end
    
    hold off;
    
end

%% 辅助函数：合并选项
function merged_options = merge_options(default_options, user_options)
    merged_options = default_options;
    if ~isempty(user_options)
        fields = fieldnames(user_options);
        for i = 1:length(fields)
            merged_options.(fields{i}) = user_options.(fields{i});
        end
    end
end

%% 辅助函数：添加统计文本
function add_statistics_text(data, ax)
    % 计算统计信息
    mean_val = mean(data.reflectance);
    std_val = std(data.reflectance);
    min_val = min(data.reflectance);
    max_val = max(data.reflectance);
    
    % 创建统计文本
    stats_text = sprintf('统计信息:\n平均值: %.2f\n标准差: %.2f\n最小值: %.2f\n最大值: %.2f', ...
        mean_val, std_val, min_val, max_val);
    
    % 添加文本框
    text(0.02, 0.98, stats_text, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'white', 'EdgeColor', 'black', ...
        'FontSize', 10, 'Parent', ax);
end