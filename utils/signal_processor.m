%% 信号处理工具函数
% 提供光谱数据的各种信号处理功能
% 包括滤波、平滑、去噪、特征提取等

%% 主函数：综合信号处理
function processed_data = signal_processor(data, processing_options)
    % 输入参数：
    %   data - 包含wavenumber和reflectance的结构体
    %   processing_options - 处理选项结构体
    % 输出参数：
    %   processed_data - 处理后的数据结构体
    
    if nargin < 2
        processing_options = struct();
    end
    
    % 默认处理选项
    default_options = struct(
        'smooth_method', 'sgolay',      % 平滑方法: 'sgolay', 'moving', 'lowess', 'gaussian'
        'smooth_window', 5,             % 平滑窗口大小
        'smooth_order', 3,              % Savitzky-Golay多项式阶数
        'detrend_method', 'none',       % 去趋势方法: 'none', 'linear', 'polynomial', 'baseline'
        'detrend_order', 2,             % 多项式去趋势阶数
        'filter_type', 'none',          % 滤波类型: 'none', 'lowpass', 'highpass', 'bandpass'
        'filter_cutoff', [0.1, 0.9],    % 滤波截止频率（归一化）
        'filter_order', 4,              % 滤波器阶数
        'denoise_method', 'none',       % 去噪方法: 'none', 'wavelet', 'median', 'wiener'
        'wavelet_type', 'db4',          % 小波类型
        'wavelet_levels', 4,            % 小波分解层数
        'outlier_method', 'none',       % 异常值处理: 'none', 'zscore', 'iqr', 'hampel'
        'outlier_threshold', 3,         % 异常值阈值
        'resample_method', 'none',      % 重采样方法: 'none', 'interp', 'decimate'
        'resample_factor', 1,           % 重采样因子
        'normalize_method', 'none',     % 归一化方法: 'none', 'minmax', 'zscore', 'robust'
        'derivative_order', 0,          % 导数阶数: 0(无), 1(一阶), 2(二阶)
        'show_plots', false             % 是否显示处理结果图
    );
    
    % 合并用户选项
    options = merge_options(default_options, processing_options);
    
    fprintf('=== 信号处理 ===\n');
    fprintf('输入数据点数: %d\n', length(data.wavenumber));
    
    % 初始化处理数据
    processed_data = data;
    processing_steps = {};
    
    %% 1. 异常值检测和处理
    if ~strcmp(options.outlier_method, 'none')
        fprintf('\n1. 异常值处理 (%s)...\n', options.outlier_method);
        [processed_data, outlier_info] = remove_outliers(processed_data, options);
        processing_steps{end+1} = sprintf('异常值处理: 移除%d个异常点', outlier_info.removed_count);
    end
    
    %% 2. 去趋势处理
    if ~strcmp(options.detrend_method, 'none')
        fprintf('\n2. 去趋势处理 (%s)...\n', options.detrend_method);
        [processed_data, trend_info] = detrend_signal(processed_data, options);
        processing_steps{end+1} = sprintf('去趋势: %s', options.detrend_method);
    end
    
    %% 3. 滤波处理
    if ~strcmp(options.filter_type, 'none')
        fprintf('\n3. 滤波处理 (%s)...\n', options.filter_type);
        [processed_data, filter_info] = filter_signal(processed_data, options);
        processing_steps{end+1} = sprintf('滤波: %s', options.filter_type);
    end
    
    %% 4. 去噪处理
    if ~strcmp(options.denoise_method, 'none')
        fprintf('\n4. 去噪处理 (%s)...\n', options.denoise_method);
        [processed_data, denoise_info] = denoise_signal(processed_data, options);
        processing_steps{end+1} = sprintf('去噪: %s', options.denoise_method);
    end
    
    %% 5. 平滑处理
    if ~strcmp(options.smooth_method, 'none')
        fprintf('\n5. 平滑处理 (%s)...\n', options.smooth_method);
        [processed_data, smooth_info] = smooth_signal(processed_data, options);
        processing_steps{end+1} = sprintf('平滑: %s', options.smooth_method);
    end
    
    %% 6. 重采样
    if ~strcmp(options.resample_method, 'none')
        fprintf('\n6. 重采样 (%s)...\n', options.resample_method);
        [processed_data, resample_info] = resample_signal(processed_data, options);
        processing_steps{end+1} = sprintf('重采样: %s', options.resample_method);
    end
    
    %% 7. 导数计算
    if options.derivative_order > 0
        fprintf('\n7. 导数计算 (%d阶)...\n', options.derivative_order);
        [processed_data, derivative_info] = compute_derivative(processed_data, options);
        processing_steps{end+1} = sprintf('%d阶导数', options.derivative_order);
    end
    
    %% 8. 归一化
    if ~strcmp(options.normalize_method, 'none')
        fprintf('\n8. 归一化 (%s)...\n', options.normalize_method);
        [processed_data, normalize_info] = normalize_signal(processed_data, options);
        processing_steps{end+1} = sprintf('归一化: %s', options.normalize_method);
    end
    
    %% 添加处理信息
    processed_data.processing_info = struct(
        'original_length', length(data.wavenumber),
        'processed_length', length(processed_data.wavenumber),
        'processing_steps', {processing_steps},
        'options_used', options,
        'processing_time', datetime('now')
    );
    
    fprintf('\n信号处理完成\n');
    fprintf('输出数据点数: %d\n', length(processed_data.wavenumber));
    
    %% 显示处理结果
    if options.show_plots
        plot_processing_results(data, processed_data);
    end
end

%% 异常值检测和移除
function [data, info] = remove_outliers(data, options)
    original_length = length(data.reflectance);
    
    switch options.outlier_method
        case 'zscore'
            % Z-score方法
            z_scores = abs(zscore(data.reflectance));
            outlier_idx = z_scores > options.outlier_threshold;
            
        case 'iqr'
            % 四分位距方法
            Q1 = quantile(data.reflectance, 0.25);
            Q3 = quantile(data.reflectance, 0.75);
            IQR = Q3 - Q1;
            lower_bound = Q1 - 1.5 * IQR;
            upper_bound = Q3 + 1.5 * IQR;
            outlier_idx = data.reflectance < lower_bound | data.reflectance > upper_bound;
            
        case 'hampel'
            % Hampel滤波器
            window_size = min(options.smooth_window, floor(length(data.reflectance)/4));
            outlier_idx = isoutlier(data.reflectance, 'hampel', 'SamplePoints', data.wavenumber, ...
                                  'ThresholdFactor', options.outlier_threshold);
            
        otherwise
            error('未知的异常值检测方法: %s', options.outlier_method);
    end
    
    % 移除异常值
    data.wavenumber = data.wavenumber(~outlier_idx);
    data.reflectance = data.reflectance(~outlier_idx);
    
    removed_count = sum(outlier_idx);
    fprintf('   移除异常值: %d个 (%.1f%%)\n', removed_count, removed_count/original_length*100);
    
    info = struct('removed_count', removed_count, 'outlier_indices', find(outlier_idx));
end

%% 去趋势处理
function [data, info] = detrend_signal(data, options)
    original_reflectance = data.reflectance;
    
    switch options.detrend_method
        case 'linear'
            % 线性去趋势
            p = polyfit(data.wavenumber, data.reflectance, 1);
            trend = polyval(p, data.wavenumber);
            data.reflectance = data.reflectance - trend;
            
        case 'polynomial'
            % 多项式去趋势
            p = polyfit(data.wavenumber, data.reflectance, options.detrend_order);
            trend = polyval(p, data.wavenumber);
            data.reflectance = data.reflectance - trend;
            
        case 'baseline'
            % 基线校正（使用最小值作为基线）
            baseline = min(data.reflectance);
            data.reflectance = data.reflectance - baseline;
            trend = baseline * ones(size(data.reflectance));
            
        otherwise
            error('未知的去趋势方法: %s', options.detrend_method);
    end
    
    trend_strength = std(trend) / std(original_reflectance);
    fprintf('   趋势强度: %.3f\n', trend_strength);
    
    info = struct('trend', trend, 'trend_strength', trend_strength);
end

%% 滤波处理
function [data, info] = filter_signal(data, options)
    % 设计滤波器
    fs = 1 / mean(diff(data.wavenumber));  % 采样频率估算
    nyquist = fs / 2;
    
    switch options.filter_type
        case 'lowpass'
            [b, a] = butter(options.filter_order, options.filter_cutoff(1));
            
        case 'highpass'
            [b, a] = butter(options.filter_order, options.filter_cutoff(1), 'high');
            
        case 'bandpass'
            [b, a] = butter(options.filter_order, options.filter_cutoff);
            
        otherwise
            error('未知的滤波类型: %s', options.filter_type);
    end
    
    % 应用滤波器
    try
        data.reflectance = filtfilt(b, a, data.reflectance);
        fprintf('   滤波成功\n');
    catch ME
        fprintf('   滤波失败，使用零相位滤波: %s\n', ME.message);
        data.reflectance = filter(b, a, data.reflectance);
    end
    
    info = struct('filter_coeffs', struct('b', b, 'a', a));
end

%% 去噪处理
function [data, info] = denoise_signal(data, options)
    original_reflectance = data.reflectance;
    
    switch options.denoise_method
        case 'wavelet'
            % 小波去噪
            try
                [data.reflectance, ~] = wdenoise(data.reflectance, options.wavelet_levels, ...
                                                'Wavelet', options.wavelet_type);
                fprintf('   小波去噪成功\n');
            catch ME
                fprintf('   小波去噪失败: %s\n', ME.message);
                % 使用简单的移动平均作为备选
                window_size = options.smooth_window;
                data.reflectance = movmean(data.reflectance, window_size);
            end
            
        case 'median'
            % 中值滤波
            window_size = options.smooth_window;
            data.reflectance = medfilt1(data.reflectance, window_size);
            fprintf('   中值滤波完成\n');
            
        case 'wiener'
            % 维纳滤波
            try
                data.reflectance = wiener2(data.reflectance, [options.smooth_window, 1]);
                fprintf('   维纳滤波完成\n');
            catch ME
                fprintf('   维纳滤波失败: %s\n', ME.message);
                % 使用高斯滤波作为备选
                sigma = options.smooth_window / 4;
                data.reflectance = imgaussfilt(data.reflectance, sigma);
            end
            
        otherwise
            error('未知的去噪方法: %s', options.denoise_method);
    end
    
    % 计算去噪效果
    noise_reduction = std(original_reflectance - data.reflectance) / std(original_reflectance);
    fprintf('   噪声降低: %.1f%%\n', noise_reduction * 100);
    
    info = struct('noise_reduction', noise_reduction);
end

%% 平滑处理
function [data, info] = smooth_signal(data, options)
    original_reflectance = data.reflectance;
    
    switch options.smooth_method
        case 'sgolay'
            % Savitzky-Golay平滑
            window_size = options.smooth_window;
            if mod(window_size, 2) == 0
                window_size = window_size + 1;  % 确保为奇数
            end
            
            if window_size > length(data.reflectance)
                window_size = length(data.reflectance);
                if mod(window_size, 2) == 0
                    window_size = window_size - 1;
                end
            end
            
            order = min(options.smooth_order, window_size - 1);
            data.reflectance = sgolayfilt(data.reflectance, order, window_size);
            
        case 'moving'
            % 移动平均
            data.reflectance = movmean(data.reflectance, options.smooth_window);
            
        case 'lowess'
            % LOWESS平滑
            span = options.smooth_window / length(data.reflectance);
            data.reflectance = smooth(data.reflectance, span, 'lowess');
            
        case 'gaussian'
            % 高斯平滑
            sigma = options.smooth_window / 4;
            data.reflectance = imgaussfilt(data.reflectance, sigma);
            
        otherwise
            error('未知的平滑方法: %s', options.smooth_method);
    end
    
    % 计算平滑效果
    smoothness = std(diff(data.reflectance)) / std(diff(original_reflectance));
    fprintf('   平滑度改善: %.1f%%\n', (1 - smoothness) * 100);
    
    info = struct('smoothness_improvement', 1 - smoothness);
end

%% 重采样
function [data, info] = resample_signal(data, options)
    original_length = length(data.wavenumber);
    
    switch options.resample_method
        case 'interp'
            % 插值重采样
            if options.resample_factor > 1
                % 上采样
                new_wavenumber = linspace(min(data.wavenumber), max(data.wavenumber), ...
                                        original_length * options.resample_factor);
                data.reflectance = interp1(data.wavenumber, data.reflectance, new_wavenumber, 'spline');
                data.wavenumber = new_wavenumber(:);
            else
                % 下采样
                downsample_factor = round(1 / options.resample_factor);
                indices = 1:downsample_factor:length(data.wavenumber);
                data.wavenumber = data.wavenumber(indices);
                data.reflectance = data.reflectance(indices);
            end
            
        case 'decimate'
            % 抽取重采样
            if options.resample_factor < 1
                decimate_factor = round(1 / options.resample_factor);
                data.reflectance = decimate(data.reflectance, decimate_factor);
                data.wavenumber = decimate(data.wavenumber, decimate_factor);
            end
            
        otherwise
            error('未知的重采样方法: %s', options.resample_method);
    end
    
    new_length = length(data.wavenumber);
    fprintf('   重采样: %d -> %d 点\n', original_length, new_length);
    
    info = struct('original_length', original_length, 'new_length', new_length);
end

%% 导数计算
function [data, info] = compute_derivative(data, options)
    switch options.derivative_order
        case 1
            % 一阶导数
            data.reflectance = gradient(data.reflectance, data.wavenumber);
            
        case 2
            % 二阶导数
            first_derivative = gradient(data.reflectance, data.wavenumber);
            data.reflectance = gradient(first_derivative, data.wavenumber);
            
        otherwise
            error('不支持的导数阶数: %d', options.derivative_order);
    end
    
    fprintf('   计算%d阶导数完成\n', options.derivative_order);
    
    info = struct('derivative_order', options.derivative_order);
end

%% 归一化
function [data, info] = normalize_signal(data, options)
    original_range = [min(data.reflectance), max(data.reflectance)];
    
    switch options.normalize_method
        case 'minmax'
            % 最小-最大归一化
            data.reflectance = (data.reflectance - min(data.reflectance)) / ...
                              (max(data.reflectance) - min(data.reflectance));
            
        case 'zscore'
            % Z-score标准化
            data.reflectance = zscore(data.reflectance);
            
        case 'robust'
            % 鲁棒归一化（使用中位数和MAD）
            median_val = median(data.reflectance);
            mad_val = mad(data.reflectance, 1);
            data.reflectance = (data.reflectance - median_val) / mad_val;
            
        otherwise
            error('未知的归一化方法: %s', options.normalize_method);
    end
    
    new_range = [min(data.reflectance), max(data.reflectance)];
    fprintf('   归一化: [%.3f, %.3f] -> [%.3f, %.3f]\n', ...
            original_range(1), original_range(2), new_range(1), new_range(2));
    
    info = struct('original_range', original_range, 'new_range', new_range);
end

%% 合并选项
function merged_options = merge_options(default_options, user_options)
    merged_options = default_options;
    
    if ~isempty(user_options)
        field_names = fieldnames(user_options);
        for i = 1:length(field_names)
            merged_options.(field_names{i}) = user_options.(field_names{i});
        end
    end
end

%% 绘制处理结果
function plot_processing_results(original_data, processed_data)
    figure('Name', '信号处理结果', 'Position', [100, 100, 1400, 800]);
    
    % 原始信号
    subplot(2, 3, 1);
    plot(original_data.wavenumber, original_data.reflectance, 'b-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率');
    title('原始信号');
    grid on;
    
    % 处理后信号
    subplot(2, 3, 2);
    plot(processed_data.wavenumber, processed_data.reflectance, 'r-', 'LineWidth', 1.5);
    xlabel('波数 (cm^{-1})');
    ylabel('反射率');
    title('处理后信号');
    grid on;
    
    % 对比图
    subplot(2, 3, 3);
    if length(original_data.wavenumber) == length(processed_data.wavenumber)
        plot(original_data.wavenumber, original_data.reflectance, 'b-', 'LineWidth', 1, 'DisplayName', '原始');
        hold on;
        plot(processed_data.wavenumber, processed_data.reflectance, 'r-', 'LineWidth', 1.5, 'DisplayName', '处理后');
        legend('show');
    else
        % 长度不同时分别绘制
        yyaxis left;
        plot(original_data.wavenumber, original_data.reflectance, 'b-', 'LineWidth', 1);
        ylabel('原始反射率', 'Color', 'b');
        
        yyaxis right;
        plot(processed_data.wavenumber, processed_data.reflectance, 'r-', 'LineWidth', 1.5);
        ylabel('处理后反射率', 'Color', 'r');
    end
    xlabel('波数 (cm^{-1})');
    title('信号对比');
    grid on;
    
    % 频谱分析
    subplot(2, 3, 4);
    if length(original_data.reflectance) > 10
        [psd_orig, f_orig] = periodogram(original_data.reflectance);
        semilogy(f_orig, psd_orig, 'b-', 'LineWidth', 1.5);
        xlabel('归一化频率');
        ylabel('功率谱密度');
        title('原始信号频谱');
        grid on;
    end
    
    subplot(2, 3, 5);
    if length(processed_data.reflectance) > 10
        [psd_proc, f_proc] = periodogram(processed_data.reflectance);
        semilogy(f_proc, psd_proc, 'r-', 'LineWidth', 1.5);
        xlabel('归一化频率');
        ylabel('功率谱密度');
        title('处理后信号频谱');
        grid on;
    end
    
    % 处理步骤信息
    subplot(2, 3, 6);
    if isfield(processed_data, 'processing_info') && isfield(processed_data.processing_info, 'processing_steps')
        steps = processed_data.processing_info.processing_steps;
        if ~isempty(steps)
            step_text = strjoin(steps, '\n');
            text(0.05, 0.95, step_text, 'FontSize', 10, 'Units', 'normalized', ...
                 'VerticalAlignment', 'top', 'Interpreter', 'none');
        end
    end
    axis off;
    title('处理步骤', 'FontSize', 12, 'FontWeight', 'bold');
    
    sgtitle('信号处理结果对比', 'FontSize', 16, 'FontWeight', 'bold');
end