%% 数学工具函数集合
% 提供各种数学计算和数值分析功能

%% 数值积分函数
function result = numerical_integration(func, a, b, method, options)
    % 数值积分
    % 输入参数：
    %   func - 函数句柄或数据点
    %   a, b - 积分区间
    %   method - 积分方法: 'trapz', 'simpson', 'quad', 'quadgk'
    %   options - 选项结构体
    
    if nargin < 4
        method = 'trapz';
    end
    
    if nargin < 5
        options = struct();
    end
    
    % 默认选项
    default_options = struct(
        'n_points', 1000,       % 积分点数
        'tolerance', 1e-6,      % 容差
        'max_iter', 1000        % 最大迭代次数
    );
    
    options = merge_options(default_options, options);
    
    switch lower(method)
        case 'trapz'
            if isa(func, 'function_handle')
                x = linspace(a, b, options.n_points);
                y = func(x);
                result = trapz(x, y);
            else
                % func是数据点
                result = trapz(func);
            end
            
        case 'simpson'
            if isa(func, 'function_handle')
                result = simpson_rule(func, a, b, options.n_points);
            else
                result = simpson_rule_data(func);
            end
            
        case 'quad'
            if isa(func, 'function_handle')
                result = quad(func, a, b, options.tolerance);
            else
                error('quad方法需要函数句柄');
            end
            
        case 'quadgk'
            if isa(func, 'function_handle')
                result = quadgk(func, a, b, 'RelTol', options.tolerance);
            else
                error('quadgk方法需要函数句柄');
            end
            
        otherwise
            error('未知的积分方法: %s', method);
    end
end

%% 辛普森积分规则
function result = simpson_rule(func, a, b, n)
    % 辛普森1/3规则
    if mod(n, 2) == 1
        n = n + 1;  % 确保为偶数
    end
    
    h = (b - a) / n;
    x = a:h:b;
    y = func(x);
    
    result = h/3 * (y(1) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end));
end

%% 数据点辛普森积分
function result = simpson_rule_data(y)
    n = length(y);
    if mod(n, 2) == 0
        n = n - 1;
        y = y(1:n);
    end
    
    result = (y(1) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end)) / 3;
end

%% 数值微分函数
function dy = numerical_derivative(y, x, method, order)
    % 数值微分
    % 输入参数：
    %   y - 函数值数组
    %   x - 自变量数组（可选，默认为等间距）
    %   method - 微分方法: 'forward', 'backward', 'central', 'gradient'
    %   order - 微分阶数（1或2）
    
    if nargin < 2 || isempty(x)
        x = 1:length(y);
    end
    
    if nargin < 3
        method = 'central';
    end
    
    if nargin < 4
        order = 1;
    end
    
    if length(y) ~= length(x)
        error('y和x的长度必须相同');
    end
    
    switch order
        case 1
            dy = first_derivative(y, x, method);
        case 2
            dy = second_derivative(y, x, method);
        otherwise
            error('只支持1阶和2阶微分');
    end
end

%% 一阶微分
function dy = first_derivative(y, x, method)
    n = length(y);
    dy = zeros(size(y));
    
    switch lower(method)
        case 'forward'
            % 前向差分
            for i = 1:n-1
                dy(i) = (y(i+1) - y(i)) / (x(i+1) - x(i));
            end
            dy(n) = dy(n-1);  % 边界处理
            
        case 'backward'
            % 后向差分
            dy(1) = dy(2);  % 边界处理
            for i = 2:n
                dy(i) = (y(i) - y(i-1)) / (x(i) - x(i-1));
            end
            
        case 'central'
            % 中心差分
            dy(1) = (y(2) - y(1)) / (x(2) - x(1));  % 前向
            for i = 2:n-1
                dy(i) = (y(i+1) - y(i-1)) / (x(i+1) - x(i-1));
            end
            dy(n) = (y(n) - y(n-1)) / (x(n) - x(n-1));  % 后向
            
        case 'gradient'
            % MATLAB内置gradient函数
            dy = gradient(y, x);
            
        otherwise
            error('未知的微分方法: %s', method);
    end
end

%% 二阶微分
function d2y = second_derivative(y, x, method)
    % 计算一阶微分
    dy = first_derivative(y, x, method);
    % 再次微分得到二阶微分
    d2y = first_derivative(dy, x, method);
end

%% 插值函数
function yi = interpolation(x, y, xi, method, extrapolation)
    % 插值函数
    % 输入参数：
    %   x, y - 原始数据点
    %   xi - 插值点
    %   method - 插值方法
    %   extrapolation - 外推方法
    
    if nargin < 4
        method = 'linear';
    end
    
    if nargin < 5
        extrapolation = 'extrap';
    end
    
    % 检查数据
    if length(x) ~= length(y)
        error('x和y的长度必须相同');
    end
    
    % 移除重复点
    [x, unique_idx] = unique(x);
    y = y(unique_idx);
    
    % 执行插值
    try
        yi = interp1(x, y, xi, method, extrapolation);
    catch ME
        warning('插值失败: %s，使用线性插值', ME.message);
        yi = interp1(x, y, xi, 'linear', extrapolation);
    end
end

%% 曲线拟合函数
function [fit_params, fit_func, gof] = curve_fitting(x, y, fit_type, options)
    % 曲线拟合
    % 输入参数：
    %   x, y - 数据点
    %   fit_type - 拟合类型
    %   options - 拟合选项
    % 输出参数：
    %   fit_params - 拟合参数
    %   fit_func - 拟合函数句柄
    %   gof - 拟合优度
    
    if nargin < 3
        fit_type = 'polynomial';
    end
    
    if nargin < 4
        options = struct();
    end
    
    % 默认选项
    default_options = struct(
        'poly_order', 2,        % 多项式阶数
        'robust', false,        % 鲁棒拟合
        'weights', [],          % 权重
        'start_point', [],      % 初始参数
        'lower_bounds', [],     % 参数下界
        'upper_bounds', []      % 参数上界
    );
    
    options = merge_options(default_options, options);
    
    switch lower(fit_type)
        case 'polynomial'
            [fit_params, gof] = polynomial_fit(x, y, options);
            fit_func = @(x_new) polyval(fit_params, x_new);
            
        case 'exponential'
            [fit_params, gof] = exponential_fit(x, y, options);
            fit_func = @(x_new) fit_params(1) * exp(fit_params(2) * x_new) + fit_params(3);
            
        case 'gaussian'
            [fit_params, gof] = gaussian_fit(x, y, options);
            fit_func = @(x_new) fit_params(1) * exp(-((x_new - fit_params(2)) / fit_params(3)).^2);
            
        case 'sinusoidal'
            [fit_params, gof] = sinusoidal_fit(x, y, options);
            fit_func = @(x_new) fit_params(1) * sin(fit_params(2) * x_new + fit_params(3)) + fit_params(4);
            
        case 'custom'
            if ~isfield(options, 'custom_func')
                error('自定义拟合需要提供custom_func');
            end
            [fit_params, gof] = custom_fit(x, y, options);
            fit_func = @(x_new) options.custom_func(fit_params, x_new);
            
        otherwise
            error('未知的拟合类型: %s', fit_type);
    end
end

%% 多项式拟合
function [params, gof] = polynomial_fit(x, y, options)
    if options.robust
        [params, S] = polyfit(x, y, options.poly_order);
        % 计算拟合优度
        y_fit = polyval(params, x);
        gof = calculate_goodness_of_fit(y, y_fit);
    else
        [params, S] = polyfit(x, y, options.poly_order);
        y_fit = polyval(params, x);
        gof = calculate_goodness_of_fit(y, y_fit);
        gof.S = S;
    end
end

%% 指数拟合
function [params, gof] = exponential_fit(x, y, options)
    % 指数函数: y = a * exp(b * x) + c
    
    % 初始参数估计
    if isempty(options.start_point)
        a0 = max(y) - min(y);
        b0 = (log(max(y)) - log(min(y))) / (max(x) - min(x));
        c0 = min(y);
        start_point = [a0, b0, c0];
    else
        start_point = options.start_point;
    end
    
    % 拟合函数
    exp_func = @(params, x) params(1) * exp(params(2) * x) + params(3);
    
    % 非线性拟合
    try
        params = lsqcurvefit(exp_func, start_point, x, y, ...
                           options.lower_bounds, options.upper_bounds);
    catch
        % 使用fminsearch作为备选
        objective = @(params) sum((y - exp_func(params, x)).^2);
        params = fminsearch(objective, start_point);
    end
    
    % 计算拟合优度
    y_fit = exp_func(params, x);
    gof = calculate_goodness_of_fit(y, y_fit);
end

%% 高斯拟合
function [params, gof] = gaussian_fit(x, y, options)
    % 高斯函数: y = a * exp(-((x - b) / c)^2)
    
    % 初始参数估计
    if isempty(options.start_point)
        [max_val, max_idx] = max(y);
        a0 = max_val;
        b0 = x(max_idx);
        c0 = (max(x) - min(x)) / 4;
        start_point = [a0, b0, c0];
    else
        start_point = options.start_point;
    end
    
    % 拟合函数
    gauss_func = @(params, x) params(1) * exp(-((x - params(2)) / params(3)).^2);
    
    % 非线性拟合
    try
        params = lsqcurvefit(gauss_func, start_point, x, y, ...
                           options.lower_bounds, options.upper_bounds);
    catch
        objective = @(params) sum((y - gauss_func(params, x)).^2);
        params = fminsearch(objective, start_point);
    end
    
    % 计算拟合优度
    y_fit = gauss_func(params, x);
    gof = calculate_goodness_of_fit(y, y_fit);
end

%% 正弦拟合
function [params, gof] = sinusoidal_fit(x, y, options)
    % 正弦函数: y = a * sin(b * x + c) + d
    
    % 初始参数估计
    if isempty(options.start_point)
        a0 = (max(y) - min(y)) / 2;
        d0 = mean(y);
        
        % 估计频率
        [psd, f] = periodogram(y - d0);
        [~, max_idx] = max(psd);
        b0 = 2 * pi * f(max_idx) * length(x) / (max(x) - min(x));
        
        c0 = 0;
        start_point = [a0, b0, c0, d0];
    else
        start_point = options.start_point;
    end
    
    % 拟合函数
    sin_func = @(params, x) params(1) * sin(params(2) * x + params(3)) + params(4);
    
    % 非线性拟合
    try
        params = lsqcurvefit(sin_func, start_point, x, y, ...
                           options.lower_bounds, options.upper_bounds);
    catch
        objective = @(params) sum((y - sin_func(params, x)).^2);
        params = fminsearch(objective, start_point);
    end
    
    % 计算拟合优度
    y_fit = sin_func(params, x);
    gof = calculate_goodness_of_fit(y, y_fit);
end

%% 自定义拟合
function [params, gof] = custom_fit(x, y, options)
    if isempty(options.start_point)
        error('自定义拟合需要提供初始参数');
    end
    
    % 非线性拟合
    try
        params = lsqcurvefit(options.custom_func, options.start_point, x, y, ...
                           options.lower_bounds, options.upper_bounds);
    catch
        objective = @(params) sum((y - options.custom_func(params, x)).^2);
        params = fminsearch(objective, options.start_point);
    end
    
    % 计算拟合优度
    y_fit = options.custom_func(params, x);
    gof = calculate_goodness_of_fit(y, y_fit);
end

%% 计算拟合优度
function gof = calculate_goodness_of_fit(y_actual, y_fit)
    % 计算各种拟合优度指标
    
    n = length(y_actual);
    
    % 残差
    residuals = y_actual - y_fit;
    
    % 均方误差
    mse = mean(residuals.^2);
    
    % 均方根误差
    rmse = sqrt(mse);
    
    % 平均绝对误差
    mae = mean(abs(residuals));
    
    % 决定系数 R²
    ss_res = sum(residuals.^2);
    ss_tot = sum((y_actual - mean(y_actual)).^2);
    r_squared = 1 - ss_res / ss_tot;
    
    % 调整决定系数
    if n > 2
        adj_r_squared = 1 - (1 - r_squared) * (n - 1) / (n - 2);
    else
        adj_r_squared = r_squared;
    end
    
    % 相关系数
    correlation = corr(y_actual, y_fit);
    
    % 构建结果结构体
    gof = struct(
        'mse', mse,
        'rmse', rmse,
        'mae', mae,
        'r_squared', r_squared,
        'adj_r_squared', adj_r_squared,
        'correlation', correlation,
        'residuals', residuals
    );
end

%% 峰值检测函数
function [peaks, peak_locs, valleys, valley_locs] = peak_detection(y, options)
    % 峰值和谷值检测
    
    if nargin < 2
        options = struct();
    end
    
    % 默认选项
    default_options = struct(
        'min_peak_height', -Inf,
        'min_peak_distance', 1,
        'min_peak_prominence', 0,
        'threshold', 0
    );
    
    options = merge_options(default_options, options);
    
    % 检测峰值
    [peaks, peak_locs] = findpeaks(y, ...
        'MinPeakHeight', options.min_peak_height, ...
        'MinPeakDistance', options.min_peak_distance, ...
        'MinPeakProminence', options.min_peak_prominence, ...
        'Threshold', options.threshold);
    
    % 检测谷值（通过检测负信号的峰值）
    [valleys_neg, valley_locs] = findpeaks(-y, ...
        'MinPeakHeight', -max(y), ...
        'MinPeakDistance', options.min_peak_distance, ...
        'MinPeakProminence', options.min_peak_prominence, ...
        'Threshold', options.threshold);
    
    valleys = -valleys_neg;
end

%% 统计分析函数
function stats = statistical_analysis(data, confidence_level)
    % 统计分析
    
    if nargin < 2
        confidence_level = 0.95;
    end
    
    n = length(data);
    
    % 基本统计量
    mean_val = mean(data);
    median_val = median(data);
    std_val = std(data);
    var_val = var(data);
    min_val = min(data);
    max_val = max(data);
    range_val = max_val - min_val;
    
    % 分位数
    q25 = quantile(data, 0.25);
    q75 = quantile(data, 0.75);
    iqr_val = q75 - q25;
    
    % 偏度和峰度
    skewness_val = skewness(data);
    kurtosis_val = kurtosis(data);
    
    % 置信区间
    alpha = 1 - confidence_level;
    t_critical = tinv(1 - alpha/2, n - 1);
    margin_error = t_critical * std_val / sqrt(n);
    ci_lower = mean_val - margin_error;
    ci_upper = mean_val + margin_error;
    
    % 异常值检测（IQR方法）
    outlier_threshold_lower = q25 - 1.5 * iqr_val;
    outlier_threshold_upper = q75 + 1.5 * iqr_val;
    outliers = data < outlier_threshold_lower | data > outlier_threshold_upper;
    num_outliers = sum(outliers);
    
    % 构建结果结构体
    stats = struct(
        'n', n,
        'mean', mean_val,
        'median', median_val,
        'std', std_val,
        'var', var_val,
        'min', min_val,
        'max', max_val,
        'range', range_val,
        'q25', q25,
        'q75', q75,
        'iqr', iqr_val,
        'skewness', skewness_val,
        'kurtosis', kurtosis_val,
        'confidence_interval', [ci_lower, ci_upper],
        'confidence_level', confidence_level,
        'outliers', outliers,
        'num_outliers', num_outliers
    );
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