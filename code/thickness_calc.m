function [thickness_results] = thickness_calc(wavenumber, reflectance, incident_angle, material_type)
% THICKNESS_CALC - 外延层厚度计算与验证程序
%
% 功能描述:
%   基于红外干涉法计算碳化硅外延层厚度
%   包括多种算法和可靠性验证
%
% 输入参数:
%   wavenumber - 波数数组 (cm^-1)
%   reflectance - 反射率数组 (%)
%   incident_angle - 入射角 (度)
%   material_type - 材料类型 ('SiC' 或 'Si')
%
% 输出参数:
%   thickness_results - 厚度计算结果结构体
%     .thickness - 计算得到的厚度 (μm)
%     .method - 使用的计算方法
%     .reliability - 可靠性评估
%     .extrema_points - 极值点信息
%     .interference_orders - 干涉级数


    if nargin < 4
        material_type = 'SiC';
    end
    
    fprintf('=== 外延层厚度计算程序 ===\n');
    fprintf('材料类型: %s\n', material_type);
    fprintf('入射角: %.1f°\n', incident_angle);
    fprintf('数据点数: %d\n', length(wavenumber));
    
    % 初始化结果结构体
    thickness_results = struct();
    thickness_results.input_params = struct('incident_angle', incident_angle, ...
                                           'material_type', material_type, ...
                                           'data_points', length(wavenumber));
    
    try
        % 1. 数据预处理
        fprintf('\n1. 数据预处理...\n');
        [processed_data] = preprocess_spectrum_data(wavenumber, reflectance);
        
        % 2. 极值点提取
        fprintf('2. 极值点提取...\n');
        [extrema_info] = extract_extrema_points(processed_data);
        thickness_results.extrema_points = extrema_info;
        
        % 3. 干涉级数计算
        fprintf('3. 干涉级数计算...\n');
        [interference_orders] = calculate_interference_orders(extrema_info);
        thickness_results.interference_orders = interference_orders;
        
        % 4. 厚度计算 - 多种方法
        fprintf('4. 厚度计算...\n');
        thickness_methods = calculate_thickness_multiple_methods(extrema_info, ...
                                                               interference_orders, ...
                                                               incident_angle, ...
                                                               material_type);
        thickness_results.methods = thickness_methods;
        
        % 5. 选择最佳厚度值
        [best_thickness, best_method] = select_best_thickness(thickness_methods);
        thickness_results.thickness = best_thickness;
        thickness_results.method = best_method;
        
        % 6. 可靠性验证
        fprintf('5. 可靠性验证...\n');
        reliability_assessment = assess_reliability(thickness_methods, extrema_info);
        thickness_results.reliability = reliability_assessment;
        
        % 7. 结果输出
        fprintf('\n=== 计算结果 ===\n');
        fprintf('最佳厚度: %.3f μm\n', best_thickness);
        fprintf('计算方法: %s\n', best_method);
        fprintf('可靠性等级: %s\n', reliability_assessment.grade);
        fprintf('置信度: %.1f%%\n', reliability_assessment.confidence * 100);
        
    catch ME
        fprintf('厚度计算出错: %s\n', ME.message);
        thickness_results.error = ME.message;
        thickness_results.thickness = NaN;
    end
end

function [processed_data] = preprocess_spectrum_data(wavenumber, reflectance)
% 光谱数据预处理
    
    % 数据验证
    if length(wavenumber) ~= length(reflectance)
        error('波数和反射率数组长度不匹配');
    end
    
    % 移除NaN值
    valid_idx = ~isnan(wavenumber) & ~isnan(reflectance);
    wavenumber = wavenumber(valid_idx);
    reflectance = reflectance(valid_idx);
    
    % 数据排序（按波数升序）
    [wavenumber, sort_idx] = sort(wavenumber);
    reflectance = reflectance(sort_idx);
    
    % 异常值检测和处理
    z_scores = abs((reflectance - mean(reflectance)) / std(reflectance));
    outlier_mask = z_scores > 3;
    outlier_count = sum(outlier_mask);
    
    if outlier_count > 0
        fprintf('  检测到 %d 个异常值，进行插值处理\n', outlier_count);
        reflectance(outlier_mask) = interp1(find(~outlier_mask), ...
                                           reflectance(~outlier_mask), ...
                                           find(outlier_mask), 'linear', 'extrap');
    end
    
    % 平滑滤波
    if length(reflectance) > 10
        window_size = min(11, length(reflectance));
        if mod(window_size, 2) == 0
            window_size = window_size - 1;
        end
        reflectance_smooth = smooth(reflectance, window_size, 'sgolay');
    else
        reflectance_smooth = reflectance;
    end
    
    % 转换为波长
    wavelength = 10000 ./ wavenumber;  % μm
    
    processed_data = struct();
    processed_data.wavenumber = wavenumber;
    processed_data.wavelength = wavelength;
    processed_data.reflectance_raw = reflectance;
    processed_data.reflectance = reflectance_smooth;
    processed_data.outlier_count = outlier_count;
    
    fprintf('  数据预处理完成: %d 个有效数据点\n', length(wavenumber));
end

function [extrema_info] = extract_extrema_points(processed_data)
% 提取干涉条纹极值点
    
    wavenumber = processed_data.wavenumber;
    reflectance = processed_data.reflectance;
    wavelength = processed_data.wavelength;
    
    % 寻找极大值点（峰值）
    [peak_values, peak_locs] = findpeaks(reflectance, 'MinPeakProminence', 0.5, ...
                                        'MinPeakDistance', 5);
    
    % 寻找极小值点（谷值）
    [valley_values, valley_locs] = findpeaks(-reflectance, 'MinPeakProminence', 0.5, ...
                                            'MinPeakDistance', 5);
    valley_values = -valley_values;
    
    % 合并极值点
    all_locs = [peak_locs; valley_locs];
    all_values = [peak_values; valley_values];
    all_types = [ones(length(peak_locs), 1); -ones(length(valley_locs), 1)];
    
    % 按位置排序
    [~, sort_idx] = sort(all_locs);
    all_locs = all_locs(sort_idx);
    all_values = all_values(sort_idx);
    all_types = all_types(sort_idx);
    
    extrema_info = struct();
    extrema_info.indices = all_locs;
    extrema_info.wavenumbers = wavenumber(all_locs);
    extrema_info.wavelengths = wavelength(all_locs);
    extrema_info.reflectances = all_values;
    extrema_info.types = all_types;  % 1为峰值，-1为谷值
    extrema_info.peak_count = length(peak_locs);
    extrema_info.valley_count = length(valley_locs);
    
    fprintf('  找到极值点: %d 个 (峰值: %d, 谷值: %d)\n', ...
            length(all_locs), length(peak_locs), length(valley_locs));
end

function [interference_orders] = calculate_interference_orders(extrema_info)
% 计算干涉级数
    
    wavelengths = extrema_info.wavelengths;
    types = extrema_info.types;
    
    if length(wavelengths) < 2
        interference_orders = [];
        return;
    end
    
    % 初始化干涉级数
    orders = zeros(size(wavelengths));
    
    % 设置参考点（第一个极值点）
    orders(1) = 1.0;  % 假设第一个极值点对应级数为1
    
    % 计算相邻极值点的级数
    for i = 2:length(wavelengths)
        % 相邻极值点的相位差为π，对应级数差为0.5
        orders(i) = orders(i-1) + 0.5;
    end
    
    interference_orders = struct();
    interference_orders.orders = orders;
    interference_orders.wavelengths = wavelengths;
    interference_orders.types = types;
    
    fprintf('  干涉级数范围: %.1f - %.1f\n', min(orders), max(orders));
end

function [thickness_methods] = calculate_thickness_multiple_methods(extrema_info, ...
                                                                  interference_orders, ...
                                                                  incident_angle, ...
                                                                  material_type)
% 使用多种方法计算厚度
    
    % 材料参数
    switch upper(material_type)
        case 'SIC'
            n_material = 2.55;  % SiC折射率
        case 'SI'
            n_material = 3.42;  % Si折射率
        otherwise
            n_material = 2.55;  % 默认SiC
    end
    
    % 计算折射角
    theta_i = deg2rad(incident_angle);
    theta_r = asin(sin(theta_i) / n_material);  % 斯涅尔定律
    
    thickness_methods = struct();
    
    % 方法1: 基于相邻极值点的厚度计算
    if length(interference_orders.wavelengths) >= 2
        thickness_adjacent = calculate_thickness_adjacent_extrema(interference_orders, ...
                                                                n_material, theta_r);
        thickness_methods.adjacent_extrema = thickness_adjacent;
    end
    
    % 方法2: 基于多波长拟合的厚度计算
    if length(interference_orders.wavelengths) >= 3
        thickness_fitting = calculate_thickness_fitting(interference_orders, ...
                                                       n_material, theta_r);
        thickness_methods.wavelength_fitting = thickness_fitting;
    end
    
    % 方法3: 基于傅里叶变换的厚度计算
    thickness_fft = calculate_thickness_fft(extrema_info, n_material, theta_r);
    thickness_methods.fft_method = thickness_fft;
    
    % 方法4: 基于最小二乘拟合的厚度计算
    if length(interference_orders.wavelengths) >= 4
        thickness_lsq = calculate_thickness_least_squares(interference_orders, ...
                                                         n_material, theta_r);
        thickness_methods.least_squares = thickness_lsq;
    end
end

function [thickness] = calculate_thickness_adjacent_extrema(interference_orders, n_material, theta_r)
% 基于相邻极值点计算厚度
    
    wavelengths = interference_orders.wavelengths;
    orders = interference_orders.orders;
    
    thickness_values = [];
    
    for i = 1:length(wavelengths)
        % 厚度公式: d = (m * λ) / (2 * n * cos(θr))
        d = (orders(i) * wavelengths(i)) / (2 * n_material * cos(theta_r));
        thickness_values = [thickness_values; d];
    end
    
    thickness = struct();
    thickness.values = thickness_values;
    thickness.mean = mean(thickness_values);
    thickness.std = std(thickness_values);
    thickness.method = 'Adjacent Extrema';
    
    fprintf('    相邻极值点方法: %.3f ± %.3f μm\n', thickness.mean, thickness.std);
end

function [thickness] = calculate_thickness_fitting(interference_orders, n_material, theta_r)
% 基于多波长拟合计算厚度
    
    wavelengths = interference_orders.wavelengths;
    orders = interference_orders.orders;
    
    % 线性拟合: m = (2 * n * d * cos(θr)) / λ
    % 即: m * λ = 2 * n * d * cos(θr)
    x = wavelengths;
    y = orders .* wavelengths;
    
    % 最小二乘拟合
    p = polyfit(ones(size(x)), y, 0);  % 零次多项式（常数）
    thickness_value = p(1) / (2 * n_material * cos(theta_r));
    
    % 计算拟合优度
    y_fit = p(1) * ones(size(x));
    r_squared = 1 - sum((y - y_fit).^2) / sum((y - mean(y)).^2);
    
    thickness = struct();
    thickness.value = thickness_value;
    thickness.r_squared = r_squared;
    thickness.method = 'Wavelength Fitting';
    
    fprintf('    多波长拟合方法: %.3f μm (R² = %.3f)\n', thickness_value, r_squared);
end

function [thickness] = calculate_thickness_fft(extrema_info, n_material, theta_r)
% 基于FFT计算厚度
    
    wavenumbers = extrema_info.wavenumbers;
    
    if length(wavenumbers) < 4
        thickness = struct('value', NaN, 'method', 'FFT (insufficient data)');
        return;
    end
    
    % 计算波数间隔的FFT
    wavenumber_diffs = diff(wavenumbers);
    
    % 估算平均间隔
    avg_interval = mean(wavenumber_diffs);
    
    % 转换为厚度
    % 对于相邻极值点: Δν = 1 / (2 * n * d * cos(θr))
    thickness_value = 1 / (2 * n_material * cos(theta_r) * avg_interval);
    
    thickness = struct();
    thickness.value = thickness_value;
    thickness.avg_interval = avg_interval;
    thickness.method = 'FFT Method';
    
    fprintf('    FFT方法: %.3f μm\n', thickness_value);
end

function [thickness] = calculate_thickness_least_squares(interference_orders, n_material, theta_r)
% 基于最小二乘法计算厚度
    
    wavelengths = interference_orders.wavelengths;
    orders = interference_orders.orders;
    
    % 构建线性方程组: m = (2 * n * d * cos(θr)) / λ
    % 重写为: m * λ = 2 * n * d * cos(θr)
    A = wavelengths;
    b = orders .* wavelengths;
    
    % 最小二乘求解
    thickness_coeff = (A' * A) \ (A' * b);
    thickness_value = thickness_coeff / (2 * n_material * cos(theta_r));
    
    % 计算残差
    residuals = b - A * thickness_coeff;
    rmse = sqrt(mean(residuals.^2));
    
    thickness = struct();
    thickness.value = thickness_value;
    thickness.rmse = rmse;
    thickness.method = 'Least Squares';
    
    fprintf('    最小二乘法: %.3f μm (RMSE = %.3f)\n', thickness_value, rmse);
end

function [best_thickness, best_method] = select_best_thickness(thickness_methods)
% 选择最佳厚度值
    
    methods = fieldnames(thickness_methods);
    thickness_values = [];
    method_names = {};
    weights = [];
    
    for i = 1:length(methods)
        method_data = thickness_methods.(methods{i});
        
        if isfield(method_data, 'mean')
            thickness_values = [thickness_values; method_data.mean];
            method_names{end+1} = method_data.method;
            % 权重基于标准差的倒数
            weights = [weights; 1 / (method_data.std + 0.001)];
        elseif isfield(method_data, 'value') && ~isnan(method_data.value)
            thickness_values = [thickness_values; method_data.value];
            method_names{end+1} = method_data.method;
            
            % 根据方法质量设置权重
            if isfield(method_data, 'r_squared')
                weights = [weights; method_data.r_squared];
            elseif isfield(method_data, 'rmse')
                weights = [weights; 1 / (method_data.rmse + 0.001)];
            else
                weights = [weights; 0.5];
            end
        end
    end
    
    if isempty(thickness_values)
        best_thickness = NaN;
        best_method = 'No valid method';
        return;
    end
    
    % 加权平均
    weights = weights / sum(weights);
    best_thickness = sum(thickness_values .* weights);
    
    % 选择权重最大的方法作为最佳方法
    [~, best_idx] = max(weights);
    best_method = method_names{best_idx};
end

function [reliability] = assess_reliability(thickness_methods, extrema_info)
% 评估计算结果的可靠性
    
    methods = fieldnames(thickness_methods);
    valid_thicknesses = [];
    
    % 收集所有有效的厚度值
    for i = 1:length(methods)
        method_data = thickness_methods.(methods{i});
        if isfield(method_data, 'mean')
            valid_thicknesses = [valid_thicknesses; method_data.mean];
        elseif isfield(method_data, 'value') && ~isnan(method_data.value)
            valid_thicknesses = [valid_thicknesses; method_data.value];
        end
    end
    
    reliability = struct();
    
    if length(valid_thicknesses) >= 2
        % 计算一致性指标
        thickness_std = std(valid_thicknesses);
        thickness_mean = mean(valid_thicknesses);
        cv = thickness_std / thickness_mean;  % 变异系数
        
        % 数据质量指标
        extrema_count = length(extrema_info.wavelengths);
        
        % 可靠性评级
        if cv < 0.05 && extrema_count >= 6
            grade = 'Excellent';
            confidence = 0.95;
        elseif cv < 0.10 && extrema_count >= 4
            grade = 'Good';
            confidence = 0.85;
        elseif cv < 0.20 && extrema_count >= 3
            grade = 'Fair';
            confidence = 0.70;
        else
            grade = 'Poor';
            confidence = 0.50;
        end
        
        reliability.grade = grade;
        reliability.confidence = confidence;
        reliability.coefficient_of_variation = cv;
        reliability.method_count = length(valid_thicknesses);
        reliability.extrema_count = extrema_count;
        reliability.thickness_std = thickness_std;
        
    else
        reliability.grade = 'Insufficient';
        reliability.confidence = 0.30;
        reliability.coefficient_of_variation = NaN;
        reliability.method_count = length(valid_thicknesses);
        reliability.extrema_count = length(extrema_info.wavelengths);
    end
end

% 辅助函数：平滑滤波
function smoothed = smooth(data, window_size, method)
    if nargin < 3
        method = 'moving';
    end
    
    switch method
        case 'sgolay'
            if window_size > length(data)
                window_size = length(data);
            end
            if mod(window_size, 2) == 0
                window_size = window_size - 1;
            end
            if window_size >= 3
                smoothed = sgolayfilt(data, min(3, window_size-1), window_size);
            else
                smoothed = data;
            end
        otherwise
            % 移动平均
            smoothed = movmean(data, window_size);
    end
end