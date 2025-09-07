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
    % 计算折射角
    theta_i = deg2rad(incident_angle);
    theta_r = asin(sin(theta_i) / n_material);  % 斯涅尔定律
    
    thickness_methods = struct();
 % 计算折射角 基于最小二乘拟合的厚度计算
    if length(interference_orders.wavelengths) >= 4
        thickness_lsq = calculate_thickness_least_squares(interference_orders, ...
                                                         n_material, theta_r);
        thickness_methods.least_squares = thickness_lsq;
    end
end

function [thickness] = calculate_thickness_adjacent_extrema(interference_orders, n_material, theta_r)

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