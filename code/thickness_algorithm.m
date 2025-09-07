function [thickness, results] = thickness_algorithm(wavenumber, reflectance, angle, n_epi)
% THICKNESS_ALGORITHM - 厚度确定算法
%
% 基于红外干涉光谱数据计算SiC外延层厚度
% 使用多种算法进行厚度计算并进行结果验证
%
% 输入参数:
%   wavenumber - 波数数组 (cm^-1)
%   reflectance - 反射率数组 (%)
%   angle - 入射角 (度)
%   n_epi - 外延层折射率 (默认为SiC折射率2.55)
%
% 输出参数:
%   thickness - 计算得到的厚度 (μm)
%   results - 详细计算结果结构体
%
% 算法原理:
%   1. 极值法：基于干涉极值位置计算厚度
%   2. 相位拟合法：拟合相位-波数关系
%   3. 傅里叶变换法：频域分析周期性
%
% 作者: CUMCU数学建模团队
% 日期: 2024

    % 输入参数验证和默认值设置
    if nargin < 3
        error('thickness_algorithm:InvalidInput', '至少需要三个输入参数');
    end
    
    if nargin < 4
        const = constants();
        n_epi = const.n_sic;  % 默认使用SiC折射率
    end
    
    if length(wavenumber) ~= length(reflectance)
        error('thickness_algorithm:SizeMismatch', '波数和反射率数组长度必须相同');
    end
    
    % 数据预处理
    [wavenumber_clean, reflectance_clean] = preprocess_data(wavenumber, reflectance);
    
    % 初始化结果结构体
    results = struct();
    results.input_angle = angle;
    results.refractive_index = n_epi;
    results.data_points = length(wavenumber_clean);
    
    % 方法1：极值法计算厚度
    try
        [thickness_extrema, extrema_info] = extrema_method(wavenumber_clean, reflectance_clean, n_epi, angle);
        results.extrema_method.thickness = thickness_extrema;
        results.extrema_method.info = extrema_info;
        results.extrema_method.success = true;
    catch ME
        results.extrema_method.thickness = NaN;
        results.extrema_method.error = ME.message;
        results.extrema_method.success = false;
        thickness_extrema = NaN;
    end
    
    % 方法2：相位拟合法
    try
        [thickness_phase, phase_info] = phase_fitting_method(wavenumber_clean, reflectance_clean, n_epi, angle);
        results.phase_fitting.thickness = thickness_phase;
        results.phase_fitting.info = phase_info;
        results.phase_fitting.success = true;
    catch ME
        results.phase_fitting.thickness = NaN;
        results.phase_fitting.error = ME.message;
        results.phase_fitting.success = false;
        thickness_phase = NaN;
    end
    
    % 方法3：傅里叶变换法
    try
        [thickness_fft, fft_info] = fft_method(wavenumber_clean, reflectance_clean, n_epi, angle);
        results.fft_method.thickness = thickness_fft;
        results.fft_method.info = fft_info;
        results.fft_method.success = true;
    catch ME
        results.fft_method.thickness = NaN;
        results.fft_method.error = ME.message;
        results.fft_method.success = false;
        thickness_fft = NaN;
    end
    
    % 综合多种方法的结果
    valid_thicknesses = [];
    method_weights = [];
    
    if ~isnan(thickness_extrema)
        valid_thicknesses = [valid_thicknesses, thickness_extrema];
        method_weights = [method_weights, 0.4];  % 极值法权重
    end
    
    if ~isnan(thickness_phase)
        valid_thicknesses = [valid_thicknesses, thickness_phase];
        method_weights = [method_weights, 0.4];  % 相位拟合法权重
    end
    
    if ~isnan(thickness_fft)
        valid_thicknesses = [valid_thicknesses, thickness_fft];
        method_weights = [method_weights, 0.2];  % FFT法权重
    end
    
    if isempty(valid_thicknesses)
        error('thickness_algorithm:NoValidResult', '所有计算方法都失败');
    end
    
    % 加权平均计算最终厚度
    method_weights = method_weights / sum(method_weights);  % 归一化权重
    thickness = sum(valid_thicknesses .* method_weights);
    
    % 计算结果统计信息
    results.final_thickness = thickness;
    results.thickness_std = std(valid_thicknesses);
    results.thickness_range = [min(valid_thicknesses), max(valid_thicknesses)];
    results.num_valid_methods = length(valid_thicknesses);
    results.method_weights = method_weights;
    
    % 输出计算结果
    fprintf('\n=== 厚度计算结果 ===\n');
    fprintf('入射角: %.1f°\n', angle);
    fprintf('折射率: %.3f\n', n_epi);
    fprintf('数据点数: %d\n', results.data_points);
    fprintf('\n各方法计算结果:\n');
    
    if results.extrema_method.success
        fprintf('  极值法: %.3f μm\n', thickness_extrema);
    end
    if results.phase_fitting.success
        fprintf('  相位拟合法: %.3f μm\n', thickness_phase);
    end
    if results.fft_method.success
        fprintf('  FFT法: %.3f μm\n', thickness_fft);
    end
    
    fprintf('\n最终结果: %.3f ± %.3f μm\n', thickness, results.thickness_std);
    
end

function [wavenumber_clean, reflectance_clean] = preprocess_data(wavenumber, reflectance)
% 数据预处理：去除异常值和平滑
    
    % 去除NaN和Inf值
    valid_idx = isfinite(wavenumber) & isfinite(reflectance);
    wavenumber_clean = wavenumber(valid_idx);
    reflectance_clean = reflectance(valid_idx);
    
    % 排序数据
    [wavenumber_clean, sort_idx] = sort(wavenumber_clean);
    reflectance_clean = reflectance_clean(sort_idx);
    
    % 简单平滑滤波
    if length(reflectance_clean) > 5
        reflectance_clean = smooth(reflectance_clean, 3);
    end
    
end

function [thickness, info] = extrema_method(wavenumber, reflectance, n_epi, angle)
% 极值法计算厚度
    
    const = constants();
    theta_rad = angle * const.deg2rad;
    
    % 计算折射角
    sin_theta_t = (const.n_air / n_epi) * sin(theta_rad);
    cos_theta_t = sqrt(1 - sin_theta_t^2);
    
    % 寻找极值点
    [~, max_locs] = findpeaks(reflectance, 'MinPeakDistance', 5);
    [~, min_locs] = findpeaks(-reflectance, 'MinPeakDistance', 5);
    
    % 合并并排序极值点
    all_locs = sort([max_locs, min_locs]);
    
    if length(all_locs) < 2
        error('extrema_method:InsufficientPeaks', '找到的极值点不足');
    end
    
    % 计算相邻极值点间的厚度
    thickness_estimates = [];
    for i = 1:length(all_locs)-1
        wn1 = wavenumber(all_locs(i));
        wn2 = wavenumber(all_locs(i+1));
        
        % 厚度公式：d = 1/(2*n*cos(θ)*Δν)
        delta_wn = abs(wn2 - wn1);
        d_est = 1 / (2 * n_epi * cos_theta_t * delta_wn);
        thickness_estimates = [thickness_estimates, d_est * 10000];  % 转换为μm
    end
    
    thickness = median(thickness_estimates);
    
    info.num_extrema = length(all_locs);
    info.thickness_estimates = thickness_estimates;
    info.std_dev = std(thickness_estimates);
    
end

function [thickness, info] = phase_fitting_method(wavenumber, reflectance, n_epi, angle)
% 相位拟合法计算厚度
    
    % 提取相位信息（简化处理）
    % 这里使用反射率的对数来近似相位
    log_R = log(max(reflectance, 0.01));  % 避免log(0)
    
    % 线性拟合
    p = polyfit(wavenumber, log_R, 1);
    
    % 从拟合斜率估计厚度
    const = constants();
    theta_rad = angle * const.deg2rad;
    sin_theta_t = (const.n_air / n_epi) * sin(theta_rad);
    cos_theta_t = sqrt(1 - sin_theta_t^2);
    
    % 厚度估计（经验公式）
    thickness = abs(p(1)) * 10000 / (4 * pi * n_epi * cos_theta_t);
    
    info.fit_coefficients = p;
    info.r_squared = corrcoef(wavenumber, log_R);
    info.r_squared = info.r_squared(1,2)^2;
    
end

function [thickness, info] = fft_method(wavenumber, reflectance, n_epi, angle)
% 傅里叶变换法计算厚度
    
    % 确保数据等间距
    wn_min = min(wavenumber);
    wn_max = max(wavenumber);
    n_points = length(wavenumber);
    wn_uniform = linspace(wn_min, wn_max, n_points);
    
    % 插值到等间距网格
    R_uniform = interp1(wavenumber, reflectance, wn_uniform, 'linear');
    
    % 去除直流分量
    R_ac = R_uniform - mean(R_uniform);
    
    % FFT分析
    Y = fft(R_ac);
    P = abs(Y).^2;
    
    % 频率轴
    fs = 1 / (wn_uniform(2) - wn_uniform(1));  % 采样频率
    f = (0:n_points-1) * fs / n_points;
    
    % 寻找主频率（排除直流分量）
    [~, max_idx] = max(P(2:floor(n_points/2)));
    max_idx = max_idx + 1;  % 补偿索引偏移
    
    dominant_freq = f(max_idx);
    
    % 从主频率计算厚度
    const = constants();
    theta_rad = angle * const.deg2rad;
    sin_theta_t = (const.n_air / n_epi) * sin(theta_rad);
    cos_theta_t = sqrt(1 - sin_theta_t^2);
    
    thickness = 1 / (2 * n_epi * cos_theta_t * dominant_freq) * 10000;  % 转换为μm
    
    info.dominant_frequency = dominant_freq;
    info.power_spectrum = P;
    info.frequency_axis = f;
    
end