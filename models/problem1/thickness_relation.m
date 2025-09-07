function [thickness, relation_func] = thickness_relation(wavenumber, reflectance, n_epi, theta_i)
% THICKNESS_RELATION - 建立外延层厚度与光学参数的关系
%
% 基于干涉理论建立厚度与反射光谱的数学关系
% 通过分析反射率极值位置确定厚度
%
% 输入参数:
%   wavenumber - 波数数组 (cm^-1)
%   reflectance - 反射率数组 (%)
%   n_epi - 外延层折射率 (无量纲)
%   theta_i - 入射角 (度)
%
% 输出参数:
%   thickness - 计算得到的厚度 (μm)
%   relation_func - 厚度关系函数句柄
%
% 数学原理:
%   当相位差 δ = 2mπ 时出现反射极大
%   当相位差 δ = (2m+1)π 时出现反射极小
%   厚度公式: d = (P_i - 0.5) * λ_i / (2 * n_epi * cos(θ_t))
%


    % 输入参数验证
    if nargin < 4
        error('thickness_relation:InvalidInput', '需要提供四个输入参数');
    end
    
    if length(wavenumber) ~= length(reflectance)
        error('thickness_relation:SizeMismatch', '波数和反射率数组长度必须相同');
    end
    
    if n_epi <= 0
        error('thickness_relation:InvalidRefractiveIndex', '折射率必须为正数');
    end
    
    % 加载物理常数
    const = constants();
    
    % 角度转换
    theta_i_rad = theta_i * const.deg2rad;
    
    % 计算外延层中的折射角
    sin_theta_t = (const.n_air / n_epi) * sin(theta_i_rad);
    if abs(sin_theta_t) > 1
        error('thickness_relation:TotalReflection', '发生全反射');
    end
    theta_t = asin(sin_theta_t);
    cos_theta_t = cos(theta_t);
    
    % 波数转换为波长 (μm)
    wavelength = 10000 ./ wavenumber;  % cm^-1 转 μm
    
    % 寻找反射率的极值点
    [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakHeight', max(reflectance)*0.1);
    [valleys, valley_locs] = findpeaks(-reflectance, 'MinPeakHeight', -max(reflectance)*0.9);
    valleys = -valleys;
    
    % 合并极值点并排序
    all_extrema_locs = [peak_locs, valley_locs];
    all_extrema_vals = [peaks, valleys];
    [sorted_locs, sort_idx] = sort(all_extrema_locs);
    sorted_vals = all_extrema_vals(sort_idx);
    
    % 确定极值类型（极大或极小）
    extrema_types = zeros(size(sorted_locs));
    for i = 1:length(sorted_locs)
        if ismember(sorted_locs(i), peak_locs)
            extrema_types(i) = 1;  % 极大
        else
            extrema_types(i) = -1; % 极小
        end
    end
    
    % 计算每个极值对应的干涉级数
    num_extrema = length(sorted_locs);
    if num_extrema < 2
        error('thickness_relation:InsufficientExtrema', '找到的极值点不足，无法计算厚度');
    end
    
    % 使用相邻极值点计算厚度
    thickness_estimates = [];
    
    for i = 1:num_extrema-1
        % 获取相邻两个极值点的波长
        lambda1 = wavelength(sorted_locs(i));
        lambda2 = wavelength(sorted_locs(i+1));
        
        % 计算干涉级数差
        if extrema_types(i) == extrema_types(i+1)
            % 同类型极值，级数差为1
            delta_m = 1;
        else
            % 不同类型极值，级数差为0.5
            delta_m = 0.5;
        end
        
        % 根据厚度公式计算
        % Δm = 2*n*d*cos(θ_t) * (1/λ2 - 1/λ1)
        thickness_est = delta_m / (2 * n_epi * cos_theta_t * abs(1/lambda2 - 1/lambda1));
        thickness_estimates = [thickness_estimates, thickness_est];
    end
    
    % 取厚度估计的平均值
    thickness = mean(thickness_estimates);
    
    % 定义厚度关系函数
    relation_func = @(d) calculate_theoretical_reflectance(d, wavenumber, n_epi, theta_i);
    
    % 输出计算结果
    fprintf('找到 %d 个极值点\n', num_extrema);
    fprintf('厚度估计值: ');
    for i = 1:length(thickness_estimates)
        fprintf('%.3f ', thickness_estimates(i));
    end
    fprintf('μm\n');
    fprintf('平均厚度: %.3f μm\n', thickness);
    fprintf('标准差: %.3f μm\n', std(thickness_estimates));
    
end

function R_theoretical = calculate_theoretical_reflectance(thickness, wavenumber, n_epi, theta_i)
% 计算给定厚度下的理论反射率
    
    const = constants();
    theta_i_rad = theta_i * const.deg2rad;
    
    % 计算每个波数对应的反射率
    R_theoretical = zeros(size(wavenumber));
    
    for i = 1:length(wavenumber)
        lambda = 10000 / wavenumber(i);  % 转换为μm
        
        % 计算相位差
        delta = phase_difference(thickness, n_epi, lambda, theta_i);
        
        % 计算菲涅尔系数
        [r_s, ~, ~, ~] = fresnel_formula(const.n_air, n_epi, theta_i_rad);
        
        % 考虑多次反射的反射率（简化模型）
        R_theoretical(i) = abs(r_s)^2 * (1 + cos(delta));
    end
    
end