function [corrected_thickness, correction_factors] = model_correction(original_thickness, interference_data, correction_method)
% MODEL_CORRECTION - 多光束干涉模型修正算法
%
% 基于多光束干涉理论修正厚度计算结果
% 消除多光束干涉对测量精度的影响
%
% 输入参数:
%   original_thickness - 原始厚度计算结果 (μm)
%   interference_data - 干涉数据结构体
%     .wavenumber - 波数数组 (cm^-1)
%     .reflectance - 反射率数组
%     .angle - 入射角 (度)
%     .material_params - 材料参数
%   correction_method - 修正方法 ('airy', 'fabry_perot', 'multi_reflection')
%
% 输出参数:
%   corrected_thickness - 修正后的厚度 (μm)
%   correction_factors - 修正因子和详细信息
%
% 修正方法:
%   1. Airy函数修正 - 基于Airy函数的多光束干涉理论
%   2. Fabry-Perot修正 - 考虑多次反射的F-P干涉仪模型
%   3. 多次反射修正 - 逐级计算多次反射贡献



    % 输入参数验证
    if nargin < 2
        error('model_correction:InvalidInput', '需要提供原始厚度和干涉数据');
    end
    
    if nargin < 3
        correction_method = 'airy'; % 默认使用Airy函数修正
    end
    
    % 参数检查
    if isempty(original_thickness) || original_thickness <= 0
        error('model_correction:InvalidThickness', '原始厚度必须为正数');
    end
    
    % 加载常数和参数
    const = constants();
    params = parameters();
    
    % 提取干涉数据
    wavenumber = interference_data.wavenumber;
    reflectance = interference_data.reflectance;
    angle = interference_data.angle;
    material_params = interference_data.material_params;
    
    fprintf('\n=== 多光束干涉模型修正 ===\n');
    fprintf('原始厚度: %.3f μm\n', original_thickness);
    fprintf('修正方法: %s\n', correction_method);
    fprintf('入射角度: %.1f°\n', angle);
    
    % 初始化修正因子结构
    correction_factors = struct();
    correction_factors.method = correction_method;
    correction_factors.original_thickness = original_thickness;
    correction_factors.angle = angle;
    correction_factors.timestamp = datetime('now');
    
    % 计算基本光学参数
    angle_rad = deg2rad(angle);
    n1 = 1.0; % 空气折射率
    n2 = material_params.n_substrate; % 衬底折射率
    
    if isfield(material_params, 'n_epilayer')
        n_epi = material_params.n_epilayer;
    else
        n_epi = n2; % 如果没有外延层信息，使用衬底折射率
    end
    
    % 计算折射角
    theta_t = asin(n1 * sin(angle_rad) / n_epi);
    
    % 根据修正方法进行计算
    switch lower(correction_method)
        case 'airy'
            [corrected_thickness, correction_factors] = airy_correction(...
                original_thickness, wavenumber, reflectance, n1, n_epi, n2, theta_t, correction_factors);
            
        case 'fabry_perot'
            [corrected_thickness, correction_factors] = fabry_perot_correction(...
                original_thickness, wavenumber, reflectance, n1, n_epi, n2, theta_t, correction_factors);
            
        case 'multi_reflection'
            [corrected_thickness, correction_factors] = multi_reflection_correction(...
                original_thickness, wavenumber, reflectance, n1, n_epi, n2, theta_t, correction_factors);
            
        otherwise
            warning('model_correction:UnknownMethod', '未知修正方法，使用默认Airy修正');
            [corrected_thickness, correction_factors] = airy_correction(...
                original_thickness, wavenumber, reflectance, n1, n_epi, n2, theta_t, correction_factors);
    end
    
    % 计算修正效果
    correction_factors.thickness_change = corrected_thickness - original_thickness;
    correction_factors.relative_change = correction_factors.thickness_change / original_thickness * 100;
    
    fprintf('修正后厚度: %.3f μm\n', corrected_thickness);
    fprintf('厚度变化: %.3f μm (%.2f%%)\n', ...
        correction_factors.thickness_change, correction_factors.relative_change);
    
end

%% Airy函数修正
function [corrected_thickness, factors] = airy_correction(d0, wavenumber, reflectance, n1, n2, n3, theta_t, factors)
    % 基于Airy函数的多光束干涉修正
    
    fprintf('使用Airy函数修正...\n');
    
    % 计算菲涅尔反射系数
    r12 = (n1 - n2) / (n1 + n2); % 空气-外延层界面
    r23 = (n2 - n3) / (n2 + n3); % 外延层-衬底界面
    
    % 计算反射率
    R12 = r12^2;
    R23 = r23^2;
    
    % 计算精细度系数
    F = 4 * R23 / (1 - R23)^2;
    factors.finesse_coefficient = F;
    factors.r12 = r12;
    factors.r23 = r23;
    factors.R12 = R12;
    factors.R23 = R23;
    
    % 波长转换
    lambda = 1 ./ wavenumber * 1e4; % 转换为μm
    
    % 计算相位
    delta = 4 * pi * n2 * d0 * cos(theta_t) ./ lambda;
    
    % Airy函数理论反射率
    R_theory = R12 + (1 - R12)^2 * R23 ./ (1 - R12 * R23 + 2 * sqrt(R12 * R23) .* cos(delta));
    
    % 计算实测与理论的差异
    R_measured = reflectance / 100; % 转换为小数
    residual = R_measured - R_theory;
    rms_error = sqrt(mean(residual.^2));
    
    factors.rms_error_original = rms_error;
    factors.R_theory = R_theory;
    factors.residual = residual;
    
    % 基于残差进行厚度修正
    % 使用最小二乘法优化厚度
    options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8);
    
    objective_function = @(d) calculate_airy_residual(d, wavenumber, R_measured, n1, n2, n3, theta_t);
    
    % 在原始厚度附近搜索最优值
    d_range = [d0 * 0.8, d0 * 1.2];
    corrected_thickness = fminbnd(objective_function, d_range(1), d_range(2), options);
    
    % 计算修正后的误差
    final_residual = objective_function(corrected_thickness);
    factors.rms_error_corrected = final_residual;
    factors.improvement_ratio = factors.rms_error_original / final_residual;
    
    fprintf('Airy修正完成，RMS误差改善: %.2f倍\n', factors.improvement_ratio);
end

%% Fabry-Perot修正
function [corrected_thickness, factors] = fabry_perot_correction(d0, wavenumber, reflectance, n1, n2, n3, theta_t, factors)
    % 基于Fabry-Perot干涉仪模型的修正
    
    fprintf('使用Fabry-Perot修正...\n');
    
    % 计算透射系数
    t12 = 2 * n1 / (n1 + n2);
    t21 = 2 * n2 / (n2 + n1);
    t23 = 2 * n2 / (n2 + n3);
    t32 = 2 * n3 / (n3 + n2);
    
    % 反射系数
    r12 = (n1 - n2) / (n1 + n2);
    r21 = -r12;
    r23 = (n2 - n3) / (n2 + n3);
    r32 = -r23;
    
    factors.transmission_coeffs = [t12, t21, t23, t32];
    factors.reflection_coeffs = [r12, r21, r23, r32];
    
    % 波长转换
    lambda = 1 ./ wavenumber * 1e4; % μm
    
    % F-P腔的自由光谱范围和精细度
    FSR = lambda.^2 ./ (2 * n2 * d0 * cos(theta_t)); % 自由光谱范围
    finesse = pi * sqrt(abs(r23)) / (1 - abs(r23)); % 精细度
    
    factors.free_spectral_range = mean(FSR);
    factors.finesse = finesse;
    
    % F-P传输函数
    delta = 4 * pi * n2 * d0 * cos(theta_t) ./ lambda;
    T_fp = (t12 * t23)^2 ./ (1 - r21 * r23 * exp(1i * delta)).^2;
    R_fp = abs(1 - T_fp); % 反射率
    
    % 优化厚度
    R_measured = reflectance / 100;
    objective_function = @(d) calculate_fp_residual(d, wavenumber, R_measured, n1, n2, n3, theta_t);
    
    options = optimset('Display', 'off', 'TolFun', 1e-8);
    d_range = [d0 * 0.8, d0 * 1.2];
    corrected_thickness = fminbnd(objective_function, d_range(1), d_range(2), options);
    
    factors.optimization_range = d_range;
    factors.rms_error_corrected = objective_function(corrected_thickness);
    
    fprintf('Fabry-Perot修正完成\n');
end

%% 多次反射修正
function [corrected_thickness, factors] = multi_reflection_correction(d0, wavenumber, reflectance, n1, n2, n3, theta_t, factors)
    % 考虑多次反射的修正方法
    
    fprintf('使用多次反射修正...\n');
    
    % 计算多次反射系数
    r12 = (n1 - n2) / (n1 + n2);
    r23 = (n2 - n3) / (n2 + n3);
    t12 = 2 * n1 / (n1 + n2);
    t21 = 2 * n2 / (n2 + n1);
    
    % 波长转换
    lambda = 1 ./ wavenumber * 1e4;
    
    % 相位因子
    beta = 2 * pi * n2 * cos(theta_t) ./ lambda;
    
    % 多次反射级数展开（前10项）
    max_order = 10;
    R_total = zeros(size(lambda));
    
    for m = 0:max_order
        % m次内部反射的贡献
        phase_factor = exp(1i * 2 * m * beta * d0);
        amplitude = t12 * t21 * r23^m * r12^m * phase_factor;
        R_contribution = abs(amplitude).^2;
        R_total = R_total + R_contribution;
        
        % 记录各阶贡献
        if m <= 3
            factors.(sprintf('reflection_order_%d', m)) = mean(R_contribution);
        end
    end
    
    % 加上直接反射
    R_total = abs(r12)^2 + R_total;
    
    factors.max_reflection_order = max_order;
    factors.total_reflectance_theory = R_total;
    
    % 优化厚度
    R_measured = reflectance / 100;
    objective_function = @(d) calculate_multi_reflection_residual(d, wavenumber, R_measured, n1, n2, n3, theta_t, max_order);
    
    options = optimset('Display', 'off', 'TolFun', 1e-8);
    d_range = [d0 * 0.8, d0 * 1.2];
    corrected_thickness = fminbnd(objective_function, d_range(1), d_range(2), options);
    
    factors.rms_error_corrected = objective_function(corrected_thickness);
    
    fprintf('多次反射修正完成\n');
end

%% 辅助函数：Airy残差计算
function rms_error = calculate_airy_residual(d, wavenumber, R_measured, n1, n2, n3, theta_t)
    r12 = (n1 - n2) / (n1 + n2);
    r23 = (n2 - n3) / (n2 + n3);
    R12 = r12^2;
    R23 = r23^2;
    
    % 确保波数为正值
    valid_idx = wavenumber > 0;
    if ~any(valid_idx)
        rms_error = inf;
        return;
    end
    
    lambda = 1e4 ./ wavenumber(valid_idx);
    delta = 4 * pi * n2 * d * cos(theta_t) ./ lambda;
    
    R_theory = R12 + (1 - R12)^2 * R23 ./ (1 - R12 * R23 + 2 * sqrt(R12 * R23) .* cos(delta));
    
    residual = R_measured(valid_idx) - R_theory;
    rms_error = sqrt(mean(residual.^2));
end

%% 辅助函数：F-P残差计算
function rms_error = calculate_fp_residual(d, wavenumber, R_measured, n1, n2, n3, theta_t)
    r23 = (n2 - n3) / (n2 + n3);
    t12 = 2 * n1 / (n1 + n2);
    t23 = 2 * n2 / (n2 + n3);
    
    % 确保波数为正值
    valid_idx = wavenumber > 0;
    if ~any(valid_idx)
        rms_error = inf;
        return;
    end
    
    lambda = 1e4 ./ wavenumber(valid_idx);
    delta = 4 * pi * n2 * d * cos(theta_t) ./ lambda;
    
    T_fp = abs(t12 * t23).^2 ./ abs(1 - r23^2 * exp(1i * delta)).^2;
    R_fp = 1 - T_fp;
    
    residual = R_measured(valid_idx) - R_fp;
    rms_error = sqrt(mean(residual.^2));
end

%% 辅助函数：多次反射残差计算
function rms_error = calculate_multi_reflection_residual(d, wavenumber, R_measured, n1, n2, n3, theta_t, max_order)
    r12 = (n1 - n2) / (n1 + n2);
    r23 = (n2 - n3) / (n2 + n3);
    t12 = 2 * n1 / (n1 + n2);
    t21 = 2 * n2 / (n2 + n1);
    
    % 确保波数为正值
    valid_idx = wavenumber > 0;
    if ~any(valid_idx)
        rms_error = inf;
        return;
    end
    
    lambda = 1e4 ./ wavenumber(valid_idx);
    beta = 2 * pi * n2 * cos(theta_t) ./ lambda;
    
    R_total = abs(r12)^2 * ones(size(lambda));
    
    for m = 1:max_order
        phase_factor = exp(1i * 2 * m * beta * d);
        amplitude = t12 * t21 * r23^m * r12^(m-1) * phase_factor;
        R_total = R_total + abs(amplitude).^2;
    end
    
    residual = R_measured(valid_idx) - R_total;
    rms_error = sqrt(mean(residual.^2));
end