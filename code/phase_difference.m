function delta = phase_difference(thickness, n_epi, lambda, theta_i)
% PHASE_DIFFERENCE - 计算外延层中光的相位差
%
% 计算光在外延层中传播产生的相位差，考虑多次反射的影响
% 基于光程差和界面反射相位变化
%
% 输入参数:
%   thickness - 外延层厚度 (μm)
%   n_epi - 外延层折射率 (无量纲)
%   lambda - 真空中波长 (μm)
%   theta_i - 入射角 (度)
%
% 输出参数:
%   delta - 相位差 (弧度)
%
% 数学公式:
%   delta = (4*π*n_epi*d*cos(theta_t))/λ + φ1 - φ2
%   其中 φ1, φ2 是界面反射引起的相位变化
%

    % 输入参数验证
    if nargin < 4
        error('phase_difference:InvalidInput', '需要提供四个输入参数');
    end
    
    if thickness <= 0
        error('phase_difference:InvalidThickness', '厚度必须为正数');
    end
    
    if n_epi <= 0
        error('phase_difference:InvalidRefractiveIndex', '折射率必须为正数');
    end
    
    if lambda <= 0
        error('phase_difference:InvalidWavelength', '波长必须为正数');
    end
    
    if theta_i < 0 || theta_i > 90
        error('phase_difference:InvalidAngle', '入射角必须在0到90度之间');
    end
    
    % 加载物理常数
    const = constants();
    
    % 角度转换为弧度
    theta_i_rad = theta_i * const.deg2rad;
    
    % 假设从空气入射到外延层
    n_air = const.n_air;
    
    % 根据斯涅尔定律计算外延层中的折射角
    sin_theta_t = (n_air / n_epi) * sin(theta_i_rad);
    
    % 检查是否发生全反射
    if abs(sin_theta_t) > 1
        error('phase_difference:TotalReflection', '发生全反射，无法计算相位差');
    end
    
    theta_t = asin(sin_theta_t);
    cos_theta_t = cos(theta_t);
    
    % 计算光程差产生的相位变化
    optical_path_phase = (4 * pi * n_epi * thickness * cos_theta_t) / lambda;
    
    % 计算界面反射相位变化
    % 空气-外延层界面
    [r_s1, r_p1, ~, ~] = fresnel_formula(n_air, n_epi, theta_i_rad);
    
    % 外延层-衬底界面（假设衬底折射率）
    n_substrate = const.n_si;  % 假设为硅衬底
    [r_s2, r_p2, ~, ~] = fresnel_formula(n_epi, n_substrate, theta_t);
    
    % 计算反射相位（取s偏振光为例，也可以考虑p偏振光）
    phi1 = angle(r_s1);  % 第一个界面反射相位
    phi2 = angle(r_s2);  % 第二个界面反射相位
    
    % 总相位差
    delta = optical_path_phase + phi1 - phi2;
    
    % 将相位差限制在[-π, π]范围内
    delta = mod(delta + pi, 2*pi) - pi;
    
    % 输出调试信息（可选）
    if nargout == 0
        fprintf('厚度: %.3f μm\n', thickness);
        fprintf('入射角: %.1f°\n', theta_i);
        fprintf('折射角: %.1f°\n', theta_t * const.rad2deg);
        fprintf('光程相位: %.3f rad\n', optical_path_phase);
        fprintf('界面相位差: %.3f rad\n', phi1 - phi2);
        fprintf('总相位差: %.3f rad\n', delta);
    end
    
end