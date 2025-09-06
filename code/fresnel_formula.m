function [r_s, r_p, t_s, t_p] = fresnel_formula(n1, n2, theta_i)
% FRESNEL_FORMULA - 计算菲涅尔反射和透射系数
%
% 基于菲涅尔公式计算光在两种介质界面的反射和透射系数
% 适用于s偏振光和p偏振光
%
% 输入参数:
%   n1 - 入射介质折射率 (无量纲)
%   n2 - 透射介质折射率 (无量纲)
%   theta_i - 入射角 (弧度)
%
% 输出参数:
%   r_s - s偏振光反射系数 (复数)
%   r_p - p偏振光反射系数 (复数)
%   t_s - s偏振光透射系数 (复数)
%   t_p - p偏振光透射系数 (复数)
%
% 数学公式:
%   r_s = (n1*cos(theta_i) - n2*cos(theta_t)) / (n1*cos(theta_i) + n2*cos(theta_t))
%   r_p = (n2*cos(theta_i) - n1*cos(theta_t)) / (n2*cos(theta_i) + n1*cos(theta_t))
%   t_s = 2*n1*cos(theta_i) / (n1*cos(theta_i) + n2*cos(theta_t))
%   t_p = 2*n1*cos(theta_i) / (n2*cos(theta_i) + n1*cos(theta_t))
%
    % 输入参数验证
    if nargin < 3
        error('fresnel_formula:InvalidInput', '需要提供三个输入参数: n1, n2, theta_i');
    end
    
    if n1 <= 0 || n2 <= 0
        error('fresnel_formula:InvalidRefractiveIndex', '折射率必须为正数');
    end
    
    if theta_i < 0 || theta_i > pi/2
        error('fresnel_formula:InvalidAngle', '入射角必须在0到π/2之间');
    end
    
    % 根据斯涅尔定律计算透射角
    sin_theta_t = (n1 / n2) * sin(theta_i);
    
    % 检查是否发生全反射
    if abs(sin_theta_t) > 1
        % 全反射情况
        cos_theta_t = 1i * sqrt(sin_theta_t^2 - 1);
        theta_t = asin(sin_theta_t);
    else
        % 正常折射情况
        cos_theta_t = sqrt(1 - sin_theta_t^2);
        theta_t = asin(sin_theta_t);
    end
    
    cos_theta_i = cos(theta_i);
    
    % 计算s偏振光的菲涅尔系数
    denominator_s = n1 * cos_theta_i + n2 * cos_theta_t;
    r_s = (n1 * cos_theta_i - n2 * cos_theta_t) / denominator_s;
    t_s = (2 * n1 * cos_theta_i) / denominator_s;
    
    % 计算p偏振光的菲涅尔系数
    denominator_p = n2 * cos_theta_i + n1 * cos_theta_t;
    r_p = (n2 * cos_theta_i - n1 * cos_theta_t) / denominator_p;
    t_p = (2 * n1 * cos_theta_i) / denominator_p;
    
    % 验证能量守恒（对于实数折射率）
    if isreal(n1) && isreal(n2) && isreal(theta_i)
        R_s = abs(r_s)^2;
        T_s = (n2 * cos_theta_t) / (n1 * cos_theta_i) * abs(t_s)^2;
        R_p = abs(r_p)^2;
        T_p = (n2 * cos_theta_t) / (n1 * cos_theta_i) * abs(t_p)^2;
        
        % 检查能量守恒（允许小的数值误差）
        energy_conservation_s = abs(R_s + real(T_s) - 1);
        energy_conservation_p = abs(R_p + real(T_p) - 1);
        
        if energy_conservation_s > 1e-10 || energy_conservation_p > 1e-10
            warning('fresnel_formula:EnergyConservation', ...
                    '能量守恒检查失败，可能存在数值误差');
        end
    end
end