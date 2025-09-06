function [model_results] = model()
% MODEL - 建立红外干涉数学模型的主程序
%
% 功能描述:
%   建立基于红外干涉法的碳化硅外延层厚度测量数学模型
%   包括菲涅尔公式、相位差计算和厚度关系推导
%
% 输出参数:
%   model_results - 模型结果结构体
%     .fresnel_coeffs - 菲涅尔系数
%     .phase_relation - 相位差关系
%     .thickness_formula - 厚度计算公式

    fprintf('=== 红外干涉法数学模型建立 ===\n');
    
    % 初始化结果结构体
    model_results = struct();
    
    % 1. 建立菲涅尔公式模型
    fprintf('1.1 建立菲涅尔公式模型\n');
    model_results.fresnel_coeffs = build_fresnel_model();
    
    % 2. 建立相位差计算模型
    fprintf('1.2 建立相位差计算模型\n');
    model_results.phase_relation = build_phase_model();
    
    % 3. 建立厚度关系模型
    fprintf('1.3 建立厚度关系模型\n');
    model_results.thickness_formula = build_thickness_model();
    
    % 4. 干涉光强度模型
    fprintf('1.4 建立干涉光强度模型\n');
    model_results.intensity_model = build_intensity_model();
    
    fprintf('数学模型建立完成\n');
end

function fresnel_coeffs = build_fresnel_model()
% 建立菲涅尔公式模型
    
    % 定义材料参数
    n_air = 1.0;        % 空气折射率
    n_sic = 2.55;       % SiC折射率
    n_substrate = 3.42; % 衬底折射率
    
    % 菲涅尔反射和透射系数公式
    fresnel_coeffs = struct();
    fresnel_coeffs.n_air = n_air;
    fresnel_coeffs.n_sic = n_sic;
    fresnel_coeffs.n_substrate = n_substrate;
    
    % s偏振反射系数公式: rs = (n1*cos(theta1) - n2*cos(theta2)) / (n1*cos(theta1) + n2*cos(theta2))
    % p偏振反射系数公式: rp = (n2*cos(theta1) - n1*cos(theta2)) / (n2*cos(theta1) + n1*cos(theta2))
    % s偏振透射系数公式: ts = 2*n1*cos(theta1) / (n1*cos(theta1) + n2*cos(theta2))
    % p偏振透射系数公式: tp = 2*n1*cos(theta1) / (n2*cos(theta1) + n1*cos(theta2))
    
    fresnel_coeffs.formulas = {
        's偏振反射系数: rs = (n1*cos(θ1) - n2*cos(θ2)) / (n1*cos(θ1) + n2*cos(θ2))';
        'p偏振反射系数: rp = (n2*cos(θ1) - n1*cos(θ2)) / (n2*cos(θ1) + n1*cos(θ2))';
        's偏振透射系数: ts = 2*n1*cos(θ1) / (n1*cos(θ1) + n2*cos(θ2))';
        'p偏振透射系数: tp = 2*n1*cos(θ1) / (n2*cos(θ1) + n1*cos(θ2))'
    };
    
    fprintf('  菲涅尔公式模型建立完成\n');
end

function phase_relation = build_phase_model()
% 建立相位差计算模型
    
    phase_relation = struct();
    
    % 相位差公式: δ = (4π * n * d * cos(θr)) / λ
    % 其中: n - 外延层折射率, d - 厚度, θr - 折射角, λ - 波长
    phase_relation.formula = 'δ = (4π * n * d * cos(θr)) / λ';
    phase_relation.description = '相位差与厚度、折射率、入射角和波长的关系';
    
    % 光程差公式
    phase_relation.optical_path_diff = '2 * n * d * cos(θr)';
    
    fprintf('  相位差计算模型建立完成\n');
end

function thickness_formula = build_thickness_model()
% 建立厚度关系模型
    
    thickness_formula = struct();
    
    % 基于干涉极值条件的厚度公式
    % 对于相邻极值点: Δφ = 2π
    % 厚度公式: d = (m * λ) / (2 * n * cos(θr))
    thickness_formula.basic_formula = 'd = (m * λ) / (2 * n * cos(θr))';
    thickness_formula.description = '基于干涉级数的厚度计算公式';
    
    % 考虑相位跳变的修正公式
    thickness_formula.corrected_formula = 'd = ((m - 0.5) * λ) / (2 * n * cos(θr))';
    
    % 多波长厚度计算公式
    thickness_formula.multi_wavelength = 'd = (λ1 * λ2) / (2 * n * (λ2 - λ1) * cos(θr))';
    
    fprintf('  厚度关系模型建立完成\n');
end

function intensity_model = build_intensity_model()
% 建立干涉光强度模型
    
    intensity_model = struct();
    
    % 干涉光强度公式: I(λ) = a(λ) + b(λ) * cos(φ(λ))
    % 其中: a(λ) - 背景强度项, b(λ) - 调制幅度项, φ(λ) - 相位项
    intensity_model.formula = 'I(λ) = a(λ) + b(λ) * cos(φ(λ))';
    intensity_model.description = '干涉光强度与波长的关系模型';
    
    % 相位项公式
    intensity_model.phase_term = 'φ(λ) = (4π * n * d * cos(θr)) / λ';
    
    % 背景强度项和调制幅度项的物理意义
    intensity_model.background = 'a(λ): 平均反射强度，与材料本身光学性质相关';
    intensity_model.modulation = 'b(λ): 干涉调制深度，与界面反射率差异相关';
    
    fprintf('  干涉光强度模型建立完成\n');
end