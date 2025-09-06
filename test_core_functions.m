% TEST_CORE_FUNCTIONS - 测试核心功能模块
%
% 不依赖外部数据文件的功能测试脚本
% 验证所有核心算法是否正常工作
%
% 作者: CUMCU数学建模团队
% 日期: 2024

clc; clear; close all;

fprintf('\n=== 核心功能测试 ===\n');

% 添加路径
addpath(genpath('models'));
addpath(genpath('utils'));
addpath(genpath('config'));

% 加载常数和参数
try
    const = constants();
    params = parameters();
    fprintf('✓ 常数和参数加载成功\n');
catch ME
    fprintf('✗ 常数和参数加载失败: %s\n', ME.message);
    return;
end

%% 测试问题1：数学模型
fprintf('\n--- 测试问题1：数学模型 ---\n');

% 测试菲涅尔公式
try
    n1 = const.n_air;
    n2 = const.n_sic;
    theta_i = deg2rad(10);
    
    [r_s, r_p, t_s, t_p] = fresnel_formula(n1, n2, theta_i);
    
    fprintf('✓ 菲涅尔公式计算成功\n');
    fprintf('  s偏振反射系数: %.4f\n', r_s);
    fprintf('  p偏振反射系数: %.4f\n', r_p);
    
catch ME
    fprintf('✗ 菲涅尔公式测试失败: %s\n', ME.message);
end

% 测试相位差计算
try
    thickness = 5.0;
    lambda = 10.0;
    angle = 10;
    
    delta = phase_difference(thickness, n2, lambda, angle);
    
    fprintf('✓ 相位差计算成功\n');
    fprintf('  相位差: %.4f rad\n', delta);
    
catch ME
    fprintf('✗ 相位差计算测试失败: %s\n', ME.message);
end

% 测试厚度关系
try
    % 创建模拟数据
    wavenumber = linspace(400, 4000, 200)';
    % 模拟干涉条纹
    phase = 4 * pi * n2 * 3.0 * cos(deg2rad(10)) ./ (1e4 ./ wavenumber);  % 假设3μm厚度
    reflectance = (40 + 30 * cos(phase))';  % 模拟反射率，确保列向量
    
    [calc_thickness, relation_func] = thickness_relation(wavenumber, reflectance, n2, 10);
    
    fprintf('✓ 厚度关系推导成功\n');
    fprintf('  计算厚度: %.3f μm\n', calc_thickness);
    
catch ME
    fprintf('✗ 厚度关系测试失败: %s\n', ME.message);
end

%% 测试问题2：厚度算法
fprintf('\n--- 测试问题2：厚度算法 ---\n');

% 创建模拟SiC数据
try
    % 模拟光谱数据
    sim_wavenumber = linspace(500, 3500, 300)';
    true_thickness = 4.5;  % 真实厚度
    
    % 基于理论公式生成模拟反射率
    phase_sim = 4 * pi * const.n_sic * true_thickness * cos(deg2rad(10)) ./ (1e4 ./ sim_wavenumber);
    r12 = (const.n_air - const.n_sic) / (const.n_air + const.n_sic);
    sim_reflectance = abs(r12)^2 * (1 + 0.8 * cos(phase_sim)) * 100;
    
    % 添加少量噪声
    sim_reflectance = sim_reflectance + 2 * randn(size(sim_reflectance));
    
    % 测试厚度算法
    [calc_thickness, thickness_results] = thickness_algorithm(sim_wavenumber, sim_reflectance, 10, const.n_sic);
    
    fprintf('✓ 厚度算法测试成功\n');
    fprintf('  真实厚度: %.3f μm\n', true_thickness);
    fprintf('  计算厚度: %.3f μm\n', calc_thickness);
    fprintf('  相对误差: %.2f%%\n', abs(calc_thickness - true_thickness) / true_thickness * 100);
    
    % 测试可靠性分析
    dummy_processed_data = struct('quality_score', 0.85, 'snr', 25, 'wavenumber', sim_wavenumber, 'reflectance', sim_reflectance, 'angle', 10, 'n_epi', const.n_sic);
    thickness_struct = struct('final_thickness', calc_thickness, 'results', thickness_results);
    [reliability_score, reliability_result] = reliability_analysis(thickness_struct, dummy_processed_data);
    
    fprintf('✓ 可靠性分析成功\n');
    fprintf('  可靠性等级: %s\n', reliability_result.reliability_grade);
    
catch ME
    fprintf('✗ 厚度算法测试失败: %s\n', ME.message);
end

%% 测试问题3：多光束干涉
fprintf('\n--- 测试问题3：多光束干涉 ---\n');

% 测试多光束干涉条件判断
try
    % 创建具有多光束干涉特征的数据
    test_wavenumber = linspace(600, 3000, 250)';
    
    % 模拟多光束干涉（高精细度）
    R = 0.3;  % 反射率
    F = 4 * R / (1 - R)^2;  % 精细度系数
    phase_mb = 4 * pi * const.n_si * 2.0 * cos(deg2rad(15)) ./ (1e4 ./ test_wavenumber);
    
    % Airy函数
    test_reflectance = R + (1 - R)^2 * R ./ (1 - R^2 + 2*R*cos(phase_mb)) * 100;
    
    % 材料参数
    material_params = struct();
    material_params.n_substrate = const.n_si;
    material_params.n_epilayer = const.n_si;
    
    [is_multi_beam, analysis_details] = multi_beam_conditions(test_wavenumber, test_reflectance, 15, material_params);
    
    fprintf('✓ 多光束干涉条件判断成功\n');
    fprintf('  是否多光束干涉: %s\n', char(string(is_multi_beam)));
    
    if is_multi_beam
        % 测试模型修正
        [initial_thickness_val, initial_results] = thickness_algorithm(test_wavenumber, test_reflectance, 15, const.n_si);
        
        interference_data = struct();
        interference_data.wavenumber = test_wavenumber;
        interference_data.reflectance = test_reflectance;
        interference_data.angle = 15;
        interference_data.material_params = material_params;
        
        [corrected_thickness, correction_factors] = model_correction(initial_thickness_val, interference_data, 'airy');
        
        fprintf('✓ 模型修正测试成功\n');
        fprintf('  初始厚度: %.3f μm\n', initial_thickness_val);
        fprintf('  修正厚度: %.3f μm\n', corrected_thickness);
        fprintf('  修正幅度: %.2f%%\n', correction_factors.relative_change);
    end
    
catch ME
    fprintf('✗ 多光束干涉测试失败: %s\n', ME.message);
end

%% 测试工具函数
fprintf('\n--- 测试工具函数 ---\n');

% 测试数据加载器（使用模拟数据）
try
    % 创建临时测试数据文件
    test_data = [sim_wavenumber, sim_reflectance];
    test_file = 'temp_test_data.csv';
    writematrix(test_data, test_file);
    
    [spectrum_data, metadata] = data_loader(test_file, 10, 'SiC');
    
    fprintf('✓ 数据加载器测试成功\n');
    fprintf('  加载数据点数: %d\n', spectrum_data.num_points);
    fprintf('  数据质量: %s\n', metadata.quality.overall_quality);
    
    % 清理临时文件
    delete(test_file);
    
catch ME
    fprintf('✗ 数据加载器测试失败: %s\n', ME.message);
    if exist('test_file', 'var') && exist(test_file, 'file')
        delete(test_file);
    end
end

% 测试绘图功能
try
    % 创建测试数据
    plot_data = struct();
    plot_data.wavenumber = sim_wavenumber;
    plot_data.reflectance = sim_reflectance;
    plot_data.material = 'SiC';
    plot_data.angle = 10;
    
    % 测试绘图（不保存）
    plot_options = struct();
    plot_options.title = '测试光谱图';
    
    fig_handle = plot_spectrum(plot_data, plot_options);
    
    if ishandle(fig_handle)
        fprintf('✓ 绘图功能测试成功\n');
        close(fig_handle);
    else
        fprintf('✗ 绘图功能测试失败\n');
    end
    
catch ME
    fprintf('✗ 绘图功能测试失败: %s\n', ME.message);
end

%% 测试总结
fprintf('\n=== 测试总结 ===\n');
fprintf('核心功能模块测试完成\n');
fprintf('如果所有测试都显示✓，说明系统功能正常\n');
fprintf('如果有✗标记，请检查对应的错误信息\n');

fprintf('\n测试完成。\n');