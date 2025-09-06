%% 综合测试脚本
% 测试红外干涉法测量碳化硅外延层厚度项目的所有功能模块

function test_results = test_all()
    % 运行所有测试
    % 输出参数：
    %   test_results - 测试结果结构体
    
    fprintf('=== 红外干涉法测量项目综合测试 ===\n');
    fprintf('开始时间: %s\n', datestr(now));
    
    % 初始化测试结果
    test_results = struct(
        'total_tests', 0,
        'passed_tests', 0,
        'failed_tests', 0,
        'test_details', {},
        'start_time', now,
        'end_time', [],
        'duration', []
    );
    
    %% 测试列表
    test_functions = {
        @test_constants,
        @test_parameters,
        @test_fresnel_formula,
        @test_phase_difference,
        @test_thickness_relation,
        @test_thickness_algorithm,
        @test_reliability_analysis,
        @test_multi_beam_conditions,
        @test_data_loader,
        @test_signal_processor,
        @test_plot_utils,
        @test_math_utils
    };
    
    test_names = {
        '物理常数测试',
        '系统参数测试',
        '菲涅尔公式测试',
        '相位差计算测试',
        '厚度关系推导测试',
        '厚度算法测试',
        '可靠性分析测试',
        '多光束干涉条件测试',
        '数据加载器测试',
        '信号处理器测试',
        '绘图工具测试',
        '数学工具测试'
    };
    
    %% 执行测试
    for i = 1:length(test_functions)
        fprintf('\n--- %s ---\n', test_names{i});
        
        try
            % 执行测试函数
            test_result = test_functions{i}();
            
            % 记录测试结果
            test_results.total_tests = test_results.total_tests + 1;
            
            if test_result.passed
                test_results.passed_tests = test_results.passed_tests + 1;
                fprintf('✓ 测试通过\n');
            else
                test_results.failed_tests = test_results.failed_tests + 1;
                fprintf('✗ 测试失败: %s\n', test_result.error_message);
            end
            
            % 保存详细结果
            test_results.test_details{i} = struct(
                'name', test_names{i},
                'passed', test_result.passed,
                'error_message', test_result.error_message,
                'execution_time', test_result.execution_time
            );
            
        catch ME
            % 测试函数本身出错
            test_results.total_tests = test_results.total_tests + 1;
            test_results.failed_tests = test_results.failed_tests + 1;
            
            fprintf('✗ 测试执行失败: %s\n', ME.message);
            
            test_results.test_details{i} = struct(
                'name', test_names{i},
                'passed', false,
                'error_message', ME.message,
                'execution_time', 0
            );
        end
    end
    
    %% 完成测试
    test_results.end_time = now;
    test_results.duration = test_results.end_time - test_results.start_time;
    
    %% 输出测试摘要
    fprintf('\n=== 测试摘要 ===\n');
    fprintf('总测试数: %d\n', test_results.total_tests);
    fprintf('通过测试: %d\n', test_results.passed_tests);
    fprintf('失败测试: %d\n', test_results.failed_tests);
    fprintf('成功率: %.1f%%\n', test_results.passed_tests / test_results.total_tests * 100);
    fprintf('总耗时: %.2f 秒\n', test_results.duration * 24 * 3600);
    fprintf('结束时间: %s\n', datestr(test_results.end_time));
    
    %% 显示失败的测试
    if test_results.failed_tests > 0
        fprintf('\n=== 失败的测试 ===\n');
        for i = 1:length(test_results.test_details)
            if ~test_results.test_details{i}.passed
                fprintf('- %s: %s\n', test_results.test_details{i}.name, ...
                       test_results.test_details{i}.error_message);
            end
        end
    end
end

%% 测试物理常数
function result = test_constants()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 测试常数加载
        constants = load_constants();
        
        % 检查必要的常数
        required_fields = {'c', 'n_air', 'n_si', 'n_sic'};
        for i = 1:length(required_fields)
            if ~isfield(constants, required_fields{i})
                error('缺少必要常数: %s', required_fields{i});
            end
        end
        
        % 检查常数值的合理性
        if constants.c <= 0
            error('光速值不合理');
        end
        
        if constants.n_air < 0.9 || constants.n_air > 1.1
            error('空气折射率值不合理');
        end
        
        if constants.n_si < 3 || constants.n_si > 4
            error('硅折射率值不合理');
        end
        
        if constants.n_sic < 2 || constants.n_sic > 3
            error('SiC折射率值不合理');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试系统参数
function result = test_parameters()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 测试参数加载
        params = load_parameters();
        
        % 检查必要的参数
        required_fields = {'measurement', 'algorithm', 'processing'};
        for i = 1:length(required_fields)
            if ~isfield(params, required_fields{i})
                error('缺少必要参数组: %s', required_fields{i});
            end
        end
        
        % 检查测量参数
        if ~isfield(params.measurement, 'angle_range')
            error('缺少角度范围参数');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试菲涅尔公式
function result = test_fresnel_formula()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 测试基本功能
        n1 = 1.0;  % 空气
        n2 = 3.4;  % 硅
        angle = 45; % 度
        
        [r_s, r_p, t_s, t_p] = fresnel_formula(n1, n2, angle);
        
        % 检查输出格式
        if ~isnumeric(r_s) || ~isnumeric(r_p) || ~isnumeric(t_s) || ~isnumeric(t_p)
            error('输出不是数值类型');
        end
        
        % 检查能量守恒（近似）
        R_s = abs(r_s)^2;
        T_s = abs(t_s)^2 * n2 / n1;
        energy_conservation_s = abs(R_s + T_s - 1);
        
        if energy_conservation_s > 0.01
            error('s偏振能量不守恒，误差: %.4f', energy_conservation_s);
        end
        
        % 测试边界条件
        [r_normal, ~, ~, ~] = fresnel_formula(n1, n2, 0);
        expected_r_normal = (n1 - n2) / (n1 + n2);
        
        if abs(r_normal - expected_r_normal) > 1e-10
            error('垂直入射反射系数计算错误');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试相位差计算
function result = test_phase_difference()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 测试基本功能
        wavenumber = 1000; % cm^-1
        thickness = 1.0;   % μm
        n_film = 2.5;
        angle = 30; % 度
        
        phase_diff = calculate_phase_difference(wavenumber, thickness, n_film, angle);
        
        % 检查输出
        if ~isnumeric(phase_diff) || ~isreal(phase_diff)
            error('相位差输出格式错误');
        end
        
        % 检查相位差的合理性
        if phase_diff < 0
            error('相位差不应为负值');
        end
        
        % 测试厚度为零的情况
        phase_diff_zero = calculate_phase_difference(wavenumber, 0, n_film, angle);
        if abs(phase_diff_zero) > 1e-10
            error('零厚度时相位差应为零');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试厚度关系推导
function result = test_thickness_relation()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 创建测试数据
        wavenumber = linspace(800, 1200, 100);
        thickness_true = 2.0; % μm
        angle = 45;
        
        % 生成理论反射率数据
        reflectance = generate_theoretical_reflectance(wavenumber, thickness_true, angle);
        
        % 测试厚度关系推导
        relation_result = derive_thickness_relation(wavenumber, reflectance, angle);
        
        % 检查输出结构
        required_fields = {'thickness_estimate', 'theory', 'analysis'};
        for i = 1:length(required_fields)
            if ~isfield(relation_result, required_fields{i})
                error('缺少输出字段: %s', required_fields{i});
            end
        end
        
        % 检查厚度估计的准确性
        thickness_error = abs(relation_result.thickness_estimate - thickness_true) / thickness_true;
        if thickness_error > 0.1  % 10%误差容限
            error('厚度估计误差过大: %.2f%%', thickness_error * 100);
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试厚度算法
function result = test_thickness_algorithm()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 创建测试数据
        wavenumber = linspace(800, 1200, 200);
        thickness_true = 1.5; % μm
        angle = 30;
        
        % 生成带噪声的测试数据
        reflectance = generate_theoretical_reflectance(wavenumber, thickness_true, angle);
        noise = 0.01 * randn(size(reflectance));
        reflectance_noisy = reflectance + noise;
        
        data = struct('wavenumber', wavenumber, 'reflectance', reflectance_noisy);
        
        % 测试厚度算法
        algorithm_result = determine_thickness_algorithm(data, angle);
        
        % 检查输出结构
        if ~isfield(algorithm_result, 'thickness_estimates')
            error('缺少厚度估计结果');
        end
        
        % 检查是否有合理的厚度估计
        estimates = algorithm_result.thickness_estimates;
        if isempty(estimates) || all(isnan(estimates))
            error('未能获得有效的厚度估计');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试可靠性分析
function result = test_reliability_analysis()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 创建测试数据
        thickness_estimates = [1.48, 1.52, 1.49, 1.51, 1.50, 1.47, 1.53];
        
        % 测试可靠性分析
        reliability_result = analyze_reliability(thickness_estimates);
        
        % 检查输出结构
        required_fields = {'statistics', 'precision', 'outliers', 'consistency', 'reliability_score'};
        for i = 1:length(required_fields)
            if ~isfield(reliability_result, required_fields{i})
                error('缺少输出字段: %s', required_fields{i});
            end
        end
        
        % 检查统计结果的合理性
        if reliability_result.statistics.mean < 0 || reliability_result.statistics.std < 0
            error('统计结果不合理');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试多光束干涉条件
function result = test_multi_beam_conditions()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 创建测试数据
        wavenumber = linspace(800, 1200, 150);
        
        % 生成具有多光束干涉特征的数据
        reflectance = 0.3 + 0.2 * cos(0.1 * wavenumber) + 0.1 * cos(0.05 * wavenumber);
        
        % 测试多光束干涉条件判断
        [is_multi_beam, analysis_result] = multi_beam_conditions(reflectance);
        
        % 检查输出格式
        if ~islogical(is_multi_beam)
            error('多光束判断结果应为逻辑值');
        end
        
        if ~isstruct(analysis_result)
            error('分析结果应为结构体');
        end
        
        % 检查分析结果结构
        if ~isfield(analysis_result, 'judgment')
            error('缺少判断结果字段');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试数据加载器
function result = test_data_loader()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 创建临时测试数据文件
        test_file = 'test_data_temp.csv';
        test_wavenumber = linspace(800, 1200, 50)';
        test_reflectance = 0.3 + 0.1 * sin(0.01 * test_wavenumber);
        
        % 写入CSV文件
        csvwrite(test_file, [test_wavenumber, test_reflectance]);
        
        % 测试数据加载
        [data, info] = data_loader(test_file);
        
        % 检查加载结果
        if ~isstruct(data) || ~isfield(data, 'wavenumber') || ~isfield(data, 'reflectance')
            error('数据结构不正确');
        end
        
        if length(data.wavenumber) ~= length(test_wavenumber)
            error('加载的数据长度不正确');
        end
        
        % 清理临时文件
        if exist(test_file, 'file')
            delete(test_file);
        end
        
        result.passed = true;
        
    catch ME
        % 清理临时文件
        if exist('test_data_temp.csv', 'file')
            delete('test_data_temp.csv');
        end
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试信号处理器
function result = test_signal_processor()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 创建测试数据
        wavenumber = linspace(800, 1200, 100)';
        reflectance = 0.3 + 0.1 * sin(0.02 * wavenumber) + 0.02 * randn(size(wavenumber));
        data = struct('wavenumber', wavenumber, 'reflectance', reflectance);
        
        % 测试信号处理
        options = struct('smooth_method', 'sgolay', 'smooth_window', 5);
        processed_data = signal_processor(data, options);
        
        % 检查处理结果
        if ~isstruct(processed_data)
            error('处理结果应为结构体');
        end
        
        if ~isfield(processed_data, 'wavenumber') || ~isfield(processed_data, 'reflectance')
            error('处理结果缺少必要字段');
        end
        
        % 检查数据长度
        if length(processed_data.wavenumber) ~= length(data.wavenumber)
            error('处理后数据长度改变');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试绘图工具
function result = test_plot_utils()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 创建测试数据
        wavenumber = linspace(800, 1200, 50)';
        reflectance = 0.3 + 0.1 * sin(0.01 * wavenumber);
        data = struct('wavenumber', wavenumber, 'reflectance', reflectance);
        
        % 测试绘图功能（不显示图形）
        options = struct('figure_size', [400, 300]);
        
        % 创建图形但不显示
        fig_handle = plot_spectral_data(data, options);
        
        % 检查图形句柄
        if ~ishandle(fig_handle)
            error('未能创建有效的图形句柄');
        end
        
        % 关闭图形
        close(fig_handle);
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 测试数学工具
function result = test_math_utils()
    start_time = tic;
    result = struct('passed', false, 'error_message', '', 'execution_time', 0);
    
    try
        % 测试数值积分
        func = @(x) x.^2;
        integral_result = numerical_integration(func, 0, 1, 'trapz');
        expected_result = 1/3;
        
        if abs(integral_result - expected_result) > 0.01
            error('数值积分结果不准确');
        end
        
        % 测试数值微分
        x = linspace(0, 2*pi, 100);
        y = sin(x);
        dy = numerical_derivative(y, x, 'central');
        expected_dy = cos(x);
        
        % 检查微分精度（中间部分）
        mid_range = 20:80;
        max_error = max(abs(dy(mid_range) - expected_dy(mid_range)));
        
        if max_error > 0.1
            error('数值微分精度不足');
        end
        
        % 测试插值
        x_orig = [1, 2, 3, 4, 5];
        y_orig = [1, 4, 9, 16, 25];
        x_interp = 2.5;
        y_interp = interpolation(x_orig, y_orig, x_interp, 'linear');
        expected_y = 6.5;
        
        if abs(y_interp - expected_y) > 0.1
            error('插值结果不准确');
        end
        
        result.passed = true;
        
    catch ME
        result.error_message = ME.message;
    end
    
    result.execution_time = toc(start_time);
end

%% 辅助函数：生成理论反射率
function reflectance = generate_theoretical_reflectance(wavenumber, thickness, angle)
    % 简化的理论反射率计算
    
    % 材料参数
    n1 = 1.0;   % 空气
    n2 = 2.5;   % 薄膜
    
    % 转换为波长
    wavelength = 1e7 ./ wavenumber;  % nm
    
    % 计算相位差
    angle_rad = angle * pi / 180;
    cos_theta = cos(angle_rad);
    
    phase_diff = 4 * pi * n2 * thickness * 1000 * cos_theta ./ wavelength;
    
    % 菲涅尔反射系数
    r = (n1 - n2) / (n1 + n2);
    
    % 简化的多光束干涉公式
    reflectance = abs(r)^2 * (1 + 0.5 * cos(phase_diff));