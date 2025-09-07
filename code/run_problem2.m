% RUN_PROBLEM2 - 运行问题二的所有计算
% 这个脚本调用问题二的所有函数并生成结果
%
% 作者: CUMCU数学建模团队
% 日期: 2024

clc; clear; close all;

fprintf('\n=== 开始运行问题二 ===\n');

% 添加路径
addpath(genpath('../utils'));
addpath(genpath('.'));

try
    % 1. 读取附件1数据
    fprintf('\n1. 读取附件1数据...\n');
    data_file = '../附件1.xlsx';
    
    % 使用excel_reader读取数据
    addpath('../utils');
    addpath('../config');
    [data, metadata] = excel_reader(data_file);
    
    % 从结构体中提取数据
    if isstruct(data) && isfield(data, 'wavenumber') && isfield(data, 'reflectance')
        wavenumber = data.wavenumber;
        reflectance = data.reflectance;
        fprintf('从结构体中提取数据成功\n');
    else
        error('数据格式不正确，无法提取波数和反射率数据');
    end
    
    fprintf('数据读取成功: %d个数据点\n', length(wavenumber));
    fprintf('波数范围: %.1f - %.1f cm^-1\n', min(wavenumber), max(wavenumber));
    fprintf('反射率范围: %.2f - %.2f%%\n', min(reflectance), max(reflectance));
    
    % 2. 设置计算参数
    fprintf('\n2. 设置计算参数...\n');
    angles = [10, 15, 20];  % 不同入射角度
    n_sic = 2.55;           % SiC折射率
    
    % 3. 对每个角度进行厚度计算
    fprintf('\n3. 进行厚度计算...\n');
    thickness_results = [];
    angle_results = [];
    
    for i = 1:length(angles)
        angle = angles(i);
        fprintf('\n计算入射角 %.1f° 的情况...\n', angle);
        
        try
            % 调用thickness_algorithm函数
            [thickness, detailed_results] = thickness_algorithm(wavenumber, reflectance, angle, n_sic);
            
            % 存储结果
            thickness_results = [thickness_results, thickness];
            angle_results = [angle_results, angle];
            
            % 存储详细结果
            problem2_results.(['angle_' num2str(angle)]) = detailed_results;
            problem2_results.(['angle_' num2str(angle)]).thickness = thickness;
            
            fprintf('角度 %.1f° 计算成功: %.3f μm\n', angle, thickness);
            
        catch ME
            fprintf('角度 %.1f° 计算失败: %s\n', angle, ME.message);
        end
    end
    
    % 4. 进行可靠性分析
    fprintf('\n4. 进行可靠性分析...\n');
    if ~isempty(thickness_results)
        % 创建测量数据结构
        measurement_data.wavenumber = wavenumber;
        measurement_data.reflectance = reflectance;
        measurement_data.filename = data_file;
        
        % 对第一个成功的结果进行可靠性分析
        first_angle = angle_results(1);
        first_results = problem2_results.(['angle_' num2str(first_angle)]);
        
        try
            [reliability_score, reliability_detail] = reliability_analysis(first_results, measurement_data);
            problem2_results.reliability_score = reliability_score;
            problem2_results.reliability_detail = reliability_detail;
            fprintf('可靠性评分: %.1f/100\n', reliability_score);
        catch ME
            fprintf('可靠性分析失败: %s\n', ME.message);
        end
    end
    
    % 5. 汇总结果
    fprintf('\n5. 汇总计算结果...\n');
    problem2_results.angles = angle_results;
    problem2_results.thickness_values = thickness_results;
    problem2_results.input_data.wavenumber = wavenumber;
    problem2_results.input_data.reflectance = reflectance;
    problem2_results.material = 'SiC';
    problem2_results.refractive_index = n_sic;
    problem2_results.timestamp = datetime('now');
    
    if ~isempty(thickness_results)
        problem2_results.mean_thickness = mean(thickness_results);
        problem2_results.std_thickness = std(thickness_results);
        problem2_results.thickness_range = [min(thickness_results), max(thickness_results)];
        
        fprintf('\n=== 问题二计算结果汇总 ===\n');
        fprintf('测试角度数量: %d\n', length(angle_results));
        fprintf('平均厚度: %.3f μm\n', problem2_results.mean_thickness);
        fprintf('厚度标准差: %.3f μm\n', problem2_results.std_thickness);
        fprintf('厚度范围: %.3f - %.3f μm\n', problem2_results.thickness_range(1), problem2_results.thickness_range(2));
        
        for i = 1:length(angle_results)
            fprintf('角度 %.1f°: %.3f μm\n', angle_results(i), thickness_results(i));
        end
    end
    
    % 6. 保存结果
    fprintf('\n6. 保存结果...\n');
    save('results/problem2_results.mat', 'problem2_results');
    fprintf('结果已保存到 results/problem2_results.mat\n');
    
    fprintf('\n=== 问题二运行完成 ===\n');
    
catch ME
    fprintf('\n错误: %s\n', ME.message);
    fprintf('错误位置: %s (第%d行)\n', ME.stack(1).name, ME.stack(1).line);
end