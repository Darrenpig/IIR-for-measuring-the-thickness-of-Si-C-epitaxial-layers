%% 红外干涉法碳化硅外延层厚度测量 - 主程序
% 作者：CUMCU团队
% 日期：2024
% 描述：基于红外干涉法的SiC外延层厚度测量系统主程序

function main()
    % 清理工作空间
    clc; clear; close all;
    
    % 添加项目路径
    addpath(genpath(pwd));
    
    % 显示欢迎信息
    fprintf('=== 红外干涉法碳化硅外延层厚度测量系统 ===\n');
    fprintf('系统初始化中...\n');
    
    % 加载配置参数
    constants = load_constants();
    params = load_parameters();
    
    % 显示菜单
    while true
        fprintf('\n请选择要执行的问题：\n');
        fprintf('1. 问题一：建立数学模型\n');
        fprintf('2. 问题二：厚度确定算法\n');
        fprintf('3. 问题三：多光束干涉分析\n');
        fprintf('4. 运行所有问题\n');
        fprintf('0. 退出\n');
        
        choice = input('请输入选择 (0-4): ');
        
        switch choice
            case 1
                run_problem1(constants, params);
            case 2
                run_problem2(constants, params);
            case 3
                run_problem3(constants, params);
            case 4
                run_all_problems(constants, params);
            case 0
                fprintf('程序退出。\n');
                break;
            otherwise
                fprintf('无效选择，请重新输入。\n');
        end
    end
end

%% 运行问题一
function run_problem1(constants, params)
    fprintf('\n=== 问题一：建立数学模型 ===\n');
    
    % 菲涅尔公式计算
    [r_s, r_p, t_s, t_p] = fresnel_formula(constants.n_air, constants.n_sic, params.theta_i);
    
    % 相位差计算
    delta = phase_difference(params.thickness_test, constants.n_sic, params.lambda, params.theta_i);
    
    % 厚度关系推导
    thickness_relation_demo(constants, params);
    
    % 保存结果
    results.r_s = r_s;
    results.r_p = r_p;
    results.delta = delta;
    save('results/problem1_results.mat', 'results');
    
    fprintf('问题一计算完成，结果已保存到 results/problem1_results.mat\n');
end

%% 运行问题二
function run_problem2(constants, params)
    fprintf('\n=== 问题二：厚度确定算法 ===\n');
    
    % 加载SiC数据
    data1 = load_excel_data('data/raw/attachment1.xlsx');
    data2 = load_excel_data('data/raw/attachment2.xlsx');
    
    % 处理SiC数据
    thickness1 = thickness_algorithm(data1.wavenumber, data1.reflectance, 10, constants.n_sic);
    thickness2 = thickness_algorithm(data2.wavenumber, data2.reflectance, 15, constants.n_sic);
    
    % 可靠性分析
    reliability = reliability_analysis(thickness1, thickness2);
    
    % 保存结果
    results.thickness_10deg = thickness1;
    results.thickness_15deg = thickness2;
    results.reliability = reliability;
    save('results/problem2_results.mat', 'results');
    
    fprintf('问题二计算完成，结果已保存到 results/problem2_results.mat\n');
end

%% 运行问题三
function run_problem3(constants, params)
    fprintf('\n=== 问题三：多光束干涉分析 ===\n');
    
    % 加载硅片数据
    si_data1 = load_excel_data('data/raw/attachment3.xlsx');
    si_data2 = load_excel_data('data/raw/attachment4.xlsx');
    
    % 多光束干涉条件判断
    is_multi_beam1 = multi_beam_conditions(si_data1.reflectance, params.multi_beam_threshold);
    is_multi_beam2 = multi_beam_conditions(si_data2.reflectance, params.multi_beam_threshold);
    
    % 硅片数据分析
    si_analysis1 = si_data_analyzer(si_data1, 10);
    si_analysis2 = si_data_analyzer(si_data2, 15);
    
    % 模型修正
    corrected_model = model_correction(si_analysis1, si_analysis2);
    
    % 保存结果
    results.is_multi_beam_10deg = is_multi_beam1;
    results.is_multi_beam_15deg = is_multi_beam2;
    results.si_analysis_10deg = si_analysis1;
    results.si_analysis_15deg = si_analysis2;
    results.corrected_model = corrected_model;
    save('results/problem3_results.mat', 'results');
    
    fprintf('问题三计算完成，结果已保存到 results/problem3_results.mat\n');
end

%% 运行所有问题
function run_all_problems(constants, params)
    fprintf('\n=== 运行所有问题 ===\n');
    run_problem1(constants, params);
    run_problem2(constants, params);
    run_problem3(constants, params);
    fprintf('\n所有问题计算完成！\n');
end