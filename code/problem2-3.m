function [reliability_score, analysis_results] = reliability_analysis(thickness_results, measurement_data)
% RELIABILITY_ANALYSIS - 厚度测量结果可靠性分析
%
% 对厚度计算结果进行全面的可靠性评估
% 包括统计分析、误差评估、置信区间计算等
%
% 输入参数:
%   thickness_results - 厚度计算结果结构体
%   measurement_data - 测量数据结构体（包含wavenumber, reflectance等）
%
% 输出参数:
%   reliability_score - 可靠性评分 (0-100)
%   analysis_results - 详细分析结果结构体
%
% 可靠性评估指标:
%   1. 方法一致性：多种方法结果的一致程度
%   2. 数据质量：信噪比、数据完整性
%   3. 拟合优度：理论模型与实测数据的符合程度
%   4. 统计置信度：置信区间和不确定度
%
% 作者: CUMCU数学建模团队
% 日期: 2024

    % 输入参数验证
    if nargin < 2
        error('reliability_analysis:InvalidInput', '需要提供厚度结果和测量数据');
    end
    
    % 初始化分析结果结构体
    analysis_results = struct();
    analysis_results.timestamp = datetime('now');
    
    % 加载参数配置
    params = parameters();
    
    % 1. 方法一致性分析
    [consistency_score, consistency_info] = analyze_method_consistency(thickness_results);
    analysis_results.consistency = consistency_info;
    
    % 2. 数据质量分析
    [quality_score, quality_info] = analyze_data_quality(measurement_data);
    analysis_results.data_quality = quality_info;
    
    % 3. 拟合优度分析
    [fitting_score, fitting_info] = analyze_fitting_quality(thickness_results, measurement_data);
    analysis_results.fitting_quality = fitting_info;
    
    % 4. 统计置信度分析
    [confidence_score, confidence_info] = analyze_statistical_confidence(thickness_results);
    analysis_results.statistical_confidence = confidence_info;
    
    % 5. 物理合理性检查
    [physics_score, physics_info] = analyze_physical_reasonableness(thickness_results);
    analysis_results.physical_reasonableness = physics_info;
    
    % 计算综合可靠性评分
    weights = [0.25, 0.20, 0.25, 0.20, 0.10];  % 各指标权重
    scores = [consistency_score, quality_score, fitting_score, confidence_score, physics_score];
    
    reliability_score = sum(scores .* weights);
    
    % 存储评分详情
    analysis_results.scores.consistency = consistency_score;
    analysis_results.scores.data_quality = quality_score;
    analysis_results.scores.fitting_quality = fitting_score;
    analysis_results.scores.statistical_confidence = confidence_score;
    analysis_results.scores.physical_reasonableness = physics_score;
    analysis_results.scores.overall = reliability_score;
    analysis_results.scores.weights = weights;
    
    % 生成可靠性等级
    if reliability_score >= 90
        reliability_grade = '优秀';
    elseif reliability_score >= 80
        reliability_grade = '良好';
    elseif reliability_score >= 70
        reliability_grade = '中等';
    elseif reliability_score >= 60
        reliability_grade = '一般';
    else
        reliability_grade = '较差';
    end
    
    analysis_results.reliability_grade = reliability_grade;
    
    % 生成建议和警告
    [recommendations, warnings] = generate_recommendations(analysis_results);
    analysis_results.recommendations = recommendations;
    analysis_results.warnings = warnings;
    
    % 输出分析报告
    print_reliability_report(reliability_score, reliability_grade, analysis_results);
    
end

function [score, info] = analyze_method_consistency(thickness_results)
% 分析多种方法结果的一致性
    
    valid_methods = {};
    valid_thicknesses = [];
    
    % 收集有效的厚度结果
    if isfield(thickness_results, 'extrema_method') && thickness_results.extrema_method.success
        valid_methods{end+1} = 'extrema';
        valid_thicknesses(end+1) = thickness_results.extrema_method.thickness;
    end
    
    if isfield(thickness_results, 'phase_fitting') && thickness_results.phase_fitting.success
        valid_methods{end+1} = 'phase_fitting';
        valid_thicknesses(end+1) = thickness_results.phase_fitting.thickness;
    end
    
    if isfield(thickness_results, 'fft_method') && thickness_results.fft_method.success
        valid_methods{end+1} = 'fft';
        valid_thicknesses(end+1) = thickness_results.fft_method.thickness;
    end
    
    info.num_valid_methods = length(valid_methods);
    info.valid_methods = valid_methods;
    info.thickness_values = valid_thicknesses;
    
    if length(valid_thicknesses) < 2
        score = 50;  % 只有一种方法，中等评分
        info.consistency_measure = NaN;
        info.relative_std = NaN;
    else
        % 计算相对标准差
        mean_thickness = mean(valid_thicknesses);
        std_thickness = std(valid_thicknesses);
        relative_std = std_thickness / mean_thickness * 100;
        
        info.mean_thickness = mean_thickness;
        info.std_thickness = std_thickness;
        info.relative_std = relative_std;
        
        % 根据相对标准差评分
        if relative_std < 1
            score = 100;
        elseif relative_std < 2
            score = 90;
        elseif relative_std < 5
            score = 80;
        elseif relative_std < 10
            score = 70;
        else
            score = 50;
        end
        
        info.consistency_measure = 100 - relative_std;
    end
    
end

function [score, info] = analyze_data_quality(measurement_data)
% 分析测量数据质量
    
    wavenumber = measurement_data.wavenumber;
    reflectance = measurement_data.reflectance;
    
    % 数据完整性检查
    total_points = length(wavenumber);
    valid_points = sum(isfinite(wavenumber) & isfinite(reflectance));
    completeness = valid_points / total_points * 100;
    
    % 信噪比估计
    signal_power = var(reflectance);
    noise_estimate = estimate_noise_level(reflectance);
    snr = 10 * log10(signal_power / noise_estimate^2);
    
    % 数据平滑度评估
    smoothness = assess_data_smoothness(reflectance);
    
    % 动态范围评估
    dynamic_range = (max(reflectance) - min(reflectance)) / max(reflectance) * 100;
    
    info.total_points = total_points;
    info.valid_points = valid_points;
    info.completeness = completeness;
    info.snr = snr;
    info.smoothness = smoothness;
    info.dynamic_range = dynamic_range;
    
    % 综合评分
    completeness_score = min(completeness, 100);
    snr_score = min(max((snr - 10) * 5, 0), 100);  % SNR > 10dB 为好
    smoothness_score = smoothness;
    range_score = min(dynamic_range * 2, 100);  % 动态范围 > 50% 为好
    
    score = mean([completeness_score, snr_score, smoothness_score, range_score]);
    
end

function noise_level = estimate_noise_level(data)
% 估计数据噪声水平
    
    % 使用高频分量估计噪声
    if length(data) > 10
        diff_data = diff(data);
        noise_level = std(diff_data) / sqrt(2);
    else
        noise_level = std(data) * 0.1;  % 粗略估计
    end
    
end

function smoothness_score = assess_data_smoothness(data)
% 评估数据平滑度
    
    if length(data) < 3
        smoothness_score = 50;
        return;
    end
    
    % 计算二阶差分
    second_diff = diff(data, 2);
    roughness = std(second_diff);
    signal_level = std(data);
    
    % 平滑度评分（越小越平滑）
    relative_roughness = roughness / signal_level;
    smoothness_score = max(0, 100 - relative_roughness * 1000);
    
end

function [score, info] = analyze_fitting_quality(thickness_results, measurement_data)
% 分析拟合质量
    
    if ~isfield(thickness_results, 'final_thickness')
        score = 0;
        info.error = '无有效厚度结果';
        return;
    end
    
    thickness = thickness_results.final_thickness;
    wavenumber = measurement_data.wavenumber;
    reflectance = measurement_data.reflectance;
    angle = measurement_data.angle;
    n_epi = measurement_data.n_epi;
    
    % 计算理论反射率
    try
        R_theory = calculate_theoretical_spectrum(thickness, wavenumber, n_epi, angle);
        
        % 计算拟合优度指标
        R_squared = calculate_r_squared(reflectance, R_theory);
        rmse = sqrt(mean((reflectance - R_theory).^2));
        mae = mean(abs(reflectance - R_theory));
        
        info.r_squared = R_squared;
        info.rmse = rmse;
        info.mae = mae;
        
        % 根据R²评分
        if R_squared > 0.95
            score = 100;
        elseif R_squared > 0.90
            score = 90;
        elseif R_squared > 0.80
            score = 80;
        elseif R_squared > 0.70
            score = 70;
        else
            score = 50;
        end
        
    catch ME
        score = 0;
        info.error = ME.message;
    end
    
end

function R_theory = calculate_theoretical_spectrum(thickness, wavenumber, n_epi, angle)
% 计算理论反射光谱
    
    const = constants();
    R_theory = zeros(size(wavenumber));
    
    for i = 1:length(wavenumber)
        lambda = 10000 / wavenumber(i);  % 转换为μm
        
        % 计算相位差
        delta = phase_difference(thickness, n_epi, lambda, angle);
        
        % 简化的反射率模型
        theta_rad = angle * const.deg2rad;
        [r_s, ~, ~, ~] = fresnel_formula(const.n_air, n_epi, theta_rad);
        
        R_theory(i) = abs(r_s)^2 * (1 + 0.5 * cos(delta));
    end
    
end

function r_squared = calculate_r_squared(y_actual, y_predicted)
% 计算决定系数R²
    
    ss_res = sum((y_actual - y_predicted).^2);
    ss_tot = sum((y_actual - mean(y_actual)).^2);
    r_squared = 1 - ss_res / ss_tot;
    
end

function [score, info] = analyze_statistical_confidence(thickness_results)
% 分析统计置信度
    
    if ~isfield(thickness_results, 'thickness_std')
        score = 50;
        info.error = '无标准差信息';
        return;
    end
    
    thickness = thickness_results.final_thickness;
    thickness_std = thickness_results.thickness_std;
    
    % 计算相对不确定度
    relative_uncertainty = thickness_std / thickness * 100;
    
    % 计算置信区间（假设正态分布）
    confidence_level = 0.95;
    z_score = 1.96;  % 95%置信区间
    margin_error = z_score * thickness_std;
    
    info.thickness = thickness;
    info.std_dev = thickness_std;
    info.relative_uncertainty = relative_uncertainty;
    info.confidence_interval = [thickness - margin_error, thickness + margin_error];
    info.margin_error = margin_error;
    
    % 根据相对不确定度评分
    if relative_uncertainty < 1
        score = 100;
    elseif relative_uncertainty < 2
        score = 90;
    elseif relative_uncertainty < 5
        score = 80;
    elseif relative_uncertainty < 10
        score = 70;
    else
        score = 50;
    end
    
end

function [score, info] = analyze_physical_reasonableness(thickness_results)
% 分析物理合理性
    
    thickness = thickness_results.final_thickness;
    
    % 厚度合理性检查
    reasonable_range = [0.1, 200];  % μm
    
    info.thickness = thickness;
    info.reasonable_range = reasonable_range;
    info.within_range = thickness >= reasonable_range(1) && thickness <= reasonable_range(2);
    
    if info.within_range
        score = 100;
    else
        score = 0;
    end
    
end

function [recommendations, warnings] = generate_recommendations(analysis_results)
% 生成建议和警告
    
    recommendations = {};
    warnings = {};
    
    % 基于各项分析结果生成建议
    if analysis_results.scores.consistency < 80
        recommendations{end+1} = '建议增加测量次数以提高方法一致性';
    end
    
    if analysis_results.scores.data_quality < 70
        recommendations{end+1} = '建议改善测量条件以提高数据质量';
        warnings{end+1} = '数据质量较低，可能影响结果可靠性';
    end
    
    if analysis_results.scores.fitting_quality < 80
        recommendations{end+1} = '建议检查理论模型的适用性';
    end
    
    if analysis_results.scores.statistical_confidence < 70
        warnings{end+1} = '统计置信度较低，结果不确定度较大';
    end
    
end

function print_reliability_report(score, grade, results)
% 打印可靠性分析报告
    
    fprintf('\n=== 可靠性分析报告 ===\n');
    fprintf('分析时间: %s\n', char(results.timestamp));
    fprintf('\n综合可靠性评分: %.1f/100 (%s)\n', score, grade);
    
    fprintf('\n各项指标评分:\n');
    fprintf('  方法一致性: %.1f/100\n', results.scores.consistency);
    fprintf('  数据质量: %.1f/100\n', results.scores.data_quality);
    fprintf('  拟合质量: %.1f/100\n', results.scores.fitting_quality);
    fprintf('  统计置信度: %.1f/100\n', results.scores.statistical_confidence);
    fprintf('  物理合理性: %.1f/100\n', results.scores.physical_reasonableness);
    
    if ~isempty(results.warnings)
        fprintf('\n警告:\n');
        for i = 1:length(results.warnings)
            fprintf('  - %s\n', results.warnings{i});
        end
    end
    
    if ~isempty(results.recommendations)
        fprintf('\n建议:\n');
        for i = 1:length(results.recommendations)
            fprintf('  - %s\n', results.recommendations{i});
        end
    end
    
    fprintf('\n');
    
end