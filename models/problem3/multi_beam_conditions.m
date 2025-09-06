function [is_multi_beam, analysis_results] = multi_beam_conditions(wavenumber, reflectance, angle, material_params)
% MULTI_BEAM_CONDITIONS - 多光束干涉条件判断
%
% 判断给定的光谱数据是否满足多光束干涉条件
% 分析干涉条件、光束强度分布和相干性要求
%
% 输入参数:
%   wavenumber - 波数数组 (cm^-1)
%   reflectance - 反射率数组
%   angle - 入射角度 (度)
%   material_params - 材料参数结构体
%     .n_substrate - 衬底折射率
%     .n_epilayer - 外延层折射率 (可选)
%     .thickness - 外延层厚度 (可选，μm)
%
% 输出参数:
%   is_multi_beam - 逻辑值，是否满足多光束干涉条件
%   analysis_results - 详细分析结果结构体
%
% 判断标准:
%   1. 反射率调制深度 > 阈值
%   2. 干涉条纹周期性
%   3. 相位相干性
%   4. 多次反射强度比
%   5. 光学厚度条件
%
% 作者: CUMCU数学建模团队
% 日期: 2024

    % 输入参数验证
    if nargin < 3
        error('multi_beam_conditions:InvalidInput', '至少需要波数、反射率和角度参数');
    end
    
    if nargin < 4
        material_params = struct();
    end
    
    % 参数检查
    if length(wavenumber) ~= length(reflectance)
        error('multi_beam_conditions:SizeMismatch', '波数和反射率数组长度不匹配');
    end
    
    if length(wavenumber) < 20
        error('multi_beam_conditions:InsufficientData', '数据点数太少，无法进行可靠分析');
    end
    
    % 加载默认参数
    const = constants();
    params = parameters();
    
    % 设置默认材料参数
    if ~isfield(material_params, 'n_substrate')
        material_params.n_substrate = const.n_si;  % 默认硅衬底
    end
    
    % 初始化分析结果
    analysis_results = struct();
    analysis_results.timestamp = datetime('now');
    analysis_results.angle = angle;
    analysis_results.material_params = material_params;
    analysis_results.data_points = length(wavenumber);
    analysis_results.wavenumber_range = [min(wavenumber), max(wavenumber)];
    
    fprintf('\n=== 多光束干涉条件分析 ===\n');
    fprintf('入射角度: %.1f°\n', angle);
    fprintf('数据点数: %d\n', length(wavenumber));
    fprintf('波数范围: %.1f - %.1f cm⁻¹\n', min(wavenumber), max(wavenumber));
    
    % 1. 反射率调制深度分析
    [modulation_analysis] = analyze_modulation_depth(wavenumber, reflectance);
    analysis_results.modulation = modulation_analysis;
    
    % 2. 干涉条纹周期性分析
    [periodicity_analysis] = analyze_interference_periodicity(wavenumber, reflectance);
    analysis_results.periodicity = periodicity_analysis;
    
    % 3. 相位相干性分析
    [coherence_analysis] = analyze_phase_coherence(wavenumber, reflectance, angle);
    analysis_results.coherence = coherence_analysis;
    
    % 4. 多次反射强度分析
    [reflection_analysis] = analyze_multiple_reflections(wavenumber, reflectance, angle, material_params);
    analysis_results.multiple_reflections = reflection_analysis;
    
    % 5. 光学厚度条件检查
    [thickness_analysis] = analyze_optical_thickness_conditions(wavenumber, reflectance, angle, material_params);
    analysis_results.optical_thickness = thickness_analysis;
    
    % 6. 综合判断
    [is_multi_beam, overall_assessment] = make_overall_assessment(...
        modulation_analysis, periodicity_analysis, coherence_analysis, ...
        reflection_analysis, thickness_analysis, params);
    
    analysis_results.overall_assessment = overall_assessment;
    analysis_results.is_multi_beam = is_multi_beam;
    
    % 生成分析报告
    generate_analysis_report(analysis_results);
    
end

function [modulation_analysis] = analyze_modulation_depth(wavenumber, reflectance)
% 分析反射率调制深度
    
    modulation_analysis = struct();
    
    % 计算调制深度
    R_max = max(reflectance);
    R_min = min(reflectance);
    R_mean = mean(reflectance);
    
    modulation_depth = (R_max - R_min) / (R_max + R_min);
    modulation_contrast = (R_max - R_min) / R_mean;
    
    % 计算局部调制深度变化
    window_size = min(50, floor(length(reflectance) / 10));
    local_modulations = [];
    
    for i = 1:window_size:length(reflectance)-window_size+1
        end_idx = min(i+window_size-1, length(reflectance));
        local_R = reflectance(i:end_idx);
        local_mod = (max(local_R) - min(local_R)) / (max(local_R) + min(local_R));
        local_modulations = [local_modulations, local_mod];
    end
    
    modulation_analysis.global_depth = modulation_depth;
    modulation_analysis.contrast = modulation_contrast;
    modulation_analysis.local_depths = local_modulations;
    modulation_analysis.mean_local_depth = mean(local_modulations);
    modulation_analysis.std_local_depth = std(local_modulations);
    
    % 判断标准
    params = parameters();
    modulation_analysis.threshold = params.multi_beam_modulation_threshold;
    modulation_analysis.meets_criterion = modulation_depth > params.multi_beam_modulation_threshold;
    
    fprintf('\n1. 调制深度分析:\n');
    fprintf('   全局调制深度: %.3f\n', modulation_depth);
    fprintf('   平均局部调制深度: %.3f ± %.3f\n', mean(local_modulations), std(local_modulations));
    fprintf('   判断阈值: %.3f\n', params.multi_beam_modulation_threshold);
    fprintf('   满足条件: %s\n', char(string(modulation_analysis.meets_criterion)));
    
end

function [periodicity_analysis] = analyze_interference_periodicity(wavenumber, reflectance)
% 分析干涉条纹周期性
    
    periodicity_analysis = struct();
    
    % 寻找极值点
    [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakDistance', 5);
    [valleys, valley_locs] = findpeaks(-reflectance, 'MinPeakDistance', 5);
    valleys = -valleys;
    
    % 计算峰值间距
    if length(peak_locs) >= 3
        peak_spacings = diff(wavenumber(peak_locs));
        mean_peak_spacing = mean(peak_spacings);
        std_peak_spacing = std(peak_spacings);
        peak_regularity = 1 - std_peak_spacing / mean_peak_spacing;
    else
        peak_spacings = [];
        mean_peak_spacing = NaN;
        std_peak_spacing = NaN;
        peak_regularity = 0;
    end
    
    % 计算谷值间距
    if length(valley_locs) >= 3
        valley_spacings = diff(wavenumber(valley_locs));
        mean_valley_spacing = mean(valley_spacings);
        std_valley_spacing = std(valley_spacings);
        valley_regularity = 1 - std_valley_spacing / mean_valley_spacing;
    else
        valley_spacings = [];
        mean_valley_spacing = NaN;
        std_valley_spacing = NaN;
        valley_regularity = 0;
    end
    
    % 傅里叶变换分析周期性
    N = length(reflectance);
    if N > 10
        % 去除直流分量
        R_ac = reflectance - mean(reflectance);
        
        % FFT分析
        fft_R = fft(R_ac);
        power_spectrum = abs(fft_R).^2;
        
        % 寻找主频率
        [~, max_freq_idx] = max(power_spectrum(2:floor(N/2)));
        max_freq_idx = max_freq_idx + 1;  % 调整索引
        
        dominant_frequency = max_freq_idx / N;
        spectral_purity = power_spectrum(max_freq_idx) / sum(power_spectrum(2:floor(N/2)));
    else
        dominant_frequency = NaN;
        spectral_purity = 0;
    end
    
    % 存储结果
    periodicity_analysis.num_peaks = length(peaks);
    periodicity_analysis.num_valleys = length(valleys);
    periodicity_analysis.peak_spacings = peak_spacings;
    periodicity_analysis.valley_spacings = valley_spacings;
    periodicity_analysis.mean_peak_spacing = mean_peak_spacing;
    periodicity_analysis.mean_valley_spacing = mean_valley_spacing;
    periodicity_analysis.peak_regularity = peak_regularity;
    periodicity_analysis.valley_regularity = valley_regularity;
    periodicity_analysis.dominant_frequency = dominant_frequency;
    periodicity_analysis.spectral_purity = spectral_purity;
    
    % 综合周期性评分
    overall_regularity = (peak_regularity + valley_regularity) / 2;
    periodicity_analysis.overall_regularity = overall_regularity;
    
    % 判断标准
    params = parameters();
    periodicity_analysis.regularity_threshold = params.multi_beam_regularity_threshold;
    periodicity_analysis.meets_criterion = overall_regularity > params.multi_beam_regularity_threshold && ...
                                          length(peaks) >= 3;
    
    fprintf('\n2. 周期性分析:\n');
    fprintf('   峰值数量: %d, 谷值数量: %d\n', length(peaks), length(valleys));
    if ~isnan(mean_peak_spacing)
        fprintf('   平均峰值间距: %.2f ± %.2f cm⁻¹\n', mean_peak_spacing, std_peak_spacing);
    end
    fprintf('   整体规律性: %.3f\n', overall_regularity);
    fprintf('   满足条件: %s\n', char(string(periodicity_analysis.meets_criterion)));
    
end

function [coherence_analysis] = analyze_phase_coherence(wavenumber, reflectance, angle)
% 分析相位相干性
    
    coherence_analysis = struct();
    
    % 计算相位信息（通过希尔伯特变换）
    analytic_signal = hilbert(reflectance - mean(reflectance));
    instantaneous_phase = angle(analytic_signal);
    instantaneous_amplitude = abs(analytic_signal);
    
    % 相位连续性分析
    phase_diff = diff(unwrap(instantaneous_phase));
    phase_smoothness = 1 / (1 + std(phase_diff));
    
    % 振幅稳定性分析
    amplitude_variation = std(instantaneous_amplitude) / mean(instantaneous_amplitude);
    amplitude_stability = 1 / (1 + amplitude_variation);
    
    % 相干长度估计
    autocorr = xcorr(reflectance - mean(reflectance), 'normalized');
    coherence_length_idx = find(autocorr(length(reflectance):end) < 0.5, 1);
    if ~isempty(coherence_length_idx)
        coherence_length = coherence_length_idx;
    else
        coherence_length = length(reflectance);
    end
    
    % 存储结果
    coherence_analysis.phase_smoothness = phase_smoothness;
    coherence_analysis.amplitude_stability = amplitude_stability;
    coherence_analysis.coherence_length = coherence_length;
    coherence_analysis.amplitude_variation = amplitude_variation;
    
    % 综合相干性评分
    overall_coherence = (phase_smoothness + amplitude_stability) / 2;
    coherence_analysis.overall_coherence = overall_coherence;
    
    % 判断标准
    params = parameters();
    coherence_analysis.coherence_threshold = params.multi_beam_coherence_threshold;
    coherence_analysis.meets_criterion = overall_coherence > params.multi_beam_coherence_threshold;
    
    fprintf('\n3. 相干性分析:\n');
    fprintf('   相位平滑度: %.3f\n', phase_smoothness);
    fprintf('   振幅稳定性: %.3f\n', amplitude_stability);
    fprintf('   整体相干性: %.3f\n', overall_coherence);
    fprintf('   满足条件: %s\n', char(string(coherence_analysis.meets_criterion)));
    
end

function [reflection_analysis] = analyze_multiple_reflections(wavenumber, reflectance, angle, material_params)
% 分析多次反射强度
    
    reflection_analysis = struct();
    
    % 计算理论反射系数
    n1 = 1;  % 空气
    n2 = material_params.n_substrate;
    
    % 菲涅尔反射系数
    angle_rad = deg2rad(angle);
    cos_theta1 = cos(angle_rad);
    sin_theta1 = sin(angle_rad);
    sin_theta2 = n1 / n2 * sin_theta1;
    
    if sin_theta2 <= 1  % 无全反射
        cos_theta2 = sqrt(1 - sin_theta2^2);
        
        % s偏振
        r_s = (n1 * cos_theta1 - n2 * cos_theta2) / (n1 * cos_theta1 + n2 * cos_theta2);
        % p偏振
        r_p = (n2 * cos_theta1 - n1 * cos_theta2) / (n2 * cos_theta1 + n1 * cos_theta2);
        
        % 平均反射系数
        r_avg = (abs(r_s)^2 + abs(r_p)^2) / 2;
    else
        r_avg = 1;  % 全反射
    end
    
    % 多次反射强度比估算
    R1 = r_avg;  % 第一次反射
    R2 = R1 * (1 - R1)^2 * R1;  % 第二次反射（简化）
    
    multiple_reflection_ratio = R2 / R1;
    
    % 实验数据中的多次反射特征
    % 通过反射率的高频振荡来估计
    R_smooth = smooth(reflectance, max(5, floor(length(reflectance)/20)));
    high_freq_component = reflectance - R_smooth;
    high_freq_amplitude = std(high_freq_component);
    
    reflection_analysis.theoretical_primary = R1;
    reflection_analysis.theoretical_secondary = R2;
    reflection_analysis.multiple_reflection_ratio = multiple_reflection_ratio;
    reflection_analysis.high_freq_amplitude = high_freq_amplitude;
    reflection_analysis.mean_reflectance = mean(reflectance);
    
    % 多次反射显著性
    reflection_significance = high_freq_amplitude / mean(reflectance);
    reflection_analysis.reflection_significance = reflection_significance;
    
    % 判断标准
    params = parameters();
    reflection_analysis.significance_threshold = params.multi_beam_reflection_threshold;
    reflection_analysis.meets_criterion = reflection_significance > params.multi_beam_reflection_threshold;
    
    fprintf('\n4. 多次反射分析:\n');
    fprintf('   理论一次反射率: %.3f\n', R1);
    fprintf('   多次反射比: %.4f\n', multiple_reflection_ratio);
    fprintf('   反射显著性: %.4f\n', reflection_significance);
    fprintf('   满足条件: %s\n', char(string(reflection_analysis.meets_criterion)));
    
end

function [thickness_analysis] = analyze_optical_thickness_conditions(wavenumber, reflectance, angle, material_params)
% 分析光学厚度条件
    
    thickness_analysis = struct();
    
    % 如果提供了厚度信息，进行详细分析
    if isfield(material_params, 'thickness') && ~isempty(material_params.thickness)
        thickness = material_params.thickness;  % μm
        n_epi = material_params.n_epilayer;
        
        % 计算光学厚度
        angle_rad = deg2rad(angle);
        cos_theta = cos(angle_rad);
        optical_thickness = 2 * n_epi * thickness * cos_theta;  % μm
        
        % 计算波长范围
        wavelength_range = 1e4 ./ wavenumber;  % μm (从cm^-1转换)
        
        % 光学厚度与波长比
        thickness_to_wavelength = optical_thickness ./ wavelength_range;
        
        % 多光束干涉条件：光学厚度 >> 波长
        min_ratio = min(thickness_to_wavelength);
        max_ratio = max(thickness_to_wavelength);
        mean_ratio = mean(thickness_to_wavelength);
        
        thickness_analysis.thickness = thickness;
        thickness_analysis.optical_thickness = optical_thickness;
        thickness_analysis.thickness_to_wavelength_ratio = [min_ratio, max_ratio, mean_ratio];
        
        % 判断标准：光学厚度应该是波长的几倍以上
        params = parameters();
        thickness_analysis.ratio_threshold = params.multi_beam_thickness_ratio;
        thickness_analysis.meets_criterion = mean_ratio > params.multi_beam_thickness_ratio;
        
        fprintf('\n5. 光学厚度分析:\n');
        fprintf('   物理厚度: %.3f μm\n', thickness);
        fprintf('   光学厚度: %.3f μm\n', optical_thickness);
        fprintf('   厚度/波长比: %.2f - %.2f (平均 %.2f)\n', min_ratio, max_ratio, mean_ratio);
        fprintf('   满足条件: %s\n', char(string(thickness_analysis.meets_criterion)));
        
    else
        % 无厚度信息，通过光谱特征估算
        [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakDistance', 5);
        
        if length(peak_locs) >= 2
            % 通过干涉条纹间距估算厚度
            peak_spacings = diff(wavenumber(peak_locs));
            mean_spacing = mean(peak_spacings);
            
            % 估算厚度（简化公式）
            n_est = material_params.n_substrate;  % 使用衬底折射率作为估计
            estimated_thickness = 1 / (2 * n_est * mean_spacing * 1e-4);  % μm
            
            thickness_analysis.estimated_thickness = estimated_thickness;
            thickness_analysis.peak_spacing = mean_spacing;
            thickness_analysis.meets_criterion = estimated_thickness > 1;  % 简单判断
            
            fprintf('\n5. 光学厚度分析（估算）:\n');
            fprintf('   估算厚度: %.3f μm\n', estimated_thickness);
            fprintf('   干涉条纹间距: %.2f cm⁻¹\n', mean_spacing);
        else
            thickness_analysis.meets_criterion = false;
            fprintf('\n5. 光学厚度分析：数据不足\n');
        end
    end
    
end

function [is_multi_beam, overall_assessment] = make_overall_assessment(...
    modulation_analysis, periodicity_analysis, coherence_analysis, ...
    reflection_analysis, thickness_analysis, params)
% 综合评估
    
    overall_assessment = struct();
    
    % 收集各项判断结果
    criteria_met = [
        modulation_analysis.meets_criterion, ...
        periodicity_analysis.meets_criterion, ...
        coherence_analysis.meets_criterion, ...
        reflection_analysis.meets_criterion, ...
        thickness_analysis.meets_criterion
    ];
    
    criteria_names = {
        '调制深度', '周期性', '相干性', '多次反射', '光学厚度'
    };
    
    num_criteria_met = sum(criteria_met);
    total_criteria = length(criteria_met);
    
    % 计算综合评分
    weights = [0.25, 0.25, 0.2, 0.15, 0.15];  % 各项权重
    
    scores = [
        modulation_analysis.global_depth / params.multi_beam_modulation_threshold, ...
        periodicity_analysis.overall_regularity / params.multi_beam_regularity_threshold, ...
        coherence_analysis.overall_coherence / params.multi_beam_coherence_threshold, ...
        reflection_analysis.reflection_significance / params.multi_beam_reflection_threshold, ...
        double(thickness_analysis.meets_criterion)
    ];
    
    % 限制评分在合理范围内
    scores = min(scores, 2);  % 最高2倍阈值
    
    weighted_score = sum(scores .* weights);
    
    overall_assessment.criteria_met = criteria_met;
    overall_assessment.criteria_names = criteria_names;
    overall_assessment.num_criteria_met = num_criteria_met;
    overall_assessment.total_criteria = total_criteria;
    overall_assessment.individual_scores = scores;
    overall_assessment.weighted_score = weighted_score;
    
    % 判断逻辑：至少满足3个主要条件，且综合评分 > 1
    is_multi_beam = (num_criteria_met >= 3) && (weighted_score > 1.0);
    
    % 置信度评估
    if weighted_score >= 1.5
        confidence_level = '高';
    elseif weighted_score >= 1.0
        confidence_level = '中';
    else
        confidence_level = '低';
    end
    
    overall_assessment.is_multi_beam = is_multi_beam;
    overall_assessment.confidence_level = confidence_level;
    
end

function generate_analysis_report(analysis_results)
% 生成分析报告
    
    fprintf('\n=== 多光束干涉条件分析报告 ===\n');
    fprintf('分析时间: %s\n', char(analysis_results.timestamp));
    
    assessment = analysis_results.overall_assessment;
    
    fprintf('\n综合评估结果:\n');
    fprintf('  满足条件数: %d/%d\n', assessment.num_criteria_met, assessment.total_criteria);
    fprintf('  综合评分: %.3f\n', assessment.weighted_score);
    fprintf('  置信水平: %s\n', assessment.confidence_level);
    
    fprintf('\n各项条件详情:\n');
    for i = 1:length(assessment.criteria_names)
        status = '✓';
        if ~assessment.criteria_met(i)
            status = '✗';
        end
        fprintf('  %s %s: %.3f\n', status, assessment.criteria_names{i}, assessment.individual_scores(i));
    end
    
    fprintf('\n最终判断: ');
    if analysis_results.is_multi_beam
        fprintf('满足多光束干涉条件\n');
    else
        fprintf('不满足多光束干涉条件\n');
    end
    
    % 建议
    fprintf('\n建议:\n');
    if analysis_results.is_multi_beam
        fprintf('  - 可以使用多光束干涉模型进行分析\n');
        fprintf('  - 建议考虑高阶反射项的影响\n');
    else
        fprintf('  - 建议使用双光束干涉模型\n');
        fprintf('  - 检查测量条件和样品质量\n');
        
        % 具体改进建议
        if ~assessment.criteria_met(1)
            fprintf('  - 提高测量精度以增加调制深度\n');
        end
        if ~assessment.criteria_met(2)
            fprintf('  - 检查样品表面质量和厚度均匀性\n');
        end
        if ~assessment.criteria_met(3)
            fprintf('  - 改善光源相干性或减少环境干扰\n');
        end
    end
    
end