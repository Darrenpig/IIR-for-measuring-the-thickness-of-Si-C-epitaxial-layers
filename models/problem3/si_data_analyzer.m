function [analysis_results, model_recommendations] = si_data_analyzer(data_files, angles)
% SI_DATA_ANALYZER - 硅片数据分析器
%
% 专门分析硅片的红外光谱数据，判断干涉模型类型
% 并提供模型选择和参数优化建议
%
% 输入参数:
%   data_files - 硅片数据文件路径cell数组 {'attachment3.xlsx', 'attachment4.xlsx'}
%   angles - 对应的入射角度数组 [10, 15] (度)
%
% 输出参数:
%   analysis_results - 详细分析结果结构体数组
%   model_recommendations - 模型推荐和参数建议
%
% 分析内容:
%   1. 硅片光谱特征分析
%   2. 干涉模型类型判断
%   3. 多光束vs双光束条件评估
%   4. 模型参数优化建议
%   5. 测量条件评估


    % 输入参数验证
    if nargin < 2
        error('si_data_analyzer:InvalidInput', '需要提供数据文件和角度信息');
    end
    
    if length(data_files) ~= length(angles)
        error('si_data_analyzer:SizeMismatch', '数据文件数量与角度数量不匹配');
    end
    
    % 加载常数和参数
    const = constants();
    params = parameters();
    
    % 初始化输出结构
    analysis_results = [];
    model_recommendations = struct();
    model_recommendations.timestamp = datetime('now');
    model_recommendations.num_datasets = length(data_files);
    model_recommendations.angles = angles;
    model_recommendations.material = 'Si';
    model_recommendations.refractive_index = const.n_si;
    
    fprintf('\n=== 硅片数据分析开始 ===\n');
    fprintf('数据集数量: %d\n', length(data_files));
    fprintf('入射角度: ');
    fprintf('%.1f° ', angles);
    fprintf('\n材料: 硅片 (n = %.3f)\n\n', const.n_si);
    
    % 逐个分析数据文件
    multi_beam_scores = [];
    model_types = {};
    quality_scores = [];
    
    for i = 1:length(data_files)
        fprintf('分析数据集 %d/%d: %s (%.1f°)\n', i, length(data_files), data_files{i}, angles(i));
        
        try
            % 读取和预处理数据
            [wavenumber, reflectance] = read_si_spectrum_data(data_files{i});
            [processed_data] = preprocess_si_data(wavenumber, reflectance, angles(i));
            
            % 硅片特征分析
            [si_features] = analyze_si_spectral_features(processed_data);
            
            % 干涉模型判断
            [model_analysis] = determine_interference_model(processed_data, const.n_si);
            
            % 多光束条件检查
            material_params = struct('n_substrate', const.n_si);
            [is_multi_beam, multi_beam_analysis] = multi_beam_conditions(...
                processed_data.wavenumber, processed_data.reflectance, angles(i), material_params);
            
            % 数据质量评估
            [quality_assessment] = assess_data_quality(processed_data, si_features);
            
            % 模型参数估算
            [parameter_estimates] = estimate_model_parameters(processed_data, model_analysis, const.n_si);
            
            % 存储单个数据集结果
            dataset_result = struct();
            dataset_result.index = i;
            dataset_result.filename = data_files{i};
            dataset_result.angle = angles(i);
            dataset_result.processed_data = processed_data;
            dataset_result.si_features = si_features;
            dataset_result.model_analysis = model_analysis;
            dataset_result.is_multi_beam = is_multi_beam;
            dataset_result.multi_beam_analysis = multi_beam_analysis;
            dataset_result.quality_assessment = quality_assessment;
            dataset_result.parameter_estimates = parameter_estimates;
            dataset_result.analysis_success = true;
            
            analysis_results = [analysis_results, dataset_result];
            
            % 收集统计信息
            multi_beam_scores = [multi_beam_scores, multi_beam_analysis.overall_assessment.weighted_score];
            model_types{end+1} = model_analysis.recommended_model;
            quality_scores = [quality_scores, quality_assessment.overall_score];
            
            fprintf('  模型类型: %s\n', model_analysis.recommended_model);
            fprintf('  多光束评分: %.3f\n', multi_beam_analysis.overall_assessment.weighted_score);
            fprintf('  数据质量: %.1f/100\n', quality_assessment.overall_score);
            
        catch ME
            fprintf('  分析失败: %s\n', ME.message);
            
            % 记录失败信息
            dataset_result = struct();
            dataset_result.index = i;
            dataset_result.filename = data_files{i};
            dataset_result.angle = angles(i);
            dataset_result.analysis_success = false;
            dataset_result.error_message = ME.message;
            
            analysis_results = [analysis_results, dataset_result];
        end
    end
    
    % 综合分析和模型推荐
    if ~isempty(multi_beam_scores)
        [final_recommendations] = generate_model_recommendations(...
            analysis_results, multi_beam_scores, model_types, quality_scores, angles);
        model_recommendations = merge_structs(model_recommendations, final_recommendations);
    end
    
    % 生成分析报告
    generate_si_analysis_report(model_recommendations, analysis_results);
    
    fprintf('\n=== 硅片数据分析完成 ===\n');
    
end

function [wavenumber, reflectance] = read_si_spectrum_data(filename)
% 读取硅片光谱数据
    
    try
        % 使用excel_reader工具函数
        data = excel_reader(filename);
        
        % 假设第一列是波数，第二列是反射率
        if size(data, 2) >= 2
            wavenumber = data(:, 1);
            reflectance = data(:, 2);
        else
            error('数据文件格式不正确');
        end
        
        % 硅片数据特定验证
        if length(wavenumber) < 20
            error('硅片数据点数太少');
        end
        
        % 检查波数范围（硅片通常在中红外区域）
        if max(wavenumber) < 500 || min(wavenumber) > 5000
            warning('波数范围可能不适合硅片分析');
        end
        
    catch ME
        error('read_si_spectrum_data:ReadError', '读取硅片数据失败: %s', ME.message);
    end
    
end

function [processed_data] = preprocess_si_data(wavenumber, reflectance, angle)
% 硅片数据预处理
    
    processed_data = struct();
    processed_data.angle = angle;
    processed_data.material = 'Si';
    
    % 1. 数据清理
    valid_idx = isfinite(wavenumber) & isfinite(reflectance) & ...
                wavenumber > 0 & reflectance >= 0 & reflectance <= 1;
    
    wavenumber_clean = wavenumber(valid_idx);
    reflectance_clean = reflectance(valid_idx);
    
    % 2. 数据排序
    [wavenumber_clean, sort_idx] = sort(wavenumber_clean);
    reflectance_clean = reflectance_clean(sort_idx);
    
    % 3. 硅片特定的波数范围限制
    const = constants();
    si_range_idx = wavenumber_clean >= const.si_wavenumber_range(1) & ...
                   wavenumber_clean <= const.si_wavenumber_range(2);
    
    wavenumber_final = wavenumber_clean(si_range_idx);
    reflectance_final = reflectance_clean(si_range_idx);
    
    % 4. 异常值处理（硅片反射率通常较高且稳定）
    [reflectance_final] = remove_si_outliers(reflectance_final);
    
    % 5. 数据平滑（保留干涉特征）
    params = parameters();
    if params.si_smoothing_window > 1 && length(reflectance_final) > params.si_smoothing_window
        reflectance_final = smooth(reflectance_final, params.si_smoothing_window);
    end
    
    processed_data.wavenumber = wavenumber_final;
    processed_data.reflectance = reflectance_final;
    
    % 记录预处理信息
    processed_data.preprocessing_info.original_points = length(wavenumber);
    processed_data.preprocessing_info.final_points = length(wavenumber_final);
    processed_data.preprocessing_info.data_retention = length(wavenumber_final) / length(wavenumber);
    
end

function [reflectance_clean] = remove_si_outliers(reflectance)
% 移除硅片数据中的异常值
    
    % 硅片反射率通常在0.3-0.7范围内
    median_R = median(reflectance);
    mad_R = mad(reflectance, 1);  % 中位数绝对偏差
    
    % 使用修正的Z-score方法
    modified_z_scores = 0.6745 * (reflectance - median_R) / mad_R;
    outlier_threshold = 3.5;
    
    outlier_idx = abs(modified_z_scores) > outlier_threshold;
    
    if sum(outlier_idx) > 0
        fprintf('    移除 %d 个硅片数据异常值\n', sum(outlier_idx));
        reflectance_clean = reflectance(~outlier_idx);
    else
        reflectance_clean = reflectance;
    end
    
end

function [si_features] = analyze_si_spectral_features(processed_data)
% 分析硅片光谱特征
    
    wavenumber = processed_data.wavenumber;
    reflectance = processed_data.reflectance;
    
    si_features = struct();
    
    % 1. 基本统计特征
    si_features.mean_reflectance = mean(reflectance);
    si_features.std_reflectance = std(reflectance);
    si_features.min_reflectance = min(reflectance);
    si_features.max_reflectance = max(reflectance);
    si_features.reflectance_range = si_features.max_reflectance - si_features.min_reflectance;
    
    % 2. 硅片特征反射率检查
    const = constants();
    expected_si_reflectance = ((const.n_si - 1) / (const.n_si + 1))^2;  % 垂直入射理论值
    si_features.expected_reflectance = expected_si_reflectance;
    si_features.reflectance_deviation = abs(si_features.mean_reflectance - expected_si_reflectance);
    
    % 3. 干涉条纹特征
    [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakDistance', 3);
    [valleys, valley_locs] = findpeaks(-reflectance, 'MinPeakDistance', 3);
    valleys = -valleys;
    
    si_features.num_peaks = length(peaks);
    si_features.num_valleys = length(valleys);
    
    if length(peaks) >= 2
        peak_spacings = diff(wavenumber(peak_locs));
        si_features.mean_peak_spacing = mean(peak_spacings);
        si_features.std_peak_spacing = std(peak_spacings);
        si_features.peak_regularity = 1 - si_features.std_peak_spacing / si_features.mean_peak_spacing;
    else
        si_features.mean_peak_spacing = NaN;
        si_features.peak_regularity = 0;
    end
    
    % 4. 表面质量指标
    % 通过反射率的平滑度评估表面质量
    reflectance_gradient = gradient(reflectance);
    si_features.surface_roughness_indicator = std(reflectance_gradient);
    
    % 5. 干涉可见度
    if si_features.num_peaks >= 1 && si_features.num_valleys >= 1
        visibility = (si_features.max_reflectance - si_features.min_reflectance) / ...
                    (si_features.max_reflectance + si_features.min_reflectance);
        si_features.interference_visibility = visibility;
    else
        si_features.interference_visibility = 0;
    end
    
    % 6. 硅片厚度估算（如果有干涉条纹）
    if ~isnan(si_features.mean_peak_spacing) && si_features.mean_peak_spacing > 0
        % 简化的厚度估算公式
        angle_rad = deg2rad(processed_data.angle);
        cos_theta = cos(angle_rad);
        estimated_thickness = 1 / (2 * const.n_si * si_features.mean_peak_spacing * 1e-4 * cos_theta);  % μm
        si_features.estimated_thickness = estimated_thickness;
    else
        si_features.estimated_thickness = NaN;
    end
    
end

function [model_analysis] = determine_interference_model(processed_data, n_si)
% 确定干涉模型类型
    
    model_analysis = struct();
    
    wavenumber = processed_data.wavenumber;
    reflectance = processed_data.reflectance;
    angle = processed_data.angle;
    
    % 1. 分析干涉条纹特征
    [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakDistance', 5);
    [valleys, valley_locs] = findpeaks(-reflectance, 'MinPeakDistance', 5);
    
    % 2. 计算调制深度
    modulation_depth = (max(reflectance) - min(reflectance)) / (max(reflectance) + min(reflectance));
    
    % 3. 分析反射率水平
    mean_reflectance = mean(reflectance);
    theoretical_si_reflectance = ((n_si - 1) / (n_si + 1))^2;
    
    % 4. 检查多次反射特征
    % 高频振荡可能表明多次反射
    R_smooth = smooth(reflectance, max(5, floor(length(reflectance)/10)));
    high_freq_component = reflectance - R_smooth;
    high_freq_amplitude = std(high_freq_component);
    
    % 5. 模型判断逻辑
    model_scores = struct();
    
    % 双光束干涉模型评分
    two_beam_score = 0;
    if length(peaks) >= 2 && modulation_depth > 0.1
        two_beam_score = two_beam_score + 30;
    end
    if abs(mean_reflectance - theoretical_si_reflectance) < 0.1
        two_beam_score = two_beam_score + 25;
    end
    if high_freq_amplitude / mean_reflectance < 0.02
        two_beam_score = two_beam_score + 25;
    end
    if length(peaks) < 10  % 不太复杂的干涉图样
        two_beam_score = two_beam_score + 20;
    end
    
    % 多光束干涉模型评分
    multi_beam_score = 0;
    if modulation_depth > 0.3
        multi_beam_score = multi_beam_score + 30;
    end
    if length(peaks) >= 5
        multi_beam_score = multi_beam_score + 25;
    end
    if high_freq_amplitude / mean_reflectance > 0.05
        multi_beam_score = multi_beam_score + 25;
    end
    if mean_reflectance > theoretical_si_reflectance * 1.2  % 多次反射增强
        multi_beam_score = multi_beam_score + 20;
    end
    
    % 简单反射模型评分（无明显干涉）
    simple_reflection_score = 0;
    if length(peaks) <= 1
        simple_reflection_score = simple_reflection_score + 40;
    end
    if modulation_depth < 0.05
        simple_reflection_score = simple_reflection_score + 30;
    end
    if abs(mean_reflectance - theoretical_si_reflectance) < 0.05
        simple_reflection_score = simple_reflection_score + 30;
    end
    
    model_scores.two_beam = two_beam_score;
    model_scores.multi_beam = multi_beam_score;
    model_scores.simple_reflection = simple_reflection_score;
    
    % 确定推荐模型
    [max_score, max_idx] = max([two_beam_score, multi_beam_score, simple_reflection_score]);
    model_names = {'双光束干涉', '多光束干涉', '简单反射'};
    recommended_model = model_names{max_idx};
    
    model_analysis.model_scores = model_scores;
    model_analysis.recommended_model = recommended_model;
    model_analysis.confidence = max_score / 100;
    model_analysis.modulation_depth = modulation_depth;
    model_analysis.num_interference_features = length(peaks) + length(valleys);
    model_analysis.high_freq_amplitude = high_freq_amplitude;
    
end

function [quality_assessment] = assess_data_quality(processed_data, si_features)
% 评估数据质量
    
    quality_assessment = struct();
    
    % 1. 数据完整性
    data_completeness = processed_data.preprocessing_info.data_retention * 100;
    
    % 2. 信噪比估算
    signal_power = var(processed_data.reflectance);
    noise_estimate = si_features.surface_roughness_indicator;
    snr_estimate = 10 * log10(signal_power / (noise_estimate^2 + eps));
    
    % 3. 干涉特征清晰度
    interference_clarity = si_features.interference_visibility * 100;
    
    % 4. 数据一致性
    const = constants();
    expected_reflectance = ((const.n_si - 1) / (const.n_si + 1))^2;
    consistency_score = max(0, 100 - abs(si_features.mean_reflectance - expected_reflectance) * 500);
    
    % 5. 测量范围适宜性
    wavenumber_range = max(processed_data.wavenumber) - min(processed_data.wavenumber);
    range_score = min(100, wavenumber_range / 1000 * 100);  % 假设1000 cm^-1为满分
    
    % 综合评分
    weights = [0.2, 0.25, 0.25, 0.15, 0.15];
    individual_scores = [data_completeness, min(snr_estimate*10, 100), interference_clarity, consistency_score, range_score];
    overall_score = sum(individual_scores .* weights);
    
    quality_assessment.data_completeness = data_completeness;
    quality_assessment.snr_estimate = snr_estimate;
    quality_assessment.interference_clarity = interference_clarity;
    quality_assessment.consistency_score = consistency_score;
    quality_assessment.range_score = range_score;
    quality_assessment.individual_scores = individual_scores;
    quality_assessment.overall_score = overall_score;
    
    % 质量等级
    if overall_score >= 80
        quality_grade = '优秀';
    elseif overall_score >= 60
        quality_grade = '良好';
    elseif overall_score >= 40
        quality_grade = '一般';
    else
        quality_grade = '较差';
    end
    
    quality_assessment.quality_grade = quality_grade;
    
end

function [parameter_estimates] = estimate_model_parameters(processed_data, model_analysis, n_si)
% 估算模型参数
    
    parameter_estimates = struct();
    parameter_estimates.model_type = model_analysis.recommended_model;
    
    wavenumber = processed_data.wavenumber;
    reflectance = processed_data.reflectance;
    angle = processed_data.angle;
    
    % 基本光学参数
    parameter_estimates.refractive_index = n_si;
    parameter_estimates.incident_angle = angle;
    
    % 根据模型类型估算参数
    switch model_analysis.recommended_model
        case '双光束干涉'
            % 估算厚度和反射系数
            [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakDistance', 3);
            if length(peak_locs) >= 2
                peak_spacing = mean(diff(wavenumber(peak_locs)));
                angle_rad = deg2rad(angle);
                cos_theta = cos(angle_rad);
                thickness_est = 1 / (2 * n_si * peak_spacing * 1e-4 * cos_theta);  % μm
                parameter_estimates.thickness = thickness_est;
            end
            
            % 反射系数估算
            R_mean = mean(reflectance);
            parameter_estimates.reflection_coefficient = sqrt(R_mean);
            
        case '多光束干涉'
            % 多光束参数估算
            R_max = max(reflectance);
            R_min = min(reflectance);
            
            % 精细度估算
            finesse_est = pi * sqrt(R_max) / (2 * (1 - sqrt(R_max)));
            parameter_estimates.finesse = finesse_est;
            
            % 反射率估算
            parameter_estimates.surface_reflectivity = sqrt(R_max);
            
        case '简单反射'
            % 简单反射参数
            parameter_estimates.surface_reflectance = mean(reflectance);
            parameter_estimates.surface_roughness = std(reflectance);
    end
    
    % 拟合优度评估
    parameter_estimates.model_fit_quality = model_analysis.confidence;
    
end

function [recommendations] = generate_model_recommendations(analysis_results, multi_beam_scores, model_types, quality_scores, angles)
% 生成模型推荐
    
    recommendations = struct();
    
    % 统计模型类型分布
    unique_models = unique(model_types);
    model_counts = zeros(size(unique_models));
    for i = 1:length(unique_models)
        model_counts(i) = sum(strcmp(model_types, unique_models{i}));
    end
    
    [~, dominant_model_idx] = max(model_counts);
    dominant_model = unique_models{dominant_model_idx};
    
    recommendations.dominant_model = dominant_model;
    recommendations.model_distribution = containers.Map(unique_models, num2cell(model_counts));
    
    % 综合评分
    recommendations.mean_multi_beam_score = mean(multi_beam_scores);
    recommendations.mean_quality_score = mean(quality_scores);
    
    % 角度依赖性分析
    if length(angles) >= 2
        angle_analysis = analyze_angle_effects(analysis_results, angles);
        recommendations.angle_analysis = angle_analysis;
    end
    
    % 最终建议
    if strcmp(dominant_model, '多光束干涉') && recommendations.mean_multi_beam_score > 1.2
        recommendations.final_recommendation = '推荐使用多光束干涉模型';
        recommendations.model_parameters = 'high_finesse';
    elseif strcmp(dominant_model, '双光束干涉')
        recommendations.final_recommendation = '推荐使用双光束干涉模型';
        recommendations.model_parameters = 'standard';
    else
        recommendations.final_recommendation = '推荐使用简单反射模型';
        recommendations.model_parameters = 'basic';
    end
    
    % 改进建议
    improvement_suggestions = {};
    if recommendations.mean_quality_score < 60
        improvement_suggestions{end+1} = '提高测量精度和数据质量';
    end
    if recommendations.mean_multi_beam_score < 0.8
        improvement_suggestions{end+1} = '检查样品表面质量';
    end
    if length(unique_models) > 1
        improvement_suggestions{end+1} = '统一测量条件以获得一致的模型类型';
    end
    
    recommendations.improvement_suggestions = improvement_suggestions;
    
end

function [angle_analysis] = analyze_angle_effects(analysis_results, angles)
% 分析角度效应
    
    angle_analysis = struct();
    
    % 提取各角度的关键参数
    multi_beam_scores = [];
    quality_scores = [];
    model_types = {};
    
    for i = 1:length(analysis_results)
        if analysis_results(i).analysis_success
            multi_beam_scores(end+1) = analysis_results(i).multi_beam_analysis.overall_assessment.weighted_score;
            quality_scores(end+1) = analysis_results(i).quality_assessment.overall_score;
            model_types{end+1} = analysis_results(i).model_analysis.recommended_model;
        end
    end
    
    angle_analysis.angles = angles;
    angle_analysis.multi_beam_scores = multi_beam_scores;
    angle_analysis.quality_scores = quality_scores;
    angle_analysis.model_types = model_types;
    
    % 角度依赖性评估
    if length(angles) >= 2 && length(multi_beam_scores) >= 2
        % 计算相关性
        corr_mb = corrcoef(angles, multi_beam_scores);
        corr_quality = corrcoef(angles, quality_scores);
        
        angle_analysis.multi_beam_angle_correlation = corr_mb(1,2);
        angle_analysis.quality_angle_correlation = corr_quality(1,2);
        
        % 最优角度推荐
        [~, best_angle_idx] = max(multi_beam_scores .* quality_scores);
        angle_analysis.recommended_angle = angles(best_angle_idx);
    end
    
end

function [merged] = merge_structs(struct1, struct2)
% 合并结构体
    merged = struct1;
    fields = fieldnames(struct2);
    for i = 1:length(fields)
        merged.(fields{i}) = struct2.(fields{i});
    end
end

function generate_si_analysis_report(recommendations, analysis_results)
% 生成硅片分析报告
    
    fprintf('\n=== 硅片数据分析报告 ===\n');
    fprintf('分析时间: %s\n', char(recommendations.timestamp));
    fprintf('材料: %s (n = %.3f)\n', recommendations.material, recommendations.refractive_index);
    
    fprintf('\n模型推荐结果:\n');
    fprintf('  主导模型: %s\n', recommendations.dominant_model);
    fprintf('  最终建议: %s\n', recommendations.final_recommendation);
    fprintf('  平均多光束评分: %.3f\n', recommendations.mean_multi_beam_score);
    fprintf('  平均数据质量: %.1f/100\n', recommendations.mean_quality_score);
    
    fprintf('\n各数据集详情:\n');
    for i = 1:length(analysis_results)
        result = analysis_results(i);
        if result.analysis_success
            fprintf('  数据集%d (%.1f°): %s, 质量=%.1f/100\n', ...
                i, result.angle, result.model_analysis.recommended_model, ...
                result.quality_assessment.overall_score);
        else
            fprintf('  数据集%d (%.1f°): 分析失败\n', i, result.angle);
        end
    end
    
    if isfield(recommendations, 'angle_analysis')
        fprintf('\n角度效应分析:\n');
        if isfield(recommendations.angle_analysis, 'recommended_angle')
            fprintf('  推荐测量角度: %.1f°\n', recommendations.angle_analysis.recommended_angle);
        end
    end
    
    if ~isempty(recommendations.improvement_suggestions)
        fprintf('\n改进建议:\n');
        for i = 1:length(recommendations.improvement_suggestions)
            fprintf('  - %s\n', recommendations.improvement_suggestions{i});
        end
    end
    
end