function [processed_data, processing_results] = sic_data_processor(data_files, angles)
% SIC_DATA_PROCESSOR - SiC外延层数据处理器
%
% 专门处理SiC外延层的红外光谱数据
% 包括数据读取、预处理、厚度计算和结果分析
%
% 输入参数:
%   data_files - 数据文件路径cell数组 {'attachment1.xlsx', 'attachment2.xlsx'}
%   angles - 对应的入射角度数组 [10, 15] (度)
%
% 输出参数:
%   processed_data - 处理后的数据结构体数组
%   processing_results - 处理结果汇总
%
% 处理流程:
%   1. 读取Excel数据文件
%   2. 数据预处理和质量检查
%   3. 厚度计算（多种方法）
%   4. 可靠性分析
%   5. 结果对比和验证
%


    % 输入参数验证
    if nargin < 2
        error('sic_data_processor:InvalidInput', '需要提供数据文件和角度信息');
    end
    
    if length(data_files) ~= length(angles)
        error('sic_data_processor:SizeMismatch', '数据文件数量与角度数量不匹配');
    end
    
    % 加载常数和参数
    const = constants();
    params = parameters();
    
    % 初始化输出结构
    processed_data = [];
    processing_results = struct();
    processing_results.timestamp = datetime('now');
    processing_results.num_datasets = length(data_files);
    processing_results.angles = angles;
    processing_results.material = 'SiC';
    processing_results.refractive_index = const.n_sic;
    
    fprintf('\n=== SiC外延层数据处理开始 ===\n');
    fprintf('处理数据集数量: %d\n', length(data_files));
    fprintf('入射角度: ');
    fprintf('%.1f° ', angles);
    fprintf('\n\n');
    
    % 逐个处理数据文件
    thickness_results = [];
    reliability_scores = [];
    
    for i = 1:length(data_files)
        fprintf('处理数据集 %d/%d: %s (%.1f°)\n', i, length(data_files), data_files{i}, angles(i));
        
        try
            % 读取数据
            [wavenumber, reflectance] = read_spectrum_data(data_files{i});
            
            % 创建测量数据结构
            measurement_data = struct();
            measurement_data.filename = data_files{i};
            measurement_data.wavenumber = wavenumber;
            measurement_data.reflectance = reflectance;
            measurement_data.angle = angles(i);
            measurement_data.n_epi = const.n_sic;
            measurement_data.material = 'SiC';
            
            % 数据预处理
            [processed_measurement] = preprocess_sic_data(measurement_data);
            
            % 厚度计算
            [thickness, thickness_detail] = thickness_algorithm(...
                processed_measurement.wavenumber, ...
                processed_measurement.reflectance, ...
                angles(i), ...
                const.n_sic);
            
            % 可靠性分析
            [reliability_score, reliability_detail] = reliability_analysis(...
                thickness_detail, processed_measurement);
            
            % 存储结果
            dataset_result = struct();
            dataset_result.index = i;
            dataset_result.filename = data_files{i};
            dataset_result.angle = angles(i);
            dataset_result.raw_data = measurement_data;
            dataset_result.processed_data = processed_measurement;
            dataset_result.thickness = thickness;
            dataset_result.thickness_detail = thickness_detail;
            dataset_result.reliability_score = reliability_score;
            dataset_result.reliability_detail = reliability_detail;
            dataset_result.processing_success = true;
            
            processed_data = [processed_data, dataset_result];
            thickness_results = [thickness_results, thickness];
            reliability_scores = [reliability_scores, reliability_score];
            
            fprintf('  厚度: %.3f μm, 可靠性: %.1f/100\n', thickness, reliability_score);
            
        catch ME
            fprintf('  处理失败: %s\n', ME.message);
            
            % 记录失败信息
            dataset_result = struct();
            dataset_result.index = i;
            dataset_result.filename = data_files{i};
            dataset_result.angle = angles(i);
            dataset_result.processing_success = false;
            dataset_result.error_message = ME.message;
            
            processed_data = [processed_data, dataset_result];
        end
    end
    
    % 汇总分析结果
    processing_results.thickness_results = thickness_results;
    processing_results.reliability_scores = reliability_scores;
    
    if ~isempty(thickness_results)
        processing_results.mean_thickness = mean(thickness_results);
        processing_results.std_thickness = std(thickness_results);
        processing_results.thickness_range = [min(thickness_results), max(thickness_results)];
        processing_results.mean_reliability = mean(reliability_scores);
        
        % 角度依赖性分析
        if length(thickness_results) >= 2
            [angle_analysis] = analyze_angle_dependence(angles, thickness_results, reliability_scores);
            processing_results.angle_analysis = angle_analysis;
        end
        
        % 最终推荐厚度
        [final_thickness, confidence] = determine_final_thickness(thickness_results, reliability_scores);
        processing_results.final_thickness = final_thickness;
        processing_results.confidence_level = confidence;
        
    else
        processing_results.mean_thickness = NaN;
        processing_results.std_thickness = NaN;
        processing_results.final_thickness = NaN;
    end
    
    % 生成处理报告
    generate_processing_report(processing_results, processed_data);
    
    fprintf('\n=== SiC外延层数据处理完成 ===\n');
    
end

function [wavenumber, reflectance] = read_spectrum_data(filename)
% 读取光谱数据文件
    
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
        
        % 数据验证
        if length(wavenumber) < 10
            error('数据点数太少');
        end
        
    catch ME
        error('read_spectrum_data:ReadError', '读取文件失败: %s', ME.message);
    end
    
end

function [processed_data] = preprocess_sic_data(measurement_data)
% SiC数据预处理
    
    processed_data = measurement_data;
    
    wavenumber = measurement_data.wavenumber;
    reflectance = measurement_data.reflectance;
    
    % 1. 数据清理
    valid_idx = isfinite(wavenumber) & isfinite(reflectance) & ...
                wavenumber > 0 & reflectance >= 0;
    
    wavenumber_clean = wavenumber(valid_idx);
    reflectance_clean = reflectance(valid_idx);
    
    % 2. 数据排序
    [wavenumber_clean, sort_idx] = sort(wavenumber_clean);
    reflectance_clean = reflectance_clean(sort_idx);
    
    % 3. 异常值检测和处理
    [reflectance_clean] = remove_outliers(reflectance_clean);
    
    % 4. 数据平滑（可选）
    params = parameters();
    if params.smoothing_window > 1 && length(reflectance_clean) > params.smoothing_window
        reflectance_clean = smooth(reflectance_clean, params.smoothing_window);
    end
    
    % 5. 波数范围限制（SiC特定）
    const = constants();
    valid_range_idx = wavenumber_clean >= const.wavenumber_range(1) & ...
                      wavenumber_clean <= const.wavenumber_range(2);
    
    processed_data.wavenumber = wavenumber_clean(valid_range_idx);
    processed_data.reflectance = reflectance_clean(valid_range_idx);
    
    % 记录预处理信息
    processed_data.preprocessing_info.original_points = length(wavenumber);
    processed_data.preprocessing_info.valid_points = length(processed_data.wavenumber);
    processed_data.preprocessing_info.data_loss_ratio = 1 - length(processed_data.wavenumber) / length(wavenumber);
    
end

function [data_clean] = remove_outliers(data)
% 移除异常值
    
    % 使用四分位数方法检测异常值
    Q1 = quantile(data, 0.25);
    Q3 = quantile(data, 0.75);
    IQR = Q3 - Q1;
    
    % 定义异常值范围
    lower_bound = Q1 - 1.5 * IQR;
    upper_bound = Q3 + 1.5 * IQR;
    
    % 移除异常值
    outlier_idx = data < lower_bound | data > upper_bound;
    
    if sum(outlier_idx) > 0
        fprintf('    检测到 %d 个异常值，已移除\n', sum(outlier_idx));
        data_clean = data(~outlier_idx);
    else
        data_clean = data;
    end
    
end

function [angle_analysis] = analyze_angle_dependence(angles, thicknesses, reliabilities)
% 分析角度依赖性
    
    angle_analysis = struct();
    angle_analysis.angles = angles;
    angle_analysis.thicknesses = thicknesses;
    angle_analysis.reliabilities = reliabilities;
    
    % 计算厚度随角度的变化
    if length(angles) >= 2
        % 线性拟合
        p = polyfit(angles, thicknesses, 1);
        angle_analysis.thickness_slope = p(1);  % μm/degree
        angle_analysis.thickness_intercept = p(2);
        
        % 计算相关系数
        corr_matrix = corrcoef(angles, thicknesses);
        angle_analysis.correlation = corr_matrix(1, 2);
        
        % 角度一致性评估
        thickness_variation = std(thicknesses) / mean(thicknesses) * 100;
        angle_analysis.thickness_variation = thickness_variation;
        
        if thickness_variation < 2
            angle_analysis.consistency_grade = '优秀';
        elseif thickness_variation < 5
            angle_analysis.consistency_grade = '良好';
        else
            angle_analysis.consistency_grade = '一般';
        end
    end
    
end

function [final_thickness, confidence] = determine_final_thickness(thicknesses, reliabilities)
% 确定最终推荐厚度
    
    if isempty(thicknesses)
        final_thickness = NaN;
        confidence = 0;
        return;
    end
    
    % 加权平均（根据可靠性评分）
    weights = reliabilities / sum(reliabilities);
    final_thickness = sum(thicknesses .* weights);
    
    % 计算置信度
    thickness_std = std(thicknesses);
    mean_reliability = mean(reliabilities);
    
    % 综合置信度评估
    consistency_factor = max(0, 100 - thickness_std / mean(thicknesses) * 100);
    confidence = (mean_reliability + consistency_factor) / 2;
    
end

function generate_processing_report(results, processed_data)
% 生成处理报告
    
    fprintf('\n=== SiC数据处理报告 ===\n');
    fprintf('处理时间: %s\n', char(results.timestamp));
    fprintf('材料类型: %s\n', results.material);
    fprintf('折射率: %.3f\n', results.refractive_index);
    
    fprintf('\n各数据集结果:\n');
    for i = 1:length(processed_data)
        data = processed_data(i);
        if data.processing_success
            fprintf('  数据集%d (%.1f°): 厚度=%.3f μm, 可靠性=%.1f/100\n', ...
                i, data.angle, data.thickness, data.reliability_score);
        else
            fprintf('  数据集%d (%.1f°): 处理失败\n', i, data.angle);
        end
    end
    
    if ~isnan(results.mean_thickness)
        fprintf('\n统计结果:\n');
        fprintf('  平均厚度: %.3f ± %.3f μm\n', results.mean_thickness, results.std_thickness);
        fprintf('  厚度范围: [%.3f, %.3f] μm\n', results.thickness_range(1), results.thickness_range(2));
        fprintf('  平均可靠性: %.1f/100\n', results.mean_reliability);
        
        if isfield(results, 'final_thickness')
            fprintf('\n最终推荐:\n');
            fprintf('  推荐厚度: %.3f μm\n', results.final_thickness);
            fprintf('  置信水平: %.1f/100\n', results.confidence_level);
        end
        
        if isfield(results, 'angle_analysis')
            fprintf('\n角度依赖性分析:\n');
            fprintf('  厚度变异系数: %.2f%%\n', results.angle_analysis.thickness_variation);
            fprintf('  一致性等级: %s\n', results.angle_analysis.consistency_grade);
        end
    end
    
end