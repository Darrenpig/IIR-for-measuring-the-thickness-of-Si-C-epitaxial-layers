function [spectrum_data, metadata] = data_loader(file_path, angle, material_type)
% DATA_LOADER - 通用数据加载器
%
% 加载和预处理光谱数据文件，支持多种格式
% 提供统一的数据接口和质量检查
%
% 输入参数:
%   file_path - 数据文件路径
%   angle - 入射角度 (度)
%   material_type - 材料类型 ('SiC' 或 'Si')
%
% 输出参数:
%   spectrum_data - 标准化的光谱数据结构体
%     .wavenumber - 波数数组 (cm^-1)
%     .reflectance - 反射率数组 (%)
%     .wavelength - 波长数组 (μm)
%     .frequency - 频率数组 (Hz)
%     .angle - 入射角 (度)
%     .material - 材料类型
%     .num_points - 数据点数
%   metadata - 元数据信息
%


    % 输入参数验证
    if nargin < 1
        error('data_loader:InvalidInput', '必须提供文件路径');
    end
    
    if nargin < 2
        angle = 0; % 默认垂直入射
    end
    
    if nargin < 3
        material_type = 'Unknown';
    end
    
    % 检查文件是否存在
    if ~exist(file_path, 'file')
        error('data_loader:FileNotFound', '文件不存在: %s', file_path);
    end
    
    fprintf('\n=== 数据加载器 ===\n');
    fprintf('文件路径: %s\n', file_path);
    fprintf('入射角度: %.1f°\n', angle);
    fprintf('材料类型: %s\n', material_type);
    
    % 获取文件信息
    [~, filename, ext] = fileparts(file_path);
    file_info = dir(file_path);
    
    % 初始化元数据
    metadata = struct();
    metadata.filename = filename;
    metadata.extension = ext;
    metadata.file_size = file_info.bytes;
    metadata.load_time = datetime('now');
    metadata.angle = angle;
    metadata.material = material_type;
    
    % 根据文件扩展名选择加载方法
    switch lower(ext)
        case {'.xlsx', '.xls'}
            [raw_data, metadata] = load_excel_data(file_path, metadata);
        case '.csv'
            [raw_data, metadata] = load_csv_data(file_path, metadata);
        case '.txt'
            [raw_data, metadata] = load_text_data(file_path, metadata);
        case '.mat'
            [raw_data, metadata] = load_mat_data(file_path, metadata);
        otherwise
            error('data_loader:UnsupportedFormat', '不支持的文件格式: %s', ext);
    end
    
    % 数据质量检查
    [raw_data, quality_info] = data_quality_check(raw_data);
    metadata.quality = quality_info;
    
    % 标准化数据结构
    spectrum_data = standardize_spectrum_data(raw_data, angle, material_type, metadata);
    
    % 数据预处理
    spectrum_data = preprocess_spectrum_data(spectrum_data);
    
    % 输出加载结果
    fprintf('数据加载完成:\n');
    fprintf('  数据点数: %d\n', spectrum_data.num_points);
    fprintf('  波数范围: %.1f - %.1f cm^-1\n', min(spectrum_data.wavenumber), max(spectrum_data.wavenumber));
    fprintf('  反射率范围: %.2f - %.2f%%\n', min(spectrum_data.reflectance), max(spectrum_data.reflectance));
    fprintf('  数据质量: %s\n', quality_info.overall_quality);
    
end

%% Excel数据加载
function [data, metadata] = load_excel_data(file_path, metadata)
    fprintf('加载Excel文件...\n');
    
    try
        % 读取数值数据
        [num_data, txt_data, raw_data] = xlsread(file_path);
        
        % 检查数据格式
        if size(num_data, 2) < 2
            error('data_loader:InvalidFormat', 'Excel文件至少需要两列数据（波数和反射率）');
        end
        
        % 提取波数和反射率
        data.wavenumber = num_data(:, 1);
        data.reflectance = num_data(:, 2);
        
        % 移除NaN值
        valid_idx = ~isnan(data.wavenumber) & ~isnan(data.reflectance);
        data.wavenumber = data.wavenumber(valid_idx);
        data.reflectance = data.reflectance(valid_idx);
        
        % 更新元数据
        metadata.original_points = size(num_data, 1);
        metadata.valid_points = length(data.wavenumber);
        metadata.header_info = txt_data(1, :); % 表头信息
        
        fprintf('Excel数据加载成功，有效数据点: %d\n', metadata.valid_points);
        
    catch ME
        error('data_loader:ExcelError', 'Excel文件读取失败: %s', ME.message);
    end
end

%% CSV数据加载
function [data, metadata] = load_csv_data(file_path, metadata)
    fprintf('加载CSV文件...\n');
    
    try
        csv_data = readmatrix(file_path);
        
        if size(csv_data, 2) < 2
            error('data_loader:InvalidFormat', 'CSV文件至少需要两列数据');
        end
        
        data.wavenumber = csv_data(:, 1);
        data.reflectance = csv_data(:, 2);
        
        % 移除NaN值
        valid_idx = ~isnan(data.wavenumber) & ~isnan(data.reflectance);
        data.wavenumber = data.wavenumber(valid_idx);
        data.reflectance = data.reflectance(valid_idx);
        
        metadata.original_points = size(csv_data, 1);
        metadata.valid_points = length(data.wavenumber);
        
        fprintf('CSV数据加载成功\n');
        
    catch ME
        error('data_loader:CSVError', 'CSV文件读取失败: %s', ME.message);
    end
end

%% 文本数据加载
function [data, metadata] = load_text_data(file_path, metadata)
    fprintf('加载文本文件...\n');
    
    try
        txt_data = readmatrix(file_path);
        
        if size(txt_data, 2) < 2
            error('data_loader:InvalidFormat', '文本文件至少需要两列数据');
        end
        
        data.wavenumber = txt_data(:, 1);
        data.reflectance = txt_data(:, 2);
        
        % 移除NaN值
        valid_idx = ~isnan(data.wavenumber) & ~isnan(data.reflectance);
        data.wavenumber = data.wavenumber(valid_idx);
        data.reflectance = data.reflectance(valid_idx);
        
        metadata.original_points = size(txt_data, 1);
        metadata.valid_points = length(data.wavenumber);
        
        fprintf('文本数据加载成功\n');
        
    catch ME
        error('data_loader:TextError', '文本文件读取失败: %s', ME.message);
    end
end

%% MAT数据加载
function [data, metadata] = load_mat_data(file_path, metadata)
    fprintf('加载MAT文件...\n');
    
    try
        mat_data = load(file_path);
        field_names = fieldnames(mat_data);
        
        % 查找波数和反射率字段
        wavenumber_fields = {'wavenumber', 'wave_number', 'k', 'frequency'};
        reflectance_fields = {'reflectance', 'reflection', 'R', 'intensity'};
        
        wavenumber_field = '';
        reflectance_field = '';
        
        for i = 1:length(wavenumber_fields)
            if ismember(wavenumber_fields{i}, field_names)
                wavenumber_field = wavenumber_fields{i};
                break;
            end
        end
        
        for i = 1:length(reflectance_fields)
            if ismember(reflectance_fields{i}, field_names)
                reflectance_field = reflectance_fields{i};
                break;
            end
        end
        
        if isempty(wavenumber_field) || isempty(reflectance_field)
            error('data_loader:InvalidFormat', 'MAT文件中找不到波数或反射率数据');
        end
        
        data.wavenumber = mat_data.(wavenumber_field);
        data.reflectance = mat_data.(reflectance_field);
        
        metadata.original_points = length(data.wavenumber);
        metadata.valid_points = length(data.wavenumber);
        metadata.mat_fields = field_names;
        
        fprintf('MAT数据加载成功\n');
        
    catch ME
        error('data_loader:MATError', 'MAT文件读取失败: %s', ME.message);
    end
end

%% 数据质量检查
function [data, quality_info] = data_quality_check(data)
    fprintf('执行数据质量检查...\n');
    
    quality_info = struct();
    quality_info.check_time = datetime('now');
    
    % 检查数据点数
    n_points = length(data.wavenumber);
    quality_info.num_points = n_points;
    
    if n_points < 10
        quality_info.point_count_status = 'Poor';
        warning('data_loader:InsufficientData', '数据点数过少: %d', n_points);
    elseif n_points < 50
        quality_info.point_count_status = 'Fair';
    else
        quality_info.point_count_status = 'Good';
    end
    
    % 检查数据范围
    wavenumber_range = [min(data.wavenumber), max(data.wavenumber)];
    reflectance_range = [min(data.reflectance), max(data.reflectance)];
    
    quality_info.wavenumber_range = wavenumber_range;
    quality_info.reflectance_range = reflectance_range;
    
    % 检查波数单调性
    is_monotonic = all(diff(data.wavenumber) > 0) || all(diff(data.wavenumber) < 0);
    quality_info.is_monotonic = is_monotonic;
    
    if ~is_monotonic
        fprintf('警告: 波数数据不单调，正在排序...\n');
        [data.wavenumber, sort_idx] = sort(data.wavenumber);
        data.reflectance = data.reflectance(sort_idx);
    end
    
    % 检查异常值
    reflectance_mean = mean(data.reflectance);
    reflectance_std = std(data.reflectance);
    outlier_threshold = 3 * reflectance_std;
    
    outlier_idx = abs(data.reflectance - reflectance_mean) > outlier_threshold;
    num_outliers = sum(outlier_idx);
    quality_info.num_outliers = num_outliers;
    quality_info.outlier_percentage = num_outliers / n_points * 100;
    
    if num_outliers > 0
        fprintf('发现 %d 个异常值 (%.1f%%)\n', num_outliers, quality_info.outlier_percentage);
    end
    
    % 检查数据间隔
    wavenumber_intervals = diff(data.wavenumber);
    interval_mean = mean(wavenumber_intervals);
    interval_std = std(wavenumber_intervals);
    quality_info.interval_uniformity = interval_std / interval_mean;
    
    % 综合质量评估
    quality_score = 0;
    
    if strcmp(quality_info.point_count_status, 'Good')
        quality_score = quality_score + 3;
    elseif strcmp(quality_info.point_count_status, 'Fair')
        quality_score = quality_score + 2;
    else
        quality_score = quality_score + 1;
    end
    
    if is_monotonic
        quality_score = quality_score + 2;
    end
    
    if quality_info.outlier_percentage < 5
        quality_score = quality_score + 2;
    elseif quality_info.outlier_percentage < 10
        quality_score = quality_score + 1;
    end
    
    if quality_info.interval_uniformity < 0.1
        quality_score = quality_score + 1;
    end
    
    % 质量等级
    if quality_score >= 7
        quality_info.overall_quality = 'Excellent';
    elseif quality_score >= 5
        quality_info.overall_quality = 'Good';
    elseif quality_score >= 3
        quality_info.overall_quality = 'Fair';
    else
        quality_info.overall_quality = 'Poor';
    end
    
    quality_info.quality_score = quality_score;
    
    fprintf('数据质量检查完成，质量等级: %s\n', quality_info.overall_quality);
end

%% 标准化光谱数据结构
function spectrum_data = standardize_spectrum_data(raw_data, angle, material_type, metadata)
    % 创建标准化的光谱数据结构
    
    spectrum_data = struct();
    
    % 基本数据
    spectrum_data.wavenumber = raw_data.wavenumber; % cm^-1
    spectrum_data.reflectance = raw_data.reflectance; % %
    spectrum_data.angle = angle; % 度
    spectrum_data.material = material_type;
    spectrum_data.num_points = length(raw_data.wavenumber);
    
    % 计算派生量
    const = constants();
    
    % 波长 (μm)
    spectrum_data.wavelength = 1 ./ spectrum_data.wavenumber * 1e4;
    
    % 频率 (Hz)
    spectrum_data.frequency = const.c * spectrum_data.wavenumber * 100; % 转换为Hz
    
    % 能量 (eV)
    spectrum_data.photon_energy = const.h * spectrum_data.frequency / const.eV;
    
    % 统计信息
    spectrum_data.stats = struct();
    spectrum_data.stats.wavenumber_range = [min(spectrum_data.wavenumber), max(spectrum_data.wavenumber)];
    spectrum_data.stats.wavelength_range = [min(spectrum_data.wavelength), max(spectrum_data.wavelength)];
    spectrum_data.stats.reflectance_mean = mean(spectrum_data.reflectance);
    spectrum_data.stats.reflectance_std = std(spectrum_data.reflectance);
    spectrum_data.stats.reflectance_range = [min(spectrum_data.reflectance), max(spectrum_data.reflectance)];
    
    % 元数据
    spectrum_data.metadata = metadata;
    spectrum_data.creation_time = datetime('now');
end

%% 数据预处理
function spectrum_data = preprocess_spectrum_data(spectrum_data)
    % 可选的数据预处理步骤
    
    % 平滑处理（如果需要）
    params = parameters();
    
    if params.enable_smoothing
        fprintf('应用数据平滑...\n');
        spectrum_data.reflectance_raw = spectrum_data.reflectance; % 保存原始数据
        spectrum_data.reflectance = smooth(spectrum_data.reflectance, params.smoothing_window);
        spectrum_data.is_smoothed = true;
    else
        spectrum_data.is_smoothed = false;
    end
    
    % 基线校正（如果需要）
    if params.enable_baseline_correction
        fprintf('应用基线校正...\n');
        if ~isfield(spectrum_data, 'reflectance_raw')
            spectrum_data.reflectance_raw = spectrum_data.reflectance;
        end
        
        % 简单的线性基线校正
        baseline = linspace(spectrum_data.reflectance(1), spectrum_data.reflectance(end), length(spectrum_data.reflectance));
        spectrum_data.reflectance = spectrum_data.reflectance - baseline' + mean(spectrum_data.reflectance);
        spectrum_data.baseline_corrected = true;
    else
        spectrum_data.baseline_corrected = false;
    end
end