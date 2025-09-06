function [data, metadata] = excel_reader(file_path, options)
% EXCEL_READER - Excel文件读取器
%
% 专门用于读取光谱数据的Excel文件
% 支持多种Excel格式和数据结构
%
% 输入参数:
%   file_path - Excel文件路径
%   options - 读取选项结构体（可选）
%     .sheet - 工作表名称或索引（默认为1）
%     .header_row - 表头行号（默认为1）
%     .data_start_row - 数据开始行号（默认为2）
%     .wavenumber_col - 波数列号（默认为1）
%     .reflectance_col - 反射率列号（默认为2）
%
% 输出参数:
%   data - 数据结构体
%     .wavenumber - 波数数组
%     .reflectance - 反射率数组
%   metadata - 元数据信息
%
% 作者: CUMCU数学建模团队
% 日期: 2024

    % 输入参数验证
    if nargin < 1
        error('excel_reader:InvalidInput', '必须提供文件路径');
    end
    
    if nargin < 2
        options = struct();
    end
    
    % 默认选项
    default_options = struct(...
        'sheet', 1, ...
        'header_row', 1, ...
        'data_start_row', 2, ...
        'wavenumber_col', 1, ...
        'reflectance_col', 2 ...
    );
    
    % 合并选项
    field_names = fieldnames(default_options);
    for i = 1:length(field_names)
        if ~isfield(options, field_names{i})
            options = setfield(options, field_names{i}, getfield(default_options, field_names{i}));
        end
    end
    
    % 检查文件是否存在
    if ~exist(file_path, 'file')
        error('excel_reader:FileNotFound', '文件不存在: %s', file_path);
    end
    
    fprintf('读取Excel文件: %s\n', file_path);
    
    try
        % 使用xlsread读取数据
        [num_data, txt_data, raw_data] = xlsread(file_path, options.sheet);
        
        % 检查数据有效性
        if isempty(num_data)
            error('excel_reader:NoData', '文件中没有数值数据');
        end
        
        if size(num_data, 2) < 2
            error('excel_reader:InsufficientColumns', '数据至少需要两列（波数和反射率）');
        end
        
        % 提取波数和反射率数据
        wavenumber = num_data(:, options.wavenumber_col);
        reflectance = num_data(:, options.reflectance_col);
        
        % 移除无效数据
        valid_idx = ~isnan(wavenumber) & ~isnan(reflectance) & isfinite(wavenumber) & isfinite(reflectance);
        wavenumber = wavenumber(valid_idx);
        reflectance = reflectance(valid_idx);
        
        % 构建数据结构
        data = struct();
        data.wavenumber = wavenumber;
        data.reflectance = reflectance;
        
        % 构建元数据
        metadata = struct();
        metadata.file_path = file_path;
        metadata.sheet_name = options.sheet;
        metadata.original_rows = size(num_data, 1);
        metadata.valid_rows = length(wavenumber);
        metadata.read_time = datetime('now');
        
        if ~isempty(txt_data)
            metadata.header_info = txt_data;
        end
        
        fprintf('成功读取 %d 行有效数据\n', length(wavenumber));
        
    catch ME
        % 如果xlsread失败，尝试使用readtable
        fprintf('xlsread失败，尝试使用readtable...\n');
        
        try
            table_data = readtable(file_path, 'Sheet', options.sheet);
            
            if width(table_data) < 2
                error('excel_reader:InsufficientColumns', '数据至少需要两列');
            end
            
            % 提取数据
            wavenumber = table2array(table_data(:, options.wavenumber_col));
            reflectance = table2array(table_data(:, options.reflectance_col));
            
            % 移除无效数据
            valid_idx = ~isnan(wavenumber) & ~isnan(reflectance) & isfinite(wavenumber) & isfinite(reflectance);
            wavenumber = wavenumber(valid_idx);
            reflectance = reflectance(valid_idx);
            
            % 构建数据结构
            data = struct();
            data.wavenumber = wavenumber;
            data.reflectance = reflectance;
            
            % 构建元数据
            metadata = struct();
            metadata.file_path = file_path;
            metadata.sheet_name = options.sheet;
            metadata.original_rows = height(table_data);
            metadata.valid_rows = length(wavenumber);
            metadata.read_time = datetime('now');
            metadata.column_names = table_data.Properties.VariableNames;
            
            fprintf('使用readtable成功读取 %d 行有效数据\n', length(wavenumber));
            
        catch ME2
            error('excel_reader:ReadError', 'Excel文件读取失败: %s', ME2.message);
        end
    end
    
    % 数据质量检查
    if length(wavenumber) < 10
        warning('excel_reader:InsufficientData', '数据点数较少: %d', length(wavenumber));
    end
    
    % 检查数据范围
    wavenumber_range = max(wavenumber) - min(wavenumber);
    if wavenumber_range < 100
        warning('excel_reader:SmallRange', '波数范围较小: %.1f cm^-1', wavenumber_range);
    end
    
    % 添加统计信息到元数据
    metadata.wavenumber_range = [min(wavenumber), max(wavenumber)];
    metadata.reflectance_range = [min(reflectance), max(reflectance)];
    metadata.wavenumber_mean = mean(wavenumber);
    metadata.reflectance_mean = mean(reflectance);
    metadata.data_quality = 'Good';  % 简单的质量评估
    
    if any(reflectance < 0)
        metadata.data_quality = 'Warning: Negative reflectance values';
        warning('excel_reader:NegativeValues', '存在负反射率值');
    end
    
    if any(reflectance > 100)
        metadata.data_quality = 'Warning: Reflectance > 100%';
        warning('excel_reader:HighReflectance', '存在大于100%的反射率值');
    end
    
end