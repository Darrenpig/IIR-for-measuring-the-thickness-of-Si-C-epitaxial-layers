function params = parameters()
% PARAMETERS - 定义红外干涉法测量系统的参数配置
% 
% 输出参数:
%   params - 包含所有系统参数的结构体
%


    % 测量角度设置
    params.incident_angles = [10, 15];  % 入射角度 (度)
    
    % 数据文件配置
    params.data_files = {
        'data/附件1.xlsx',  % SiC 10度数据
'data/附件2.xlsx',  % SiC 15度数据
'data/附件3.xlsx',  % Si 10度数据
'data/附件4.xlsx'   % Si 15度数据
    };
    
    % 数据处理参数
    params.enable_smoothing = false;           % 是否启用数据平滑
    params.smoothing_window = 5;               % 平滑窗口大小
    params.enable_baseline_correction = false; % 是否启用基线校正
    params.peak_detection_threshold = 0.01; % 峰值检测阈值
    params.noise_filter_cutoff = 0.1;   % 噪声滤波截止频率
    
    % 多光束干涉参数
    params.multi_beam_modulation_threshold = 0.15;  % 多光束干涉调制深度阈值
    params.multi_beam_regularity_threshold = 0.8;   % 规律性阈值
    params.multi_beam_coherence_threshold = 0.8;    % 相干性阈值
    params.multi_beam_reflection_threshold = 0.3;   % 多次反射阈值
    params.finesse_threshold = 2.0;                 % 精细度阈值
    params.coherence_threshold = 0.8;               % 相干性阈值
    
    % 厚度计算参数
    params.thickness_range = [1, 100];  % 厚度搜索范围 (μm)
    params.thickness_step = 0.01;       % 厚度计算步长 (μm)
    
    % 多光束干涉判断参数
    params.multi_beam_threshold = 0.05; % 多光束干涉判断阈值
    params.interference_order_max = 20; % 最大干涉级数
    
    % 可靠性分析参数
    params.confidence_level = 0.95;     % 置信水平
    params.error_tolerance = 0.02;      % 误差容限
    
    % 绘图参数
    params.plot_line_width = 1.5;      % 绘图线宽
    params.plot_font_size = 12;        % 绘图字体大小
    params.plot_colors = {'b', 'r', 'g', 'k', 'm', 'c'}; % 绘图颜色
    
    % 输出格式参数
    params.output_precision = 4;        % 输出精度（小数位数）
    params.save_format = 'mat';         % 保存格式
    
end