function constants = load_constants()
% LOAD_CONSTANTS - 加载物理常数和系统参数
% 
% 输出参数:
%   constants - 包含所有物理常数的结构体
%


    % 物理常数
    constants.c = 2.998e8;              % 光速 (m/s)
    constants.h = 6.626e-34;            % 普朗克常数 (J·s)
    constants.k_B = 1.381e-23;          % 玻尔兹曼常数 (J/K)
    
    % 材料折射率 (在红外波段的典型值)
    constants.n_SiC = 2.65;             % SiC折射率
    constants.n_Si = 3.42;              % Si折射率
    constants.n_air = 1.0;              % 空气折射率
    
    % 测量参数
    constants.wavelength_range = [8, 14]; % 波长范围 (μm)
    constants.temperature = 300;         % 室温 (K)
    
    % 系统参数
    constants.measurement_precision = 1e-3; % 测量精度
    constants.angle_tolerance = 0.1;     % 角度容差 (度)
    
end