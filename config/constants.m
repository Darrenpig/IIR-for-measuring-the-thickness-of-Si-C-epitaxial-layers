function constants = constants()
% CONSTANTS - 定义红外干涉法测量中使用的物理常数
% 
% 输出参数:
%   constants - 包含所有物理常数的结构体
%


    % 基本物理常数
    constants.c = 2.998e8;              % 光速 (m/s)
    constants.h = 6.626e-34;            % 普朗克常数 (J·s)
    constants.pi = pi;                  % 圆周率
    
    % 材料光学参数
    constants.n_air = 1.0;              % 空气折射率
    constants.n_sic = 2.55;             % SiC外延层折射率 (红外波段)
    constants.n_si = 3.42;              % Si衬底折射率 (红外波段)
    
    % 测量参数
    constants.wavelength_range = [8, 12]; % 红外波长范围 (μm)
    constants.wavenumber_range = [833, 1250]; % 波数范围 (cm^-1)
    
    % 角度转换
    constants.deg2rad = pi/180;         % 度转弧度
    constants.rad2deg = 180/pi;         % 弧度转度
    
    % 单位转换
    constants.um2m = 1e-6;              % 微米转米
    constants.nm2m = 1e-9;              % 纳米转米
    constants.cm2m = 1e-2;              % 厘米转米
    
    % 数值计算参数
    constants.tolerance = 1e-10;        % 数值计算容差
    constants.max_iterations = 1000;    % 最大迭代次数
    
end