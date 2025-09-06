function smoothed_data = smooth(data, window_size)
% SMOOTH - 简单的数据平滑函数
%
% 使用移动平均法对数据进行平滑处理
% 替代MATLAB Signal Processing Toolbox中的smooth函数
%
% 输入参数:
%   data - 输入数据向量
%   window_size - 平滑窗口大小（默认为5）
%
% 输出参数:
%   smoothed_data - 平滑后的数据
%
% 作者: CUMCU数学建模团队
% 日期: 2024

    if nargin < 2
        window_size = 5;
    end
    
    % 确保window_size为奇数
    if mod(window_size, 2) == 0
        window_size = window_size + 1;
    end
    
    % 确保数据为列向量
    data = data(:);
    n = length(data);
    
    % 初始化输出
    smoothed_data = zeros(size(data));
    
    % 计算半窗口大小
    half_window = floor(window_size / 2);
    
    % 对每个点进行平滑
    for i = 1:n
        % 确定窗口范围
        start_idx = max(1, i - half_window);
        end_idx = min(n, i + half_window);
        
        % 计算移动平均
        smoothed_data(i) = mean(data(start_idx:end_idx));
    end
    
end