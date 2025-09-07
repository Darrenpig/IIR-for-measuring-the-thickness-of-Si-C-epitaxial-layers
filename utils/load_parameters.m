function params = load_parameters()
% LOAD_PARAMETERS - 加载系统参数配置
% 
% 输出参数:
%   params - 包含所有系统参数的结构体
%


    % 调用配置文件中的参数函数
    params = parameters();
    
end