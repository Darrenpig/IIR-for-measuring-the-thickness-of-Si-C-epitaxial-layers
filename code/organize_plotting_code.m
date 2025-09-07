function organize_plotting_code()
% ORGANIZE_PLOTTING_CODE - MATLAB绘图代码整理脚本
%
% 该脚本根据整理报告自动执行代码重组：
% 1. 创建新的目录结构
% 2. 移动使用中的代码到core目录
% 3. 归档未使用的代码到archive目录
% 4. 生成整理日志


    fprintf('\n=== MATLAB绘图代码整理工具 ===\n');
    
    % 获取当前目录
    current_dir = pwd;
    fprintf('当前工作目录: %s\n', current_dir);
    
    % 自动执行整理（非交互式）
    fprintf('\n开始自动执行代码整理...\n');
    execute_organization(true);
end

function execute_organization(do_execute)
    % 定义文件分类
    core_files = {
        'plot_problem1_results.m', '问题一结果绘图函数';
        'plot_problem2_results.m', '问题二结果绘图函数';
        'plot_problem3_results.m', '问题三结果绘图函数';
        'plot_problem3_advanced.m', '问题三高级绘图技术';
    };
    
    utils_files = {
        '../utils/plot_utils.m', '绘图工具函数集合';
        '../utils/plot_spectrum.m', '光谱数据绘图函数';
    };
    
    archive_files = {
        'plot_all_with_color51.m', 'color_variants', '使用color51配色重绘所有图表';
        'plot_with_color51.m', 'color_variants', '使用color51配色的基础绘图';
        'redraw_figure2_figure3.m', 'specific_purpose', '重绘特定图表';
        'compare_problem3_plots.m', 'specific_purpose', '问题三绘图对比';
        'image_analysis.m', 'development', '图像分析相关';
    };
    
    if do_execute
        fprintf('\n=== 开始执行代码整理 ===\n');
        
        % 创建目录结构
        create_directories();
        
        % 移动核心文件
        move_core_files(core_files);
        
        % 复制工具文件
        copy_utils_files(utils_files);
        
        % 归档未使用文件
        archive_unused_files(archive_files);
        
        % 生成整理日志
        generate_log(core_files, utils_files, archive_files);
        
        fprintf('\n=== 代码整理完成 ===\n');
        fprintf('请查看生成的日志文件: plotting_organization_log.txt\n');
        
    else
        fprintf('\n=== 预览模式 - 将要执行的操作 ===\n');
        preview_operations(core_files, utils_files, archive_files);
    end
end

function create_directories()
    fprintf('\n1. 创建目录结构...\n');
    
    dirs_to_create = {
        'plotting';
        'plotting/core';
        'plotting/utils';
        'archive';
        'archive/plotting';
        'archive/plotting/color_variants';
        'archive/plotting/specific_purpose';
        'archive/plotting/development';
    };
    
    for i = 1:length(dirs_to_create)
        dir_path = dirs_to_create{i};
        if ~exist(dir_path, 'dir')
            mkdir(dir_path);
            fprintf('  创建目录: %s\n', dir_path);
        else
            fprintf('  目录已存在: %s\n', dir_path);
        end
    end
end

function move_core_files(core_files)
    fprintf('\n2. 移动核心绘图文件...\n');
    
    for i = 1:size(core_files, 1)
        source_file = core_files{i, 1};
        description = core_files{i, 2};
        dest_file = fullfile('plotting', 'core', source_file);
        
        if exist(source_file, 'file')
            try
                movefile(source_file, dest_file);
                fprintf('  移动: %s -> %s (%s)\n', source_file, dest_file, description);
            catch ME
                fprintf('  错误: 无法移动 %s - %s\n', source_file, ME.message);
            end
        else
            fprintf('  警告: 文件不存在 %s\n', source_file);
        end
    end
end

function copy_utils_files(utils_files)
    fprintf('\n3. 复制工具文件...\n');
    
    for i = 1:size(utils_files, 1)
        source_file = utils_files{i, 1};
        description = utils_files{i, 2};
        [~, filename, ext] = fileparts(source_file);
        dest_file = fullfile('plotting', 'utils', [filename, ext]);
        
        if exist(source_file, 'file')
            try
                copyfile(source_file, dest_file);
                fprintf('  复制: %s -> %s (%s)\n', source_file, dest_file, description);
            catch ME
                fprintf('  错误: 无法复制 %s - %s\n', source_file, ME.message);
            end
        else
            fprintf('  警告: 文件不存在 %s\n', source_file);
        end
    end
end

function archive_unused_files(archive_files)
    fprintf('\n4. 归档未使用文件...\n');
    
    for i = 1:size(archive_files, 1)
        source_file = archive_files{i, 1};
        category = archive_files{i, 2};
        description = archive_files{i, 3};
        dest_file = fullfile('archive', 'plotting', category, source_file);
        
        if exist(source_file, 'file')
            try
                movefile(source_file, dest_file);
                fprintf('  归档: %s -> %s (%s)\n', source_file, dest_file, description);
            catch ME
                fprintf('  错误: 无法归档 %s - %s\n', source_file, ME.message);
            end
        else
            fprintf('  警告: 文件不存在 %s\n', source_file);
        end
    end
end

function generate_log(core_files, utils_files, archive_files)
    fprintf('\n5. 生成整理日志...\n');
    
    log_file = 'plotting_organization_log.txt';
    fid = fopen(log_file, 'w');
    
    if fid == -1
        fprintf('  错误: 无法创建日志文件\n');
        return;
    end
    
    fprintf(fid, 'MATLAB绘图代码整理日志\n');
    fprintf(fid, '整理时间: %s\n', datestr(now));
    fprintf(fid, '========================================\n\n');
    
    fprintf(fid, '核心绘图文件 (移动到 plotting/core/):\n');
    for i = 1:size(core_files, 1)
        fprintf(fid, '  - %s: %s\n', core_files{i, 1}, core_files{i, 2});
    end
    
    fprintf(fid, '\n工具文件 (复制到 plotting/utils/):\n');
    for i = 1:size(utils_files, 1)
        fprintf(fid, '  - %s: %s\n', utils_files{i, 1}, utils_files{i, 2});
    end
    
    fprintf(fid, '\n归档文件 (移动到 archive/plotting/):\n');
    for i = 1:size(archive_files, 1)
        fprintf(fid, '  - %s -> %s/: %s\n', archive_files{i, 1}, archive_files{i, 2}, archive_files{i, 3});
    end
    
    fprintf(fid, '\n新目录结构:\n');
    fprintf(fid, 'plotting/\n');
    fprintf(fid, '├── core/\n');
    fprintf(fid, '│   ├── plot_problem1_results.m\n');
    fprintf(fid, '│   ├── plot_problem2_results.m\n');
    fprintf(fid, '│   ├── plot_problem3_results.m\n');
    fprintf(fid, '│   └── plot_problem3_advanced.m\n');
    fprintf(fid, '└── utils/\n');
    fprintf(fid, '    ├── plot_utils.m\n');
    fprintf(fid, '    └── plot_spectrum.m\n');
    fprintf(fid, '\n');
    fprintf(fid, 'archive/plotting/\n');
    fprintf(fid, '├── color_variants/\n');
    fprintf(fid, '├── specific_purpose/\n');
    fprintf(fid, '└── development/\n');
    
    fclose(fid);
    fprintf('  日志文件已生成: %s\n', log_file);
end

function preview_operations(core_files, utils_files, archive_files)
    fprintf('\n将要创建的目录:\n');
    fprintf('  plotting/core/\n');
    fprintf('  plotting/utils/\n');
    fprintf('  archive/plotting/color_variants/\n');
    fprintf('  archive/plotting/specific_purpose/\n');
    fprintf('  archive/plotting/development/\n');
    
    fprintf('\n将要移动的核心文件:\n');
    for i = 1:size(core_files, 1)
        fprintf('  %s -> plotting/core/\n', core_files{i, 1});
    end
    
    fprintf('\n将要复制的工具文件:\n');
    for i = 1:size(utils_files, 1)
        fprintf('  %s -> plotting/utils/\n', utils_files{i, 1});
    end
    
    fprintf('\n将要归档的文件:\n');
    for i = 1:size(archive_files, 1)
        fprintf('  %s -> archive/plotting/%s/\n', archive_files{i, 1}, archive_files{i, 2});
    end
    
    fprintf('\n如需执行，请重新运行并选择选项1。\n');
end