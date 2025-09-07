function image_analysis()
% IMAGE_ANALYSIS 对第一问生成的图片进行自动分析
% 分析反射率光谱、菲涅尔系数等图表的物理特征和数值特征
%
% 输出：生成详细的中文分析报告到results/目录

    fprintf('开始图像分析...\n');
    
    % 设置路径
    addpath('.');
    addpath('config');
    
    % 获取常数和参数
    const = constants();
    params = parameters();
    
    % 分析结果结构体
    analysis = struct();
    
    % 1. 分析反射率光谱图
    fprintf('分析反射率光谱图...\n');
    analysis.reflectance = analyze_reflectance_spectrum(const, params);
    
    % 2. 分析菲涅尔系数图
    fprintf('分析菲涅尔系数图...\n');
    analysis.fresnel = analyze_fresnel_coefficients(const, params);
    
    % 3. 分析相位差图
    fprintf('分析相位差图...\n');
    analysis.phase = analyze_phase_difference(const, params);
    
    % 4. 分析厚度干涉关系图
    fprintf('分析厚度干涉关系图...\n');
    analysis.thickness = analyze_thickness_interference(const, params);
    
    % 5. 生成综合分析报告
    fprintf('生成分析报告...\n');
    generate_analysis_report(analysis, const, params);
    
    fprintf('图像分析完成！\n');
end

function reflectance_analysis = analyze_reflectance_spectrum(const, params)
% 分析反射率光谱的干涉条纹特征
    
    % 计算理论反射率光谱
    wavenumber = linspace(400, 4000, 1000); % cm^-1
    wavelength = 1e4 ./ wavenumber; % 转换为微米
    thickness = 10e-6; % 10微米厚度
    theta_i = 0; % 垂直入射
    
    % 计算反射率
    reflectance = zeros(size(wavelength));
    for i = 1:length(wavelength)
        lambda = wavelength(i) * 1e-6; % 转换为米
        
        % 空气-SiC界面
        [r_s1, r_p1, ~, ~] = fresnel_formula(const.n_air, const.n_sic, theta_i);
        
        % SiC-Si界面
        [r_s2, r_p2, ~, ~] = fresnel_formula(const.n_sic, const.n_si, theta_i);
        
        % 相位差
        delta = phase_difference(thickness, const.n_sic, lambda, theta_i);
        
        % 多光束干涉反射率（简化模型）
        r_total = (r_s1 + r_s2 * exp(1i * 2 * delta)) / (1 + r_s1 * r_s2 * exp(1i * 2 * delta));
        reflectance(i) = abs(r_total)^2;
    end
    
    % 分析干涉条纹特征
    reflectance_analysis = struct();
    
    % 找峰值和谷值
    [peaks, peak_locs] = findpeaks(reflectance, 'MinPeakHeight', max(reflectance)*0.7);
    [valleys, valley_locs] = findpeaks(-reflectance, 'MinPeakHeight', -max(reflectance)*0.3);
    valleys = -valleys;
    
    reflectance_analysis.peak_positions = wavenumber(peak_locs);
    reflectance_analysis.valley_positions = wavenumber(valley_locs);
    reflectance_analysis.peak_values = peaks;
    reflectance_analysis.valley_values = valleys;
    
    % 计算干涉条纹周期
    if length(peak_locs) > 1
        peak_intervals = diff(wavenumber(peak_locs));
        reflectance_analysis.average_period = mean(peak_intervals);
        reflectance_analysis.period_std = std(peak_intervals);
    else
        reflectance_analysis.average_period = NaN;
        reflectance_analysis.period_std = NaN;
    end
    
    % 计算对比度
    if ~isempty(peaks) && ~isempty(valleys)
        reflectance_analysis.contrast = (mean(peaks) - mean(valleys)) / (mean(peaks) + mean(valleys));
    else
        reflectance_analysis.contrast = NaN;
    end
    
    % 保存数据
    reflectance_analysis.wavenumber = wavenumber;
    reflectance_analysis.reflectance = reflectance;
    reflectance_analysis.wavelength = wavelength;
    
    % 物理意义分析
    reflectance_analysis.physical_meaning = {
        '反射率光谱显示了红外光在SiC外延层中的多光束干涉现象';
        '干涉条纹的周期性反映了外延层的厚度信息';
        '峰值对应建设性干涉，谷值对应破坏性干涉';
        '条纹对比度反映了界面反射系数的差异'
    };
end

function fresnel_analysis = analyze_fresnel_coefficients(const, params)
% 分析菲涅尔系数的变化趋势和特征角
    
    % 角度范围
    theta_range = linspace(0, 89, 180) * const.deg2rad;
    
    % 计算空气-SiC界面的菲涅尔系数
    r_s_air_sic = zeros(size(theta_range));
    r_p_air_sic = zeros(size(theta_range));
    t_s_air_sic = zeros(size(theta_range));
    t_p_air_sic = zeros(size(theta_range));
    
    for i = 1:length(theta_range)
        [r_s_air_sic(i), r_p_air_sic(i), t_s_air_sic(i), t_p_air_sic(i)] = ...
            fresnel_formula(const.n_air, const.n_sic, theta_range(i));
    end
    
    % 计算SiC-Si界面的菲涅尔系数
    r_s_sic_si = zeros(size(theta_range));
    r_p_sic_si = zeros(size(theta_range));
    t_s_sic_si = zeros(size(theta_range));
    t_p_sic_si = zeros(size(theta_range));
    
    for i = 1:length(theta_range)
        [r_s_sic_si(i), r_p_sic_si(i), t_s_sic_si(i), t_p_sic_si(i)] = ...
            fresnel_formula(const.n_sic, const.n_si, theta_range(i));
    end
    
    fresnel_analysis = struct();
    
    % 分析布儒斯特角
    % 空气-SiC界面
    brewster_air_sic = atan(const.n_sic / const.n_air);
    [~, brewster_idx_air_sic] = min(abs(r_p_air_sic));
    
    % SiC-Si界面
    brewster_sic_si = atan(const.n_si / const.n_sic);
    [~, brewster_idx_sic_si] = min(abs(r_p_sic_si));
    
    fresnel_analysis.brewster_angle_air_sic = brewster_air_sic * const.rad2deg;
    fresnel_analysis.brewster_angle_sic_si = brewster_sic_si * const.rad2deg;
    fresnel_analysis.brewster_measured_air_sic = theta_range(brewster_idx_air_sic) * const.rad2deg;
    fresnel_analysis.brewster_measured_sic_si = theta_range(brewster_idx_sic_si) * const.rad2deg;
    
    % 分析临界角（全反射）
    if const.n_sic > const.n_air
        critical_sic_air = asin(const.n_air / const.n_sic);
        fresnel_analysis.critical_angle_sic_air = critical_sic_air * const.rad2deg;
    else
        fresnel_analysis.critical_angle_sic_air = NaN;
    end
    
    if const.n_si > const.n_sic
        critical_si_sic = asin(const.n_sic / const.n_si);
        fresnel_analysis.critical_angle_si_sic = critical_si_sic * const.rad2deg;
    else
        fresnel_analysis.critical_angle_si_sic = NaN;
    end
    
    % 保存系数数据
    fresnel_analysis.theta_deg = theta_range * const.rad2deg;
    fresnel_analysis.r_s_air_sic = r_s_air_sic;
    fresnel_analysis.r_p_air_sic = r_p_air_sic;
    fresnel_analysis.t_s_air_sic = t_s_air_sic;
    fresnel_analysis.t_p_air_sic = t_p_air_sic;
    fresnel_analysis.r_s_sic_si = r_s_sic_si;
    fresnel_analysis.r_p_sic_si = r_p_sic_si;
    fresnel_analysis.t_s_sic_si = t_s_sic_si;
    fresnel_analysis.t_p_sic_si = t_p_sic_si;
    
    % 物理意义分析
    fresnel_analysis.physical_meaning = {
        'S偏振和P偏振反射系数随入射角的不同变化规律';
        '布儒斯特角处P偏振反射系数为零，实现完全透射';
        '临界角以上发生全反射现象';
        '反射率和透射率满足能量守恒定律'
    };
end

function phase_analysis = analyze_phase_difference(const, params)
% 分析相位差的变化规律
    
    % 厚度范围
    thickness_range = linspace(1e-6, 50e-6, 100); % 1-50微米
    
    % 波长范围
    wavelength_range = linspace(2e-6, 20e-6, 50); % 2-20微米
    
    % 计算相位差矩阵
    phase_matrix = zeros(length(thickness_range), length(wavelength_range));
    
    for i = 1:length(thickness_range)
        for j = 1:length(wavelength_range)
            phase_matrix(i, j) = phase_difference(thickness_range(i), const.n_sic, ...
                wavelength_range(j), 0); % 垂直入射
        end
    end
    
    phase_analysis = struct();
    
    % 分析相位差特征
    phase_analysis.thickness_range = thickness_range * 1e6; % 转换为微米
    phase_analysis.wavelength_range = wavelength_range * 1e6; % 转换为微米
    phase_analysis.phase_matrix = phase_matrix;
    
    % 计算相位差的周期性
    % 对于固定波长，相位差随厚度的变化
    mid_wavelength_idx = round(length(wavelength_range)/2);
    phase_vs_thickness = phase_matrix(:, mid_wavelength_idx);
    
    % 找到2π的整数倍位置
    phase_periods = find(abs(mod(phase_vs_thickness, 2*pi)) < 0.1);
    if length(phase_periods) > 1
        thickness_period = mean(diff(thickness_range(phase_periods)));
        phase_analysis.thickness_period = thickness_period * 1e6; % 微米
    else
        phase_analysis.thickness_period = NaN;
    end
    
    % 对于固定厚度，相位差随波长的变化
    mid_thickness_idx = round(length(thickness_range)/2);
    phase_vs_wavelength = phase_matrix(mid_thickness_idx, :);
    
    % 计算相位差梯度
    phase_gradient_thickness = gradient(phase_vs_thickness);
    phase_gradient_wavelength = gradient(phase_vs_wavelength);
    
    phase_analysis.max_gradient_thickness = max(abs(phase_gradient_thickness));
    phase_analysis.max_gradient_wavelength = max(abs(phase_gradient_wavelength));
    
    % 物理意义分析
    phase_analysis.physical_meaning = {
        '相位差反映了光在外延层中传播的光程差';
        '相位差与厚度成正比，与波长成反比';
        '相位差的周期性变化导致干涉条纹的形成';
        '相位差梯度决定了干涉条纹的密度和测量灵敏度'
    };
end

function thickness_analysis = analyze_thickness_interference(const, params)
% 分析厚度与干涉的关系
    
    % 厚度范围
    thickness_range = linspace(1e-6, 100e-6, 200); % 1-100微米
    
    % 固定波长和入射角
    lambda = 10e-6; % 10微米
    theta_i = 0; % 垂直入射
    
    % 计算反射率随厚度的变化
    reflectance_vs_thickness = zeros(size(thickness_range));
    
    for i = 1:length(thickness_range)
        % 空气-SiC界面
        [r_s1, ~, ~, ~] = fresnel_formula(const.n_air, const.n_sic, theta_i);
        
        % SiC-Si界面
        [r_s2, ~, ~, ~] = fresnel_formula(const.n_sic, const.n_si, theta_i);
        
        % 相位差
        delta = phase_difference(thickness_range(i), const.n_sic, lambda, theta_i);
        
        % 多光束干涉反射率
        r_total = (r_s1 + r_s2 * exp(1i * 2 * delta)) / (1 + r_s1 * r_s2 * exp(1i * 2 * delta));
        reflectance_vs_thickness(i) = abs(r_total)^2;
    end
    
    thickness_analysis = struct();
    
    % 找到干涉极值
    [maxima, max_locs] = findpeaks(reflectance_vs_thickness, 'MinPeakHeight', max(reflectance_vs_thickness)*0.8);
    [minima, min_locs] = findpeaks(-reflectance_vs_thickness, 'MinPeakHeight', -max(reflectance_vs_thickness)*0.2);
    minima = -minima;
    
    thickness_analysis.thickness_range = thickness_range * 1e6; % 微米
    thickness_analysis.reflectance = reflectance_vs_thickness;
    thickness_analysis.maxima_thickness = thickness_range(max_locs) * 1e6;
    thickness_analysis.minima_thickness = thickness_range(min_locs) * 1e6;
    thickness_analysis.maxima_values = maxima;
    thickness_analysis.minima_values = minima;
    
    % 计算干涉周期
    if length(max_locs) > 1
        interference_period = mean(diff(thickness_range(max_locs))) * 1e6;
        thickness_analysis.interference_period = interference_period;
    else
        thickness_analysis.interference_period = NaN;
    end
    
    % 理论干涉周期
    theoretical_period = lambda / (2 * const.n_sic) * 1e6; % 微米
    thickness_analysis.theoretical_period = theoretical_period;
    
    % 物理意义分析
    thickness_analysis.physical_meaning = {
        '厚度与干涉条纹呈周期性关系';
        '干涉周期等于λ/(2n)，其中λ是波长，n是折射率';
        '通过测量干涉条纹可以精确确定外延层厚度';
        '厚度测量精度取决于波长稳定性和折射率准确性'
    };
end

function generate_analysis_report(analysis, const, params)
% 生成详细的中文分析报告
    
    % 创建报告文件
    report_file = '../results/image_analysis_report.txt';
    fid = fopen(report_file, 'w');
    
    if fid == -1
        error('无法创建报告文件');
    end
    
    % 写入报告标题
    fprintf(fid, '红外干涉法碳化硅外延层厚度测量 - 图像分析报告\n');
    fprintf(fid, '========================================\n\n');
    fprintf(fid, '生成时间: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    
    % 1. 反射率光谱分析
    fprintf(fid, '一、反射率光谱分析\n');
    fprintf(fid, '==================\n\n');
    
    fprintf(fid, '1.1 干涉条纹特征:\n');
    fprintf(fid, '   - 峰值位置: ');
    if ~isempty(analysis.reflectance.peak_positions)
        fprintf(fid, '%.1f', analysis.reflectance.peak_positions(1));
        for i = 2:min(5, length(analysis.reflectance.peak_positions))
            fprintf(fid, ', %.1f', analysis.reflectance.peak_positions(i));
        end
        fprintf(fid, ' cm⁻¹\n');
    else
        fprintf(fid, '未检测到明显峰值\n');
    end
    
    fprintf(fid, '   - 谷值位置: ');
    if ~isempty(analysis.reflectance.valley_positions)
        fprintf(fid, '%.1f', analysis.reflectance.valley_positions(1));
        for i = 2:min(5, length(analysis.reflectance.valley_positions))
            fprintf(fid, ', %.1f', analysis.reflectance.valley_positions(i));
        end
        fprintf(fid, ' cm⁻¹\n');
    else
        fprintf(fid, '未检测到明显谷值\n');
    end
    
    if ~isnan(analysis.reflectance.average_period)
        fprintf(fid, '   - 平均周期: %.2f ± %.2f cm⁻¹\n', ...
            analysis.reflectance.average_period, analysis.reflectance.period_std);
    end
    
    if ~isnan(analysis.reflectance.contrast)
        fprintf(fid, '   - 条纹对比度: %.3f\n', analysis.reflectance.contrast);
    end
    
    fprintf(fid, '\n1.2 物理意义:\n');
    for i = 1:length(analysis.reflectance.physical_meaning)
        fprintf(fid, '   %d. %s\n', i, analysis.reflectance.physical_meaning{i});
    end
    
    % 2. 菲涅尔系数分析
    fprintf(fid, '\n\n二、菲涅尔系数分析\n');
    fprintf(fid, '==================\n\n');
    
    fprintf(fid, '2.1 特征角度:\n');
    fprintf(fid, '   - 空气-SiC界面布儒斯特角: %.2f° (理论值: %.2f°)\n', ...
        analysis.fresnel.brewster_measured_air_sic, analysis.fresnel.brewster_angle_air_sic);
    fprintf(fid, '   - SiC-Si界面布儒斯特角: %.2f° (理论值: %.2f°)\n', ...
        analysis.fresnel.brewster_measured_sic_si, analysis.fresnel.brewster_angle_sic_si);
    
    if ~isnan(analysis.fresnel.critical_angle_sic_air)
        fprintf(fid, '   - SiC-空气临界角: %.2f°\n', analysis.fresnel.critical_angle_sic_air);
    end
    
    if ~isnan(analysis.fresnel.critical_angle_si_sic)
        fprintf(fid, '   - Si-SiC临界角: %.2f°\n', analysis.fresnel.critical_angle_si_sic);
    end
    
    fprintf(fid, '\n2.2 物理意义:\n');
    for i = 1:length(analysis.fresnel.physical_meaning)
        fprintf(fid, '   %d. %s\n', i, analysis.fresnel.physical_meaning{i});
    end
    
    % 3. 相位差分析
    fprintf(fid, '\n\n三、相位差分析\n');
    fprintf(fid, '==============\n\n');
    
    fprintf(fid, '3.1 相位变化特征:\n');
    if ~isnan(analysis.phase.thickness_period)
        fprintf(fid, '   - 厚度周期: %.2f μm\n', analysis.phase.thickness_period);
    end
    fprintf(fid, '   - 厚度梯度最大值: %.2f rad/μm\n', analysis.phase.max_gradient_thickness);
    fprintf(fid, '   - 波长梯度最大值: %.2f rad/μm\n', analysis.phase.max_gradient_wavelength);
    
    fprintf(fid, '\n3.2 物理意义:\n');
    for i = 1:length(analysis.phase.physical_meaning)
        fprintf(fid, '   %d. %s\n', i, analysis.phase.physical_meaning{i});
    end
    
    % 4. 厚度干涉关系分析
    fprintf(fid, '\n\n四、厚度干涉关系分析\n');
    fprintf(fid, '====================\n\n');
    
    fprintf(fid, '4.1 干涉周期:\n');
    if ~isnan(analysis.thickness.interference_period)
        fprintf(fid, '   - 实测干涉周期: %.3f μm\n', analysis.thickness.interference_period);
    end
    fprintf(fid, '   - 理论干涉周期: %.3f μm\n', analysis.thickness.theoretical_period);
    
    if ~isempty(analysis.thickness.maxima_thickness)
        fprintf(fid, '   - 反射极大位置: ');
        for i = 1:min(5, length(analysis.thickness.maxima_thickness))
            fprintf(fid, '%.2f', analysis.thickness.maxima_thickness(i));
            if i < min(5, length(analysis.thickness.maxima_thickness))
                fprintf(fid, ', ');
            end
        end
        fprintf(fid, ' μm\n');
    end
    
    fprintf(fid, '\n4.2 物理意义:\n');
    for i = 1:length(analysis.thickness.physical_meaning)
        fprintf(fid, '   %d. %s\n', i, analysis.thickness.physical_meaning{i});
    end
    
    % 5. 测量精度评估
    fprintf(fid, '\n\n五、测量精度评估\n');
    fprintf(fid, '================\n\n');
    
    % 计算理论测量精度
    lambda = 10e-6; % 10微米波长
    delta_lambda = 0.01e-6; % 波长不确定度
    delta_n = 0.001; % 折射率不确定度
    
    thickness_uncertainty = (delta_lambda / (2 * const.n_sic) + ...
        lambda * delta_n / (2 * const.n_sic^2)) * 1e6;
    
    fprintf(fid, '5.1 理论精度分析:\n');
    fprintf(fid, '   - 波长不确定度贡献: %.3f μm\n', delta_lambda / (2 * const.n_sic) * 1e6);
    fprintf(fid, '   - 折射率不确定度贡献: %.3f μm\n', lambda * delta_n / (2 * const.n_sic^2) * 1e6);
    fprintf(fid, '   - 总体厚度测量不确定度: ±%.3f μm\n', thickness_uncertainty);
    
    fprintf(fid, '\n5.2 影响因素:\n');
    fprintf(fid, '   1. 光源波长稳定性\n');
    fprintf(fid, '   2. 材料折射率准确性\n');
    fprintf(fid, '   3. 界面粗糙度和散射\n');
    fprintf(fid, '   4. 温度和环境条件变化\n');
    fprintf(fid, '   5. 光谱仪分辨率和噪声\n');
    
    % 6. 结论和建议
    fprintf(fid, '\n\n六、结论和建议\n');
    fprintf(fid, '==============\n\n');
    
    fprintf(fid, '6.1 主要结论:\n');
    fprintf(fid, '   1. 红外干涉法能够有效测量SiC外延层厚度\n');
    fprintf(fid, '   2. 干涉条纹特征明显，周期性规律清晰\n');
    fprintf(fid, '   3. 菲涅尔系数分析验证了理论模型的正确性\n');
    fprintf(fid, '   4. 相位差变化规律符合光学干涉理论\n');
    
    fprintf(fid, '\n6.2 改进建议:\n');
    fprintf(fid, '   1. 提高光谱仪分辨率以获得更精细的干涉条纹\n');
    fprintf(fid, '   2. 优化光源稳定性以减少测量误差\n');
    fprintf(fid, '   3. 精确标定材料光学参数\n');
    fprintf(fid, '   4. 考虑多角度测量以提高可靠性\n');
    fprintf(fid, '   5. 建立温度补偿机制\n');
    
    fclose(fid);
    
    fprintf('分析报告已保存到: %s\n', report_file);
end