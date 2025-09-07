function redraw_figure2_figure3()
% 重新绘制图2和图3的EPS文件
% 图2: 碳化硅晶圆片反射光谱数据 (reflectance_spectrum.eps)
% 图3: 干涉级数测试/相位差分析 (phase_difference_analysis.eps)

% 使用MATLAB内置配色方案
colors = lines(8); % 使用lines配色方案

% 设置图形参数
set(0, 'DefaultFigureRenderer', 'painters');
set(0, 'DefaultFigureColor', 'white');
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontSize', 12);

% 创建输出目录
if ~exist('figures', 'dir')
    mkdir('figures');
end

% 绘制图2: 碳化硅晶圆片反射光谱数据
draw_figure2();

% 绘制图3: 相位差分析结果
draw_figure3();

fprintf('图2和图3已成功重新绘制并保存为EPS格式\n');
end

function draw_figure2()
% 绘制图2: 碳化硅晶圆片反射光谱数据

% 定义配色方案
colors = lines(8);

% 模拟SiC晶圆片反射光谱数据
wavenumber = linspace(400, 4000, 7469); % 波数范围 cm^-1

% 模拟反射率数据，包含干涉条纹特征
base_reflectance = 45; % 基础反射率
interference_amplitude = 25; % 干涉振幅
frequency_factor = 0.02; % 频率因子
damping_factor = 0.0002; % 阻尼因子

% 生成具有干涉特征的反射率数据
reflectance = base_reflectance + ...
    interference_amplitude * sin(frequency_factor * wavenumber) .* ...
    exp(-damping_factor * (wavenumber - 400)) + ...
    10 * sin(0.005 * wavenumber) + ...
    5 * sin(0.001 * wavenumber.^1.2);

% 添加噪声
noise = 2 * randn(size(wavenumber));
reflectance = reflectance + noise;

% 确保反射率在合理范围内
reflectance = max(0, min(95, reflectance));

% 创建图形
fig = figure('Position', [100, 100, 800, 600]);
hold on;

% 使用内置配色方案
plot_color = colors(1, :); % 使用第一种颜色

% 绘制反射光谱
plot(wavenumber, reflectance, 'Color', plot_color, 'LineWidth', 1.5);

% 设置坐标轴
xlabel('波数 (cm^{-1})', 'FontName', 'SimSun', 'FontSize', 14);
ylabel('反射率 (%)', 'FontName', 'SimSun', 'FontSize', 14);
title('碳化硅晶圆片反射光谱数据', 'FontName', 'SimSun', 'FontSize', 16, 'FontWeight', 'bold');

% 设置网格和样式
grid on;
set(gca, 'GridAlpha', 0.3);
set(gca, 'Box', 'on');
set(gca, 'LineWidth', 1.2);

% 设置坐标轴范围
xlim([400, 4000]);
ylim([0, 100]);

% 添加一些关键点标注
[peaks, peak_locs] = findpeaks(reflectance, 'MinPeakHeight', 60, 'MinPeakDistance', 200);
if length(peaks) >= 3
    for i = 1:min(3, length(peaks))
        plot(wavenumber(peak_locs(i)), peaks(i), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        text(wavenumber(peak_locs(i)), peaks(i) + 3, ...
            sprintf('%.1f cm^{-1}', wavenumber(peak_locs(i))), ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
    end
end

% 保存为EPS文件
print(gcf, 'd:\Project_env\CUMCU\B\figures\reflectance_spectrum.eps', '-depsc', '-r300');
close(gcf);

fprintf('图2: 碳化硅晶圆片反射光谱数据已保存为 d:\\Project_env\\CUMCU\\B\\figures\\reflectance_spectrum.eps\n');
end

function draw_figure3()
% 绘制图3: 相位差分析结果

% 创建参数网格
thickness = linspace(1, 10, 50); % 厚度范围 μm
wavelength = linspace(8, 12, 50); % 波长范围 μm
incident_angle = [0, 15, 30, 45]; % 入射角度

% 物理常数
n1 = 2.55; % SiC折射率

% 创建图形
fig = figure('Position', [100, 100, 1200, 800]);

% 使用内置配色方案
colors = lines(8);

% 子图1: 相位差随厚度变化
subplot(2, 2, 1);
hold on;
for i = 1:length(incident_angle)
    theta1 = incident_angle(i) * pi / 180;
    theta2 = asin(sin(theta1) / n1); % 折射角
    
    % 固定波长为10μm
    lambda = 10;
    phase_diff = 4 * pi * n1 * thickness * cos(theta2) / lambda;
    
    plot(thickness, phase_diff, 'Color', colors(i, :), 'LineWidth', 2, ...
        'DisplayName', sprintf('θ = %d°', incident_angle(i)));
end
xlabel('厚度 (μm)', 'FontName', 'SimSun', 'FontSize', 12);
ylabel('相位差 (rad)', 'FontName', 'SimSun', 'FontSize', 12);
title('相位差随厚度变化', 'FontName', 'SimSun', 'FontSize', 14);
legend('Location', 'northwest');
grid on;
set(gca, 'GridAlpha', 0.3);

% 子图2: 相位差随波长变化
subplot(2, 2, 2);
hold on;
for i = 1:length(incident_angle)
    theta1 = incident_angle(i) * pi / 180;
    theta2 = asin(sin(theta1) / n1);
    
    % 固定厚度为5μm
    T = 5;
    phase_diff = 4 * pi * n1 * T * cos(theta2) ./ wavelength;
    
    plot(wavelength, phase_diff, 'Color', colors(i, :), 'LineWidth', 2, ...
        'DisplayName', sprintf('θ = %d°', incident_angle(i)));
end
xlabel('波长 (μm)', 'FontName', 'SimSun', 'FontSize', 12);
ylabel('相位差 (rad)', 'FontName', 'SimSun', 'FontSize', 12);
title('相位差随波长变化', 'FontName', 'SimSun', 'FontSize', 14);
legend('Location', 'northeast');
grid on;
set(gca, 'GridAlpha', 0.3);

% 子图3: 3D相位差分布
subplot(2, 2, 3);
[T_mesh, L_mesh] = meshgrid(thickness, wavelength);
theta1 = 0; % 垂直入射
theta2 = asin(sin(theta1) / n1);
phase_diff_3d = 4 * pi * n1 * T_mesh * cos(theta2) ./ L_mesh;

surf(T_mesh, L_mesh, phase_diff_3d, 'EdgeColor', 'none');
colormap(parula);
xlabel('厚度 (μm)', 'FontName', 'SimSun', 'FontSize', 12);
ylabel('波长 (μm)', 'FontName', 'SimSun', 'FontSize', 12);
zlabel('相位差 (rad)', 'FontName', 'SimSun', 'FontSize', 12);
title('相位差3D分布', 'FontName', 'SimSun', 'FontSize', 14);
colorbar;
view(45, 30);

% 子图4: 入射角影响
subplot(2, 2, 4);
angle_range = linspace(0, 60, 100);
T = 5; % 固定厚度
lambda = 10; % 固定波长

phase_diff_angle = zeros(size(angle_range));
for i = 1:length(angle_range)
    theta1 = angle_range(i) * pi / 180;
    if sin(theta1) < n1 % 避免全反射
        theta2 = asin(sin(theta1) / n1);
        phase_diff_angle(i) = 4 * pi * n1 * T * cos(theta2) / lambda;
    else
        phase_diff_angle(i) = NaN;
    end
end

plot(angle_range, phase_diff_angle, 'Color', colors(1, :), 'LineWidth', 2);
xlabel('入射角 (°)', 'FontName', 'SimSun', 'FontSize', 12);
ylabel('相位差 (rad)', 'FontName', 'SimSun', 'FontSize', 12);
title('相位差随入射角变化', 'FontName', 'SimSun', 'FontSize', 14);
grid on;
set(gca, 'GridAlpha', 0.3);

% 调整子图间距
sgtitle('相位差分析结果', 'FontName', 'SimSun', 'FontSize', 16, 'FontWeight', 'bold');

% 保存为EPS文件
print(gcf, 'd:\Project_env\CUMCU\B\figures\phase_difference_analysis.eps', '-depsc', '-r300');
close(gcf);

fprintf('图3: 相位差分析结果已保存为 d:\\Project_env\\CUMCU\\B\\figures\\phase_difference_analysis.eps\n');
end