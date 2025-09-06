#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
data_process.py - 数据预处理与极值点提取程序

功能描述:
    对红外干涉光谱数据进行预处理，包括数据清洗、滤波、极值点提取等
    为后续厚度计算提供高质量的特征数据

主要功能:
    1. Excel数据读取与格式转换
    2. 数据质量检查与异常值处理
    3. 信号滤波与噪声去除
    4. 干涉条纹极值点识别
    5. 特征参数提取

作者: CUMCU数学建模团队
日期: 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import find_peaks, savgol_filter
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings('ignore')

class DataProcessor:
    """数据预处理类"""
    
    def __init__(self):
        self.raw_data = None
        self.processed_data = None
        self.extrema_points = None
        self.quality_metrics = {}
        
    def load_excel_data(self, file_path, sheet_name=0):
        """
        加载Excel数据文件
        
        参数:
            file_path: Excel文件路径
            sheet_name: 工作表名称或索引
            
        返回:
            success: 是否成功加载
        """
        try:
            print(f"正在加载数据文件: {file_path}")
            
            # 读取Excel文件
            df = pd.read_excel(file_path, sheet_name=sheet_name)
            
            # 假设第一列为波数，第二列为反射率
            if df.shape[1] >= 2:
                wavenumber = df.iloc[:, 0].values  # cm^-1
                reflectance = df.iloc[:, 1].values  # %
                
                # 数据验证
                if len(wavenumber) == 0 or len(reflectance) == 0:
                    raise ValueError("数据为空")
                    
                self.raw_data = {
                    'wavenumber': wavenumber,
                    'reflectance': reflectance,
                    'wavelength': 10000.0 / wavenumber  # 转换为波长 (μm)
                }
                
                print(f"数据加载成功: {len(wavenumber)} 个数据点")
                print(f"波数范围: {wavenumber.min():.1f} - {wavenumber.max():.1f} cm^-1")
                print(f"反射率范围: {reflectance.min():.2f} - {reflectance.max():.2f} %")
                
                return True
            else:
                raise ValueError("数据格式不正确，需要至少两列数据")
                
        except Exception as e:
            print(f"数据加载失败: {str(e)}")
            return False
    
    def quality_check(self, outlier_threshold=3.0):
        """
        数据质量检查
        
        参数:
            outlier_threshold: 异常值检测阈值（标准差倍数）
            
        返回:
            quality_report: 质量检查报告
        """
        if self.raw_data is None:
            return None
            
        print("\n=== 数据质量检查 ===")
        
        wavenumber = self.raw_data['wavenumber']
        reflectance = self.raw_data['reflectance']
        
        # 1. 检查缺失值
        missing_count = np.sum(np.isnan(wavenumber)) + np.sum(np.isnan(reflectance))
        
        # 2. 检查异常值
        z_scores = np.abs((reflectance - np.mean(reflectance)) / np.std(reflectance))
        outliers = np.sum(z_scores > outlier_threshold)
        
        # 3. 检查数据连续性
        wavenumber_diff = np.diff(wavenumber)
        avg_step = np.mean(wavenumber_diff)
        irregular_steps = np.sum(np.abs(wavenumber_diff - avg_step) > 2 * np.std(wavenumber_diff))
        
        # 4. 信噪比估计
        # 使用高频成分估计噪声水平
        high_freq = signal.detrend(reflectance)
        noise_level = np.std(high_freq)
        signal_level = np.std(reflectance)
        snr = 20 * np.log10(signal_level / noise_level) if noise_level > 0 else float('inf')
        
        self.quality_metrics = {
            'total_points': len(wavenumber),
            'missing_values': missing_count,
            'outliers': outliers,
            'outlier_ratio': outliers / len(reflectance) * 100,
            'irregular_steps': irregular_steps,
            'snr_db': snr,
            'noise_level': noise_level
        }
        
        # 质量等级评估
        if outliers / len(reflectance) < 0.05 and snr > 20:
            quality_grade = "Excellent"
        elif outliers / len(reflectance) < 0.10 and snr > 15:
            quality_grade = "Good"
        elif outliers / len(reflectance) < 0.20 and snr > 10:
            quality_grade = "Fair"
        else:
            quality_grade = "Poor"
            
        self.quality_metrics['quality_grade'] = quality_grade
        
        print(f"数据点数: {self.quality_metrics['total_points']}")
        print(f"缺失值: {missing_count}")
        print(f"异常值: {outliers} ({self.quality_metrics['outlier_ratio']:.1f}%)")
        print(f"信噪比: {snr:.1f} dB")
        print(f"质量等级: {quality_grade}")
        
        return self.quality_metrics
    
    def preprocess_data(self, filter_window=11, filter_order=3, remove_outliers=True):
        """
        数据预处理
        
        参数:
            filter_window: 滤波窗口大小
            filter_order: 滤波器阶数
            remove_outliers: 是否移除异常值
            
        返回:
            success: 预处理是否成功
        """
        if self.raw_data is None:
            return False
            
        print("\n=== 数据预处理 ===")
        
        wavenumber = self.raw_data['wavenumber'].copy()
        reflectance = self.raw_data['reflectance'].copy()
        
        # 1. 异常值处理
        if remove_outliers:
            z_scores = np.abs((reflectance - np.mean(reflectance)) / np.std(reflectance))
            valid_mask = z_scores <= 3.0
            wavenumber = wavenumber[valid_mask]
            reflectance = reflectance[valid_mask]
            print(f"移除异常值: {np.sum(~valid_mask)} 个点")
        
        # 2. 数据插值（确保等间距）
        if len(wavenumber) > 100:  # 确保有足够的数据点
            wavenumber_interp = np.linspace(wavenumber.min(), wavenumber.max(), len(wavenumber))
            f_interp = interp1d(wavenumber, reflectance, kind='cubic', fill_value='extrapolate')
            reflectance_interp = f_interp(wavenumber_interp)
            wavenumber = wavenumber_interp
            reflectance = reflectance_interp
            print("数据插值完成")
        
        # 3. Savitzky-Golay滤波
        if len(reflectance) > filter_window:
            reflectance_filtered = savgol_filter(reflectance, filter_window, filter_order)
            print(f"Savitzky-Golay滤波完成 (窗口={filter_window}, 阶数={filter_order})")
        else:
            reflectance_filtered = reflectance
            print("数据点不足，跳过滤波")
        
        # 4. 基线校正
        # 使用多项式拟合去除基线漂移
        baseline_poly = np.polyfit(wavenumber, reflectance_filtered, deg=3)
        baseline = np.polyval(baseline_poly, wavenumber)
        reflectance_corrected = reflectance_filtered - baseline + np.mean(reflectance_filtered)
        
        self.processed_data = {
            'wavenumber': wavenumber,
            'reflectance': reflectance_corrected,
            'reflectance_raw': reflectance,
            'reflectance_filtered': reflectance_filtered,
            'baseline': baseline,
            'wavelength': 10000.0 / wavenumber
        }
        
        print("数据预处理完成")
        return True
    
    def find_extrema(self, prominence=0.5, distance=10):
        """
        寻找干涉条纹的极值点
        
        参数:
            prominence: 峰值显著性阈值
            distance: 相邻峰值最小距离
            
        返回:
            extrema_info: 极值点信息
        """
        if self.processed_data is None:
            return None
            
        print("\n=== 极值点提取 ===")
        
        reflectance = self.processed_data['reflectance']
        wavenumber = self.processed_data['wavenumber']
        
        # 寻找极大值点（峰值）
        peaks, peak_properties = find_peaks(reflectance, 
                                          prominence=prominence, 
                                          distance=distance)
        
        # 寻找极小值点（谷值）
        valleys, valley_properties = find_peaks(-reflectance, 
                                              prominence=prominence, 
                                              distance=distance)
        
        # 合并极值点并排序
        all_extrema = np.concatenate([peaks, valleys])
        extrema_types = np.concatenate([np.ones(len(peaks)), -np.ones(len(valleys))])
        
        # 按波数排序
        sort_indices = np.argsort(wavenumber[all_extrema])
        all_extrema = all_extrema[sort_indices]
        extrema_types = extrema_types[sort_indices]
        
        self.extrema_points = {
            'indices': all_extrema,
            'wavenumbers': wavenumber[all_extrema],
            'wavelengths': 10000.0 / wavenumber[all_extrema],
            'reflectances': reflectance[all_extrema],
            'types': extrema_types,  # 1为峰值，-1为谷值
            'peak_indices': peaks,
            'valley_indices': valleys
        }
        
        print(f"找到极值点: {len(all_extrema)} 个 (峰值: {len(peaks)}, 谷值: {len(valleys)})")
        print(f"波数范围: {wavenumber[all_extrema].min():.1f} - {wavenumber[all_extrema].max():.1f} cm^-1")
        
        return self.extrema_points
    
    def calculate_interference_order(self):
        """
        计算干涉级数
        
        返回:
            interference_orders: 干涉级数数组
        """
        if self.extrema_points is None:
            return None
            
        wavelengths = self.extrema_points['wavelengths']
        
        # 选择参考波长（通常选择第一个极值点）
        if len(wavelengths) < 2:
            return None
            
        lambda_ref = wavelengths[0]
        
        # 计算干涉级数
        # 对于相邻极值点，相位差为π，对应级数差为0.5
        orders = np.zeros(len(wavelengths))
        
        for i in range(1, len(wavelengths)):
            # 累积级数差
            delta_order = 0.5  # 相邻极值点的级数差
            orders[i] = orders[i-1] + delta_order
            
        return orders
    
    def export_results(self, output_file='processed_data.xlsx'):
        """
        导出处理结果
        
        参数:
            output_file: 输出文件名
        """
        if self.processed_data is None or self.extrema_points is None:
            print("没有可导出的数据")
            return False
            
        try:
            # 创建Excel写入器
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                
                # 导出处理后的完整数据
                processed_df = pd.DataFrame({
                    '波数 (cm^-1)': self.processed_data['wavenumber'],
                    '波长 (μm)': self.processed_data['wavelength'],
                    '原始反射率 (%)': self.processed_data['reflectance_raw'],
                    '滤波反射率 (%)': self.processed_data['reflectance_filtered'],
                    '校正反射率 (%)': self.processed_data['reflectance']
                })
                processed_df.to_excel(writer, sheet_name='处理数据', index=False)
                
                # 导出极值点数据
                extrema_df = pd.DataFrame({
                    '序号': range(1, len(self.extrema_points['wavenumbers']) + 1),
                    '波数 (cm^-1)': self.extrema_points['wavenumbers'],
                    '波长 (μm)': self.extrema_points['wavelengths'],
                    '反射率 (%)': self.extrema_points['reflectances'],
                    '类型': ['峰值' if t > 0 else '谷值' for t in self.extrema_points['types']]
                })
                extrema_df.to_excel(writer, sheet_name='极值点', index=False)
                
                # 导出质量指标
                if self.quality_metrics:
                    quality_df = pd.DataFrame(list(self.quality_metrics.items()), 
                                            columns=['指标', '数值'])
                    quality_df.to_excel(writer, sheet_name='质量指标', index=False)
            
            print(f"结果已导出至: {output_file}")
            return True
            
        except Exception as e:
            print(f"导出失败: {str(e)}")
            return False
    
    def plot_results(self, save_figure=True, figure_name='data_processing_results.png'):
        """
        绘制处理结果
        
        参数:
            save_figure: 是否保存图片
            figure_name: 图片文件名
        """
        if self.processed_data is None:
            return
            
        plt.figure(figsize=(15, 10))
        
        # 子图1: 原始数据 vs 处理数据
        plt.subplot(2, 2, 1)
        plt.plot(self.processed_data['wavenumber'], self.processed_data['reflectance_raw'], 
                'b-', alpha=0.7, label='原始数据')
        plt.plot(self.processed_data['wavenumber'], self.processed_data['reflectance'], 
                'r-', linewidth=2, label='处理数据')
        plt.xlabel('波数 (cm⁻¹)')
        plt.ylabel('反射率 (%)')
        plt.title('数据预处理对比')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 子图2: 极值点标记
        if self.extrema_points is not None:
            plt.subplot(2, 2, 2)
            plt.plot(self.processed_data['wavenumber'], self.processed_data['reflectance'], 
                    'b-', linewidth=1, label='处理数据')
            
            # 标记峰值和谷值
            peak_mask = self.extrema_points['types'] > 0
            valley_mask = self.extrema_points['types'] < 0
            
            plt.plot(self.extrema_points['wavenumbers'][peak_mask], 
                    self.extrema_points['reflectances'][peak_mask], 
                    'ro', markersize=8, label=f'峰值 ({np.sum(peak_mask)}个)')
            plt.plot(self.extrema_points['wavenumbers'][valley_mask], 
                    self.extrema_points['reflectances'][valley_mask], 
                    'go', markersize=8, label=f'谷值 ({np.sum(valley_mask)}个)')
            
            plt.xlabel('波数 (cm⁻¹)')
            plt.ylabel('反射率 (%)')
            plt.title('极值点识别')
            plt.legend()
            plt.grid(True, alpha=0.3)
        
        # 子图3: 波长域显示
        plt.subplot(2, 2, 3)
        plt.plot(self.processed_data['wavelength'], self.processed_data['reflectance'], 
                'b-', linewidth=2)
        plt.xlabel('波长 (μm)')
        plt.ylabel('反射率 (%)')
        plt.title('波长域光谱')
        plt.grid(True, alpha=0.3)
        
        # 子图4: 质量指标
        if self.quality_metrics:
            plt.subplot(2, 2, 4)
            metrics_text = f"""数据质量报告:
            
总数据点: {self.quality_metrics.get('total_points', 'N/A')}
异常值: {self.quality_metrics.get('outliers', 'N/A')} ({self.quality_metrics.get('outlier_ratio', 0):.1f}%)
信噪比: {self.quality_metrics.get('snr_db', 0):.1f} dB
质量等级: {self.quality_metrics.get('quality_grade', 'N/A')}
            """
            plt.text(0.1, 0.5, metrics_text, fontsize=12, verticalalignment='center',
                    transform=plt.gca().transAxes, bbox=dict(boxstyle='round', facecolor='lightgray'))
            plt.axis('off')
            plt.title('数据质量指标')
        
        plt.tight_layout()
        
        if save_figure:
            plt.savefig(figure_name, dpi=300, bbox_inches='tight')
            print(f"图片已保存: {figure_name}")
        
        plt.show()

def main():
    """主函数 - 数据处理流程演示"""
    print("=== 红外干涉光谱数据预处理系统 ===")
    
    # 创建数据处理器
    processor = DataProcessor()
    
    # 示例：处理附件1数据
    file_path = "../附件1.xlsx"  # 相对于code目录的路径
    
    if processor.load_excel_data(file_path):
        # 数据质量检查
        processor.quality_check()
        
        # 数据预处理
        if processor.preprocess_data():
            # 极值点提取
            processor.find_extrema()
            
            # 计算干涉级数
            orders = processor.calculate_interference_order()
            if orders is not None:
                print(f"\n干涉级数计算完成，范围: {orders.min():.1f} - {orders.max():.1f}")
            
            # 导出结果
            processor.export_results('processed_attachment1.xlsx')
            
            # 绘制结果
            processor.plot_results(save_figure=True, figure_name='attachment1_processing.png')
            
            print("\n数据预处理完成！")
        else:
            print("数据预处理失败")
    else:
        print("数据加载失败")

if __name__ == "__main__":
    main()