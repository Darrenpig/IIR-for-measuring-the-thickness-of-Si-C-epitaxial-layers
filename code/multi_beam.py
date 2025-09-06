#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
multi_beam.py - 多光束干涉分析与修正程序

功能描述:
    分析多光束干涉现象对厚度测量的影响
    提供多光束干涉条件判断和模型修正算法

主要功能:
    1. 多光束干涉条件判断
    2. 干涉条纹特征分析
    3. 多光束干涉模型建立
    4. 厚度测量修正算法
    5. 影响评估与消除

作者: CUMCU数学建模团队
日期: 2024
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import curve_fit, minimize
from scipy.fft import fft, fftfreq
import warnings
warnings.filterwarnings('ignore')

class MultiBeamAnalyzer:
    """多光束干涉分析器"""
    
    def __init__(self):
        self.data = None
        self.analysis_results = {}
        self.multi_beam_conditions = {}
        self.correction_factors = {}
        
    def load_data(self, wavenumber, reflectance, material_type='Si', incident_angle=10.0):
        """
        加载光谱数据
        
        参数:
            wavenumber: 波数数组 (cm^-1)
            reflectance: 反射率数组 (%)
            material_type: 材料类型 ('Si' 或 'SiC')
            incident_angle: 入射角 (度)
        """
        self.data = {
            'wavenumber': np.array(wavenumber),
            'reflectance': np.array(reflectance),
            'wavelength': 10000.0 / np.array(wavenumber),  # 转换为波长 (μm)
            'material_type': material_type,
            'incident_angle': incident_angle
        }
        
        # 材料参数
        if material_type.upper() == 'SI':
            self.data['refractive_index'] = 3.42
        elif material_type.upper() == 'SIC':
            self.data['refractive_index'] = 2.55
        else:
            self.data['refractive_index'] = 3.42  # 默认硅
            
        print(f"数据加载完成: {len(wavenumber)} 个数据点")
        print(f"材料: {material_type} (n = {self.data['refractive_index']})")
        print(f"入射角: {incident_angle}°")
        
    def analyze_multi_beam_conditions(self):
        """
        分析多光束干涉条件
        
        返回:
            conditions: 多光束干涉条件分析结果
        """
        if self.data is None:
            return None
            
        print("\n=== 多光束干涉条件分析 ===")
        
        wavenumber = self.data['wavenumber']
        reflectance = self.data['reflectance']
        
        conditions = {}
        
        # 1. 反射率调制深度分析
        modulation_analysis = self._analyze_modulation_depth(reflectance)
        conditions['modulation'] = modulation_analysis
        
        # 2. 干涉条纹周期性分析
        periodicity_analysis = self._analyze_interference_periodicity(wavenumber, reflectance)
        conditions['periodicity'] = periodicity_analysis
        
        # 3. 相位相干性分析
        coherence_analysis = self._analyze_phase_coherence(wavenumber, reflectance)
        conditions['coherence'] = coherence_analysis
        
        # 4. 多次反射强度分析
        reflection_analysis = self._analyze_multiple_reflections(reflectance)
        conditions['reflection'] = reflection_analysis
        
        # 5. 综合判断
        overall_assessment = self._assess_multi_beam_presence(conditions)
        conditions['overall'] = overall_assessment
        
        self.multi_beam_conditions = conditions
        
        print(f"多光束干涉判断: {'存在' if overall_assessment['present'] else '不存在'}")
        print(f"置信度: {overall_assessment['confidence']:.1%}")
        
        return conditions
    
    def _analyze_modulation_depth(self, reflectance):
        """
        分析反射率调制深度
        
        参数:
            reflectance: 反射率数组
            
        返回:
            modulation_info: 调制深度分析结果
        """
        # 计算全局调制深度
        max_ref = np.max(reflectance)
        min_ref = np.min(reflectance)
        global_modulation = (max_ref - min_ref) / (max_ref + min_ref)
        
        # 计算局部调制深度
        window_size = len(reflectance) // 10
        local_modulations = []
        
        for i in range(0, len(reflectance) - window_size, window_size // 2):
            window_data = reflectance[i:i + window_size]
            local_max = np.max(window_data)
            local_min = np.min(window_data)
            if local_max + local_min > 0:
                local_mod = (local_max - local_min) / (local_max + local_min)
                local_modulations.append(local_mod)
        
        avg_local_modulation = np.mean(local_modulations) if local_modulations else 0
        modulation_contrast = np.std(local_modulations) if len(local_modulations) > 1 else 0
        
        # 判断是否满足多光束干涉的调制深度条件
        # 多光束干涉通常具有更高的调制深度和更好的对比度
        threshold = 0.3  # 调制深度阈值
        meets_condition = global_modulation > threshold and avg_local_modulation > threshold * 0.8
        
        modulation_info = {
            'global_modulation': global_modulation,
            'average_local_modulation': avg_local_modulation,
            'modulation_contrast': modulation_contrast,
            'threshold': threshold,
            'meets_condition': meets_condition
        }
        
        print(f"  调制深度分析: 全局={global_modulation:.3f}, 局部平均={avg_local_modulation:.3f}")
        print(f"  满足条件: {'是' if meets_condition else '否'}")
        
        return modulation_info
    
    def _analyze_interference_periodicity(self, wavenumber, reflectance):
        """
        分析干涉条纹周期性
        
        参数:
            wavenumber: 波数数组
            reflectance: 反射率数组
            
        返回:
            periodicity_info: 周期性分析结果
        """
        # 寻找极值点
        peaks, _ = signal.find_peaks(reflectance, prominence=0.5, distance=5)
        valleys, _ = signal.find_peaks(-reflectance, prominence=0.5, distance=5)
        
        # 计算峰值间距
        if len(peaks) > 1:
            peak_intervals = np.diff(wavenumber[peaks])
            avg_peak_interval = np.mean(peak_intervals)
            peak_regularity = 1 - np.std(peak_intervals) / avg_peak_interval if avg_peak_interval > 0 else 0
        else:
            peak_intervals = []
            avg_peak_interval = 0
            peak_regularity = 0
        
        # 计算谷值间距
        if len(valleys) > 1:
            valley_intervals = np.diff(wavenumber[valleys])
            avg_valley_interval = np.mean(valley_intervals)
            valley_regularity = 1 - np.std(valley_intervals) / avg_valley_interval if avg_valley_interval > 0 else 0
        else:
            valley_intervals = []
            avg_valley_interval = 0
            valley_regularity = 0
        
        # FFT分析周期性
        if len(reflectance) > 10:
            # 去除直流分量
            reflectance_ac = reflectance - np.mean(reflectance)
            
            # FFT
            fft_result = fft(reflectance_ac)
            freqs = fftfreq(len(reflectance_ac), d=np.mean(np.diff(wavenumber)))
            
            # 找到主频率
            positive_freqs = freqs[freqs > 0]
            positive_fft = np.abs(fft_result[freqs > 0])
            
            if len(positive_fft) > 0:
                main_freq_idx = np.argmax(positive_fft)
                main_frequency = positive_freqs[main_freq_idx]
                spectral_purity = positive_fft[main_freq_idx] / np.sum(positive_fft)
            else:
                main_frequency = 0
                spectral_purity = 0
        else:
            main_frequency = 0
            spectral_purity = 0
        
        # 综合周期性评分
        overall_regularity = (peak_regularity + valley_regularity) / 2
        
        # 判断是否满足周期性条件
        regularity_threshold = 0.7
        peak_count_threshold = 3
        meets_condition = (overall_regularity > regularity_threshold and 
                          len(peaks) >= peak_count_threshold)
        
        periodicity_info = {
            'peak_count': len(peaks),
            'valley_count': len(valleys),
            'peak_intervals': peak_intervals,
            'valley_intervals': valley_intervals,
            'avg_peak_interval': avg_peak_interval,
            'avg_valley_interval': avg_valley_interval,
            'peak_regularity': peak_regularity,
            'valley_regularity': valley_regularity,
            'overall_regularity': overall_regularity,
            'main_frequency': main_frequency,
            'spectral_purity': spectral_purity,
            'meets_condition': meets_condition
        }
        
        print(f"  周期性分析: 峰值={len(peaks)}个, 谷值={len(valleys)}个")
        print(f"  规律性: {overall_regularity:.3f}, 满足条件: {'是' if meets_condition else '否'}")
        
        return periodicity_info
    
    def _analyze_phase_coherence(self, wavenumber, reflectance):
        """
        分析相位相干性
        
        参数:
            wavenumber: 波数数组
            reflectance: 反射率数组
            
        返回:
            coherence_info: 相干性分析结果
        """
        # 计算瞬时相位
        analytic_signal = signal.hilbert(reflectance - np.mean(reflectance))
        instantaneous_phase = np.unwrap(np.angle(analytic_signal))
        
        # 计算相位导数（瞬时频率）
        if len(instantaneous_phase) > 1:
            phase_derivative = np.gradient(instantaneous_phase, wavenumber)
            phase_stability = 1 / (1 + np.std(phase_derivative))
        else:
            phase_stability = 0
        
        # 计算相位相关性
        if len(instantaneous_phase) > 10:
            # 分段计算相位相关性
            segment_size = len(instantaneous_phase) // 5
            correlations = []
            
            for i in range(4):
                seg1 = instantaneous_phase[i*segment_size:(i+1)*segment_size]
                seg2 = instantaneous_phase[(i+1)*segment_size:(i+2)*segment_size]
                if len(seg1) > 0 and len(seg2) > 0:
                    corr = np.corrcoef(seg1, seg2)[0, 1]
                    if not np.isnan(corr):
                        correlations.append(abs(corr))
            
            avg_correlation = np.mean(correlations) if correlations else 0
        else:
            avg_correlation = 0
        
        # 相干长度估计
        autocorr = np.correlate(reflectance, reflectance, mode='full')
        autocorr = autocorr[autocorr.size // 2:]
        autocorr = autocorr / autocorr[0]  # 归一化
        
        # 找到相关性下降到1/e的位置
        coherence_threshold = 1/np.e
        coherence_length_idx = np.where(autocorr < coherence_threshold)[0]
        if len(coherence_length_idx) > 0:
            coherence_length = coherence_length_idx[0]
        else:
            coherence_length = len(autocorr)
        
        # 判断是否满足相干性条件
        coherence_threshold_val = 0.8
        meets_condition = (phase_stability > 0.5 and 
                          avg_correlation > coherence_threshold_val and
                          coherence_length > 10)
        
        coherence_info = {
            'phase_stability': phase_stability,
            'average_correlation': avg_correlation,
            'coherence_length': coherence_length,
            'instantaneous_phase': instantaneous_phase,
            'meets_condition': meets_condition
        }
        
        print(f"  相干性分析: 相位稳定性={phase_stability:.3f}, 平均相关性={avg_correlation:.3f}")
        print(f"  满足条件: {'是' if meets_condition else '否'}")
        
        return coherence_info
    
    def _analyze_multiple_reflections(self, reflectance):
        """
        分析多次反射强度
        
        参数:
            reflectance: 反射率数组
            
        返回:
            reflection_info: 多次反射分析结果
        """
        # 计算反射率的高阶统计量
        mean_reflectance = np.mean(reflectance)
        std_reflectance = np.std(reflectance)
        skewness = self._calculate_skewness(reflectance)
        kurtosis = self._calculate_kurtosis(reflectance)
        
        # 分析反射率分布
        # 多光束干涉通常导致更复杂的反射率分布
        hist, bin_edges = np.histogram(reflectance, bins=20)
        hist_normalized = hist / np.sum(hist)
        
        # 计算分布的复杂度（熵）
        entropy = -np.sum(hist_normalized * np.log(hist_normalized + 1e-10))
        
        # 峰值尖锐度分析
        peaks, peak_properties = signal.find_peaks(reflectance, prominence=0.5)
        if len(peaks) > 0 and 'prominences' in peak_properties:
            avg_prominence = np.mean(peak_properties['prominences'])
            max_prominence = np.max(peak_properties['prominences'])
        else:
            avg_prominence = 0
            max_prominence = 0
        
        # 判断多次反射强度
        # 多光束干涉通常具有更高的峰值尖锐度和更复杂的分布
        reflection_threshold = 0.3
        meets_condition = (avg_prominence > reflection_threshold and 
                          entropy > 2.0 and 
                          abs(kurtosis) > 0.5)
        
        reflection_info = {
            'mean_reflectance': mean_reflectance,
            'std_reflectance': std_reflectance,
            'skewness': skewness,
            'kurtosis': kurtosis,
            'entropy': entropy,
            'avg_prominence': avg_prominence,
            'max_prominence': max_prominence,
            'peak_count': len(peaks),
            'meets_condition': meets_condition
        }
        
        print(f"  多次反射分析: 平均突出度={avg_prominence:.3f}, 熵={entropy:.3f}")
        print(f"  满足条件: {'是' if meets_condition else '否'}")
        
        return reflection_info
    
    def _assess_multi_beam_presence(self, conditions):
        """
        综合评估多光束干涉的存在性
        
        参数:
            conditions: 各项条件分析结果
            
        返回:
            assessment: 综合评估结果
        """
        # 各项条件的权重
        weights = {
            'modulation': 0.3,
            'periodicity': 0.3,
            'coherence': 0.2,
            'reflection': 0.2
        }
        
        # 计算加权得分
        total_score = 0
        for condition_name, weight in weights.items():
            if condition_name in conditions:
                condition_met = conditions[condition_name]['meets_condition']
                total_score += weight if condition_met else 0
        
        # 判断多光束干涉是否存在
        threshold = 0.6  # 阈值
        present = total_score > threshold
        confidence = total_score
        
        # 生成详细报告
        report = {
            'modulation_met': conditions.get('modulation', {}).get('meets_condition', False),
            'periodicity_met': conditions.get('periodicity', {}).get('meets_condition', False),
            'coherence_met': conditions.get('coherence', {}).get('meets_condition', False),
            'reflection_met': conditions.get('reflection', {}).get('meets_condition', False)
        }
        
        assessment = {
            'present': present,
            'confidence': confidence,
            'total_score': total_score,
            'threshold': threshold,
            'detailed_report': report
        }
        
        return assessment
    
    def build_multi_beam_model(self, thickness_estimate=None):
        """
        建立多光束干涉模型
        
        参数:
            thickness_estimate: 厚度初始估计值 (μm)
            
        返回:
            model: 多光束干涉模型
        """
        if self.data is None or not self.multi_beam_conditions:
            return None
            
        print("\n=== 建立多光束干涉模型 ===")
        
        wavenumber = self.data['wavenumber']
        reflectance = self.data['reflectance']
        wavelength = self.data['wavelength']
        n = self.data['refractive_index']
        theta = np.radians(self.data['incident_angle'])
        
        # 计算折射角
        theta_r = np.arcsin(np.sin(theta) / n)
        
        # 如果没有提供厚度估计，使用简单方法估算
        if thickness_estimate is None:
            thickness_estimate = self._estimate_thickness_simple(wavenumber, reflectance, n, theta_r)
        
        # 多光束干涉模型参数
        model_params = {
            'thickness': thickness_estimate,
            'refractive_index': n,
            'incident_angle': self.data['incident_angle'],
            'refraction_angle': np.degrees(theta_r)
        }
        
        # 拟合多光束干涉模型
        try:
            fitted_params = self._fit_multi_beam_model(wavenumber, reflectance, model_params)
            model_params.update(fitted_params)
        except Exception as e:
            print(f"  模型拟合警告: {str(e)}")
        
        # 计算理论光谱
        theoretical_spectrum = self._calculate_multi_beam_spectrum(wavenumber, model_params)
        
        # 计算拟合优度
        r_squared = self._calculate_r_squared(reflectance, theoretical_spectrum)
        
        model = {
            'parameters': model_params,
            'theoretical_spectrum': theoretical_spectrum,
            'r_squared': r_squared,
            'wavenumber': wavenumber,
            'experimental_spectrum': reflectance
        }
        
        print(f"  模型建立完成")
        print(f"  拟合优度 R² = {r_squared:.3f}")
        print(f"  修正厚度: {model_params['thickness']:.3f} μm")
        
        return model
    
    def correct_thickness_measurement(self, original_thickness, correction_method='analytical'):
        """
        修正厚度测量结果
        
        参数:
            original_thickness: 原始厚度测量值 (μm)
            correction_method: 修正方法 ('analytical' 或 'numerical')
            
        返回:
            corrected_results: 修正后的结果
        """
        if not self.multi_beam_conditions.get('overall', {}).get('present', False):
            print("未检测到多光束干涉，无需修正")
            return {
                'original_thickness': original_thickness,
                'corrected_thickness': original_thickness,
                'correction_factor': 1.0,
                'correction_applied': False
            }
        
        print("\n=== 厚度测量修正 ===")
        
        if correction_method == 'analytical':
            corrected_thickness, correction_factor = self._analytical_correction(original_thickness)
        else:
            corrected_thickness, correction_factor = self._numerical_correction(original_thickness)
        
        corrected_results = {
            'original_thickness': original_thickness,
            'corrected_thickness': corrected_thickness,
            'correction_factor': correction_factor,
            'correction_applied': True,
            'method': correction_method
        }
        
        print(f"  原始厚度: {original_thickness:.3f} μm")
        print(f"  修正厚度: {corrected_thickness:.3f} μm")
        print(f"  修正因子: {correction_factor:.3f}")
        
        return corrected_results
    
    def _estimate_thickness_simple(self, wavenumber, reflectance, n, theta_r):
        """
        简单厚度估算
        """
        # 寻找相邻极值点
        peaks, _ = signal.find_peaks(reflectance, prominence=0.5)
        if len(peaks) < 2:
            return 1.0  # 默认值
        
        # 使用前两个峰值估算厚度
        lambda1 = 10000.0 / wavenumber[peaks[0]]
        lambda2 = 10000.0 / wavenumber[peaks[1]]
        
        # 简单公式: d = λ1*λ2 / (2*n*cos(θr)*(λ1-λ2))
        if abs(lambda1 - lambda2) > 0.001:
            thickness = (lambda1 * lambda2) / (2 * n * np.cos(theta_r) * abs(lambda1 - lambda2))
        else:
            thickness = 1.0
        
        return thickness
    
    def _fit_multi_beam_model(self, wavenumber, reflectance, initial_params):
        """
        拟合多光束干涉模型
        """
        def model_function(wn, thickness, r1, r2, phase_offset):
            wavelength = 10000.0 / wn
            n = initial_params['refractive_index']
            theta_r = np.radians(initial_params['refraction_angle'])
            
            # 多光束干涉公式（简化版）
            delta = 4 * np.pi * n * thickness * np.cos(theta_r) / wavelength + phase_offset
            
            # Airy函数近似
            finesse = 4 * r1 * r2 / (1 - r1 * r2)**2
            intensity = (1 - r1) / (1 + finesse * np.sin(delta/2)**2)
            
            return intensity * 100  # 转换为百分比
        
        # 初始参数猜测
        p0 = [
            initial_params['thickness'],  # 厚度
            0.3,  # r1
            0.3,  # r2
            0.0   # 相位偏移
        ]
        
        # 参数边界
        bounds = (
            [0.1, 0.01, 0.01, -np.pi],
            [10.0, 0.9, 0.9, np.pi]
        )
        
        try:
            popt, _ = curve_fit(model_function, wavenumber, reflectance, 
                              p0=p0, bounds=bounds, maxfev=1000)
            
            fitted_params = {
                'thickness': popt[0],
                'r1': popt[1],
                'r2': popt[2],
                'phase_offset': popt[3]
            }
        except:
            # 如果拟合失败，返回初始参数
            fitted_params = {
                'r1': 0.3,
                'r2': 0.3,
                'phase_offset': 0.0
            }
        
        return fitted_params
    
    def _calculate_multi_beam_spectrum(self, wavenumber, params):
        """
        计算多光束干涉理论光谱
        """
        wavelength = 10000.0 / wavenumber
        n = params['refractive_index']
        thickness = params['thickness']
        theta_r = np.radians(params['refraction_angle'])
        r1 = params.get('r1', 0.3)
        r2 = params.get('r2', 0.3)
        phase_offset = params.get('phase_offset', 0.0)
        
        # 相位差
        delta = 4 * np.pi * n * thickness * np.cos(theta_r) / wavelength + phase_offset
        
        # 多光束干涉强度（Airy函数）
        finesse = 4 * r1 * r2 / (1 - r1 * r2)**2
        intensity = (1 - r1) / (1 + finesse * np.sin(delta/2)**2)
        
        return intensity * 100
    
    def _analytical_correction(self, original_thickness):
        """
        解析修正方法
        """
        # 基于多光束干涉理论的解析修正
        conditions = self.multi_beam_conditions
        
        # 修正因子基于各项条件的强度
        modulation_factor = conditions.get('modulation', {}).get('global_modulation', 0)
        periodicity_factor = conditions.get('periodicity', {}).get('overall_regularity', 0)
        
        # 经验修正公式
        correction_factor = 1.0 + 0.1 * modulation_factor + 0.05 * periodicity_factor
        corrected_thickness = original_thickness / correction_factor
        
        return corrected_thickness, correction_factor
    
    def _numerical_correction(self, original_thickness):
        """
        数值修正方法
        """
        # 使用数值优化进行修正
        def objective_function(thickness):
            # 计算理论光谱
            params = {
                'thickness': thickness[0],
                'refractive_index': self.data['refractive_index'],
                'refraction_angle': np.arcsin(np.sin(np.radians(self.data['incident_angle'])) / self.data['refractive_index']) * 180 / np.pi,
                'r1': 0.3,
                'r2': 0.3,
                'phase_offset': 0.0
            }
            
            theoretical = self._calculate_multi_beam_spectrum(self.data['wavenumber'], params)
            
            # 计算与实验数据的差异
            error = np.sum((self.data['reflectance'] - theoretical)**2)
            return error
        
        # 优化
        result = minimize(objective_function, [original_thickness], 
                         bounds=[(0.1, 10.0)], method='L-BFGS-B')
        
        if result.success:
            corrected_thickness = result.x[0]
            correction_factor = original_thickness / corrected_thickness
        else:
            corrected_thickness = original_thickness
            correction_factor = 1.0
        
        return corrected_thickness, correction_factor
    
    def _calculate_skewness(self, data):
        """计算偏度"""
        mean_val = np.mean(data)
        std_val = np.std(data)
        if std_val == 0:
            return 0
        return np.mean(((data - mean_val) / std_val) ** 3)
    
    def _calculate_kurtosis(self, data):
        """计算峰度"""
        mean_val = np.mean(data)
        std_val = np.std(data)
        if std_val == 0:
            return 0
        return np.mean(((data - mean_val) / std_val) ** 4) - 3
    
    def _calculate_r_squared(self, observed, predicted):
        """计算决定系数"""
        ss_res = np.sum((observed - predicted) ** 2)
        ss_tot = np.sum((observed - np.mean(observed)) ** 2)
        if ss_tot == 0:
            return 0
        return 1 - (ss_res / ss_tot)
    
    def plot_analysis_results(self, save_figure=True, figure_name='multi_beam_analysis.png'):
        """
        绘制分析结果
        
        参数:
            save_figure: 是否保存图片
            figure_name: 图片文件名
        """
        if self.data is None:
            return
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        wavenumber = self.data['wavenumber']
        reflectance = self.data['reflectance']
        wavelength = self.data['wavelength']
        
        # 子图1: 原始光谱
        axes[0, 0].plot(wavenumber, reflectance, 'b-', linewidth=2)
        axes[0, 0].set_xlabel('波数 (cm⁻¹)')
        axes[0, 0].set_ylabel('反射率 (%)')
        axes[0, 0].set_title('原始光谱数据')
        axes[0, 0].grid(True, alpha=0.3)
        
        # 子图2: 波长域光谱
        axes[0, 1].plot(wavelength, reflectance, 'r-', linewidth=2)
        axes[0, 1].set_xlabel('波长 (μm)')
        axes[0, 1].set_ylabel('反射率 (%)')
        axes[0, 1].set_title('波长域光谱')
        axes[0, 1].grid(True, alpha=0.3)
        
        # 子图3: 极值点标记
        axes[0, 2].plot(wavenumber, reflectance, 'b-', linewidth=1, alpha=0.7)
        peaks, _ = signal.find_peaks(reflectance, prominence=0.5)
        valleys, _ = signal.find_peaks(-reflectance, prominence=0.5)
        axes[0, 2].plot(wavenumber[peaks], reflectance[peaks], 'ro', markersize=8, label=f'峰值 ({len(peaks)}个)')
        axes[0, 2].plot(wavenumber[valleys], reflectance[valleys], 'go', markersize=8, label=f'谷值 ({len(valleys)}个)')
        axes[0, 2].set_xlabel('波数 (cm⁻¹)')
        axes[0, 2].set_ylabel('反射率 (%)')
        axes[0, 2].set_title('极值点识别')
        axes[0, 2].legend()
        axes[0, 2].grid(True, alpha=0.3)
        
        # 子图4: FFT频谱分析
        reflectance_ac = reflectance - np.mean(reflectance)
        fft_result = np.abs(fft(reflectance_ac))
        freqs = fftfreq(len(reflectance_ac), d=np.mean(np.diff(wavenumber)))
        positive_mask = freqs > 0
        axes[1, 0].plot(freqs[positive_mask], fft_result[positive_mask], 'g-', linewidth=2)
        axes[1, 0].set_xlabel('频率 (cm)')
        axes[1, 0].set_ylabel('幅度')
        axes[1, 0].set_title('FFT频谱分析')
        axes[1, 0].grid(True, alpha=0.3)
        
        # 子图5: 多光束干涉条件评估
        if self.multi_beam_conditions:
            conditions = self.multi_beam_conditions
            condition_names = ['调制深度', '周期性', '相干性', '多次反射']
            condition_values = [
                conditions.get('modulation', {}).get('meets_condition', False),
                conditions.get('periodicity', {}).get('meets_condition', False),
                conditions.get('coherence', {}).get('meets_condition', False),
                conditions.get('reflection', {}).get('meets_condition', False)
            ]
            
            colors = ['green' if val else 'red' for val in condition_values]
            bars = axes[1, 1].bar(condition_names, [1 if val else 0 for val in condition_values], color=colors, alpha=0.7)
            axes[1, 1].set_ylabel('满足条件')
            axes[1, 1].set_title('多光束干涉条件评估')
            axes[1, 1].set_ylim(0, 1.2)
            
            # 添加文本标签
            for bar, val in zip(bars, condition_values):
                height = bar.get_height()
                axes[1, 1].text(bar.get_x() + bar.get_width()/2., height + 0.05,
                               '满足' if val else '不满足',
                               ha='center', va='bottom')
        
        # 子图6: 综合评估结果
        if self.multi_beam_conditions and 'overall' in self.multi_beam_conditions:
            overall = self.multi_beam_conditions['overall']
            confidence = overall['confidence']
            present = overall['present']
            
            # 绘制置信度条形图
            axes[1, 2].bar(['多光束干涉'], [confidence], 
                          color='green' if present else 'red', alpha=0.7)
            axes[1, 2].axhline(y=overall['threshold'], color='black', linestyle='--', label='阈值')
            axes[1, 2].set_ylabel('置信度')
            axes[1, 2].set_title(f"综合判断: {'存在' if present else '不存在'}")
            axes[1, 2].set_ylim(0, 1)
            axes[1, 2].legend()
            
            # 添加置信度数值
            axes[1, 2].text(0, confidence + 0.05, f'{confidence:.1%}', 
                           ha='center', va='bottom', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        
        if save_figure:
            plt.savefig(figure_name, dpi=300, bbox_inches='tight')
            print(f"分析结果图已保存: {figure_name}")
        
        plt.show()
    
    def export_results(self, output_file='multi_beam_analysis.xlsx'):
        """
        导出分析结果
        
        参数:
            output_file: 输出文件名
        """
        if not self.multi_beam_conditions:
            print("没有分析结果可导出")
            return False
        
        try:
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                # 导出条件分析结果
                conditions_data = []
                for condition_name, condition_data in self.multi_beam_conditions.items():
                    if isinstance(condition_data, dict) and 'meets_condition' in condition_data:
                        conditions_data.append({
                            '条件': condition_name,
                            '满足': '是' if condition_data['meets_condition'] else '否',
                            '详细信息': str(condition_data)
                        })
                
                if conditions_data:
                    conditions_df = pd.DataFrame(conditions_data)
                    conditions_df.to_excel(writer, sheet_name='条件分析', index=False)
                
                # 导出综合评估结果
                if 'overall' in self.multi_beam_conditions:
                    overall = self.multi_beam_conditions['overall']
                    assessment_df = pd.DataFrame([
                        {'项目': '多光束干涉存在性', '结果': '是' if overall['present'] else '否'},
                        {'项目': '置信度', '结果': f"{overall['confidence']:.1%}"},
                        {'项目': '总得分', '结果': f"{overall['total_score']:.3f}"},
                        {'项目': '判断阈值', '结果': f"{overall['threshold']:.3f}"}
                    ])
                    assessment_df.to_excel(writer, sheet_name='综合评估', index=False)
            
            print(f"分析结果已导出至: {output_file}")
            return True
            
        except Exception as e:
            print(f"导出失败: {str(e)}")
            return False

def main():
    """主函数 - 多光束干涉分析演示"""
    print("=== 多光束干涉分析系统 ===")
    
    # 创建分析器
    analyzer = MultiBeamAnalyzer()
    
    # 示例：分析附件3数据（硅片）
    try:
        # 这里需要实际的数据加载
        # 示例数据
        wavenumber = np.linspace(400, 4000, 1000)
        # 模拟多光束干涉光谱
        wavelength = 10000.0 / wavenumber
        thickness = 2.0  # μm
        n_si = 3.42
        theta_r = np.radians(10.0 / n_si)
        
        # 多光束干涉模拟
        delta = 4 * np.pi * n_si * thickness * np.cos(theta_r) / wavelength
        r1, r2 = 0.3, 0.3
        finesse = 4 * r1 * r2 / (1 - r1 * r2)**2
        reflectance = (1 - r1) / (1 + finesse * np.sin(delta/2)**2) * 100
        
        # 添加噪声
        reflectance += np.random.normal(0, 0.5, len(reflectance))
        
        # 加载数据
        analyzer.load_data(wavenumber, reflectance, 'Si', 10.0)
        
        # 分析多光束干涉条件
        conditions = analyzer.analyze_multi_beam_conditions()
        
        if conditions and conditions.get('overall', {}).get('present', False):
            print("\n检测到多光束干涉现象")
            
            # 建立多光束干涉模型
            model = analyzer.build_multi_beam_model(thickness_estimate=2.0)
            
            # 厚度修正
            corrected_results = analyzer.correct_thickness_measurement(2.1, 'analytical')
            
        else:
            print("\n未检测到显著的多光束干涉现象")
        
        # 绘制分析结果
        analyzer.plot_analysis_results(save_figure=True, figure_name='multi_beam_demo.png')
        
        # 导出结果
        analyzer.export_results('multi_beam_demo_results.xlsx')
        
        print("\n多光束干涉分析完成！")
        
    except Exception as e:
        print(f"分析过程出错: {str(e)}")

if __name__ == "__main__":
    main()