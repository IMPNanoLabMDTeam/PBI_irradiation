#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PBI辐射计算器 - CLI版本
基于SRIM数据和Waligorski-Zhang模型的PBI材料辐射损伤计算
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os
from datetime import datetime
import argparse

# 设置matplotlib字体和参数，以确保图形可正确生成
plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial', 'sans-serif']
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.autolayout'] = True
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['font.size'] = 10

class SRIM数据处理器:
    """
    SRIM数据处理类
    
    功能说明：
    - 处理来自SRIM计算的阻止本领数据
    - 进行密度修正以适应实际PBI材料参数
    - 计算离子在PBI中的能量损失
    
    数据来源：SRIM-2013计算结果，86Kr在PBI中1.56GeV时的阻止本领
    """
    
    def __init__(self):
        self.srim阻止本领_keV_um = 5.165E+03
        self.srim密度 = 1.5
        
    def 计算能量损失(self, 初始能量_MeV_amu, pbi厚度_um, pbi密度_g_cm3=1.947):
        密度修正因子 = pbi密度_g_cm3 / self.srim密度
        修正后阻止本领 = self.srim阻止本领_keV_um * 密度修正因子
        总能量损失_keV = 修正后阻止本领 * pbi厚度_um
        能量损失_MeV_amu = 总能量损失_keV / (1000.0 * 86.0)
        最终能量 = max(0, 初始能量_MeV_amu - 能量损失_MeV_amu)
        
        return {
            '初始能量': 初始能量_MeV_amu,
            '最终能量': 最终能量,
            '能量损失': 能量损失_MeV_amu,
            '损失百分比': (能量损失_MeV_amu / 初始能量_MeV_amu) * 100,
            '密度修正因子': 密度修正因子,
            '修正后阻止本领': 修正后阻止本领
        }

class WZ模型计算器:
    """
    Waligorski-Zhang径向剂量分布模型计算器
    """
    
    def __init__(self):
        self.电子质量_keV = 511.0
        self.水常数_keV_per_mm = 8.5
        self.射程参数_g_per_cm2_keV__alpha = 6e-6
        self.阿伏加德罗常数 = 6.022e23
        self.eV转K = 11605
        self.pbi分子量 = 308.0
        
    def 计算β和γ(self, 能量_MeV_per_amu):
        γ = 1 + 能量_MeV_per_amu / 931.5
        β = math.sqrt(1 - 1/γ**2)
        return β, γ
    
    def 计算有效电荷(self, Z, β):
        return Z * (1 - math.exp(-125 * β * Z**(-2/3)))
    
    def 计算射程参数α(self, β):
        return 1.079 if β < 0.03 else 1.667
    
    def 计算径向剂量(self, 半径_nm, Z, 能量_MeV_per_amu, 密度_g_per_cm3, 电离能_eV=12.0):
        β, γ = self.计算β和γ(能量_MeV_per_amu)
        Z星 = self.计算有效电荷(Z, β)
        α = self.计算射程参数α(β)
        
        电离能_keV = 电离能_eV / 1000.0
        θ = self.射程参数_g_per_cm2_keV__alpha * (电离能_keV ** α)
        W = 2 * self.电子质量_keV * β**2 * γ**2
        T = self.射程参数_g_per_cm2_keV__alpha * (W ** α)
        
        半径_cm = 半径_nm * 1e-7
        t_g_per_cm2 = 半径_cm * 密度_g_per_cm3
        
        if t_g_per_cm2 <= 0 or (T + θ) <= 0 or α <= 0:
            return 0
        
        水常数_keV_per_cm = self.水常数_keV_per_mm * 1e1
        常数因子_keV_cm_per_g = 水常数_keV_per_cm * Z星**2 / (2 * np.pi * α * β**2 * t_g_per_cm2)
        分数 = (t_g_per_cm2 + θ) / (T + θ)
        
        if 分数 >= 1:
            return 0
        
        幂次项 = (1 - 分数) ** (1/α)
        剂量_keV_cm3_per_g2 = 常数因子_keV_cm_per_g * 幂次项 / (t_g_per_cm2 + θ)
        剂量_keV_per_cm3 = 剂量_keV_cm3_per_g2 * (密度_g_per_cm3)**2
        剂量_keV_per_nm3 = 剂量_keV_per_cm3 * 1e-21

        return 剂量_keV_per_nm3
    
    def 剂量转温度(self, 剂量_keV_per_nm3, 密度_g_cm3):
        if 剂量_keV_per_nm3 <= 0:
            return 0
        
        pbi分子量 = 308.0
        每分子原子数 = 36
        
        分子数密度_per_nm3 = (密度_g_cm3 * 1e-21 * self.阿伏加德罗常数) / pbi分子量
        原子数密度_per_nm3 = 分子数密度_per_nm3 * 每分子原子数
        每个原子能量_eV = 剂量_keV_per_nm3 * 1e3 / 原子数密度_per_nm3
        温度_K = 每个原子能量_eV * self.eV转K
            
        return 温度_K
    
    def 计算径向能量分布(self, 半径数组_nm, Z, 能量_MeV_amu, 密度_g_cm3, 电离能_eV=12.0):
        剂量数组 = np.array([
            self.计算径向剂量(r, Z, 能量_MeV_amu, 密度_g_cm3, 电离能_eV)
            for r in 半径数组_nm
        ])
        
        pbi分子量 = 308.0
        每分子原子数 = 36
        
        分子数密度_nm3 = (密度_g_cm3 * self.阿伏加德罗常数 * 1e-21) / pbi分子量
        原子数密度_nm3 = 分子数密度_nm3 * 每分子原子数
        
        单分子能量数组_eV = (剂量数组 * 1000.0) / 分子数密度_nm3
        单原子能量数组_eV = (剂量数组 * 1000.0) / 原子数密度_nm3
        
        累积能量 = np.zeros_like(半径数组_nm)
        for i in range(1, len(半径数组_nm)):
            dr = 半径数组_nm[i] - 半径数组_nm[i-1]
            圆环面积 = 2 * np.pi * 半径数组_nm[i] * dr
            圆环能量 = 剂量数组[i] * 圆环面积 * 1.0
            累积能量[i] = 累积能量[i-1] + 圆环能量
        
        return {
            '半径': 半径数组_nm,
            '剂量密度': 剂量数组,
            '单分子能量_eV': 单分子能量数组_eV,
            '单原子能量_eV': 单原子能量数组_eV,
            '累积能量': 累积能量,
            '分子数密度_nm3': 分子数密度_nm3,
            '原子数密度_nm3': 原子数密度_nm3,
            '每分子原子数': 每分子原子数,
            '总能量_keV': 累积能量[-1] if len(累积能量) > 0 else 0
        }

def plot_results(figure, 半径, 剂量, 归一化剂量, 温度, 能量信息, g因子, 输入参数, 分析信息):
    """绘制计算结果图表"""
    figure.clear()
    
    子图1 = figure.add_subplot(221)
    子图2 = figure.add_subplot(222)
    子图3 = figure.add_subplot(223)
    子图4 = figure.add_subplot(224)
    
    # 子图1：剂量分布-对数
    子图1.loglog(半径, 剂量, 'b-', label='Original Dose')
    子图1.loglog(半径, 归一化剂量, 'r--', label=f'Normalized (g={g因子})')
    子图1.set_xlabel('Radial Distance (nm)')
    子图1.set_ylabel('Dose (keV/nm^3)')
    子图1.set_title('Dose Distribution (Log)')
    子图1.grid(True, alpha=0.3)
    子图1.legend()
    
    # 子图2：剂量分布-线性
    子图2.plot(半径, 剂量, 'b-', label='Original Dose')
    子图2.plot(半径, 归一化剂量, 'r--', label=f'Normalized (g={g因子})')
    子图2.set_xlabel('Radial Distance (nm)')
    子图2.set_ylabel('Dose (keV/nm^3)')
    子图2.set_title('Dose Distribution (Linear)')
    子图2.grid(True, alpha=0.3)
    子图2.legend()

    # 子图3：温度分布-对数
    子图3.loglog(半径, 温度, 'g-', label='Temperature')
    子图3.set_xlabel('Radial Distance (nm)')
    子图3.set_ylabel('Temperature (K)')
    子图3.set_title('Temperature Distribution (Log)')
    子图3.grid(True, alpha=0.3)
    子图3.legend()

    # 子图4：单原子能量分布
    单原子能量 = 分析信息['单原子能量数组_eV']
    子图4.loglog(半径, 单原子能量, 'm-', label='Energy per PBI Atom')
    子图4.set_xlabel('Radial Distance (nm)')
    子图4.set_ylabel('Energy per Atom (eV)')
    子图4.set_title('Single Atom Energy Distribution')
    子图4.grid(True, alpha=0.3)
    子图4.legend()

    总标题 = f'86Kr Radiation in PBI: E_in={能量信息["初始能量"]:.1f} MeV/u, E_final={能量信息["最终能量"]:.1f} MeV/u'
    figure.suptitle(总标题, fontsize=14, fontweight='bold')
    figure.tight_layout(rect=[0, 0, 1, 0.96])

def save_data_file(文件名, 计算结果):
    """保存详细的数据文件"""
    with open(文件名, 'w', encoding='utf-8') as f:
        f.write(f"# PBI辐射计算结果 - CLI版本\n")
        f.write(f"# 生成时间: {datetime.now()}\n\n")
        
        f.write(f"# 输入参数:\n")
        for 键, 值 in 计算结果['参数'].items():
            f.write(f"# {键}: {值}\n")
        f.write(f"\n")
        
        能量信息 = 计算结果['能量信息']
        f.write(f"# SRIM能量损失分析:\n")
        for 键, 值 in 能量信息.items():
            f.write(f"# {键}: {值}\n")
        f.write(f"\n")
        
        分析信息 = 计算结果['分析信息']
        f.write(f"# 分析结果:\n")
        f.write(f"# PBI分子数密度: {分析信息['分子数密度_nm3']:.6e} 分子/nm³\n")
        f.write(f"# 径迹中心温度: {分析信息['径迹中心温度_K']:,.0f} K\n")
        f.write(f"# 总沉积能量: {分析信息['总能量_keV']:.6e} keV\n\n")
        
        f.write("径向距离_nm\t原始剂量_keV_nm3\t归一化剂量_keV_nm3\t温度_K\t单原子能量_eV\n")
        
        半径 = 计算结果['半径']
        原始剂量 = 计算结果['原始剂量']
        归一化剂量 = 计算结果['归一化剂量']
        温度 = 计算结果['温度']
        单原子能量 = 计算结果['分析信息']['单原子能量数组_eV']
        
        for r, d_o, d_n, t, e_a in zip(半径, 原始剂量, 归一化剂量, 温度, 单原子能量):
            f.write(f"{r:.6e}\t{d_o:.6e}\t{d_n:.6e}\t{t:.6e}\t{e_a:.6e}\n")

def main():
    """主函数：解析命令行参数并执行计算"""
    parser = argparse.ArgumentParser(description="PBI辐射计算器 - CLI版本")
    parser.add_argument('--initial-energy', type=float, default=25.0, help="初始能量 (MeV/u)")
    parser.add_argument('--thickness', type=float, default=40.0, help="PBI厚度 (μm)")
    parser.add_argument('--density', type=float, default=1.947, help="PBI密度 (g/cm³)")
    parser.add_argument('--ionization-energy', type=float, default=12.0, help="平均电离能 (eV)")
    parser.add_argument('--g-factor', type=float, default=0.17, help="归一化因子")
    parser.add_argument('--max-radius', type=float, default=100.0, help="最大半径 (nm)")
    parser.add_argument('--save', action='store_true', help="保存计算结果的图形和数据文件")

    args = parser.parse_args()

    # 实例化计算器
    srim_processor = SRIM数据处理器()
    wz_calculator = WZ模型计算器()

    print("开始计算...")
    
    # 执行计算
    能量信息 = srim_processor.计算能量损失(args.initial_energy, args.thickness, args.density)
    最终能量 = 能量信息['最终能量']
    阻止本领 = 能量信息['修正后阻止本领'] / 1000
    半径数组 = np.logspace(0, np.log10(args.max_radius), 1500)
    
    能量分布信息 = wz_calculator.计算径向能量分布(半径数组, 36, 最终能量, args.density, args.ionization_energy)
    剂量数组 = 能量分布信息['剂量密度']
    归一化剂量 = 剂量数组 * args.g_factor / 能量分布信息['累积能量'][-1] * 阻止本领
    
    # debug 累计能量
    print(f"累计能量: {能量分布信息['累积能量']}")
    print(f"阻止本领: {阻止本领}")
    print(f"归一化剂量: {归一化剂量}")
    温度数组 = np.array([wz_calculator.剂量转温度(d, args.density) for d in 归一化剂量])
    
    分析信息 = {
        '分子数密度_nm3': 能量分布信息['分子数密度_nm3'],
        '原子数密度_nm3': 能量分布信息['原子数密度_nm3'],
        '每分子原子数': 能量分布信息['每分子原子数'],
        '径迹中心温度_K': 温度数组[0] if len(温度数组) > 0 else 0,
        '最高温度_K': np.max(温度数组) if len(温度数组) > 0 else 0,
        '平均温度_K': np.mean(温度数组) if len(温度数组) > 0 else 0,
        '单分子能量数组_eV': 能量分布信息['单分子能量_eV'],
        '单原子能量数组_eV': 能量分布信息['单原子能量_eV'],
        '累积能量': 能量分布信息['累积能量'],
        '总能量_keV': 能量分布信息['总能量_keV']
    }

    # 打印关键结果到控制台
    print("\n=== 计算结果摘要 ===")
    print(f"径迹中心温度: {分析信息['径迹中心温度_K']:,.0f} K")
    print(f"最高温度: {分析信息['最高温度_K']:,.0f} K")
    print(f"PBI原子数密度: {分析信息['原子数密度_nm3']:.2e} 原子/nm³")
    print(f"总沉积能量: {分析信息['总能量_keV']:.2e} keV")
    print(f"能量损失百分比: {能量信息['损失百分比']:.2f}%")
    
    # 保存结果
    if args.save:
        print("\n正在保存结果...")
        计算结果 = {
            '半径': 半径数组, '原始剂量': 剂量数组, '归一化剂量': 归一化剂量,
            '温度': 温度数组, '能量信息': 能量信息, '分析信息': 分析信息,
            '参数': vars(args)
        }

        输出目录 = "WZ_CLI_Results"
        os.makedirs(输出目录, exist_ok=True)
        时间戳 = datetime.now().strftime("%Y%m%d_%H%M%S")
        基础名称 = f"WZ_Kr86_{args.initial_energy:.0f}MeVamu_{时间戳}"
        
        # 保存图形
        fig = plt.figure(figsize=(16, 12))
        plot_results(fig, 半径数组, 剂量数组, 归一化剂量, 温度数组, 能量信息, args.g_factor, vars(args), 分析信息)
        图形文件 = os.path.join(输出目录, f"{基础名称}_图形.png")
        fig.savefig(图形文件, dpi=300, bbox_inches='tight')
        plt.close(fig) # 关闭图形避免显示
        print(f"图形已保存到: {图形文件}")

        # 保存数据
        数据文件 = os.path.join(输出目录, f"{基础名称}_数据.txt")
        save_data_file(数据文件, 计算结果)
        print(f"数据已保存到: {数据文件}")

if __name__ == "__main__":
    main() 