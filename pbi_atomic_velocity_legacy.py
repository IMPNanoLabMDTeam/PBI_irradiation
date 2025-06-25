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
import sys

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
    - 计算离子在PBI中的能量损失
    
    数据来源：SRIM-2013计算结果，86Kr在PBI中1.56GeV时的阻止本领
    """
    
    def __init__(self):
        self.srim阻止本领_keV_per_um = 4.156E+03
        self.pbi密度_g_per_cm3 = 1.5  # 实际密度等于SRIM密度，无需修正
        
    def 计算能量损失(self, 初始能量_MeV_per_amu, pbi厚度_um):
        总能量损失_keV = self.srim阻止本领_keV_per_um * pbi厚度_um
        能量损失_MeV_per_amu = 总能量损失_keV / (1000.0 * 86.0)
        最终能量 = max(0, 初始能量_MeV_per_amu - 能量损失_MeV_per_amu)
        
        return {
            '初始能量': 初始能量_MeV_per_amu,
            '最终能量': 最终能量,
            '能量损失': 能量损失_MeV_per_amu,
            '损失百分比': (能量损失_MeV_per_amu / 初始能量_MeV_per_amu) * 100,
            '阻止本领': self.srim阻止本领_keV_per_um
        }

class WZ模型计算器:
    """
    Waligorski-Zhang径向剂量分布模型计算器
    使用数据文件中的实际原子数密度
    """
    
    def __init__(self):
        self.电子质量_keV_per_c2 = 511.0
        self.水常数_keV_per_mm = 8.5
        self.射程参数_g_per_cm2_keV__alpha = 6e-6
        self.阿伏加德罗常数_per_mol = 6.022e23
        self.eV转K = 11605
        
        # 盒子尺寸（单位：Å）
        xlo, xhi = 299.862329034454, 484.85957681257116
        ylo, yhi = 299.862329034454, 484.85957681257116
        zlo, zhi = 86.0314725411649, 146.2737064987652
        # 体积(Å³)
        volume_A3 = (xhi-xlo) * (yhi-ylo) * (zhi-zlo)
        # 体积(nm³)
        volume_nm3 = volume_A3 * 1e-3
        self.总原子数 = 208224
        self.实际原子数密度_per_nm3 = self.总原子数 / volume_nm3
        self.平均原子质量_g_per_mol = 8.457464    # 总质量/总原子数
        
        # 数据文件中的原子组成
        self.原子组成 = {
            'C': 112320,  # 碳原子数
            'N': 17280,   # 氮原子数  
            'H': 72576,   # 氢原子数
            'O': 6048     # 氧原子数
        }
        # 记录体积，便于输出
        self.盒子体积_nm3 = volume_nm3
        
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
        θ = self.射程参数_g_per_cm2_keV__alpha * (电离能_keV ** α)   # g_per_cm2
        W = 2 * self.电子质量_keV_per_c2 * β**2 * γ   # keV
        T = self.射程参数_g_per_cm2_keV__alpha * (W ** α)   # g_per_cm2
        
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
    
    def 剂量转温度(self, 剂量_keV_per_nm3, 密度_g_per_cm3):
        if 剂量_keV_per_nm3 <= 0:
            return 0
        
        self.pbi分子量_g_per_mol = 308.0
        每分子原子数 = 36
        
        分子数密度_per_nm3 = (密度_g_per_cm3 * 1e-21 * self.阿伏加德罗常数_per_mol) / self.pbi分子量_g_per_mol
        原子数密度_per_nm3 = 分子数密度_per_nm3 * 每分子原子数
        每个原子能量_eV = 剂量_keV_per_nm3 * 1e3 / 原子数密度_per_nm3
        温度_K = 每个原子能量_eV * self.eV转K
            
        return 温度_K
    
    def 计算径向能量分布(self, 半径数组_nm, Z, 能量_MeV_per_amu, 密度_g_per_cm3, 电离能_eV=12.0):
        """
        计算径向能量分布，使用实际原子数密度
        """
        剂量数组 = np.array([
            self.计算径向剂量(r, Z, 能量_MeV_per_amu, 密度_g_per_cm3, 电离能_eV)
            for r in 半径数组_nm
        ])
        
        # 使用实际原子数密度计算单原子能量
        单原子能量数组_eV = (剂量数组 * 1000.0) / self.实际原子数密度_per_nm3
        
        累积能量 = np.zeros_like(半径数组_nm)
        for i in range(1, len(半径数组_nm)):
            dr = 半径数组_nm[i] - 半径数组_nm[i-1]
            圆环面积 = 2 * np.pi * 半径数组_nm[i] * dr
            圆环能量 = 剂量数组[i] * 圆环面积 * 1.0
            累积能量[i] = 累积能量[i-1] + 圆环能量
        
        return {
            '半径': 半径数组_nm,
            '剂量密度': 剂量数组,
            '单原子能量_eV': 单原子能量数组_eV,
            '累积能量': 累积能量,
            '原子数密度_nm3': self.实际原子数密度_per_nm3,
            '平均原子质量_g_per_mol': self.平均原子质量_g_per_mol,
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
    parser.add_argument('--density', type=float, default=1.5, help="PBI密度 (g/cm³) - 根据SRIM数据")
    parser.add_argument('--ionization-energy', type=float, default=12.0, help="平均电离能 (eV)")
    parser.add_argument('--g-factor', type=float, default=0.17, help="归一化因子")
    parser.add_argument('--max-radius', type=float, default=10.0, help="最大半径 (nm)")
    parser.add_argument('--save', action='store_true', help="保存计算结果的图形和数据文件")

    args = parser.parse_args()

    # 实例化计算器
    srim_processor = SRIM数据处理器()
    wz_calculator = WZ模型计算器()

    print("开始计算...")
    
    # 执行计算
    能量信息 = srim_processor.计算能量损失(args.initial_energy, args.thickness)
    最终能量 = 能量信息['最终能量']
    阻止本领 = 能量信息['阻止本领'] / 1000
    半径数组 = np.logspace(0, np.log10(args.max_radius), 150)
    
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

def plot_dose_temperature_distributions(半径数组, 剂量数组, 归一化剂量, 温度数组, 能量信息, 输入参数, 分析信息, 输出目录):
    """
    专门绘制沉积能量和温度分布图的函数
    包含对数和线性坐标版本，确保图形清晰可读
    """
    # 设置matplotlib参数以确保高质量输出
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.titlesize': 18,
        'lines.linewidth': 2.5,
        'axes.linewidth': 1.2,
        'grid.linewidth': 0.8,
        'grid.alpha': 0.3
    })
    
    # 创建时间戳
    时间戳 = datetime.now().strftime("%Y%m%d_%H%M%S")
    基础名称 = f"PBI_Kr86_E{输入参数['initial_energy']:.0f}MeV_d{输入参数['density']:.2f}gcm3_{时间戳}"
    
    # 1. 沉积能量分布图 (对数坐标)
    fig1, ax1 = plt.subplots(figsize=(10, 8))
    
    # 绘制原始剂量和归一化剂量
    line1 = ax1.loglog(半径数组, 剂量数组, 'b-', linewidth=2.5, label='Original Dose')
    line2 = ax1.loglog(半径数组, 归一化剂量, 'r--', linewidth=2.5, label=f'Normalized Dose (g={输入参数["g_factor"]:.2f})')
    
    ax1.set_xlabel('Radial Distance r (nm)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Dose (keV/nm$^3$)', fontsize=14, fontweight='bold')
    ax1.set_title(f'$^{{86}}$Kr Ion Dose Distribution in PBI\nE$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u\nPBI Density = {输入参数["density"]:.3f} g/cm$^3$', 
                  fontsize=16, fontweight='bold', pad=20)
    
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(loc='upper right', framealpha=0.9)
    
    # 添加关键数值标注
    ax1.text(0.02, 0.95, f'Track Center Dose:\n{剂量数组[0]:.2e} keV/nm$^3$\nMax Dose:\n{np.max(剂量数组):.2e} keV/nm$^3$\nTotal Energy:\n{分析信息["总能量_keV"]:.2e} keV', 
             transform=ax1.transAxes, verticalalignment='top', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, pad=0.5))
    
    plt.tight_layout()
    图形文件1 = os.path.join(输出目录, f"{基础名称}_Dose_LogScale.png")
    plt.savefig(图形文件1, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 2. 沉积能量分布图 (线性坐标)
    fig2, ax2 = plt.subplots(figsize=(10, 8))
    
    # 只显示前100nm的数据以获得更好的可视效果
    mask = 半径数组 <= 10
    
    ax2.plot(半径数组[mask], 剂量数组[mask], 'b-', linewidth=2.5, label='Original Dose')
    ax2.plot(半径数组[mask], 归一化剂量[mask], 'r--', linewidth=2.5, label=f'Normalized Dose (g={输入参数["g_factor"]:.2f})')
    
    ax2.set_xlabel('Radial Distance r (nm)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Dose (keV/nm$^3$)', fontsize=14, fontweight='bold')
    ax2.set_title(f'$^{{86}}$Kr Ion Dose Distribution in PBI (Linear Scale)\nE$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u\nPBI Density = {输入参数["density"]:.3f} g/cm$^3$', 
                  fontsize=16, fontweight='bold', pad=20)
    
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right', framealpha=0.9)
    
    # 添加关键数值标注
    ax2.text(0.97, 0.95, f'Track Center Dose:\n{剂量数组[0]:.2e} keV/nm$^3$\nDose at 10 nm:\n{剂量数组[np.argmin(np.abs(半径数组-10))]:.2e} keV/nm$^3$\nDose at 50 nm:\n{剂量数组[np.argmin(np.abs(半径数组-50))]:.2e} keV/nm$^3$', 
             transform=ax2.transAxes, verticalalignment='top', horizontalalignment='right', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9, pad=0.5))
    
    plt.tight_layout()
    图形文件2 = os.path.join(输出目录, f"{基础名称}_Dose_LinearScale.png")
    plt.savefig(图形文件2, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 3. 温度分布图 (对数坐标)
    fig3, ax3 = plt.subplots(figsize=(10, 8))
    
    ax3.loglog(半径数组, 温度数组, 'g-', linewidth=2.5, label='Temperature')
    
    ax3.set_xlabel('Radial Distance r (nm)', fontsize=14, fontweight='bold')
    ax3.set_ylabel('Temperature T (K)', fontsize=14, fontweight='bold')
    ax3.set_title(f'$^{{86}}$Kr Ion Induced Temperature Distribution in PBI\nE$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u\nPBI Density = {输入参数["density"]:.3f} g/cm$^3$', 
                  fontsize=16, fontweight='bold', pad=20)
    
    ax3.grid(True, which="both", alpha=0.3)
    ax3.legend(loc='upper right', framealpha=0.9)
    
    # 添加温度标准线
    ax3.axhline(y=300, color='orange', linestyle=':', alpha=0.7, linewidth=1.5)
    ax3.axhline(y=1000, color='red', linestyle=':', alpha=0.7, linewidth=1.5)
    ax3.axhline(y=5000, color='darkred', linestyle=':', alpha=0.7, linewidth=1.5)
    
    # 添加温度标准线标签（避免与legend重叠）
    ax3.text(80, 320, 'Room T (300 K)', fontsize=9, color='orange', fontweight='bold')
    ax3.text(80, 1100, 'High T (1000 K)', fontsize=9, color='red', fontweight='bold')
    ax3.text(80, 5500, 'Very High T (5000 K)', fontsize=9, color='darkred', fontweight='bold')
    
    # 添加关键数值标注
    ax3.text(0.02, 0.95, f'Track Center T:\n{分析信息["径迹中心温度_K"]:,.0f} K\nMax Temperature:\n{分析信息["最高温度_K"]:,.0f} K\nAvg Temperature:\n{分析信息["平均温度_K"]:,.0f} K', 
             transform=ax3.transAxes, verticalalignment='top', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9, pad=0.5))
    
    plt.tight_layout()
    图形文件3 = os.path.join(输出目录, f"{基础名称}_Temperature_LogScale.png")
    plt.savefig(图形文件3, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 4. 温度分布图 (线性坐标)
    fig4, ax4 = plt.subplots(figsize=(10, 8))
    
    # 只显示前50nm的数据
    mask = 半径数组 <= 10
    
    ax4.plot(半径数组[mask], 温度数组[mask], 'g-', linewidth=2.5, label='Temperature')
    
    ax4.set_xlabel('Radial Distance r (nm)', fontsize=14, fontweight='bold')
    ax4.set_ylabel('Temperature T (K)', fontsize=14, fontweight='bold')
    ax4.set_title(f'$^{{86}}$Kr Ion Induced Temperature Distribution in PBI (Linear Scale)\nE$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u\nPBI Density = {输入参数["density"]:.3f} g/cm$^3$', 
                  fontsize=16, fontweight='bold', pad=20)
    
    ax4.grid(True, alpha=0.3)
    ax4.legend(loc='upper right', framealpha=0.9)
    
    # 添加温度标准线
    ax4.axhline(y=300, color='orange', linestyle=':', alpha=0.7, linewidth=1.5)
    ax4.axhline(y=1000, color='red', linestyle=':', alpha=0.7, linewidth=1.5)
    if 分析信息["最高温度_K"] > 5000:
        ax4.axhline(y=5000, color='darkred', linestyle=':', alpha=0.7, linewidth=1.5)
    
    # 添加温度标准线标签（避免与数据重叠）
    ax4.text(42, 320, 'Room T (300 K)', fontsize=9, color='orange', fontweight='bold')
    ax4.text(42, 1100, 'High T (1000 K)', fontsize=9, color='red', fontweight='bold')
    if 分析信息["最高温度_K"] > 5000:
        ax4.text(42, 5500, 'Very High T (5000 K)', fontsize=9, color='darkred', fontweight='bold')
    
    # 添加关键数值标注
    ax4.text(0.97, 0.95, f'Track Center T:\n{分析信息["径迹中心温度_K"]:,.0f} K\nT at 10 nm:\n{温度数组[np.argmin(np.abs(半径数组-10))]:,.0f} K\nT at 50 nm:\n{温度数组[np.argmin(np.abs(半径数组-50))]:,.0f} K', 
             transform=ax4.transAxes, verticalalignment='top', horizontalalignment='right', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.9, pad=0.5))
    
    plt.tight_layout()
    图形文件4 = os.path.join(输出目录, f"{基础名称}_Temperature_LinearScale.png")
    plt.savefig(图形文件4, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 5. 综合对比图
    fig5, ((ax5a, ax5b), (ax5c, ax5d)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 子图5a: 剂量分布对数
    ax5a.loglog(半径数组, 剂量数组, 'b-', linewidth=2, label='Original')
    ax5a.loglog(半径数组, 归一化剂量, 'r--', linewidth=2, label='Normalized')
    ax5a.set_xlabel('r (nm)')
    ax5a.set_ylabel('Dose (keV/nm$^3$)')
    ax5a.set_title('Dose Distribution (Log)')
    ax5a.grid(True, which="both", alpha=0.3)
    ax5a.legend()
    
    # 子图5b: 剂量分布线性
    mask = 半径数组 <= 10
    ax5b.plot(半径数组[mask], 剂量数组[mask], 'b-', linewidth=2, label='Original')
    ax5b.plot(半径数组[mask], 归一化剂量[mask], 'r--', linewidth=2, label='Normalized')
    ax5b.set_xlabel('r (nm)')
    ax5b.set_ylabel('Dose (keV/nm$^3$)')
    ax5b.set_title('Dose Distribution (Linear, r ≤ 50 nm)')
    ax5b.grid(True, alpha=0.3)
    ax5b.legend()
    
    # 子图5c: 温度分布对数
    ax5c.loglog(半径数组, 温度数组, 'g-', linewidth=2, label='Temperature')
    ax5c.set_xlabel('r (nm)')
    ax5c.set_ylabel('T (K)')
    ax5c.set_title('Temperature Distribution (Log)')
    ax5c.grid(True, which="both", alpha=0.3)
    ax5c.legend()
    
    # 子图5d: 温度分布线性
    ax5d.plot(半径数组[mask], 温度数组[mask], 'g-', linewidth=2, label='Temperature')
    ax5d.set_xlabel('r (nm)')
    ax5d.set_ylabel('T (K)')
    ax5d.set_title('Temperature Distribution (Linear, r ≤ 50 nm)')
    ax5d.grid(True, alpha=0.3)
    ax5d.legend()
    
    fig5.suptitle(f'$^{{86}}$Kr Ion Energy Deposition and Temperature in PBI\nE$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u, ρ = {输入参数["density"]:.3f} g/cm$^3$', 
                  fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    图形文件5 = os.path.join(输出目录, f"{基础名称}_Combined_Plot.png")
    plt.savefig(图形文件5, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 恢复matplotlib参数
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 10,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14,
        'lines.linewidth': 1.5,
        'axes.linewidth': 0.8,
        'grid.linewidth': 0.5,
        'grid.alpha': 0.3
    })
    
    return [图形文件1, 图形文件2, 图形文件3, 图形文件4, 图形文件5]

def generate_detailed_plots():
    """
    生成详细图表的主函数
    使用默认参数进行计算并生成图表
    """
    # 默认计算参数（根据SRIM数据文件中的PBI参数）
    默认参数 = {
        'initial_energy': 25.0,
        'thickness': 40.0,
        'density': 1.5,  # 根据SRIM文件中的PBI密度
        'ionization_energy': 12.0,
        'g_factor': 0.17,
        'max_radius': 10.0
    }
    
    # 实例化计算器
    srim_processor = SRIM数据处理器()
    wz_calculator = WZ模型计算器()
    
    print("开始计算详细图表...")
    
    # 执行计算
    能量信息 = srim_processor.计算能量损失(默认参数['initial_energy'], 默认参数['thickness'])
    最终能量 = 能量信息['最终能量']
    阻止本领 = 能量信息['阻止本领'] / 1000
    半径数组 = np.logspace(0, np.log10(默认参数['max_radius']), 150)
    
    能量分布信息 = wz_calculator.计算径向能量分布(半径数组, 36, 最终能量, 默认参数['density'], 默认参数['ionization_energy'])
    剂量数组 = 能量分布信息['剂量密度']
    归一化剂量 = 剂量数组 * 默认参数['g_factor'] / 能量分布信息['累积能量'][-1] * 阻止本领
    温度数组 = np.array([wz_calculator.剂量转温度(d, 默认参数['density']) for d in 归一化剂量])
    
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
    
    # 创建输出目录
    输出目录 = r"D:\Desktop\MY Workspace\PBI_irradiation\WZ_AtomicVelocity_Results"
    os.makedirs(输出目录, exist_ok=True)
    
    # 生成图表
    图形文件列表 = plot_dose_temperature_distributions(半径数组, 剂量数组, 归一化剂量, 温度数组, 能量信息, 默认参数, 分析信息, 输出目录)
    
    # 保存详细数据文件
    计算结果 = {
        '半径': 半径数组, '原始剂量': 剂量数组, '归一化剂量': 归一化剂量,
        '温度': 温度数组, '能量信息': 能量信息, '分析信息': 分析信息,
        '参数': 默认参数
    }
    
    时间戳 = datetime.now().strftime("%Y%m%d_%H%M%S")
    数据文件 = os.path.join(输出目录, f"PBI_Kr86_E{默认参数['initial_energy']:.0f}MeV_d{默认参数['density']:.2f}gcm3_{时间戳}_DetailedData.txt")
    save_data_file(数据文件, 计算结果)
    
    # 打印结果
    print("\n=== 详细图表生成完成 ===")
    print(f"径迹中心温度: {分析信息['径迹中心温度_K']:,.0f} K")
    print(f"最高温度: {分析信息['最高温度_K']:,.0f} K")
    print(f"PBI原子数密度: {分析信息['原子数密度_nm3']:.2e} 原子/nm³")
    print(f"总沉积能量: {分析信息['总能量_keV']:.2e} keV")
    print(f"能量损失百分比: {能量信息['损失百分比']:.2f}%")
    
    print("\n生成的图形文件:")
    for 文件 in 图形文件列表:
        print(f"  - {文件}")
    
    print(f"\n生成的数据文件:")
    print(f"  - {数据文件}")
    
    return 图形文件列表

class 原子速度计算器:
    """
    基于等动能原理计算不同原子的速度分布
    考虑初始温度300K的热运动能量
    根据动能相等原理：1/2 * m * v^2 = constant
    因此 v ∝ 1/√m
    """
    
    def __init__(self):
        # 原子质量 (g/mol)，来自PBI模型数据文件
        self.原子质量_g_per_mol = {
            'C': 12.011150,  # 碳原子质量 (g/mol)
            'H': 1.007970,   # 氢原子质量 (g/mol)
            'O': 15.999400,  # 氧原子质量 (g/mol)
            'N': 14.006700   # 氮原子质量 (g/mol)
        }
        
        # 物理常数（带单位标注）
        self.阿伏加德罗常数_per_mol = 6.022e23      # 个/mol
        self.J_per_eV = 1.60218e-19                  # J/eV
        self.kB_eV_per_K = 8.617e-5                  # 玻尔兹曼常数 eV/K
        
        # 单位转换常数
        self.m_per_s_to_angstrom_per_fs = 1e-5       # 1 m/s = 1e-5 Å/fs
        
        # 初始温度
        self.初始温度_K = 300.0                       # 环境温度 300K
    
    def 计算初始热能量(self, 原子类型='C'):
        """
        计算300K时原子的初始热运动能量
        使用3kT/2 (三维平均动能)
        """
        初始热能量_eV = 1.5 * self.kB_eV_per_K * self.初始温度_K
        return 初始热能量_eV
    
    def 计算原子速度(self, 能量_eV, 原子类型='C'):
        """
        根据原子能量计算原子速度
        
        参数:
        - 能量_eV: 原子的动能 (eV)
        - 原子类型: 'C', 'H', 'O', 'N'
        
        返回:
        - 速度 (Å/fs) - Angstroms per femtosecond
        
        计算步骤:
        1. 能量转换: eV → J
        2. 质量转换: g/mol → kg/个
        3. 速度计算: v = sqrt(2*KE/m)
        4. 单位转换: m/s → Å/fs
        """
        if 能量_eV <= 0:
            return 0
        
        # 步骤1: 将能量转换为动能
        动能_J = 能量_eV * self.J_per_eV  # J
        
        # 步骤2: 原子质量转换
        质量_g_per_mol = self.原子质量_g_per_mol[原子类型]  # g/mol
        质量_kg_per_atom = (质量_g_per_mol / 1000.0) / self.阿伏加德罗常数_per_mol  # kg/原子
        
        # 步骤3: 根据动能公式计算速度: KE = 1/2 * m * v^2
        速度_m_per_s = math.sqrt(2 * 动能_J / 质量_kg_per_atom)  # m/s
        
        # 步骤4: 单位转换: m/s → Å/fs
        速度_angstrom_per_fs = 速度_m_per_s * self.m_per_s_to_angstrom_per_fs  # Å/fs
        
        return 速度_angstrom_per_fs
    
    def 计算所有原子速度分布(self, 半径数组_nm, 单原子沉积能量数组_eV):
        """
        计算所有四种原子的速度分布
        包括初速度、总速度和速度增加值
        
        参数:
        - 半径数组_nm: 径向距离数组 (nm)
        - 单原子沉积能量数组_eV: 对应的单原子沉积能量数组 (eV)
        
        返回:
        - 速度分布字典，包含每种原子的详细信息
        """
        速度分布 = {}
        
        for 原子类型 in ['C', 'H', 'O', 'N']:
            # 计算初始热能量和初速度
            初始热能量_eV = self.计算初始热能量(原子类型)
            初速度_angstrom_per_fs = self.计算原子速度(初始热能量_eV, 原子类型)
            
            # 计算各径向位置的总能量和总速度
            总能量数组_eV = 初始热能量_eV + 单原子沉积能量数组_eV
            总速度数组_angstrom_per_fs = np.array([
                self.计算原子速度(总能量, 原子类型) 
                for 总能量 in 总能量数组_eV
            ])
            
            # 计算速度增加值
            速度增加数组_angstrom_per_fs = 总速度数组_angstrom_per_fs - 初速度_angstrom_per_fs
            
            速度分布[原子类型] = {
                '半径_nm': 半径数组_nm,
                '初始热能量_eV': 初始热能量_eV,
                '初速度_angstrom_per_fs': 初速度_angstrom_per_fs,
                '沉积能量数组_eV': 单原子沉积能量数组_eV,
                '总能量数组_eV': 总能量数组_eV,
                '总速度数组_angstrom_per_fs': 总速度数组_angstrom_per_fs,
                '速度增加数组_angstrom_per_fs': 速度增加数组_angstrom_per_fs,
                '质量_g_per_mol': self.原子质量_g_per_mol[原子类型],
                '质量_kg_per_atom': (self.原子质量_g_per_mol[原子类型] / 1000.0) / self.阿伏加德罗常数_per_mol,
                '径迹中心总速度_angstrom_per_fs': 总速度数组_angstrom_per_fs[0] if len(总速度数组_angstrom_per_fs) > 0 else 0,
                '径迹中心速度增加_angstrom_per_fs': 速度增加数组_angstrom_per_fs[0] if len(速度增加数组_angstrom_per_fs) > 0 else 0,
                '最大总速度_angstrom_per_fs': np.max(总速度数组_angstrom_per_fs) if len(总速度数组_angstrom_per_fs) > 0 else 0,
                '最大速度增加_angstrom_per_fs': np.max(速度增加数组_angstrom_per_fs) if len(速度增加数组_angstrom_per_fs) > 0 else 0
            }
        
        return 速度分布

def generate_atomic_velocity_analysis():
    """
    生成原子速度分析的主函数
    使用数据文件中的实际原子数密度计算
    """
    # 默认计算参数
    默认参数 = {
        'initial_energy': 25.0,     # MeV/u
        'thickness': 40.0,          # μm
        'density': 1.5,             # g/cm³
        'ionization_energy': 12.0,  # eV
        'g_factor': 0.17,           # 无量纲
        'max_radius': 50.0          # nm
    }
    
    # 实例化计算器
    srim_processor = SRIM数据处理器()
    wz_calculator = WZ模型计算器()
    velocity_calculator = 原子速度计算器()
    
    print("开始计算原子速度分布...")
    print(f"使用数据文件: PBI_5_864_CHON-2019-innerwall(2).data")
    
    # 打印原子数密度信息
    print(f"\n=== 原子数密度计算 ===")
    print(f"数据文件信息:")
    print(f"  总原子数: {wz_calculator.总原子数}")
    print(f"  盒子体积: {wz_calculator.盒子体积_nm3:.3f} nm³")
    print(f"  原子数密度 (方法1-直接体积): {wz_calculator.实际原子数密度_per_nm3:.6f} atoms/nm³")
    print(f"  平均原子质量: {wz_calculator.平均原子质量_g_per_mol:.6f} g/mol")
    
    # 原子组成信息
    print(f"\n原子组成:")
    for 元素, 数量 in wz_calculator.原子组成.items():
        百分比 = (数量 / wz_calculator.总原子数) * 100
        print(f"  {元素}: {数量} atoms ({百分比:.1f}%)")
    
    # === 完全复制原脚本的计算部分 ===
    能量信息 = srim_processor.计算能量损失(默认参数['initial_energy'], 默认参数['thickness'])
    最终能量 = 能量信息['最终能量']
    阻止本领 = 能量信息['阻止本领'] / 1000  # keV/μm → MeV/μm
    半径数组 = np.logspace(0, np.log10(默认参数['max_radius']), 150)
    
    能量分布信息 = wz_calculator.计算径向能量分布(半径数组, 36, 最终能量, 默认参数['density'], 默认参数['ionization_energy'])
    剂量数组 = 能量分布信息['剂量密度']
    归一化剂量 = 剂量数组 * 默认参数['g_factor'] / 能量分布信息['累积能量'][-1] * 阻止本领
    
    # 使用实际原子数密度计算单原子能量
    单原子能量数组_eV = (归一化剂量 * 1000.0) / wz_calculator.实际原子数密度_per_nm3
    
    # === 完全复制pbi_cli.py中的温度计算方法 ===
    温度数组 = np.array([wz_calculator.剂量转温度(d, 默认参数['density']) for d in 归一化剂量])
    
    分析信息 = {
        #'分子数密度_nm3': 能量分布信息['分子数密度_nm3'],
        '原子数密度_nm3': 能量分布信息['原子数密度_nm3'],
        #'每分子原子数': 能量分布信息['每分子原子数'],
        '径迹中心温度_K': 温度数组[0] if len(温度数组) > 0 else 0,
        '最高温度_K': np.max(温度数组) if len(温度数组) > 0 else 0,
        '平均温度_K': np.mean(温度数组) if len(温度数组) > 0 else 0,
        #'单分子能量数组_eV': 能量分布信息['单分子能量_eV'],
        '单原子能量数组_eV': 能量分布信息['单原子能量_eV'],
        '累积能量': 能量分布信息['累积能量'],
        '总能量_keV': 能量分布信息['总能量_keV']
    }
    
    # === 新增：计算原子速度分布 ===
    速度分布 = velocity_calculator.计算所有原子速度分布(半径数组, 单原子能量数组_eV)
    
    # 打印关键结果到控制台
    print("\n=== 温度计算结果摘要 ===")
    print(f"径迹中心温度: {分析信息['径迹中心温度_K']:,.0f} K")
    print(f"最高温度: {分析信息['最高温度_K']:,.0f} K")
    print(f"平均温度: {分析信息['平均温度_K']:,.0f} K")
    
    print("\n=== 原子速度计算结果摘要 ===")
    print(f"总沉积能量: {能量分布信息['总能量_keV']:.2e} keV")
    print(f"能量损失百分比: {能量信息['损失百分比']:.2f}%")
    print(f"初始环境温度: {velocity_calculator.初始温度_K:.0f} K")
    
    print("\n=== 各原子初始热运动参数 (300K) ===")
    for 原子类型 in ['C', 'H', 'O', 'N']:
        初始热能量 = 速度分布[原子类型]['初始热能量_eV']
        初速度 = 速度分布[原子类型]['初速度_angstrom_per_fs']
        质量_g_mol = 速度分布[原子类型]['质量_g_per_mol']
        print(f"{原子类型}原子: 初始热能量 = {初始热能量:.4f} eV, 初速度 = {初速度:.4f} Å/fs (质量: {质量_g_mol:.6f} g/mol)")
    
    print("\n=== 径迹中心原子速度参数 ===")
    for 原子类型 in ['C', 'H', 'O', 'N']:
        沉积能量 = 速度分布[原子类型]['沉积能量数组_eV'][0]
        总能量 = 速度分布[原子类型]['总能量数组_eV'][0]
        总速度 = 速度分布[原子类型]['径迹中心总速度_angstrom_per_fs']
        速度增加 = 速度分布[原子类型]['径迹中心速度增加_angstrom_per_fs']
        初速度 = 速度分布[原子类型]['初速度_angstrom_per_fs']
        print(f"{原子类型}原子: 沉积能量 = {沉积能量:.2e} eV, 总能量 = {总能量:.2e} eV")
        print(f"       初速度 = {初速度:.4f} Å/fs, 总速度 = {总速度:.4f} Å/fs, 速度增加 = {速度增加:.4f} Å/fs")
    
    print("\n=== 理论速度比验证 (相对于C原子, 基于总速度) ===")
    碳质量 = 速度分布['C']['质量_g_per_mol']
    print("基于 v ∝ 1/√m 关系:")
    for 原子类型 in ['H', 'O', 'N']:
        质量 = 速度分布[原子类型]['质量_g_per_mol']
        理论比 = math.sqrt(碳质量 / 质量)
        实际比_总速度 = 速度分布[原子类型]['径迹中心总速度_angstrom_per_fs'] / 速度分布['C']['径迹中心总速度_angstrom_per_fs']
        实际比_初速度 = 速度分布[原子类型]['初速度_angstrom_per_fs'] / 速度分布['C']['初速度_angstrom_per_fs']
        print(f"v_{原子类型}/v_C: 理论={理论比:.3f}, 实际(总速度)={实际比_总速度:.3f}, 实际(初速度)={实际比_初速度:.3f}")
    
    # 创建图形样式
    颜色样式 = {'C': 'blue', 'H': 'red', 'O': 'green', 'N': 'orange'}
    线型样式 = {'C': '-', 'H': '--', 'O': '-.', 'N': ':'}
    标记样式 = {'C': 'o', 'H': 's', 'O': '^', 'N': 'd'}
    
    输出目录 = r"D:\Desktop\MY Workspace\PBI_irradiation\WZ_AtomicVelocity_Results"
    os.makedirs(输出目录, exist_ok=True)
    时间戳 = datetime.now().strftime("%Y%m%d_%H%M%S")
    基础文件名 = f"PBI_AtomicVelocity_Updated_{时间戳}"
    
    # 设置matplotlib参数
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.titlesize': 18,
        'lines.linewidth': 2.5,
        'axes.linewidth': 1.2,
        'grid.linewidth': 0.8,
        'grid.alpha': 0.3
    })
    
    # 1. 各原子速度分布图 (对数坐标) - 总速度
    plt.figure(figsize=(12, 8))
    
    for 原子类型 in ['C', 'H', 'O', 'N']:
        总速度 = 速度分布[原子类型]['总速度数组_angstrom_per_fs']
        质量 = 速度分布[原子类型]['质量_g_per_mol']
        初速度 = 速度分布[原子类型]['初速度_angstrom_per_fs']
        plt.loglog(半径数组, 总速度, color=颜色样式[原子类型], 
                  linestyle=线型样式[原子类型], linewidth=2.5, 
                  label=f'{原子类型} atom (m = {质量:.3f} g/mol, v₀ = {初速度:.3f} Å/fs)')
    
    plt.xlabel('Radial Distance r (nm)', fontsize=14, fontweight='bold')
    plt.ylabel('Total Atomic Velocity v_total (Å/fs)', fontsize=14, fontweight='bold')
    plt.title('Total Atomic Velocity Distribution in PBI under $^{86}$Kr Irradiation\n' + 
              f'E$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u\n' + 
              f'Initial Temperature = 300 K, Atom Density = {wz_calculator.实际原子数密度_per_nm3:.1f} atoms/nm³', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(loc='upper right', framealpha=0.9)
    
    # 添加关键数值标注
    标注文本 = "Track Center Total Velocities:\n"
    for 原子类型 in ['C', 'H', 'O', 'N']:
        中心总速度 = 速度分布[原子类型]['径迹中心总速度_angstrom_per_fs']
        标注文本 += f"{原子类型}: {中心总速度:.3f} Å/fs\n"
    
    plt.text(0.02, 0.95, 标注文本.strip(), transform=plt.gca().transAxes, 
             verticalalignment='top', fontsize=11,
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9, pad=0.5))
    
    plt.tight_layout()
    图形文件1 = os.path.join(输出目录, f"{基础文件名}_TotalVelocity_LogScale.png")
    plt.savefig(图形文件1, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 2. 各原子速度增加分布图 (对数坐标)
    plt.figure(figsize=(12, 8))
    
    for 原子类型 in ['C', 'H', 'O', 'N']:
        速度增加 = 速度分布[原子类型]['速度增加数组_angstrom_per_fs']
        质量 = 速度分布[原子类型]['质量_g_per_mol']
        plt.loglog(半径数组, 速度增加, color=颜色样式[原子类型], 
                  linestyle=线型样式[原子类型], linewidth=2.5, 
                  label=f'{原子类型} atom (m = {质量:.3f} g/mol)')
    
    plt.xlabel('Radial Distance r (nm)', fontsize=14, fontweight='bold')
    plt.ylabel('Velocity Increase Δv (Å/fs)', fontsize=14, fontweight='bold')
    plt.title('Atomic Velocity Increase Distribution in PBI under $^{86}$Kr Irradiation\n' + 
              f'E$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u\n' + 
              f'Initial Temperature = 300 K', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.grid(True, which="both", alpha=0.3)
    plt.legend(loc='upper right', framealpha=0.9)
    
    # 添加关键数值标注
    标注文本 = "Track Center Velocity Increases:\n"
    for 原子类型 in ['C', 'H', 'O', 'N']:
        中心速度增加 = 速度分布[原子类型]['径迹中心速度增加_angstrom_per_fs']
        标注文本 += f"{原子类型}: {中心速度增加:.3f} Å/fs\n"
    
    plt.text(0.02, 0.95, 标注文本.strip(), transform=plt.gca().transAxes, 
             verticalalignment='top', fontsize=11,
             bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.9, pad=0.5))
    
    plt.tight_layout()
    图形文件2 = os.path.join(输出目录, f"{基础文件名}_VelocityIncrease_LogScale.png")
    plt.savefig(图形文件2, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 3. 各原子总速度分布图 (线性坐标, r ≤ 10 nm)
    plt.figure(figsize=(12, 8))
    
    mask = 半径数组 <= 10
    
    for 原子类型 in ['C', 'H', 'O', 'N']:
        总速度 = 速度分布[原子类型]['总速度数组_angstrom_per_fs']
        质量 = 速度分布[原子类型]['质量_g_per_mol']
        初速度 = 速度分布[原子类型]['初速度_angstrom_per_fs']
        plt.plot(半径数组[mask], 总速度[mask], color=颜色样式[原子类型], 
                linestyle=线型样式[原子类型], linewidth=2.5, 
                marker=标记样式[原子类型], markersize=4, markevery=15,
                label=f'{原子类型} atom (m = {质量:.3f} g/mol, v₀ = {初速度:.3f} Å/fs)')
    
    plt.xlabel('Radial Distance r (nm)', fontsize=14, fontweight='bold')
    plt.ylabel('Total Atomic Velocity v_total (Å/fs)', fontsize=14, fontweight='bold')
    plt.title('Total Atomic Velocity Distribution in PBI (Linear Scale, r ≤ 10 nm)\n' + 
              f'E$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u\n' + 
              f'Initial Temperature = 300 K', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right', framealpha=0.9)
    
    # 添加关键数值标注
    标注文本 = "At Track Center (r = 1 nm):\n"
    for 原子类型 in ['C', 'H', 'O', 'N']:
        中心总速度 = 速度分布[原子类型]['径迹中心总速度_angstrom_per_fs']
        标注文本 += f"{原子类型}: {中心总速度:.3f} Å/fs\n"
    
    标注文本 += f"\nAt r = 10 nm:\n"
    idx_10nm = np.argmin(np.abs(半径数组 - 10))
    for 原子类型 in ['C', 'H', 'O', 'N']:
        速度_10nm = 速度分布[原子类型]['总速度数组_angstrom_per_fs'][idx_10nm]
        标注文本 += f"{原子类型}: {速度_10nm:.3f} Å/fs\n"
    
    plt.text(0.97, 0.95, 标注文本.strip(), transform=plt.gca().transAxes, 
             verticalalignment='top', horizontalalignment='right', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.9, pad=0.5))
    
    plt.tight_layout()
    图形文件3 = os.path.join(输出目录, f"{基础文件名}_TotalVelocity_LinearScale.png")
    plt.savefig(图形文件3, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 4. 各原子速度增加分布图 (线性坐标, r ≤ 10 nm)
    plt.figure(figsize=(12, 8))
    
    mask = 半径数组 <= 10
    
    for 原子类型 in ['C', 'H', 'O', 'N']:
        速度增加 = 速度分布[原子类型]['速度增加数组_angstrom_per_fs']
        质量 = 速度分布[原子类型]['质量_g_per_mol']
        plt.plot(半径数组[mask], 速度增加[mask], color=颜色样式[原子类型], 
                linestyle=线型样式[原子类型], linewidth=2.5, 
                marker=标记样式[原子类型], markersize=4, markevery=15,
                label=f'{原子类型} atom (m = {质量:.3f} g/mol)')
    
    plt.xlabel('Radial Distance r (nm)', fontsize=14, fontweight='bold')
    plt.ylabel('Velocity Increase Δv (Å/fs)', fontsize=14, fontweight='bold')
    plt.title('Atomic Velocity Increase Distribution in PBI (Linear Scale, r ≤ 10 nm)\n' + 
              f'E$_{{initial}}$ = {能量信息["初始能量"]:.1f} MeV/u, E$_{{final}}$ = {能量信息["最终能量"]:.1f} MeV/u\n' + 
              f'Initial Temperature = 300 K', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper right', framealpha=0.9)
    
    # 添加关键数值标注
    标注文本 = "At Track Center (r = 1 nm):\n"
    for 原子类型 in ['C', 'H', 'O', 'N']:
        中心速度增加 = 速度分布[原子类型]['径迹中心速度增加_angstrom_per_fs']
        标注文本 += f"{原子类型}: {中心速度增加:.3f} Å/fs\n"
    
    标注文本 += f"\nAt r = 10 nm:\n"
    idx_10nm = np.argmin(np.abs(半径数组 - 10))
    for 原子类型 in ['C', 'H', 'O', 'N']:
        速度增加_10nm = 速度分布[原子类型]['速度增加数组_angstrom_per_fs'][idx_10nm]
        标注文本 += f"{原子类型}: {速度增加_10nm:.3f} Å/fs\n"
    
    plt.text(0.97, 0.95, 标注文本.strip(), transform=plt.gca().transAxes, 
             verticalalignment='top', horizontalalignment='right', fontsize=10,
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9, pad=0.5))
    
    plt.tight_layout()
    图形文件4 = os.path.join(输出目录, f"{基础文件名}_VelocityIncrease_LinearScale.png")
    plt.savefig(图形文件4, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    # 保存数据文件
    数据文件 = os.path.join(输出目录, f"{基础文件名}_Data.txt")
    with open(数据文件, 'w', encoding='utf-8') as f:
        f.write("# PBI原子速度分布计算结果 (修正版本)\n")
        f.write(f"# 生成时间: {datetime.now()}\n")
        f.write("# 基于等动能原理: 1/2 * m * v^2 = constant, v ∝ 1/√m\n")
        f.write("# 考虑初始环境温度300K的热运动能量\n")
        f.write("# 原子质量来源: PBI模型数据文件 (g/mol)\n")
        f.write("# 数据文件: PBI_5_864_CHON-2019-innerwall(2).data\n")
        f.write("# 密度修正已移除: 实际密度 = SRIM密度 = 1.5 g/cm³\n\n")
        
        f.write("# 原子数密度计算:\n")
        f.write(f"# 原子数密度: {wz_calculator.实际原子数密度_per_nm3:.6f} atoms/nm³\n")
        f.write(f"# 盒子体积: {wz_calculator.盒子体积_nm3:.3f} nm³\n")
        f.write(f"# 总原子数: {wz_calculator.总原子数}\n")
        f.write(f"# 平均原子质量: {wz_calculator.平均原子质量_g_per_mol:.6f} g/mol\n\n")
        
        f.write("# 原子组成:\n")
        for 元素, 数量 in wz_calculator.原子组成.items():
            百分比 = (数量 / wz_calculator.总原子数) * 100
            f.write(f"# {元素}: {数量} atoms ({百分比:.1f}%)\n")
        f.write("\n")
        
        f.write("# SRIM能量损失分析:\n")
        for 键, 值 in 能量信息.items():
            f.write(f"# {键}: {值}\n")
        f.write("\n")
        
        f.write("# 温度分析结果:\n")
        f.write(f"# 初始环境温度: {velocity_calculator.初始温度_K:.0f} K\n")
        f.write(f"# 径迹中心温度: {分析信息['径迹中心温度_K']:,.0f} K\n")
        f.write(f"# 最高温度: {分析信息['最高温度_K']:,.0f} K\n")
        f.write(f"# 平均温度: {分析信息['平均温度_K']:,.0f} K\n")
        f.write(f"# PBI原子数密度: {分析信息['原子数密度_nm3']:.6e} 原子/nm³\n")
        f.write("\n")
        
        f.write("# 原子初始热运动参数 (300K):\n")
        for 原子类型 in ['C', 'H', 'O', 'N']:
            质量_g_mol = 速度分布[原子类型]['质量_g_per_mol']
            质量_kg = 速度分布[原子类型]['质量_kg_per_atom']
            初始热能量 = 速度分布[原子类型]['初始热能量_eV']
            初速度 = 速度分布[原子类型]['初速度_angstrom_per_fs']
            f.write(f"# {原子类型}原子: 质量 = {质量_g_mol:.5f} g/mol = {质量_kg:.6e} kg/个\n")
            f.write(f"#       初始热能量 = {初始热能量:.4f} eV, 初速度 = {初速度:.4f} Å/fs\n")
        f.write("\n")
        
        f.write("# 径迹中心原子速度参数:\n")
        for 原子类型 in ['C', 'H', 'O', 'N']:
            沉积能量 = 速度分布[原子类型]['沉积能量数组_eV'][0]
            总能量 = 速度分布[原子类型]['总能量数组_eV'][0]
            总速度 = 速度分布[原子类型]['径迹中心总速度_angstrom_per_fs']
            速度增加 = 速度分布[原子类型]['径迹中心速度增加_angstrom_per_fs']
            f.write(f"# {原子类型}原子: 沉积能量 = {沉积能量:.6e} eV, 总能量 = {总能量:.6e} eV\n")
            f.write(f"#       总速度 = {总速度:.4f} Å/fs, 速度增加 = {速度增加:.4f} Å/fs\n")
        f.write("\n")
        
        f.write("# 数据列说明:\n")
        f.write("# 径向距离_nm: 从径迹中心的径向距离 (nm)\n")
        f.write("# 单原子沉积能量_eV: 单个原子获得的沉积能量 (eV)\n")
        f.write("# 温度_K: 对应径向距离处的温度 (K)\n")
        f.write("# C初速度_Å_per_fs: 碳原子初速度 (Å/fs)\n")
        f.write("# C总速度_Å_per_fs: 碳原子总速度 (Å/fs)\n")
        f.write("# C速度增加_Å_per_fs: 碳原子速度增加 (Å/fs)\n")
        f.write("# H初速度_Å_per_fs: 氢原子初速度 (Å/fs)\n")
        f.write("# H总速度_Å_per_fs: 氢原子总速度 (Å/fs)\n")
        f.write("# H速度增加_Å_per_fs: 氢原子速度增加 (Å/fs)\n")
        f.write("# O初速度_Å_per_fs: 氧原子初速度 (Å/fs)\n")
        f.write("# O总速度_Å_per_fs: 氧原子总速度 (Å/fs)\n")
        f.write("# O速度增加_Å_per_fs: 氧原子速度增加 (Å/fs)\n")
        f.write("# N初速度_Å_per_fs: 氮原子初速度 (Å/fs)\n")
        f.write("# N总速度_Å_per_fs: 氮原子总速度 (Å/fs)\n")
        f.write("# N速度增加_Å_per_fs: 氮原子速度增加 (Å/fs)\n\n")
        
        f.write("径向距离_nm\t单原子沉积能量_eV\t温度_K\t")
        for 原子类型 in ['C', 'H', 'O', 'N']:
            f.write(f"{原子类型}初速度_Å_per_fs\t{原子类型}总速度_Å_per_fs\t{原子类型}速度增加_Å_per_fs\t")
        f.write("\n")
        
        for i, r in enumerate(半径数组):
            f.write(f"{r:.6e}\t{单原子能量数组_eV[i]:.6e}\t{温度数组[i]:.6e}\t")
            # 写入各原子的初速度、总速度和速度增加
            for 原子类型 in ['C', 'H', 'O', 'N']:
                初速度 = 速度分布[原子类型]['初速度_angstrom_per_fs']
                总速度 = 速度分布[原子类型]['总速度数组_angstrom_per_fs'][i]
                速度增加 = 速度分布[原子类型]['速度增加数组_angstrom_per_fs'][i]
                f.write(f"{初速度:.6e}\t{总速度:.6e}\t{速度增加:.6e}\t")
            f.write("\n")
    
    # 恢复matplotlib参数
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 10,
        'axes.titlesize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14,
        'lines.linewidth': 1.5,
        'axes.linewidth': 0.8,
        'grid.linewidth': 0.5,
        'grid.alpha': 0.3
    })
    
    print(f"\n=== 文件保存完成 ===")
    print("生成的图形文件:")
    print(f"  - {图形文件1} (总速度分布图 - 对数坐标)")
    print(f"  - {图形文件2} (速度增加分布图 - 对数坐标)")
    print(f"  - {图形文件3} (总速度分布图 - 线性坐标)")
    print(f"  - {图形文件4} (速度增加分布图 - 线性坐标)")
    print(f"\n生成的数据文件:")
    print(f"  - {数据文件}")
    
    return [图形文件1, 图形文件2, 图形文件3, 图形文件4]

# 注释掉的多余功能函数
"""
def plot_results(figure, 半径, 剂量, 归一化剂量, 温度, 能量信息, g因子, 输入参数, 分析信息):
    # 绘制计算结果图表
    pass

def save_data_file(文件名, 计算结果):
    # 保存详细的数据文件
    pass

def main():
    # 主函数：解析命令行参数并执行计算
    pass

def plot_dose_temperature_distributions(半径数组, 剂量数组, 归一化剂量, 温度数组, 能量信息, 输入参数, 分析信息, 输出目录):
    # 专门绘制沉积能量和温度分布图的函数
    pass

def generate_detailed_plots():
    # 生成详细图表的主函数
    pass
"""

# 如果直接运行此脚本，生成原子速度分析
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] == '--velocity'):
        # 生成原子速度分析
        generate_atomic_velocity_analysis()
    else:
        print("使用方法: python pbi_atomic_velocity.py [--velocity]")
        print("默认或使用 --velocity 参数运行原子速度分析") 