#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PBI速度更新器 - 基于LAMMPS数据文件的原子速度更新
基于SRIM数据和Waligorski-Zhang模型计算沉积能量，然后更新LAMMPS数据文件中的原子速度
使用随机分布确保形成热峰而非冲击波
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os
import random
from datetime import datetime
import argparse
import sys

# 设置matplotlib字体和参数
plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial', 'sans-serif']
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.autolayout'] = True
plt.rcParams['mathtext.fontset'] = 'dejavusans'
plt.rcParams['font.size'] = 10

class SRIM数据处理器:
    """
    SRIM数据处理类
    计算离子在PBI中的能量损失
    """
    
    def __init__(self):
        self.srim阻止本领_keV_per_um = 4.156E+03
        self.pbi密度_g_per_cm3 = 1.5
        
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
        
        # 数据文件中的原子组成
        self.原子组成 = {
            'C': 112320,  # 碳原子数
            'N': 17280,   # 氮原子数  
            'H': 72576,   # 氢原子数
            'O': 6048     # 氧原子数
        }
        # 记录体积和中心位置
        self.盒子体积_nm3 = volume_nm3
        self.盒子中心_angstrom = ((xhi+xlo)/2, (yhi+ylo)/2, (zhi+zlo)/2)
        
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
    
    def 计算径向能量分布(self, 半径数组_nm, Z, 能量_MeV_per_amu, 密度_g_per_cm3, 电离能_eV=12.0):
        """计算径向能量分布"""
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
            '总能量_keV': 累积能量[-1] if len(累积能量) > 0 else 0
        }

class LAMMPS数据文件处理器:
    """
    处理LAMMPS数据文件的读取和写入
    """
    
    def __init__(self):
        self.原子质量_g_per_mol = {
            'C': 12.011150,
            'H': 1.007970,
            'O': 15.999400,
            'N': 14.006700
        }
        # 根据数据文件中的Masses部分创建原子类型映射
        self.原子类型映射 = {
            1: 'C',   # 12.01115
            2: 'N',   # 14.0067
            3: 'N',   # 14.0067
            4: 'C',   # 12.01115
            5: 'C',   # 12.01115
            6: 'H',   # 1.00797
            7: 'H',   # 1.00797
            8: 'O',   # 15.9994
            9: 'N',   # 14.0067
            10: 'C',  # 12.01115
            11: 'O',  # 15.9994
            12: 'O',  # 15.9994
            13: 'H'   # 1.00797
        }
        self.阿伏加德罗常数_per_mol = 6.022e23
        self.J_per_eV = 1.60218e-19
        self.m_per_s_to_angstrom_per_fs = 1e-5
        
    def 读取原子数据(self, 文件路径):
        """读取LAMMPS数据文件中的原子坐标和速度"""
        原子数据 = {}
        速度数据 = {}
        
        with open(文件路径, 'r') as f:
            lines = f.readlines()
        
        # 找到Atoms部分
        atoms_start = None
        for i, line in enumerate(lines):
            if line.strip() == "Atoms # full":
                atoms_start = i + 2  # 跳过空行
                break
        
        if atoms_start is None:
            raise ValueError("找不到Atoms部分")
        
        # 读取原子数据直到遇到空行
        for i in range(atoms_start, len(lines)):
            line = lines[i].strip()
            if not line:  # 空行表示结束
                break
            
            parts = line.split()
            if len(parts) >= 8:
                atom_id = int(parts[0])
                atom_type = int(parts[2])
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                原子数据[atom_id] = {
                    'type': atom_type,
                    'element': self.原子类型映射[atom_type],
                    'x': x, 'y': y, 'z': z
                }
        
        # 找到Velocities部分
        velocities_start = None
        for i, line in enumerate(lines):
            if line.strip() == "Velocities":
                velocities_start = i + 2  # 跳过空行
                break
        
        if velocities_start is None:
            raise ValueError("找不到Velocities部分")
        
        # 读取速度数据
        for i in range(velocities_start, len(lines)):
            line = lines[i].strip()
            if not line:  # 空行表示结束
                break
                
            parts = line.split()
            if len(parts) >= 4:
                atom_id = int(parts[0])
                vx, vy, vz = float(parts[1]), float(parts[2]), float(parts[3])
                速度数据[atom_id] = {'vx': vx, 'vy': vy, 'vz': vz}
        
        print(f"成功读取 {len(原子数据)} 个原子的坐标和 {len(速度数据)} 个原子的速度")
        return 原子数据, 速度数据
    
    def 计算原子到中心距离(self, 原子数据, 中心位置):
        """计算每个原子到盒子中心的距离（单位：nm）"""
        距离数据 = {}
        cx, cy, cz = 中心位置
        
        for atom_id, data in 原子数据.items():
            dx = data['x'] - cx
            dy = data['y'] - cy  
            dz = data['z'] - cz
            距离_angstrom = math.sqrt(dx*dx + dy*dy + dz*dz)
            距离_nm = 距离_angstrom * 0.1  # Å转nm
            距离数据[atom_id] = 距离_nm
            
        return 距离数据
    
    def 速度转能量(self, vx, vy, vz, 原子类型):
        """将速度转换为动能（eV）"""
        速度_angstrom_per_fs = math.sqrt(vx*vx + vy*vy + vz*vz)
        速度_m_per_s = 速度_angstrom_per_fs / self.m_per_s_to_angstrom_per_fs
        
        质量_g_per_mol = self.原子质量_g_per_mol[原子类型]
        质量_kg_per_atom = (质量_g_per_mol / 1000.0) / self.阿伏加德罗常数_per_mol
        
        动能_J = 0.5 * 质量_kg_per_atom * 速度_m_per_s * 速度_m_per_s
        动能_eV = 动能_J / self.J_per_eV
        
        return 动能_eV
    
    def 能量转速度(self, 能量_eV, 原子类型):
        """将动能转换为速度大小（Å/fs）"""
        if 能量_eV <= 0:
            return 0
            
        动能_J = 能量_eV * self.J_per_eV
        质量_g_per_mol = self.原子质量_g_per_mol[原子类型]
        质量_kg_per_atom = (质量_g_per_mol / 1000.0) / self.阿伏加德罗常数_per_mol
        
        速度_m_per_s = math.sqrt(2 * 动能_J / 质量_kg_per_atom)
        速度_angstrom_per_fs = 速度_m_per_s * self.m_per_s_to_angstrom_per_fs
        
        return 速度_angstrom_per_fs
    
    def 生成随机方向(self):
        """生成随机的3D单位方向向量（用于热运动）"""
        # 使用球坐标系生成均匀分布的方向
        phi = random.uniform(0, 2 * math.pi)  # 方位角
        cos_theta = random.uniform(-1, 1)     # cos(极角)的均匀分布
        sin_theta = math.sqrt(1 - cos_theta * cos_theta)
        
        return (
            sin_theta * math.cos(phi),
            sin_theta * math.sin(phi), 
            cos_theta
        )
    
    def 更新速度数据文件(self, 输入文件, 输出文件, 新速度数据):
        """更新LAMMPS数据文件中的Velocities部分"""
        with open(输入文件, 'r') as f:
            lines = f.readlines()
        
        # 找到Velocities部分
        velocities_start = None
        velocities_end = None
        
        for i, line in enumerate(lines):
            if line.strip() == "Velocities":
                velocities_start = i + 1  # Velocities行的下一行
                break
        
        if velocities_start is None:
            raise ValueError("找不到Velocities部分")
        
        # 找到Velocities部分的结束位置
        for i in range(velocities_start + 1, len(lines)):
            line = lines[i].strip()
            if not line:  # 空行
                velocities_end = i
                break
        
        if velocities_end is None:
            velocities_end = len(lines)  # 文件结束
        
        # 准备新的速度行
        新速度行 = []
        新速度行.append("\n")  # 空行
        
        # 按原子ID排序
        sorted_atoms = sorted(新速度数据.keys())
        
        for atom_id in sorted_atoms:
            vx, vy, vz = 新速度数据[atom_id]
            # 使用与原文件相同的格式（非科学记数法）
            新速度行.append(f"{atom_id} {vx:.15f} {vy:.15f} {vz:.15f}\n")
        
        # 构建新文件内容
        新文件内容 = (
            lines[:velocities_start + 1] +  # Velocities之前的所有行包括Velocities行
            新速度行 +                      # 新的速度数据
            lines[velocities_end:]          # Velocities之后的所有行
        )
        
        # 写入新文件
        with open(输出文件, 'w') as f:
            f.writelines(新文件内容)
        
        print(f"成功更新 {len(新速度数据)} 个原子的速度数据到文件: {输出文件}")

def 主程序():
    """主程序：处理整个速度更新流程"""
    
    # 文件路径
    脚本目录 = os.path.dirname(os.path.abspath(__file__))
    原始数据文件 = os.path.join(脚本目录, "PBI_5_864_CHON-2019-innerwall.data")
    输出数据文件 = os.path.join(脚本目录, "PBI_5_864_CHON-2019-innerwall_irradiated.data")
    
    # 计算参数
    计算参数 = {
        'initial_energy': 25.0,
        'thickness': 40.0,
        'density': 1.5,
        'ionization_energy': 12.0,
        'g_factor': 0.70,
        'max_radius': 50.0
    }
    
    print("=== PBI原子速度更新器 ===")
    print(f"输入文件: {原始数据文件}")
    print(f"输出文件: {输出数据文件}")
    
    # 初始化计算器
    srim_processor = SRIM数据处理器()
    wz_calculator = WZ模型计算器()
    lammps_processor = LAMMPS数据文件处理器()
    
    # 1. 读取原始LAMMPS数据文件
    print("\n1. 读取原始LAMMPS数据文件...")
    原子数据, 原始速度数据 = lammps_processor.读取原子数据(原始数据文件)
    
    # 2. 计算每个原子到中心的距离
    print("\n2. 计算原子到盒子中心的距离...")
    距离数据 = lammps_processor.计算原子到中心距离(原子数据, wz_calculator.盒子中心_angstrom)
    
    # 3. 计算WZ模型的能量沉积分布
    print("\n3. 计算WZ模型能量沉积分布...")
    能量信息 = srim_processor.计算能量损失(计算参数['initial_energy'], 计算参数['thickness'])
    最终能量 = 能量信息['最终能量']
    阻止本领 = 能量信息['阻止本领'] / 1000  # keV/μm → MeV/μm
    
    # 创建径向距离数组
    半径数组 = np.logspace(-2, np.log10(计算参数['max_radius']), 150)  # 从0.01nm开始到50nm
    
    能量分布信息 = wz_calculator.计算径向能量分布(
        半径数组, 36, 最终能量, 计算参数['density'], 计算参数['ionization_energy']
    )
    剂量数组 = 能量分布信息['剂量密度']
    归一化剂量 = 剂量数组 * 计算参数['g_factor'] / 能量分布信息['累积能量'][-1] * 阻止本领
    
    # 使用实际原子数密度计算单原子沉积能量
    单原子沉积能量数组_eV = (归一化剂量 * 1000.0) / wz_calculator.实际原子数密度_per_nm3
    
    # 4. 为每个原子计算初始能量和新速度
    print("\n4. 计算每个原子的新速度...")
    新速度数据 = {}
    
    # 创建插值函数用于计算任意距离处的沉积能量
    from scipy import interpolate
    能量插值函数 = interpolate.interp1d(半径数组, 单原子沉积能量数组_eV, 
                                   kind='linear', bounds_error=False, fill_value=0)
    
    统计信息 = {'处理数量': 0, '超出范围': 0}
    
    for atom_id in 原子数据.keys():
        if atom_id not in 原始速度数据:
            continue
            
        # 获取原子信息
        原子类型 = 原子数据[atom_id]['element']
        距离_nm = 距离数据[atom_id]
        原始vx = 原始速度数据[atom_id]['vx']
        原始vy = 原始速度数据[atom_id]['vy']
        原始vz = 原始速度数据[atom_id]['vz']
        
        # 计算初始动能
        初始动能_eV = lammps_processor.速度转能量(原始vx, 原始vy, 原始vz, 原子类型)
        
        # 计算沉积能量
        if 距离_nm <= 计算参数['max_radius']:
            沉积能量_eV = float(能量插值函数(距离_nm))
        else:
            沉积能量_eV = 0
            统计信息['超出范围'] += 1
        
        # 计算总能量和新速度大小
        总能量_eV = 初始动能_eV + 沉积能量_eV
        新速度大小 = lammps_processor.能量转速度(总能量_eV, 原子类型)
        
        # 生成随机方向（热运动特征）
        方向x, 方向y, 方向z = lammps_processor.生成随机方向()
        
        # 计算新的速度分量
        新vx = 新速度大小 * 方向x
        新vy = 新速度大小 * 方向y
        新vz = 新速度大小 * 方向z
        
        新速度数据[atom_id] = (新vx, 新vy, 新vz)
        统计信息['处理数量'] += 1
    
    print(f"处理了 {统计信息['处理数量']} 个原子")
    print(f"超出计算范围的原子: {统计信息['超出范围']} 个")
    
    # 5. 更新输出文件
    print("\n5. 更新LAMMPS数据文件...")
    lammps_processor.更新速度数据文件(原始数据文件, 输出数据文件, 新速度数据)
    
    # 6. 输出统计信息
    print("\n=== 更新完成 ===")
    print(f"径迹中心能量沉积: {单原子沉积能量数组_eV[0]:.6f} eV")
    print(f"50nm处能量沉积: {单原子沉积能量数组_eV[-1]:.6e} eV")
    print(f"总沉积能量: {能量分布信息['总能量_keV']:.2e} keV")
    print(f"能量损失百分比: {能量信息['损失百分比']:.2f}%")
    print("速度分布已设置为随机方向，形成热峰特征")

if __name__ == "__main__":
    try:
        主程序()
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()