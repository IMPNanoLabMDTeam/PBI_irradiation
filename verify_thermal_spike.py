#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
验证热峰特征的脚本
检查 irradiated 数据文件中的速度分布是否符合热峰特征
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math

def 读取原子数据(文件路径):
    """读取LAMMPS数据文件中的原子坐标和速度"""
    原子数据 = {}
    速度数据 = {}
    
    with open(文件路径, 'r') as f:
        lines = f.readlines()
    
    # 找到Atoms部分
    atoms_start = None
    atoms_end = None
    velocities_start = None
    velocities_end = None
    
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atoms_start = i + 1
        elif line.strip() == "Velocities":
            velocities_start = i + 1
            atoms_end = i
    
    if atoms_start is None or velocities_start is None:
        raise ValueError("找不到Atoms或Velocities部分")
    
    # 读取原子数据
    for i in range(atoms_start + 1, atoms_end):
        line = lines[i].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 7:
            atom_id = int(parts[0])
            atom_type = int(parts[1])
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            原子数据[atom_id] = {
                'type': atom_type,
                'x': x, 'y': y, 'z': z
            }
    
    # 读取速度数据
    for i in range(velocities_start + 1, len(lines)):
        line = lines[i].strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) >= 4:
            atom_id = int(parts[0])
            vx = float(parts[1])
            vy = float(parts[2])
            vz = float(parts[3])
            速度数据[atom_id] = {
                'vx': vx, 'vy': vy, 'vz': vz
            }
    
    return 原子数据, 速度数据

def 计算径向距离(原子数据):
    """计算每个原子到盒子中心的距离"""
    距离数据 = {}
    
    # 计算盒子中心（假设盒子是立方体）
    x_coords = [atom['x'] for atom in 原子数据.values()]
    y_coords = [atom['y'] for atom in 原子数据.values()]
    z_coords = [atom['z'] for atom in 原子数据.values()]
    
    中心x = (max(x_coords) + min(x_coords)) / 2
    中心y = (max(y_coords) + min(y_coords)) / 2
    中心z = (max(z_coords) + min(z_coords)) / 2
    
    for atom_id, atom in 原子数据.items():
        dx = atom['x'] - 中心x
        dy = atom['y'] - 中心y
        dz = atom['z'] - 中心z
        距离 = math.sqrt(dx*dx + dy*dy + dz*dz)
        距离数据[atom_id] = 距离
    
    return 距离数据

def 分析速度分布(原子数据, 速度数据, 距离数据):
    """分析速度分布的热峰特征"""
    
    # 按距离分组
    距离组 = defaultdict(list)
    for atom_id in 原子数据.keys():
        if atom_id in 速度数据 and atom_id in 距离数据:
            距离 = 距离数据[atom_id]
            vx = 速度数据[atom_id]['vx']
            vy = 速度数据[atom_id]['vy']
            vz = 速度数据[atom_id]['vz']
            速度大小 = math.sqrt(vx*vx + vy*vy + vz*vz)
            
            # 按1nm间隔分组
            组号 = int(距离 / 10.0)  # 10 Å = 1 nm
            距离组[组号].append(速度大小)
    
    # 计算每个距离组的平均速度
    距离数组 = []
    平均速度数组 = []
    速度标准差数组 = []
    
    for 组号 in sorted(距离组.keys()):
        距离_nm = 组号 * 1.0  # 转换为nm
        速度列表 = 距离组[组号]
        
        if len(速度列表) > 0:
            距离数组.append(距离_nm)
            平均速度数组.append(np.mean(速度列表))
            速度标准差数组.append(np.std(速度列表))
    
    return 距离数组, 平均速度数组, 速度标准差数组

def 绘制热峰图(距离数组, 平均速度数组, 速度标准差数组):
    """绘制热峰分布图"""
    plt.figure(figsize=(12, 8))
    
    # 主图：平均速度 vs 距离
    plt.subplot(2, 2, 1)
    plt.errorbar(距离数组, 平均速度数组, yerr=速度标准差数组, 
                fmt='o-', capsize=3, capthick=1, markersize=4)
    plt.xlabel('距离 (nm)')
    plt.ylabel('平均速度 (Å/fs)')
    plt.title('热峰速度分布')
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    # 速度分布直方图
    plt.subplot(2, 2, 2)
    plt.hist(平均速度数组, bins=20, alpha=0.7, color='blue', edgecolor='black')
    plt.xlabel('平均速度 (Å/fs)')
    plt.ylabel('频次')
    plt.title('速度分布直方图')
    plt.grid(True, alpha=0.3)
    
    # 距离分布
    plt.subplot(2, 2, 3)
    plt.hist(距离数组, bins=20, alpha=0.7, color='green', edgecolor='black')
    plt.xlabel('距离 (nm)')
    plt.ylabel('频次')
    plt.title('距离分布')
    plt.grid(True, alpha=0.3)
    
    # 速度 vs 距离散点图
    plt.subplot(2, 2, 4)
    plt.scatter(距离数组, 平均速度数组, alpha=0.6, s=30)
    plt.xlabel('距离 (nm)')
    plt.ylabel('平均速度 (Å/fs)')
    plt.title('速度 vs 距离散点图')
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig('thermal_spike_verification.png', dpi=300, bbox_inches='tight')
    plt.show()

def 主程序():
    """主程序"""
    文件路径 = "PBI_5_864_CHON-2019-innerwall_irradiated.data"
    
    print("=== 热峰特征验证 ===")
    print(f"分析文件: {文件路径}")
    
    # 读取数据
    print("\n1. 读取原子和速度数据...")
    原子数据, 速度数据 = 读取原子数据(文件路径)
    print(f"读取了 {len(原子数据)} 个原子")
    print(f"读取了 {len(速度数据)} 个速度数据")
    
    # 计算距离
    print("\n2. 计算原子到中心的距离...")
    距离数据 = 计算径向距离(原子数据)
    
    # 分析速度分布
    print("\n3. 分析速度分布...")
    距离数组, 平均速度数组, 速度标准差数组 = 分析速度分布(原子数据, 速度数据, 距离数据)
    
    # 输出统计信息
    print(f"\n=== 热峰特征统计 ===")
    print(f"分析的距离组数: {len(距离数组)}")
    print(f"最近距离: {min(距离数组):.2f} nm")
    print(f"最远距离: {max(距离数组):.2f} nm")
    print(f"中心区域平均速度 (<1nm): {平均速度数组[0]:.4f} Å/fs")
    print(f"边缘区域平均速度 (>10nm): {平均速度数组[-1]:.4f} Å/fs")
    print(f"速度比 (中心/边缘): {平均速度数组[0]/平均速度数组[-1]:.2f}")
    
    # 绘制图形
    print("\n4. 生成验证图形...")
    绘制热峰图(距离数组, 平均速度数组, 速度标准差数组)
    
    print("\n验证完成！图形已保存为 thermal_spike_verification.png")

if __name__ == "__main__":
    try:
        主程序()
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc() 