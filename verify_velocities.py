#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
验证PBI速度更新的效果
检查热峰特征和速度分布
"""

import numpy as np
import matplotlib.pyplot as plt
import math

def read_velocities_and_positions(filename):
    """读取LAMMPS数据文件中的原子位置和速度"""
    atoms = {}
    velocities = {}
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 读取原子坐标
    atoms_start = None
    for i, line in enumerate(lines):
        if line.strip() == "Atoms # full":
            atoms_start = i + 2
            break
    
    if atoms_start:
        for i in range(atoms_start, len(lines)):
            line = lines[i].strip()
            if not line:
                break
            parts = line.split()
            if len(parts) >= 8:
                atom_id = int(parts[0])
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                atoms[atom_id] = (x, y, z)
    
    # 读取速度数据
    velocities_start = None
    for i, line in enumerate(lines):
        if line.strip() == "Velocities":
            velocities_start = i + 1
            break
    
    if velocities_start:
        # 跳过空行直到找到第一个有效的速度数据
        for i in range(velocities_start, len(lines)):
            line = lines[i].strip()
            if line:  # 找到非空行
                parts = line.split()
                if len(parts) >= 4 and parts[0].isdigit():
                    # 从这里开始读取速度数据
                    for j in range(i, len(lines)):
                        line_j = lines[j].strip()
                        if not line_j:  # 遇到空行结束
                            break
                        parts_j = line_j.split()
                        if len(parts_j) >= 4:
                            atom_id = int(parts_j[0])
                            vx, vy, vz = float(parts_j[1]), float(parts_j[2]), float(parts_j[3])
                            velocities[atom_id] = (vx, vy, vz)
                    break
    
    return atoms, velocities

def calculate_distance_from_center(atoms, center):
    """计算原子到中心的距离"""
    distances = {}
    cx, cy, cz = center
    
    for atom_id, (x, y, z) in atoms.items():
        dx, dy, dz = x - cx, y - cy, z - cz
        distance = math.sqrt(dx*dx + dy*dy + dz*dz) * 0.1  # Å到nm
        distances[atom_id] = distance
    
    return distances

def main():
    print("=== PBI速度验证分析 ===")
    
    # 文件路径
    original_file = "PBI_5_864_CHON-2019-innerwall.data"
    irradiated_file = "PBI_5_864_CHON-2019-innerwall_irradiated.data"
    
    # 盒子中心（从WZ计算器中获取）
    xlo, xhi = 299.862329034454, 484.85957681257116
    ylo, yhi = 299.862329034454, 484.85957681257116
    zlo, zhi = 86.0314725411649, 146.2737064987652
    center = ((xhi+xlo)/2, (yhi+ylo)/2, (zhi+zlo)/2)
    
    print("1. 读取原始文件...")
    orig_atoms, orig_velocities = read_velocities_and_positions(original_file)
    print(f"原始文件: {len(orig_atoms)} 个原子坐标, {len(orig_velocities)} 个速度")
    
    print("2. 读取辐射后文件...")
    irrad_atoms, irrad_velocities = read_velocities_and_positions(irradiated_file)
    print(f"辐射后文件: {len(irrad_atoms)} 个原子坐标, {len(irrad_velocities)} 个速度")
    
    print("3. 计算距离...")
    distances = calculate_distance_from_center(orig_atoms, center)
    
    print("4. 分析速度变化...")
    
    # 计算速度大小
    orig_speeds = {}
    irrad_speeds = {}
    
    for atom_id in orig_velocities:
        if atom_id in irrad_velocities:
            vx, vy, vz = orig_velocities[atom_id]
            orig_speeds[atom_id] = math.sqrt(vx*vx + vy*vy + vz*vz)
            
            vx, vy, vz = irrad_velocities[atom_id]
            irrad_speeds[atom_id] = math.sqrt(vx*vx + vy*vy + vz*vz)
    
    print(f"成功比较 {len(orig_speeds)} 个原子的速度")
    
    # 统计分析
    orig_values = list(orig_speeds.values())
    irrad_values = list(irrad_speeds.values())
    
    print("\n=== 速度统计 ===")
    print(f"原始速度:")
    print(f"  平均: {np.mean(orig_values):.6f} Å/fs")
    print(f"  最大: {np.max(orig_values):.6f} Å/fs")
    print(f"  最小: {np.min(orig_values):.6f} Å/fs")
    print(f"  标准差: {np.std(orig_values):.6f} Å/fs")
    
    print(f"\n辐射后速度:")
    print(f"  平均: {np.mean(irrad_values):.6f} Å/fs")
    print(f"  最大: {np.max(irrad_values):.6f} Å/fs")
    print(f"  最小: {np.min(irrad_values):.6f} Å/fs")
    print(f"  标准差: {np.std(irrad_values):.6f} Å/fs")
    
    print(f"\n速度变化:")
    print(f"  平均速度增加: {np.mean(irrad_values) - np.mean(orig_values):.6f} Å/fs")
    print(f"  速度增加倍数: {np.mean(irrad_values) / np.mean(orig_values):.2f}")
    
    # 径向分布分析
    print("\n5. 径向分布分析...")
    
    # 按距离分组
    distance_bins = [0, 1, 2, 5, 10, 20, 50]
    bin_counts = [0] * (len(distance_bins) - 1)
    bin_orig_speeds = [[] for _ in range(len(distance_bins) - 1)]
    bin_irrad_speeds = [[] for _ in range(len(distance_bins) - 1)]
    
    for atom_id in orig_speeds:
        if atom_id in distances:
            dist = distances[atom_id]
            for i in range(len(distance_bins) - 1):
                if distance_bins[i] <= dist < distance_bins[i+1]:
                    bin_counts[i] += 1
                    bin_orig_speeds[i].append(orig_speeds[atom_id])
                    bin_irrad_speeds[i].append(irrad_speeds[atom_id])
                    break
    
    print("\n=== 径向速度分布 ===")
    print("距离范围(nm)     原子数    原始平均速度    辐射后平均速度   速度增加")
    for i in range(len(distance_bins) - 1):
        if bin_counts[i] > 0:
            orig_avg = np.mean(bin_orig_speeds[i])
            irrad_avg = np.mean(bin_irrad_speeds[i])
            speed_increase = irrad_avg - orig_avg
            print(f"{distance_bins[i]:4.0f}-{distance_bins[i+1]:4.0f}        {bin_counts[i]:6d}    "
                  f"{orig_avg:.6f}      {irrad_avg:.6f}     {speed_increase:.6f}")
    
    # 热峰验证
    center_atoms = [atom_id for atom_id, dist in distances.items() if dist < 1.0]
    edge_atoms = [atom_id for atom_id, dist in distances.items() if dist > 10.0]
    
    center_irrad_speeds = [irrad_speeds[atom_id] for atom_id in center_atoms if atom_id in irrad_speeds]
    edge_irrad_speeds = [irrad_speeds[atom_id] for atom_id in edge_atoms if atom_id in irrad_speeds]
    
    if center_irrad_speeds and edge_irrad_speeds:
        center_avg = np.mean(center_irrad_speeds)
        edge_avg = np.mean(edge_irrad_speeds)
        thermal_ratio = center_avg / edge_avg
        
        print(f"\n=== 热峰特征验证 ===")
        print(f"中心区域(<1nm)平均速度: {center_avg:.6f} Å/fs ({len(center_atoms)} 个原子)")
        print(f"边缘区域(>10nm)平均速度: {edge_avg:.6f} Å/fs ({len(edge_atoms)} 个原子)")
        print(f"中心/边缘速度比: {thermal_ratio:.2f}")
        
        # 速度增加量比较
        center_orig_speeds = [orig_speeds[atom_id] for atom_id in center_atoms if atom_id in orig_speeds]
        edge_orig_speeds = [orig_speeds[atom_id] for atom_id in edge_atoms if atom_id in orig_speeds]
        
        if center_orig_speeds and edge_orig_speeds:
            center_orig_avg = np.mean(center_orig_speeds)
            edge_orig_avg = np.mean(edge_orig_speeds)
            center_increase = center_avg - center_orig_avg
            edge_increase = edge_avg - edge_orig_avg
            increase_ratio = center_increase / edge_increase if edge_increase > 0 else float('inf')
            
            print(f"中心区域速度增加: {center_increase:.6f} Å/fs")
            print(f"边缘区域速度增加: {edge_increase:.6f} Å/fs")
            print(f"中心/边缘速度增加比: {increase_ratio:.1f}")
        
        if thermal_ratio > 2.0:
            print("✓ 热峰特征明显：中心区域速度显著高于边缘")
        else:
            print("⚠ 热峰特征不够明显")
    
    print(f"\n=== 随机分布验证 ===")
    # 检查速度方向的随机性
    sample_atoms = list(irrad_velocities.keys())[:1000]  # 取样1000个原子
    vx_values = [irrad_velocities[atom_id][0] for atom_id in sample_atoms]
    vy_values = [irrad_velocities[atom_id][1] for atom_id in sample_atoms]
    vz_values = [irrad_velocities[atom_id][2] for atom_id in sample_atoms]
    
    print(f"样本原子数: {len(sample_atoms)}")
    print(f"vx 平均值: {np.mean(vx_values):.6f} Å/fs (接近0表示随机)")
    print(f"vy 平均值: {np.mean(vy_values):.6f} Å/fs (接近0表示随机)")
    print(f"vz 平均值: {np.mean(vz_values):.6f} Å/fs (接近0表示随机)")
    
    # 检查方向分布
    total_var = np.var(vx_values) + np.var(vy_values) + np.var(vz_values)
    print(f"总方向方差: {total_var:.6f} (高值表示良好的随机分布)")
    
    if abs(np.mean(vx_values)) < 0.001 and abs(np.mean(vy_values)) < 0.001 and abs(np.mean(vz_values)) < 0.001:
        print("✓ 速度方向随机分布良好：各分量平均值接近零")
    else:
        print("⚠ 速度方向可能存在偏向性")
    
    print("\n=== 验证完成 ===")

if __name__ == "__main__":
    main() 