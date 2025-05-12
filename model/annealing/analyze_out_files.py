import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker  # 添加导入ticker模块用于整数坐标
import json  # Add import for JSON parsing

def read_out_file(file_path, timestep=1.0):
    """Read LAMMPS output file and extract step, temperature, and density.
    
    Args:
        file_path: Path to the output file
        timestep: Timestep in femtoseconds (fs) for this force field (default=1.0 for CVFF)
    """
    data = {'step': [], 'time': [], 'temp': [], 'density': []}
    
    with open(file_path, 'r') as f:
        # Skip header line
        next(f)
        for line in f:
            # Split line into columns
            columns = line.strip().split()
            if len(columns) >= 3:
                try:
                    step = int(columns[0])
                    temp = float(columns[1])
                    density = float(columns[2])
                    
                    # Convert step to physical time (fs)
                    time = step * timestep
                    
                    data['step'].append(step)
                    data['time'].append(time)
                    data['temp'].append(temp)
                    data['density'].append(density)
                except (ValueError, IndexError):
                    continue
    
    return data

def ensure_output_dir(output_dir):
    """Create output directory if it doesn't exist."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir

def read_config_file(config_path):
    """Read configuration file that defines output files and their properties.
    
    Args:
        config_path: Path to the JSON configuration file
        
    Returns:
        Dictionary with configuration settings
    """
    try:
        with open(config_path, 'r') as f:
            config = json.load(f)
        return config
    except FileNotFoundError:
        print(f"Configuration file not found: {config_path}")
        return None
    except json.JSONDecodeError:
        print(f"Error parsing JSON in configuration file: {config_path}")
        return None

def combine_continuous_data(files_data):
    """Combine data from multiple files as continuous simulation.
    
    This version uses 'time' as the continuous axis instead of raw steps.
    """
    combined_data = {'step': [], 'time': [], 'temp': [], 'density': []}
    last_time = 0
    
    # 定义力场的排序顺序
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 首先按力场类型，然后按迭代数排序
    def sort_key(item):
        label, _ = item
        parts = label.split('_')
        force_field = parts[-2] if len(parts) > 2 else "cvff"
        iter_num = int(parts[-1])
        return (force_field_order.get(force_field, 999), iter_num)
    
    # 使用自定义排序函数排序
    sorted_files = sorted(files_data.items(), key=sort_key)
    
    for label, data in sorted_files:
        iter_num = int(label.split('_')[-1])
        force_field = label.split('_')[-2] if len(label.split('_')) > 2 else "cvff"
        
        # For steps, we don't adjust them - keep raw steps for reference
        combined_data['step'].extend(data['step'])
        
        # For time, we make it continuous
        adjusted_times = [time + last_time for time in data['time']]
        combined_data['time'].extend(adjusted_times)
        
        # Combine other data
        combined_data['temp'].extend(data['temp'])
        combined_data['density'].extend(data['density'])
        
        # Update last time for next iteration
        if adjusted_times:
            last_time = adjusted_times[-1]
    
    return combined_data

def plot_temperature(files_data, output_dir='plots'):
    """Plot temperature vs time for all files."""
    output_path = os.path.join(output_dir, 'temperature_plot.png')
    plt.figure(figsize=(12, 8))
    
    # Get combined data for continuous simulation
    combined_data = combine_continuous_data(files_data)
    
    # Plot combined temperature data
    plt.plot(combined_data['time'], combined_data['temp'], label='Continuous Simulation', color='blue')
    
    # 定义力场的排序顺序
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 自定义排序函数
    def sort_key(item):
        label, _ = item
        parts = label.split('_')
        force_field = parts[-2] if len(parts) > 2 else "cvff"
        iter_num = int(parts[-1])
        return (force_field_order.get(force_field, 999), iter_num)
    
    # Add vertical lines to separate iterations or force fields
    last_time = 0
    sorted_files = sorted(files_data.items(), key=sort_key)
    for label, data in sorted_files:
        if last_time > 0:
            plt.axvline(x=last_time, linestyle='--', color='gray', alpha=0.7)
            plt.text(last_time, plt.ylim()[1]*0.95, 
                    f"{label.split('_')[-2]}_{label.split('_')[-1]}", 
                    horizontalalignment='right', verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', alpha=0.2))
        last_time += data['time'][-1]
    
    plt.xlabel('Simulation Time (fs)', fontsize=14)
    plt.ylabel('Temperature (K)', fontsize=14)
    plt.title('Temperature Evolution - Continuous Simulation', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def plot_density(files_data, output_dir='plots'):
    """Plot density vs time for all files."""
    output_path = os.path.join(output_dir, 'density_plot.png')
    plt.figure(figsize=(12, 8))
    
    # Get combined data for continuous simulation
    combined_data = combine_continuous_data(files_data)
    
    # Plot combined density data
    plt.plot(combined_data['time'], combined_data['density'], label='Continuous Simulation', color='green')
    
    # 定义力场的排序顺序
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 自定义排序函数
    def sort_key(item):
        label, _ = item
        parts = label.split('_')
        force_field = parts[-2] if len(parts) > 2 else "cvff"
        iter_num = int(parts[-1])
        return (force_field_order.get(force_field, 999), iter_num)
    
    # Add vertical lines to separate iterations or force fields
    last_time = 0
    sorted_files = sorted(files_data.items(), key=sort_key)
    for label, data in sorted_files:
        if last_time > 0:
            plt.axvline(x=last_time, linestyle='--', color='gray', alpha=0.7)
            plt.text(last_time, plt.ylim()[1]*0.95, 
                    f"{label.split('_')[-2]}_{label.split('_')[-1]}", 
                    horizontalalignment='right', verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', alpha=0.2))
        last_time += data['time'][-1]
    
    plt.xlabel('Simulation Time (fs)', fontsize=14)
    plt.ylabel('Density (g/cm³)', fontsize=14)
    plt.title('Density Evolution - Continuous Simulation', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def plot_phase_diagram(files_data, output_dir='plots'):
    """Plot phase diagram (temperature vs density)."""
    output_path = os.path.join(output_dir, 'phase_diagram.png')
    plt.figure(figsize=(12, 8))
    
    # Get combined data for continuous simulation
    combined_data = combine_continuous_data(files_data)
    
    # Define colormap based on force field type
    cmap = plt.cm.viridis
    
    # Plot temperature vs density with color representing time
    scatter = plt.scatter(combined_data['density'], combined_data['temp'], 
               c=combined_data['time'], cmap=cmap, 
               alpha=0.7, s=5)
    
    # Add colorbar to show time progression
    cbar = plt.colorbar(scatter)
    cbar.set_label('Simulation Time (fs)', fontsize=12)
    
    plt.xlabel('Density (g/cm³)', fontsize=14)
    plt.ylabel('Temperature (K)', fontsize=14)
    plt.title('Phase Diagram - Temperature vs Density', fontsize=16)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def plot_temperature_heatmap(files_data, output_dir='plots'):
    """Create a heatmap of temperature evolution over time."""
    output_path = os.path.join(output_dir, 'temperature_heatmap.png')
    
    # Get combined data
    combined_data = combine_continuous_data(files_data)
    
    # Create a 2D array representing temperature over time
    # First bin the data to create a more manageable array
    num_bins = 100  # Number of time bins
    max_temp = max(combined_data['temp'])
    min_temp = min(combined_data['temp'])
    temp_bins = 50  # Number of temperature bins
    
    # Create histogram2d
    h, xedges, yedges = np.histogram2d(
        combined_data['time'], 
        combined_data['temp'],
        bins=[num_bins, temp_bins],
        range=[[min(combined_data['time']), max(combined_data['time'])], 
               [min_temp, max_temp]]
    )
    
    # Plot heatmap
    plt.figure(figsize=(14, 8))
    plt.imshow(h.T, aspect='auto', origin='lower', 
              extent=[min(combined_data['time']), max(combined_data['time']), min_temp, max_temp],
              cmap='plasma', interpolation='bilinear')
    
    plt.colorbar(label='Frequency')
    
    # 定义力场的排序顺序
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 自定义排序函数
    def sort_key(item):
        label, _ = item
        parts = label.split('_')
        force_field = parts[-2] if len(parts) > 2 else "cvff"
        iter_num = int(parts[-1])
        return (force_field_order.get(force_field, 999), iter_num)
    
    # Add vertical lines to separate iterations or force fields
    last_time = 0
    sorted_files = sorted(files_data.items(), key=sort_key)
    for label, data in sorted_files:
        if last_time > 0:
            plt.axvline(x=last_time, linestyle='--', color='white', alpha=0.5)
            force_field = label.split('_')[-2] if len(label.split('_')) > 2 else "cvff"
            iter_num = label.split('_')[-1]
            plt.text(last_time, max_temp*0.95, f"{force_field}_{iter_num}", 
                    color='white', horizontalalignment='right', verticalalignment='top')
        last_time += data['time'][-1]
    
    plt.ylabel('Temperature (K)', fontsize=14)
    plt.xlabel('Simulation Time (fs)', fontsize=14)
    plt.title('Temperature Distribution over Time', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def plot_combined_stats(files_data, output_dir='plots'):
    """Create a combined plot showing temperature and density statistics."""
    output_path = os.path.join(output_dir, 'combined_stats.png')
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), sharex=True)
    
    # Get combined data for continuous simulation
    combined_data = combine_continuous_data(files_data)
    
    # Plot temperature data on top subplot
    ax1.plot(combined_data['time'], combined_data['temp'], label='Temperature', color='red')
    ax1.set_ylabel('Temperature (K)', fontsize=14)
    ax1.set_title('Temperature Evolution - Continuous Simulation', fontsize=16)
    ax1.grid(True, alpha=0.3)
    
    # 定义力场的排序顺序
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 自定义排序函数
    def sort_key(item):
        label, _ = item
        parts = label.split('_')
        force_field = parts[-2] if len(parts) > 2 else "cvff"
        iter_num = int(parts[-1])
        return (force_field_order.get(force_field, 999), iter_num)
    
    # Add vertical lines to separate iterations or force fields
    last_time = 0
    sorted_files = sorted(files_data.items(), key=sort_key)
    for label, data in sorted_files:
        if last_time > 0:
            ax1.axvline(x=last_time, linestyle='--', color='gray', alpha=0.7)
            force_field = label.split('_')[-2] if len(label.split('_')) > 2 else "cvff"
            iter_num = label.split('_')[-1]
            ax1.text(last_time, ax1.get_ylim()[1]*0.95, f"{force_field}_{iter_num}", 
                    horizontalalignment='right', verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', alpha=0.2))
        last_time += data['time'][-1]
    
    # Plot density data on bottom subplot
    ax2.plot(combined_data['time'], combined_data['density'], label='Density', color='green')
    
    # Add vertical lines again for the density plot
    last_time = 0
    for label, data in sorted_files:
        if last_time > 0:
            ax2.axvline(x=last_time, linestyle='--', color='gray', alpha=0.7)
            force_field = label.split('_')[-2] if len(label.split('_')) > 2 else "cvff"
            iter_num = label.split('_')[-1]
            ax2.text(last_time, ax2.get_ylim()[1]*0.95, f"{force_field}_{iter_num}", 
                    horizontalalignment='right', verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', alpha=0.2))
        last_time += data['time'][-1]
    
    ax2.set_xlabel('Simulation Time (fs)', fontsize=14)
    ax2.set_ylabel('Density (g/cm³)', fontsize=14)
    ax2.set_title('Density Evolution - Continuous Simulation', fontsize=16)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def plot_convergence(files_data, output_dir='plots'):
    """Plot convergence of temperature and density over iterations and force fields."""
    output_path = os.path.join(output_dir, 'convergence_plot.png')
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 8))
    
    # Calculate statistics for each iteration and force field
    iterations = []
    labels = []
    avg_temps = []
    final_temps = []
    avg_densities = []
    final_densities = []
    
    # 定义力场的排序顺序，用于可能的额外文件
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 手动定义文件顺序
    desired_order = [
        'PBI_5_75_cvff_0', 
        'PBI_5_75_cvff_1', 
        'PBI_5_75_cvff_2', 
        'PBI_5_75_CHON-2019-innerwall_0'
    ]
    
    # 确保所有文件都存在
    available_files = list(files_data.keys())
    ordered_files = []
    
    # 按照期望顺序排列文件
    for file_name in desired_order:
        if file_name in available_files:
            ordered_files.append((file_name, files_data[file_name]))
    
    # 如果有额外的文件，按照常规排序添加
    def sort_key(item):
        label = item
        parts = label.split('_')
        force_field = parts[-2] if len(parts) > 2 else "cvff"
        iter_num = int(parts[-1])
        return (force_field_order.get(force_field, 999), iter_num)
    
    for file_name in sorted(available_files, key=sort_key):
        if file_name not in desired_order:
            ordered_files.append((file_name, files_data[file_name]))
    
    # 处理每个文件
    x_positions = list(range(len(ordered_files)))
    
    for i, (label, data) in enumerate(ordered_files):
        force_field = label.split('_')[-2] if len(label.split('_')) > 2 else "cvff"
        iter_num = label.split('_')[-1]
        
        # 使用预定义的位置
        x_pos = x_positions[i]
        iterations.append(x_pos)
        labels.append(f"{force_field}_{iter_num}")
        
        # Calculate average and final values
        avg_temps.append(np.mean(data['temp']))
        final_temps.append(data['temp'][-1])
        avg_densities.append(np.mean(data['density']))
        final_densities.append(data['density'][-1])
    
    # Plot temperature convergence
    ax1.plot(iterations, avg_temps, 'o-', label='Average Temperature', color='red')
    ax1.plot(iterations, final_temps, 's-', label='Final Temperature', color='darkred')
    ax1.set_xticks(iterations)
    ax1.set_xticklabels(labels, rotation=45, ha='right')
    ax1.set_xlabel('Iteration and Force Field', fontsize=14)
    ax1.set_ylabel('Temperature (K)', fontsize=14)
    ax1.set_title('Temperature Convergence', fontsize=16)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot density convergence
    ax2.plot(iterations, avg_densities, 'o-', label='Average Density', color='green')
    ax2.plot(iterations, final_densities, 's-', label='Final Density', color='darkgreen')
    ax2.set_xticks(iterations)
    ax2.set_xticklabels(labels, rotation=45, ha='right')
    ax2.set_xlabel('Iteration and Force Field', fontsize=14)
    ax2.set_ylabel('Density (g/cm³)', fontsize=14)
    ax2.set_title('Density Convergence', fontsize=16)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def plot_iteration_comparison(files_data, output_dir='plots'):
    """比较不同迭代和力场中的温度和密度分布"""
    output_path = os.path.join(output_dir, 'iteration_comparison.png')
    
    # 创建图表
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 定义力场的排序顺序
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 自定义排序函数
    def sort_key(item):
        label, _ = item
        parts = label.split('_')
        force_field = parts[-2] if len(parts) > 2 else "cvff"
        iter_num = int(parts[-1])
        return (force_field_order.get(force_field, 999), iter_num)
    
    # 整理每次迭代的数据
    sorted_files = sorted(files_data.items(), key=sort_key)
    
    # 为每次迭代绘制箱线图数据
    temp_data = []
    density_data = []
    iter_labels = []
    
    for label, data in sorted_files:
        force_field = label.split('_')[-2] if len(label.split('_')) > 2 else "cvff"
        iter_num = label.split('_')[-1]
        iter_labels.append(f"{force_field}_{iter_num}")
        temp_data.append(data['temp'])
        density_data.append(data['density'])
    
    # 绘制温度箱线图
    ax1.boxplot(temp_data, labels=iter_labels, patch_artist=True,
               boxprops=dict(facecolor='lightblue', color='blue'),
               whiskerprops=dict(color='blue'),
               capprops=dict(color='blue'),
               medianprops=dict(color='darkblue'))
    
    ax1.set_xticklabels(iter_labels, rotation=45, ha='right')
    ax1.set_xlabel('Iteration and Force Field', fontsize=14)
    ax1.set_ylabel('Temperature (K)', fontsize=14)
    ax1.set_title('Temperature Distribution by Iteration and Force Field', fontsize=16)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # 绘制密度箱线图
    ax2.boxplot(density_data, labels=iter_labels, patch_artist=True,
               boxprops=dict(facecolor='lightgreen', color='green'),
               whiskerprops=dict(color='green'),
               capprops=dict(color='green'),
               medianprops=dict(color='darkgreen'))
    
    ax2.set_xticklabels(iter_labels, rotation=45, ha='right')
    ax2.set_xlabel('Iteration and Force Field', fontsize=14)
    ax2.set_ylabel('Density (g/cm³)', fontsize=14)
    ax2.set_title('Density Distribution by Iteration and Force Field', fontsize=16)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def plot_force_field_comparison(files_data, output_dir='plots'):
    """对比不同力场的温度和密度数据"""
    output_path = os.path.join(output_dir, 'force_field_comparison.png')
    
    # 创建图表
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 按力场类型分组数据
    force_field_data = {}
    
    for label, data in files_data.items():
        force_field = label.split('_')[-2] if len(label.split('_')) > 2 else "cvff"
        
        if force_field not in force_field_data:
            force_field_data[force_field] = {'temp': [], 'density': []}
        
        force_field_data[force_field]['temp'].extend(data['temp'])
        force_field_data[force_field]['density'].extend(data['density'])
    
    # 准备箱线图数据
    temp_data = []
    density_data = []
    ff_labels = []
    
    # 定义力场的排序顺序
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 按指定顺序排序力场
    sorted_ff = sorted(force_field_data.items(), key=lambda x: force_field_order.get(x[0], 999))
    
    for ff, data in sorted_ff:
        temp_data.append(data['temp'])
        density_data.append(data['density'])
        ff_labels.append(ff)
    
    # 绘制温度箱线图
    ax1.boxplot(temp_data, labels=ff_labels, patch_artist=True,
               boxprops=dict(facecolor='lightblue', color='blue'),
               whiskerprops=dict(color='blue'),
               capprops=dict(color='blue'),
               medianprops=dict(color='darkblue'))
    
    ax1.set_xlabel('Force Field', fontsize=14)
    ax1.set_ylabel('Temperature (K)', fontsize=14)
    ax1.set_title('Temperature Distribution by Force Field', fontsize=16)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # 绘制密度箱线图
    ax2.boxplot(density_data, labels=ff_labels, patch_artist=True,
               boxprops=dict(facecolor='lightgreen', color='green'),
               whiskerprops=dict(color='green'),
               capprops=dict(color='green'),
               medianprops=dict(color='darkgreen'))
    
    ax2.set_xlabel('Force Field', fontsize=14)
    ax2.set_ylabel('Density (g/cm³)', fontsize=14)
    ax2.set_title('Density Distribution by Force Field', fontsize=16)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def plot_time_comparison(files_data, output_dir='plots'):
    """比较不同力场的模拟时间与物理特性的关系"""
    output_path = os.path.join(output_dir, 'time_comparison.png')
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # 获取combined数据
    combined_data = combine_continuous_data(files_data)
    
    # 计算每个物理时间点对应的温度变化率
    times = np.array(combined_data['time'])
    temps = np.array(combined_data['temp'])
    
    # 计算差分（温度变化率）
    temp_rates = np.gradient(temps, times)
    
    # 绘制温度变化率
    ax.plot(times, temp_rates, color='purple', label='Temperature Change Rate')
    
    # 定义力场的排序顺序
    force_field_order = {
        'cvff': 0,               # CVFF排在前面
        'CHON-2019-innerwall': 1 # ReaxFF排在后面
    }
    
    # 自定义排序函数
    def sort_key(item):
        label, _ = item
        parts = label.split('_')
        force_field = parts[-2] if len(parts) > 2 else "cvff"
        iter_num = int(parts[-1])
        return (force_field_order.get(force_field, 999), iter_num)
    
    # 添加垂直线分隔不同的迭代和力场
    last_time = 0
    sorted_files = sorted(files_data.items(), key=sort_key)
    for label, data in sorted_files:
        if last_time > 0:
            ax.axvline(x=last_time, linestyle='--', color='gray', alpha=0.7)
            force_field = label.split('_')[-2] if len(label.split('_')) > 2 else "cvff"
            iter_num = label.split('_')[-1]
            ax.text(last_time, ax.get_ylim()[1]*0.95, f"{force_field}_{iter_num}", 
                    horizontalalignment='right', verticalalignment='top',
                    bbox=dict(boxstyle='round,pad=0.3', alpha=0.2))
        last_time += data['time'][-1]
    
    ax.set_xlabel('Simulation Time (fs)', fontsize=14)
    ax.set_ylabel('Temperature Change Rate (K/fs)', fontsize=14)
    ax.set_title('Temperature Change Rate Over Simulation Time', fontsize=16)
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return output_path

def main():
    # Get script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Path to configuration file
    config_path = os.path.join(script_dir, 'analysis_config.json')
    
    # Read configuration
    config = read_config_file(config_path)
    if not config:
        print("Creating a default configuration file. Please edit it and run the script again.")
        default_config = {
            "output_dir": "plots",
            "force_fields": {
                "cvff": {"timestep": 1.0},
                "CHON-2019-innerwall": {"timestep": 0.25}
            },
            "files": [
                {"path": "PBI_5_75_cvff_0.out", "force_field": "cvff"},
                {"path": "PBI_5_75_cvff_1.out", "force_field": "cvff"},
                {"path": "PBI_5_75_cvff_2.out", "force_field": "cvff"},
                {"path": "PBI_5_75_CHON-2019-innerwall_0.out", "force_field": "CHON-2019-innerwall"}
            ]
        }
        with open(config_path, 'w') as f:
            json.dump(default_config, f, indent=4)
        return
    
    # Create output directory as specified in configuration
    output_dir = os.path.join(script_dir, config.get("output_dir", "plots"))
    ensure_output_dir(output_dir)
    
    # Get timestep mapping from configuration
    force_fields = config.get("force_fields", {})
    timestep_map = {ff: data.get("timestep", 1.0) for ff, data in force_fields.items()}
    
    # Read data from all files specified in configuration
    files_data = {}
    for file_info in config.get("files", []):
        file_path = file_info.get("path")
        ff_type = file_info.get("force_field")
        
        # Resolve relative paths
        if not os.path.isabs(file_path):
            file_path = os.path.join(script_dir, file_path)
        
        timestep = timestep_map.get(ff_type, 1.0)  # Default to 1.0 fs/step
        
        if os.path.exists(file_path):
            label = os.path.basename(file_path).replace('.out', '')
            files_data[label] = read_out_file(file_path, timestep)
            print(f"读取文件: {file_path} (时间步长: {timestep} fs/step)")
        else:
            print(f"警告: 文件不存在 {file_path}")
    
    # 如果没有找到任何数据，提前退出
    if not files_data:
        print("错误: 未找到任何输出文件数据，请检查文件路径")
        return
    
    # Generate plots
    plot_paths = []
    
    # Standard plots with time on x-axis
    plot_paths.append(plot_temperature(files_data, output_dir))
    plot_paths.append(plot_density(files_data, output_dir))
    plot_paths.append(plot_combined_stats(files_data, output_dir))
    
    # Additional plots for comprehensive analysis
    plot_paths.append(plot_phase_diagram(files_data, output_dir))
    plot_paths.append(plot_convergence(files_data, output_dir))
    plot_paths.append(plot_iteration_comparison(files_data, output_dir))
    plot_paths.append(plot_force_field_comparison(files_data, output_dir))
    plot_paths.append(plot_time_comparison(files_data, output_dir))
    
    try:
        plot_paths.append(plot_temperature_heatmap(files_data, output_dir))
    except Exception as e:
        print(f"Failed to create temperature heatmap: {e}")
        print("Continuing with other plots.")
    
    print(f"Analysis complete. Generated plots in: {output_dir}")
    for path in plot_paths:
        print(f"- {os.path.basename(path)}")

if __name__ == "__main__":
    main() 