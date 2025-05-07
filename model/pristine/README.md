# LAMMPS Data File Cleaner

此脚本从LAMMPS数据文件中删除力场和拓扑信息，同时保留必要的结构信息。

## 删除内容

脚本从LAMMPS数据文件中删除以下部分：
- Pair Coeffs（对势系数）
- Bond Coeffs（键系数）
- Angle Coeffs（角度系数）
- Dihedral Coeffs（二面角系数）
- Improper Coeffs（平面外系数）
- Bonds（键）
- Angles（角度）
- Dihedrals（二面角）
- Impropers（平面外角）
- 所有有关bonds、angles、dihedrals、impropers的类型和计数行

## 保留内容

脚本保留：
- 头部信息
- Masses（原子质量）
- Atoms（原子坐标）

## 使用方法

```
python clean_data.py <输出文件名> [输入文件1] [输入文件2] ...
```

### 示例：

处理特定文件：
```
python clean_data.py structure_clean.data PBI_5_msi2lmp.data
```

处理多个文件并输出到相同的文件名：
```
python clean_data.py structure_clean.data PBI_5_msi2lmp.data PBI_10_msi2lmp.data
```

如果未指定输入文件，它将处理目录中第一个不以"clean_"开头的.data文件。

## 附加功能

- 在输出文件中添加注释，指示它已被清理
- 允许自定义输出文件名

---

# LAMMPS Data File Cleaner

This script removes force field and topology information from LAMMPS data files while preserving necessary structural information.

## What does it remove?

The script removes the following sections from LAMMPS data files:
- Pair Coeffs
- Bond Coeffs
- Angle Coeffs
- Dihedral Coeffs
- Improper Coeffs
- Bonds
- Angles
- Dihedrals
- Impropers
- All type and count lines related to bonds, angles, dihedrals, impropers

## What does it keep?

The script preserves:
- Header information
- Masses (atom masses)
- Atoms (atom coordinates)

## Usage

```
python clean_data.py <output_filename> [input_file1] [input_file2] ...
```

### Examples:

Process a specific file:
```
python clean_data.py structure_clean.data PBI_5_msi2lmp.data
```

Process multiple files and output to the same filename:
```
python clean_data.py structure_clean.data PBI_5_msi2lmp.data PBI_10_msi2lmp.data
```

If no input files are specified, it will process the first .data file in the directory that doesn't start with 'clean_'.

## Additional Features

- Adds a comment to the output file indicating it has been cleaned
- Allows custom output filename 