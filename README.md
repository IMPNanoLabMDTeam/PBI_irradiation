# PBI Irradiation

本项目是一个科学计算项目，使用Python进行辐照相关的数据分析和可视化。

## 环境要求

- Python >= 3.9
- [uv](https://docs.astral.sh/uv/) - 现代Python包和项目管理工具

## 快速开始

### 1. 安装 uv

如果还没有安装uv，请先安装：

```bash
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# 或者使用 Homebrew (macOS)
brew install uv

# Windows
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

### 2. 克隆项目并复现环境

在新的本地目录中复现项目环境：

```bash
# 克隆项目
git clone <repository-url> PBI_irradiation
cd PBI_irradiation

# 使用 uv 自动创建虚拟环境并安装所有依赖
uv sync
```

就这么简单！`uv sync` 命令会：
- 自动检测并安装正确的Python版本（3.9）
- 创建虚拟环境
- 根据 `uv.lock` 文件安装精确版本的所有依赖包
- 确保环境完全一致

### 3. 激活环境

```bash
# 激活虚拟环境
source .venv/bin/activate

# 或者直接使用 uv 运行命令（推荐）
uv run python main.py
uv run jupyter lab
```

## 项目依赖

本项目主要包含以下核心依赖：

- **numpy**: 科学计算基础库
- **matplotlib**: 数据可视化
- **ipython**: 交互式Python环境
- **ipykernel**: Jupyter内核
- **jupyter**: Jupyter Notebook/Lab环境

完整的依赖列表和版本信息记录在 `uv.lock` 文件中，确保环境的可重现性。

## 常用命令

```bash
# 查看已安装的包
uv pip list

# 添加新的依赖
uv add <package-name>

# 移除依赖
uv remove <package-name>

# 更新依赖
uv sync --upgrade

# 运行Python脚本
uv run python script.py

# 启动Jupyter Lab
uv run jupyter lab

# 启动Jupyter Notebook
uv run jupyter notebook

# 查看项目信息
uv show
```

## 开发工作流

1. **添加新依赖**：
   ```bash
   uv add package-name
   ```

2. **同步环境**（在团队协作时）：
   ```bash
   uv sync
   ```

3. **运行代码**：
   ```bash
   uv run python main.py
   # 或
   uv run jupyter lab
   ```

## 环境说明

- **Python版本**: 3.9+（由`.python-version`文件指定）
- **包管理**: 使用uv进行现代化的依赖管理
- **锁文件**: `uv.lock` 确保所有依赖版本的精确重现
- **配置文件**: `pyproject.toml` 定义项目元数据和依赖

## 故障排除

### Python版本问题
如果遇到Python版本不匹配的问题：
```bash
# uv会自动安装并管理Python版本
uv python install 3.9
uv sync
```

### 依赖冲突
如果遇到依赖冲突：
```bash
# 清理环境重新同步
rm -rf .venv
uv sync
```

### 更新依赖
定期更新依赖到最新兼容版本：
```bash
uv sync --upgrade
```

## 为什么选择 uv？

- **快速**: 比pip快10-100倍
- **可靠**: 自动解决依赖冲突
- **简单**: 零配置，开箱即用
- **现代**: 支持最新的Python打包标准
- **兼容**: 与现有的Python生态系统完全兼容

## 贡献

欢迎提交Issue和Pull Request！

在贡献代码前，请确保：
1. 使用 `uv sync` 同步开发环境
2. 运行测试确保代码质量
3. 更新相关文档
