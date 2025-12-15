# GEO 数据挖掘流水线
# GEO Data Mining Pipeline

一个自动化的元数据挖掘和过滤流水线，用于从 NCBI GEO 数据库中系统地查询、解析和过滤特定研究目标的样本数据，并关联到 SRA 原始测序数据。

An automated metadata mining and filtering pipeline for systematically querying, parsing, and filtering samples from the NCBI GEO database based on specific research objectives, with linkage to SRA raw sequencing data.

## 📋 项目概述 | Project Overview

该流水线是一个通用的 GEO 数据挖掘工具，可根据不同的研究需求进行配置和定制。支持的测序技术包括：

This pipeline is a general-purpose GEO data mining tool that can be configured and customized for different research needs. Supported sequencing technologies include:

- 单细胞RNA测序 (scRNA-seq)
- 单核RNA测序 (snRNA-seq)  
- 空间转录组 (Spatial Transcriptomics)
- ATAC-seq
- 其他高通量测序技术

### 应用示例 | Application Examples

本项目已成功应用于以下研究场景：

This project has been successfully applied to the following research scenarios:

1. **肺部转移瘤数据挖掘** - 识别其他原发部位转移到肺部的肿瘤样本
2. **小鼠骨髓B细胞发育研究** - 挖掘B细胞发育相关的单细胞数据
3. **其他定制化研究** - 可根据研究需求配置搜索策略和过滤规则

## 🔄 流水线架构 | Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────┐
│  阶段 1: 广泛搜索 (Stage 1: Broad Search)                    │
│  - 使用复杂布尔逻辑搜索 GEO 数据库                           │
│  - 获取潜在相关的 GSE 数据集列表                             │
└──────────────────┬──────────────────────────────────────────┘
                   ↓
┌─────────────────────────────────────────────────────────────┐
│  阶段 2: 深度解析与精准过滤 (Stage 2: Deep Parsing)          │
│  - 下载完整 SOFT 元数据文件                                  │
│  - 样本级 (GSM) 启发式规则过滤                               │
│  - 根据研究目标应用定制化过滤规则                            │
└──────────────────┬──────────────────────────────────────────┘
                   ↓
┌─────────────────────────────────────────────────────────────┐
│  阶段 3: 关联原始数据 (Stage 3: Link to SRA)                │
│  - 将符合条件的 GSM 链接到 SRA                               │
│  - 获取 SRX (实验) 和 SRR (运行) 编号                        │
└──────────────────┬──────────────────────────────────────────┘
                   ↓
┌─────────────────────────────────────────────────────────────┐
│  阶段 4: 数据下载 (Stage 4: Data Download)                  │
│  - 使用 SRA Toolkit 下载 FASTQ 文件                          │
│  - 批量下载脚本                                              │
└─────────────────────────────────────────────────────────────┘
```

## 📁 项目结构 | Project Structure

```
GeneExpressionOmnibus_E-utilities/
├── scripts/          # Python脚本和Shell脚本
├── config/           # 配置文件
├── results/          # 挖掘结果CSV文件
├── docs/             # 项目文档和报告
├── logs/             # 运行日志文件
├── data/             # 数据文件（SRR列表等）
├── GEO_Cache/        # GEO数据缓存
├── README.md         # 本文件
├── requirements.txt  # Python依赖包列表
└── DIRECTORY_TREE.txt # 详细目录树
```

详细目录结构请查看 `DIRECTORY_TREE.txt` 文件。

For detailed directory structure, see `DIRECTORY_TREE.txt`.

## 🛠️ 安装与配置 | Installation & Configuration

### 1. 环境要求 | Requirements

- Python 3.8+
- NCBI SRA Toolkit (用于数据下载)

### 2. 安装 Python 依赖 | Install Python Dependencies

```bash
# 克隆或下载项目
cd GeneExpressionOmnibus_E-utilities

# 安装 Python 包
pip install -r requirements.txt
```

### 3. 安装 SRA Toolkit (可选，仅用于下载数据) | Install SRA Toolkit (Optional)

SRA Toolkit 用于下载原始 FASTQ 数据。

SRA Toolkit is required for downloading raw FASTQ data.

**MacOS:**
```bash
brew install sratoolkit
```

**Linux:**
```bash
# 下载预编译二进制文件
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
# 添加到 PATH
export PATH=$PATH:$PWD/sratoolkit.3.0.0-ubuntu64/bin
```

**详细安装指南:** https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

### 4. 配置 | Configuration

**重要：** 编辑配置文件，设置您的邮箱地址：

**IMPORTANT:** Edit the configuration file and set your email address:

```python
# config/config.py 或 config/config_v2.py
ENTREZ_EMAIL = "your.email@example.com"  # 替换为您的邮箱
```

NCBI 要求提供邮箱地址以追踪 API 使用情况。

NCBI requires an email address for API usage tracking.

**可选：** 如果您有 NCBI API Key，可以提高请求速率：

**Optional:** If you have an NCBI API Key, you can increase the request rate:

```python
ENTREZ_API_KEY = "your_api_key_here"
```

获取 API Key: https://www.ncbi.nlm.nih.gov/account/settings/

**注意：** 项目包含多个配置文件，可根据研究需求选择或创建新的配置：
- `config/config.py` - V1版本配置（肺部转移瘤示例）
- `config/config_v2.py` - V2改进版配置（肺部转移瘤示例）
- `config/config_mouse_bcell_taok.py` - 小鼠B细胞发育示例配置

## 🚀 使用方法 | Usage

### 运行挖掘流水线 | Run Mining Pipeline

**示例1：肺部转移瘤数据挖掘（V1版本）**
```bash
python scripts/geo_lung_metastasis_miner.py
```

**示例2：肺部转移瘤数据挖掘（V2版本，推荐）**
```bash
python scripts/run_v2_mining.py
```

**示例3：小鼠骨髓B细胞发育数据挖掘**
```bash
python scripts/run_mouse_bcell_taok_mining.py
```

**自定义挖掘任务：**

1. 创建新的配置文件（参考 `config/config_mouse_bcell_taok.py`）
2. 创建新的挖掘脚本（参考 `scripts/geo_mouse_bcell_taok_miner.py`）
3. 运行自定义脚本

流水线将执行以下操作：

The pipeline will:

1. **搜索 GEO 数据库** - 使用预定义的查询策略搜索相关数据集
2. **解析元数据** - 下载并深度分析每个样本的元数据
3. **应用过滤规则** - 使用启发式规则识别符合条件的样本
4. **关联 SRA** - 获取原始数据的下载链接
5. **生成结果文件**：
   - `results/` 目录下的详细结果表格（CSV格式）
   - `data/` 目录下的SRR编号列表
   - `logs/` 目录下的执行日志

### 输出文件说明 | Output Files

所有输出文件位于 `results/` 和 `data/` 目录中。

All output files are located in `results/` and `data/` directories.

结果表格通常包含以下列：

Result tables typically include the following columns:

| 列名 | 说明 |
|------|------|
| GSE | GEO 数据集编号 |
| GSM | GEO 样本编号 |
| SRX | SRA 实验编号 |
| Library_Strategy | 测序技术 (RNA-Seq, ATAC-seq 等) |
| Title | 样本标题 |
| Characteristics | 样本特征描述 |
| Filter_Reason | 过滤器判定理由 |
| Source_Name | 样本来源名称 |
| SRR_List | SRA 运行编号列表 |
| SRR_Count | SRR 数量 |
| Confidence | 置信度评分（V2版本） |

### 下载原始数据 | Download Raw Data

⚠️ **重要提示：** 在下载数据之前，请务必手动审核 CSV 结果文件，确认样本符合您的研究需求。

⚠️ **IMPORTANT:** Before downloading data, manually review the CSV results file to verify samples meet your research criteria.

使用提供的下载脚本：

Use the provided download script:

```bash
# 使脚本可执行
chmod +x scripts/download_sra_data.sh

# 运行下载脚本
./scripts/download_sra_data.sh
```

或手动使用 SRA Toolkit：

Or manually use SRA Toolkit:

```bash
# 下载单个 SRR
prefetch SRR12345678
fasterq-dump SRR12345678

# 批量下载
cat data/SRR_accession_list.txt | xargs -n 1 prefetch
cat data/SRR_accession_list.txt | xargs -n 1 fasterq-dump
```

## 🔧 自定义配置 | Customization

### 创建新的挖掘任务

1. **创建配置文件** - 在 `config/` 目录下创建新的配置文件
2. **定义搜索策略** - 设置搜索关键词和查询逻辑
3. **定义过滤规则** - 设置样本过滤条件和置信度评分
4. **创建挖掘脚本** - 参考现有脚本创建新的挖掘器

### 修改搜索策略 | Modify Search Strategy

编辑配置文件中的搜索参数：

Edit search parameters in configuration files:

```python
# 技术术语
TECH_TERMS = '("scRNA-seq" OR "single cell RNA-seq" OR ...)'

# 生物学术语
BIOLOGY_TERMS = '("your_keywords" AND ...)'

# 基础过滤
BASE_FILTERS = '"Organism"[Organism] AND ...'
```

### 调整过滤规则 | Adjust Filtering Rules

编辑挖掘脚本中的过滤方法：

Edit filtering methods in mining scripts:

```python
def enhanced_filter(self, gsm_metadata, gse_id):
    # 在此处自定义过滤逻辑
    # Customize filtering logic here
    ...
```

## 📊 过滤逻辑说明 | Filtering Logic

流水线使用可配置的启发式规则来识别符合条件的样本。过滤规则可根据研究需求进行定制。

The pipeline uses configurable heuristic rules to identify eligible samples. Filtering rules can be customized according to research needs.

### 通用过滤原则 | General Filtering Principles

1. **物种筛选** - 根据研究目标筛选特定物种
2. **组织类型** - 根据研究目标筛选特定组织
3. **技术平台** - 筛选特定的测序技术
4. **样本质量** - 排除低质量或不符合条件的样本

### 置信度评分系统 | Confidence Scoring System

V2版本及后续版本支持置信度评分，帮助评估样本的相关性：

V2 and later versions support confidence scoring to help assess sample relevance:

- **高置信度 (≥0.8)** - 样本高度符合研究目标，可直接使用
- **中等置信度 (0.5-0.8)** - 样本可能相关，需要人工复核
- **低置信度 (<0.5)** - 样本相关性较低，建议排除

## ⚠️ 重要注意事项 | Important Notes

### 1. 手动复核的必要性

**生物医学元数据存在固有的复杂性和不一致性。** 自动化脚本的结果必须经过人工审核：

- 阅读相关研究的论文和描述
- 验证样本确实符合研究目标
- 确认样本特征符合研究需求
- 检查是否有排除标准

### 2. API 使用限制

- NCBI 限制每秒 3 次请求（无 API Key）
- 有 API Key 可提升到每秒 10 次
- 脚本已内置延迟机制，请勿修改

### 3. 数据量考虑

- 单个 FASTQ 文件可能有几 GB 到几十 GB
- 下载前检查可用磁盘空间
- 考虑使用云存储或高性能计算集群

### 4. 缓存机制

- GEO 元数据会缓存在 `GEO_Cache/` 目录
- 重新运行时会使用缓存，加快速度
- 如需强制重新下载，删除缓存目录

### 5. 项目结构说明

- **scripts/**: 所有可执行脚本
- **config/**: 配置文件（修改配置请编辑此目录下的文件）
- **results/**: 挖掘结果CSV文件
- **docs/**: 项目文档、分析报告和使用指南
- **logs/**: 运行日志文件
- **data/**: 数据文件（SRR列表等）
- **GEO_Cache/**: GEO数据缓存（自动生成，可忽略）

## 🐛 故障排除 | Troubleshooting

### 问题 1: "Please configure your email in config.py"

**解决方案：** 编辑配置文件，设置 `ENTREZ_EMAIL` 为您的邮箱。

### 问题 2: "Error during GEO search"

**可能原因：**
- 网络连接问题
- NCBI 服务器暂时不可用

**解决方案：** 等待几分钟后重试，或检查网络连接。

### 问题 3: "No relevant samples found"

**可能原因：**
- 搜索条件过于严格
- 过滤规则过于保守

**解决方案：** 调整配置文件中的搜索参数和过滤规则。

### 问题 4: SRA Toolkit 下载失败

**解决方案：**
```bash
# 配置 SRA Toolkit
vdb-config --interactive

# 测试连接
prefetch --help
```

## 📚 参考资源 | References

- **NCBI GEO:** https://www.ncbi.nlm.nih.gov/geo/
- **NCBI E-utilities:** https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **SRA Toolkit:** https://github.com/ncbi/sra-tools
- **BioPython Documentation:** https://biopython.org/wiki/Documentation
- **GEOparse Documentation:** https://geoparse.readthedocs.io/

## 📄 许可证 | License

本项目采用 MIT 许可证。

This project is licensed under the MIT License.

```
MIT License

Copyright (c) 2025 GeneExpressionOmnibus_E-utilities Contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

**使用 NCBI 数据请遵守其使用条款。**

**Please comply with NCBI terms of use when using their data.**

## 🤝 贡献 | Contributing

欢迎提交问题报告和改进建议。如果您有更好的过滤规则、新的应用场景或发现了 bug，请提交 Issue 或 Pull Request。

Issues and improvement suggestions are welcome. If you have better filtering rules, new application scenarios, or find bugs, please submit an Issue or Pull Request.

## 📧 联系方式 | Contact

如有问题或建议，请通过 GitHub Issues 联系。

For questions or suggestions, please contact via GitHub Issues.

---

## 📖 相关文档 | Related Documentation

项目包含详细的文档和指南，位于 `docs/` 目录：

The project includes detailed documentation and guides in the `docs/` directory:

- `QUICKSTART.md` - 快速开始指南
- `HOW_TO_REVIEW_V2.md` - V2版本结果审核指南
- `REVIEW_QUICKSTART.md` - 结果审核快速指南
- `IMPLEMENTATION_SUMMARY.md` - 实现总结
- `V1_VS_V2_COMPARISON.md` - V1与V2版本对比
- `MOUSE_BCELL_TAOK_MINING_RESULTS.md` - 小鼠B细胞发育数据挖掘结果说明
- `RECOMMENDED_GEO_DATASETS.md` - 推荐GEO数据集
- 其他分析报告和文档

## 🧪 测试脚本 | Test Scripts

项目包含多个测试和辅助脚本，位于 `scripts/` 目录：

The project includes several test and utility scripts in the `scripts/` directory:

- `test_installation.py` - 测试安装和配置
- `smoke_test.py` - 冒烟测试
- `check_gsm.py` - 检查GSM样本信息
- `check_missing_dataset.py` - 检查缺失的数据集
- `review_v2_results.py` - 审核V2结果
- `example_usage.py` - 使用示例

**运行测试脚本：**

```bash
# 测试安装和配置
python scripts/test_installation.py

# 运行冒烟测试
python scripts/smoke_test.py

# 检查特定GSM样本
python scripts/check_gsm.py

# 审核V2结果
python scripts/review_v2_results.py
```

---

**最后更新 | Last Updated:** 2025-12-15
