# 小鼠骨髓单细胞B细胞发育和TAOK基因数据挖掘指南
# Mouse Bone Marrow Single-Cell B Cell Development and TAOK Gene Mining Guide

## 📋 概述 | Overview

本指南介绍如何使用数据挖掘脚本在GEO数据库中搜索小鼠骨髓单细胞RNA测序数据，重点关注B细胞发育和TAOK基因的影响。

This guide explains how to use the data mining script to search the GEO database for mouse bone marrow single-cell RNA-seq data, focusing on B cell development and TAOK gene effects.

## 🎯 搜索目标 | Search Objectives

1. **小鼠骨髓单细胞RNA测序数据** | Mouse bone marrow single-cell RNA-seq data
2. **B细胞发育相关** | B cell development related
3. **TAOK基因相关** | TAOK gene related (optional, may need expression analysis)

## 🔍 搜索关键词策略 | Search Keyword Strategy

脚本使用多阶段搜索策略，关键词遵循"标准+英文+数据类型导向"的原则：

The script uses a multi-stage search strategy with keywords following the "standard + English + data type oriented" principle:

### 阶段1：基础搜索 | Stage 1: Base Search
- `mouse bone marrow single-cell RNA-seq`
- `mouse bone marrow scRNA-seq`

### 阶段2：B细胞相关搜索 | Stage 2: B Cell Related Search
- `mouse bone marrow scRNA-seq B cell`
- `mouse bone marrow scRNA-seq B-cell`

### 阶段3：B细胞发育相关搜索 | Stage 3: B Cell Development Search
- `mouse bone marrow scRNA-seq B cell development`
- `mouse bone marrow scRNA-seq B-cell development`

### 阶段4：TAOK基因相关搜索（可选）| Stage 4: TAOK Gene Search (Optional)
- `Taok mouse B cell scRNA-seq`
- `TAOK kinase mouse immune single-cell`

**注意**：GEO中很多数据集不会在标题中直接提到具体基因名，建议先找数据，再在数据中分析TAOK表达。

**Note**: Many GEO datasets don't mention specific gene names in titles. It's recommended to find data first, then analyze TAOK expression in the data.

## 🚀 使用方法 | Usage

### 1. 运行挖掘脚本 | Run Mining Script

```bash
# 方法1：直接运行挖掘脚本
python scripts/geo_mouse_bcell_taok_miner.py

# 方法2：使用运行脚本（推荐）
python scripts/run_mouse_bcell_taok_mining.py
```

### 2. 脚本执行流程 | Script Execution Flow

1. **多阶段搜索** - 执行4个搜索阶段，逐步缩小范围
2. **数据解析** - 下载并解析每个GSE数据集的元数据
3. **智能过滤** - 使用置信度评分系统过滤样本
4. **结果保存** - 按置信度分组保存结果

### 3. 输出文件 | Output Files

脚本会在以下位置生成结果文件：

The script will generate result files in the following locations:

#### 主要结果文件 | Main Result Files

- `results/GEO_Mouse_Bone_Marrow_Bcell_TAOK_Mining_Results.csv` - 完整结果表格
- `results/Results_Mouse_Bcell_High_Confidence.csv` - 高置信度样本（≥0.8）
- `results/Results_Mouse_Bcell_Needs_Review.csv` - 需要复核的样本（0.5-0.8）
- `results/Results_Mouse_Bcell_Low_Confidence.csv` - 低置信度样本（<0.5）

#### 数据文件 | Data Files

- `data/SRR_accession_list_mouse_bcell.txt` - SRR编号列表（用于下载）

#### 日志文件 | Log Files

- `logs/mouse_bcell_taok_mining_YYYYMMDD_HHMMSS.log` - 详细执行日志

### 4. 结果表格字段说明 | Result Table Columns

| 字段名 | 说明 |
|--------|------|
| GSE | GEO数据集编号 |
| GSM | GEO样本编号 |
| SRX | SRA实验编号 |
| SRR_List | SRA运行编号列表 |
| SRR_Count | SRR数量 |
| Library_Strategy | 测序技术 |
| Title | 样本标题 |
| Characteristics | 样本特征描述 |
| Source_Name | 样本来源名称 |
| Organism | 物种信息 |
| B_Cell_Stages | 检测到的B细胞发育阶段 |
| Has_TAOK | 元数据中是否提到TAOK |
| Filter_Reason | 过滤判定理由 |
| Confidence | 置信度评分（0-1） |

## 🔧 配置说明 | Configuration

配置文件位于：`config/config_mouse_bcell_taok.py`

### 主要配置项 | Main Configuration Items

1. **搜索查询** - 可以修改搜索关键词
2. **过滤规则** - 可以调整关键词列表和置信度阈值
3. **输出设置** - 可以修改输出文件名和路径

### 修改搜索策略 | Modify Search Strategy

编辑 `config/config_mouse_bcell_taok.py`：

```python
# 修改搜索查询
SEARCH_QUERY_BASE = '''
    ({TECH_TERMS}) AND
    ("mouse" OR "Mus musculus") AND
    ("bone marrow" OR "bone-marrow" OR "BM") AND
    ...
'''
```

### 调整置信度阈值 | Adjust Confidence Thresholds

```python
CONFIDENCE_THRESHOLDS = {
    'high': 0.8,        # 高置信度
    'medium': 0.5,      # 中等置信度
    'low': 0.3          # 低置信度
}
```

## 📊 过滤逻辑说明 | Filtering Logic

### ✅ 必须满足的条件 | Required Conditions

1. **小鼠样本** - 必须包含 "mouse" 或 "Mus musculus"
2. **骨髓组织** - 必须包含 "bone marrow" 相关关键词
3. **单细胞技术** - 必须包含单细胞RNA测序相关关键词

### ❌ 排除条件 | Exclusion Criteria

- 人类样本
- 细胞系（除非是原代细胞）
- 肿瘤样本（除非明确是B细胞发育研究）

### 🎯 加分项 | Bonus Points

- **B细胞关键词** - 提高置信度 +0.2
- **B细胞发育阶段** - 提高置信度 +0.2
- **TAOK基因关键词** - 提高置信度 +0.3

## ⚠️ 重要提示 | Important Notes

### 1. TAOK基因分析

**重要**：GEO数据集的元数据中通常不会直接提到具体基因名。脚本会标记元数据中提到TAOK的样本，但大多数情况下需要：

1. 先找到符合条件的数据集
2. 下载原始数据或表达矩阵
3. 在表达数据中分析TAOK基因的表达

### 2. 手动复核

自动化脚本的结果必须经过人工审核：

- 阅读相关研究的论文和描述
- 验证样本确实来自小鼠骨髓
- 确认B细胞发育相关信息
- 检查是否有排除标准

### 3. API使用限制

- NCBI限制每秒3次请求（无API Key）
- 有API Key可提升到每秒10次
- 脚本已内置延迟机制

### 4. 数据下载

找到相关数据后，可以使用SRR列表下载原始数据：

```bash
# 使用提供的下载脚本
chmod +x scripts/download_sra_data.sh
./scripts/download_sra_data.sh

# 或手动下载
cat data/SRR_accession_list_mouse_bcell.txt | xargs -n 1 prefetch
cat data/SRR_accession_list_mouse_bcell.txt | xargs -n 1 fasterq-dump
```

## 🐛 故障排除 | Troubleshooting

### 问题1：导入错误

**解决方案**：确保在项目根目录运行脚本，或使用绝对路径。

### 问题2：搜索结果为空

**可能原因**：
- 搜索条件过于严格
- 网络连接问题
- NCBI服务器暂时不可用

**解决方案**：
- 调整配置文件中的搜索关键词
- 检查网络连接
- 等待几分钟后重试

### 问题3：过滤结果太少

**解决方案**：
- 降低置信度阈值
- 检查排除关键词是否过于严格
- 调整过滤规则

## 📚 相关资源 | References

- **NCBI GEO**: https://www.ncbi.nlm.nih.gov/geo/
- **NCBI E-utilities**: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **SRA Toolkit**: https://github.com/ncbi/sra-tools

## 📝 更新日志 | Changelog

- **2025-12-10**: 初始版本创建

---

**最后更新 | Last Updated**: 2025-12-10

