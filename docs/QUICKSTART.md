# 快速开始指南 | Quick Start Guide

这是一个5分钟快速上手指南，帮助您立即开始使用 GEO 肺部转移瘤数据挖掘流水线。

This is a 5-minute quick start guide to help you immediately start using the GEO Lung Metastasis Mining Pipeline.

## 📦 快速安装 | Quick Installation

```bash
# 1. 安装 Python 依赖
pip install -r requirements.txt

# 2. 配置邮箱地址（必需！）
# 编辑 config.py，将 ENTREZ_EMAIL 改为您的邮箱
nano config.py
# 或者使用您喜欢的编辑器
```

在 `config.py` 中修改：
```python
ENTREZ_EMAIL = "your.email@example.com"  # 改为您的真实邮箱
```

## 🚀 立即运行 | Run Immediately

```bash
# 运行数据挖掘流水线
python geo_lung_metastasis_miner.py
```

就这么简单！流水线会自动：
1. 搜索 GEO 数据库
2. 解析样本元数据
3. 应用过滤规则识别肺部转移瘤
4. 关联 SRA 原始数据
5. 生成结果文件

## 📊 查看结果 | View Results

运行完成后，您会得到：

1. **GEO_Lung_Metastasis_Mining_Results.csv** - 详细结果表格
   ```bash
   # 用 Excel 或其他工具打开
   open GEO_Lung_Metastasis_Mining_Results.csv
   ```

2. **SRR_accession_list.txt** - 可下载的 SRR 列表
   ```bash
   # 查看 SRR 数量
   wc -l SRR_accession_list.txt
   ```

3. **geo_mining_*.log** - 执行日志
   ```bash
   # 查看日志
   tail -f geo_mining_*.log
   ```

## ⚠️ 重要提示 | Important

### 在下载数据之前

**务必手动审核 CSV 结果文件！** 自动过滤可能不完美，请确认样本确实符合您的研究需求。

**Manually review the CSV results file!** Automatic filtering may not be perfect, please verify samples truly meet your research criteria.

打开 CSV 文件，检查以下列：
- `Title` - 样本标题
- `Characteristics` - 样本特征
- `Filter_Reason` - 过滤判定理由

## 💾 下载原始数据 | Download Raw Data

### 前置条件：安装 SRA Toolkit

**MacOS:**
```bash
brew install sratoolkit
```

**Linux:**
```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzf sratoolkit.current-ubuntu64.tar.gz
export PATH=$PATH:$PWD/sratoolkit.*/bin
```

### 运行下载脚本

```bash
# 确认样本无误后，运行下载脚本
./download_sra_data.sh
```

下载的数据会保存在 `./SRA_Data/` 目录下。

## 🎯 典型工作流程 | Typical Workflow

```bash
# 1. 配置邮箱（一次性）
nano config.py

# 2. 运行挖掘流水线
python geo_lung_metastasis_miner.py

# 3. 审核结果
open GEO_Lung_Metastasis_Mining_Results.csv

# 4. 如果结果满意，下载数据
./download_sra_data.sh
```

## 🔧 自定义搜索 | Customize Search

如果您想修改搜索条件，编辑 `config.py`：

```python
# 例如：只搜索单细胞数据
TECH_TERMS = '("scRNA-seq" OR "single cell RNA-seq")'

# 或者：添加特定的癌症类型
BIOLOGY_TERMS = '(("lung") AND ("metastasis") AND ("breast cancer"))'
```

然后重新运行：
```bash
python geo_lung_metastasis_miner.py
```

## 📝 查看高级用法 | Advanced Usage

查看示例脚本了解更多功能：
```bash
python example_usage.py
```

或阅读完整文档：
```bash
cat README.md
```

## 🆘 常见问题 | Common Issues

### 问题 1: "Please configure your email"
**解决：** 编辑 `config.py`，设置 `ENTREZ_EMAIL`

### 问题 2: "No module named 'Bio'"
**解决：** 运行 `pip install -r requirements.txt`

### 问题 3: "prefetch: command not found"
**解决：** 安装 SRA Toolkit（见上文）

### 问题 4: 下载速度很慢
**解决：** 
- 获取 NCBI API Key 并配置在 `config.py`
- 调整 `download_sra_data.sh` 中的 `PARALLEL_JOBS`

### 问题 5: 磁盘空间不足
**解决：** 
- 检查可用空间：`df -h`
- 考虑只下载部分样本
- 启用 FASTQ 压缩（默认已启用）

## 📈 预期运行时间 | Expected Runtime

- **搜索阶段：** 1-5 分钟
- **解析阶段：** 取决于数据集数量（每个 GSE 约 10-30 秒）
- **SRA 关联：** 1-3 分钟
- **数据下载：** 每个样本 10-60 分钟（取决于数据大小和网络速度）

## 🎓 学习资源 | Learning Resources

- [NCBI GEO 数据库](https://www.ncbi.nlm.nih.gov/geo/)
- [E-utilities API 文档](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [SRA Toolkit 使用指南](https://github.com/ncbi/sra-tools/wiki)

## 📧 获取帮助 | Get Help

遇到问题？查看详细文档：
```bash
cat README.md
```

或提交 Issue 到项目仓库。

---

**祝您数据挖掘顺利！| Happy Data Mining!** 🎉

