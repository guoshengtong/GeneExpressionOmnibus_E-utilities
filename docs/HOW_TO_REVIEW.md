# 如何进行手工复核 | How to Perform Manual Review

**快速指南 - 5分钟开始复核**

---

## 🎯 复核目标

确认挖掘到的39个样本中，哪些**真正是其他部位癌症转移到肺部**的样本。

---

## 📋 快速复核流程（推荐）

### 第1步：使用复核工具快速分析（2分钟）

```bash
# 运行交互式复核工具
python3 review_helper.py
```

在菜单中选择：
1. 先选择 `1` - 查看结果摘要
2. 再选择 `4` - 分析需要注意的样本
3. 然后选择 `6` - 导出筛选后的结果

这将自动：
- ✅ 排除 GSE266330（7个样本 - 肺癌骨转移）
- ⚠️ 询问是否排除 GSE179373（6个样本 - 脑转移）
- ✅ 生成新的筛选后 CSV 和 SRR 列表

### 第2步：打开结果文件检查（5分钟）

```bash
# 打开筛选后的结果
open GEO_Results_Filtered.csv
```

重点检查以下样本：

#### ✅ 可能符合要求（建议保留）

**GSE270231** - 7个样本
- 描述：骨肉瘤肺转移
- 原发部位：Osteosarcoma（骨肉瘤）
- 转移位置：Lung
- ✅ **建议：保留**

**GSE193594/GSE193534** - 4个样本
- 描述：患者肺转移样本
- 标题明确："Lung metastasis"
- ✅ **建议：保留**

#### ⚠️ 需要查看GEO页面（10分钟）

**GSE109281** - 14个样本
- 描述：标注为"metastasis"的肺部样本
- 问题：缺少明确的原发部位
- ⚠️ **需要查看GEO页面确认**

查看方法：
```bash
# 在浏览器中打开
open "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109281"
```

**GSE156405** - 1个样本
- 描述：简单标注"Lung metastasis"
- ⚠️ **需要查看原始信息**

### 第3步：批量查看GEO页面（可选，10分钟）

使用复核工具自动打开浏览器：

```bash
python3 review_helper.py
# 选择 5 - 在浏览器中打开 GEO 页面
# 选择 2 - 需要复核的 GSE (前5个)
```

这将自动打开需要复核的数据集页面。

---

## 🔍 手动复核方法（详细版）

### 方法1：查看CSV文件

打开 `GEO_Lung_Metastasis_Mining_Results.csv`，对每一行检查：

#### 判断标准

| 列名 | 要找什么 | 好的例子 | 坏的例子 |
|------|---------|---------|---------|
| Title | 明确的转移描述 | "Lung metastasis" | "Lung tumor"（不明确） |
| Characteristics | 原发部位信息 | "primary: breast" | 只有"tissue: lung" |
| Source_Name | 组织来源 | "Lung metastasis" | "Bone"（不是肺） |
| Filter_Reason | 高置信度证据 | "Detected X cancer with metastasis" | "Possible..." |

#### 排除标准

❌ **必须排除：**
- Source_Name 是 "Bone", "Brain" 等（不是肺组织）
- 描述为 "lung cancer metastasis to [other organ]"（肺癌向外转移）
- 没有任何转移相关的描述

⚠️ **需要进一步确认：**
- 只有 "lung tumor" 而没有其他信息
- 没有明确的原发部位
- 描述模糊或矛盾

### 方法2：查看GEO数据库页面

对于每个需要确认的GSE，访问：
```
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=[GSE编号]
```

在页面上查找：
1. **Series Summary** - 研究的总体描述
2. **Overall Design** - 实验设计说明
3. **Citation** - 如果有论文，可以看标题和摘要
4. **Sample Characteristics** - 每个样本的详细信息

---

## 📊 当前结果分析

根据初步分析：

### 统计概览
- 总样本数：39个
- 明确应该排除：7个（GSE266330）
- 可能不符合（脑转移）：6个（GSE179373）
- 高置信度样本：~10-15个
- 需要确认的样本：~15-20个

### 建议的复核优先级

#### 🚀 第一优先级（立即处理）
1. ✅ 排除 GSE266330（肺癌骨转移）
2. ⚠️ 决定是否保留 GSE179373（脑转移）

#### 📝 第二优先级（10分钟）
3. 查看 GSE270231 的 GEO 页面（骨肉瘤肺转移）
4. 查看 GSE109281 的 GEO 页面（确认原发部位）

#### 🔬 第三优先级（如果需要）
5. 查看 GSE156405 和其他小数据集

---

## 🛠️ 实用工具使用

### 工具1：复核辅助脚本

```bash
python3 review_helper.py
```

功能：
- ✅ 显示结果摘要
- ✅ 分析问题样本
- ✅ 自动打开GEO页面
- ✅ 导出筛选后结果

### 工具2：命令行快速查看

```bash
# 查看所有GSE列表
cut -d',' -f1 GEO_Lung_Metastasis_Mining_Results.csv | sort | uniq -c

# 查看特定GSE的样本
grep "GSE270231" GEO_Lung_Metastasis_Mining_Results.csv

# 统计每个GSE的样本数
cut -d',' -f1 GEO_Lung_Metastasis_Mining_Results.csv | tail -n +2 | sort | uniq -c
```

---

## ✅ 复核检查清单

使用这个清单跟踪进度：

```
□ 第1步：运行复核工具，导出筛选后结果
  □ 排除 GSE266330
  □ 决定 GSE179373 的处理

□ 第2步：查看重点数据集（优先级排序）
  □ GSE270231 - 骨肉瘤肺转移
  □ GSE109281 - 确认原发部位
  □ GSE193594/GSE193534 - 患者肺转移
  □ GSE156405 - 确认详情

□ 第3步：创建最终确认的样本列表
  □ 记录保留的样本
  □ 记录排除的样本及原因

□ 第4步：准备下载
  □ 确认 SRR 列表
  □ 评估存储空间需求
  □ 安装 SRA Toolkit（如需要）
```

---

## 📖 示例：如何判断一个样本

### 示例1：明确符合的样本

```
GSE270231, GSM8337654
Title: "NCHS-009 Single-cell RNAseq of osteosarcoma patient metastasis. Hilar lobe"
Characteristics: "tissue: lung tumor"
Source_Name: "lung tumor"
```

**判断：** ✅ **保留**
- 明确提到 "osteosarcoma"（骨肉瘤）
- 明确提到 "metastasis"
- 组织是肺部肿瘤
- 符合"其他部位癌症转移到肺"的要求

### 示例2：应该排除的样本

```
GSE266330, GSM8245049
Title: "LC_01_1"
Characteristics: "tissue: Bone; group: lung cancer bone metastasis"
Source_Name: "Bone"
```

**判断：** ❌ **排除**
- "lung cancer bone metastasis" = 肺癌转移到骨
- 组织来源是骨头，不是肺
- 不符合要求

### 示例3：需要进一步确认的样本

```
GSE109281, GSM2936442
Title: "Tumor_Lungs_Mixed_Mixed_1-1"
Characteristics: "clonality: polyclonal; site: metastasis; cell type: mixed"
Source_Name: "tumor"
```

**判断：** ⚠️ **需要查看GEO页面**
- 标注为 "metastasis"，位置在肺
- 但是缺少原发部位信息
- 需要查看研究描述确认

---

## 💡 快速决策指南

### 5分钟快速复核（最小化方案）

如果时间有限，至少做这些：

```bash
# 1. 排除明确不符合的样本
python3 review_helper.py
# 选择 6 - 导出筛选后结果
# 这会自动排除 GSE266330

# 2. 使用筛选后的结果
# 文件：GEO_Results_Filtered.csv
# 文件：SRR_accession_list_filtered.txt
```

### 30分钟标准复核（推荐）

```bash
# 1. 运行工具分析（5分钟）
python3 review_helper.py
# 选择 1, 4, 6

# 2. 查看重点数据集（15分钟）
# 访问 GSE270231, GSE109281, GSE193594 的 GEO 页面

# 3. 做出决策（10分钟）
# 记录保留和排除的样本
```

### 1小时深度复核（完整版）

- 查看所有数据集的GEO页面
- 阅读相关论文摘要
- 详细记录每个样本的判断依据

---

## 🎬 现在开始！

**最简单的开始方式：**

```bash
# 1. 运行这个命令
python3 review_helper.py

# 2. 在菜单中按顺序选择：
#    1 -> 查看摘要
#    4 -> 分析问题样本
#    6 -> 导出筛选后结果

# 3. 完成！你现在有了初步筛选后的结果
```

---

## 📞 需要帮助？

如果您在复核过程中遇到问题，或者需要：
- 解释某个样本的元数据
- 帮助查找原始论文
- 其他技术支持

随时告诉我！

---

**祝复核顺利！** 🎉


