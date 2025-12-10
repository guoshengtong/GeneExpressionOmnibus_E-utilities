# 如何手工复核 V2 结果 - 详细指南

## 📋 复核概览

**V2结果统计：**
- 总样本数：182个
- V1已有：39个
- V2新增：143个
  - 高置信度（≥0.8）：21个
  - 中等置信度（0.6-0.8）：122个

**复核目标：**
1. 确认高置信度样本的准确性
2. 筛选出真正符合要求的样本
3. 排除误报（如：肺癌的转移，而非转移到肺）

---

## 🚀 快速开始：三步复核法

### 第1步：使用交互式工具快速浏览（5分钟）

```bash
cd /Users/khugjil-devstation/Projects/GeneExpressionOmnibus_E-utilities
python3 review_v2_results.py
```

进入交互界面后，执行：
```
复核> quick_review
```

这会显示所有高置信度新增样本（21个），按数据集分组。

### 第2步：查看高置信度样本CSV（10分钟）

```bash
open V2_New_Samples_High_Confidence.csv
```

在Excel/Numbers中打开，重点查看这些列：
- `GSE`, `GSM` - 数据集和样本编号
- `Title` - 样本标题（最重要！）
- `Source_Name` - 来源名称
- `Confidence` - 置信度
- `Filter_Reason` - 判定理由

### 第3步：访问GEO页面确认（20-30分钟）

对于高置信度样本，逐个访问GEO页面确认。

---

## 🔍 详细复核流程

### 方法1：使用交互式工具（推荐）

#### 启动工具
```bash
python3 review_v2_results.py
```

#### 常用命令

**1. 查看摘要**
```
复核> summary
```
显示整体统计、置信度分布、数据集排名等。

**2. 快速复核高置信度新增样本**
```
复核> quick_review
```
专门显示V2新增的高置信度样本，按数据集分组。

**3. 过滤并查看新增样本**
```
复核> filter_new
```
只显示V2新增的143个样本。

**4. 过滤高置信度样本**
```
复核> filter_conf 0.8
```
只显示置信度≥0.8的样本。

**5. 查看样本列表**
```
复核> list 20
```
显示前20个样本的概要信息。

**6. 查看样本详情**
```
复核> detail 1
```
显示第1个样本的完整信息，包括GEO/SRA链接。

**7. 在浏览器中打开**
```
复核> open 1 geo     # 打开GSE页面
复核> open 1 gsm     # 打开GSM样本页面
复核> open 1 sra     # 打开SRA页面
```

**8. 按数据集过滤**
```
复核> filter_gse GSE234187
```
查看特定数据集的所有样本。

**9. 导出过滤结果**
```
复核> export my_filtered_results.csv
```
将当前过滤的样本导出到CSV文件。

**10. 重置过滤**
```
复核> reset
```
清除所有过滤条件。

#### 示例工作流程

```bash
# 1. 启动工具
python3 review_v2_results.py

# 2. 查看摘要
复核> summary

# 3. 快速浏览高置信度新增样本
复核> quick_review

# 4. 过滤出新增的高置信度样本
复核> filter_new
复核> filter_conf 0.8

# 5. 查看前20个
复核> list 20

# 6. 查看第1个样本详情
复核> detail 1

# 7. 在浏览器中打开确认
复核> open 1 geo

# 8. 导出这些高置信度样本
复核> export high_conf_new_samples.csv

# 9. 退出
复核> quit
```

---

### 方法2：使用CSV文件手工复核

#### 第1步：打开高置信度新增样本

```bash
open V2_New_Samples_High_Confidence.csv
```

#### 第2步：创建复核表格

在Excel中添加新列：
- `Manual_Review` - 手工复核结果（Accept/Reject/Uncertain）
- `Review_Notes` - 复核备注
- `Cancer_Origin` - 原发癌症类型
- `Is_Lung_Met` - 是否真的是肺转移（Yes/No）

#### 第3步：逐个复核

对于每个样本：

1. **查看标题和来源**
   - `Title`列：看是否明确提到"lung metastasis"
   - `Source_Name`列：看来源组织/细胞系

2. **访问GEO页面确认**
   ```
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE编号
   ```
   
   重点查看：
   - Overall design（研究设计）
   - Sample characteristics（样本特征）
   - 是否提到原发肿瘤类型
   - 是否明确是肺转移

3. **标记复核结果**
   - ✅ **Accept** - 明确是其他癌症的肺转移
   - ❌ **Reject** - 是肺癌、细胞系、或方向错误
   - ⚠️ **Uncertain** - 需要进一步查论文

---

## ⚠️ 复核中的常见陷阱

### 陷阱1：方向错误

❌ **错误示例：**
- "Lung cancer bone metastasis" - 这是**肺癌转移到骨**，不是转移到肺
- "Lung cancer lymph node metastasis" - 肺癌转移到淋巴结

✅ **正确示例：**
- "Osteosarcoma lung metastasis" - 骨肉瘤转移到肺 ✓
- "Breast cancer lung metastasis" - 乳腺癌转移到肺 ✓

**判断方法：**
- 关键词在前 = 原发肿瘤
- 关键词在后 = 转移部位
- "Cancer A **to/in** B" → A是原发，B是转移部位

### 陷阱2：细胞系 vs 患者组织

很多高置信度样本是**细胞系**而非患者组织：

❌ **细胞系示例：**
- MDA-MB-231 LM2 - 乳腺癌肺转移细胞系
- 786-O-LM - 肾癌肺转移细胞系
- PDX (patient-derived xenograft) - 患者来源异种移植

**是否接受取决于研究需求：**
- 如果研究机制 → 可能接受细胞系
- 如果研究临床样本 → 只接受患者组织

**识别方法：**
- 看`Source_Name`列
- 包含"cell line", "culture", "PDX", "xenograft" → 细胞系/模型
- 包含"tissue", "patient", "tumor" → 患者样本

### 陷阱3：原发肺癌的其他部位转移

❌ **不符合的示例：**
- "Brain metastasis from lung adenocarcinoma" - 肺癌的脑转移
- "Lung cancer liver metastasis" - 肺癌肝转移
- "NSCLC bone metastasis" - 非小细胞肺癌骨转移

即使包含"lung"和"metastasis"，但方向错误。

### 陷阱4：胸腔积液样本

⚠️ **需要仔细确认：**
- "lung; metastatic site to pleural effusion"
- "pleural effusion from lung adenocarcinoma"

**问题：** 胸腔积液样本可能是：
1. 原发肺癌的积液
2. 其他癌症肺转移的积液 ✓

需要查看完整的GEO页面或论文确认。

---

## 📊 数据集特别注意事项

### 高价值数据集（建议保留）

**1. GSE234187** ⭐⭐⭐
- 样本：GSM7453693
- 描述：OS tissue of lung metastatic patient
- 置信度：1.0
- 状态：✅ **确认接受** - 骨肉瘤肺转移，患者组织

**2. GSE112852** ⭐⭐
- 样本：13个（部分）
- 描述：包含明确的 "metastatic site: lung" 标注
- 建议：查看具体样本，保留标注为肺的

**3. GSE112856** ⭐
- 样本：26个（最多）
- 建议：需要详细复核，可能有大量有用样本

### 需要排除的数据集（建议删除）

**1. GSE274484 & GSE274485** ❌
- 样本：20个
- 描述："Lung Cancer Lymph Node Metastasis"
- 问题：**肺癌转移到淋巴结**，不是转移到肺
- 建议：**全部排除**

**2. GSE179572** ❌
- 样本：3个
- 描述："brain metastasis, primary tumor: lung adenocarcinoma"
- 问题：**肺癌的脑转移**，不是肺转移
- 建议：**全部排除**

**3. GSE266330** ⚠️
- 样本：7个（V1和V2都有）
- 描述："lung cancer bone metastasis"
- 问题：可能是**肺癌骨转移**
- 建议：详细复核，可能需要排除

### 细胞系数据集（根据需求决定）

**乳腺癌肺转移细胞系：**
- GSE185647 (4个样本) - LM2细胞系
- GSE129646 (2个样本) - LM2变体
- GSE129645 (5个样本) - 乳腺癌相关

**肾癌肺转移细胞系：**
- GSE272563 (2个样本) - 786-O-LM

**决策依据：**
- 接受细胞系？ → 如果研究分子机制
- 仅要患者样本？ → 排除所有细胞系

---

## 🎯 推荐复核优先级

### 第一优先级：高置信度新增样本（21个）

**预计时间：30分钟**

1. 使用工具快速浏览
   ```bash
   python3 review_v2_results.py
   复核> quick_review
   ```

2. 打开CSV详细查看
   ```bash
   open V2_New_Samples_High_Confidence.csv
   ```

3. 重点确认以下数据集：
   - GSE234187 (目标数据集) ✓
   - GSE112852 (肺转移明确)
   - GSE112856 (26个样本)
   - 乳腺癌相关数据集

### 第二优先级：快速排除明显错误（约30个）

**预计时间：15分钟**

使用工具过滤并快速标记：
```bash
复核> filter_gse GSE274484
# 确认全是淋巴结转移 → 标记删除

复核> filter_gse GSE274485
# 确认全是淋巴结转移 → 标记删除

复核> filter_gse GSE179572
# 确认是脑转移 → 标记删除
```

### 第三优先级：中等置信度样本（122个）

**预计时间：1-2小时**

1. 按数据集分组复核
   ```bash
   复核> filter_new
   复核> filter_conf 0.6 0.79
   复核> list 50
   ```

2. 对于每个数据集：
   - 访问GEO主页
   - 查看Overall design
   - 快速判断接受/拒绝

3. 不确定的标记为Uncertain，留待后续查论文

---

## 📝 创建最终筛选结果

### 方法1：使用工具导出

```bash
python3 review_v2_results.py

# 示例：只导出高置信度新增样本
复核> filter_new
复核> filter_conf 0.8
复核> export final_high_confidence.csv

# 导出所有V2新增
复核> reset
复核> filter_new
复核> export final_v2_new_all.csv
```

### 方法2：在Excel中手工筛选

1. 打开 `GEO_Lung_Metastasis_Mining_Results_V2.csv`
2. 添加一列 `Final_Decision` (Accept/Reject)
3. 使用筛选功能：
   - Found_In_V1 = False (新增样本)
   - Confidence >= 0.8 (高置信度)
4. 逐个标记 Final_Decision
5. 筛选出 Accept 的样本
6. 另存为 `Final_Accepted_Samples.csv`

### 方法3：使用Python脚本批量处理

创建一个简单的筛选脚本：

```python
import pandas as pd

# 读取V2结果
df = pd.read_csv('GEO_Lung_Metastasis_Mining_Results_V2.csv')

# 定义排除列表
exclude_gse = ['GSE274484', 'GSE274485', 'GSE179572']  # 肺癌转移出去
exclude_keywords = ['lymph node metastasis', 'brain metastasis']

# 应用过滤
df_filtered = df[~df['GSE'].isin(exclude_gse)]
for keyword in exclude_keywords:
    df_filtered = df_filtered[~df_filtered['Title'].str.contains(keyword, case=False, na=False)]

# 可选：只保留高置信度
df_high_conf = df_filtered[df_filtered['Confidence'] >= 0.7]

# 保存
df_high_conf.to_csv('Final_Filtered_Results.csv', index=False)
print(f"筛选后剩余 {len(df_high_conf)} 个样本")
```

---

## 🔗 复核时的有用链接

对于每个样本，可以访问：

**1. GEO数据集页面**
```
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE编号
```
查看：
- Title and Summary
- Overall design
- Citation (关联论文)

**2. GEO样本页面**
```
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM编号
```
查看：
- Sample characteristics
- Extract protocol
- Data processing

**3. SRA实验页面**
```
https://www.ncbi.nlm.nih.gov/sra/?term=SRX编号
```
查看：
- SRR列表
- 测序平台
- 数据大小

**4. PubMed论文**
在GEO页面底部找到关联的论文，查看：
- Abstract
- Methods - Sample collection
- Supplementary tables

---

## 📋 复核检查清单

为每个样本确认：

- [ ] **标题清晰** - Title明确描述了样本来源
- [ ] **方向正确** - 是转移**到**肺，不是肺癌转移**出去**
- [ ] **原发确认** - 能识别原发癌症类型（非肺癌）
- [ ] **样本类型** - 明确是患者组织还是细胞系
- [ ] **测序策略** - RNA-Seq/ATAC-seq/ChIP-Seq等符合需求
- [ ] **数据可用** - SRX可以链接到SRR
- [ ] **符合排除标准** - 不是细胞系/PDX（如果排除的话）

---

## 💡 快速决策指南

### ✅ 立即接受的特征

- Title包含：
  - "[Cancer type] lung metastasis"
  - "lung metastatic [cancer type]"
  - "metastatic lung tumor from [origin]"
  
- Source明确标注：
  - "lung tissue"
  - "lung tumor"
  - "metastatic lung"

- Characteristics包含：
  - "metastatic site: lung"
  - "tissue: lung metastasis"
  - "origin: [non-lung cancer]"

### ❌ 立即拒绝的特征

- Title包含：
  - "lung cancer [organ] metastasis"
  - "[organ] metastasis from lung"
  - "NSCLC/SCLC/adenocarcinoma [organ] met"

- Source包含：
  - "lymph node"
  - "brain metastasis"
  - Primary tumor: lung

### ⚠️ 需要详细确认的特征

- 包含"pleural effusion"（胸腔积液）
- 包含"cell line"（细胞系）
- 包含"PDX"或"xenograft"
- Title模糊不清
- 只提到"lung"和"metastasis"但无方向

---

## 🎓 复核技巧

### 技巧1：批量处理同一数据集

如果一个GSE有多个样本，先查看GSE主页的Overall design，一次性判断整个数据集是否符合，然后批量处理所有样本。

### 技巧2：使用搜索功能

在GEO页面按Ctrl+F搜索关键词：
- "primary"
- "origin"
- "metastasis"
- "lung"
- 癌症类型（"osteosarcoma", "breast"等）

### 技巧3：查看论文摘要

不确定时，访问关联的PubMed论文，摘要通常会明确说明：
- 研究的癌症类型
- 样本来源
- 是否包含肺转移

### 技巧4：记录不确定的样本

创建一个"Uncertain"列表，包含：
- GSE/GSM编号
- 不确定的原因
- 需要进一步查证的内容

后续有时间再详细复核。

---

## 📊 复核进度追踪

建议创建一个简单的追踪表格：

| 日期 | 复核内容 | 样本数 | 接受 | 拒绝 | 待定 | 备注 |
|------|---------|--------|------|------|------|------|
| 11/29 | 高置信度新增 | 21 | ? | ? | ? | 第一优先级 |
| 11/29 | 明显错误排除 | 30 | 0 | 30 | 0 | 肺癌转移出去 |
| 待定 | 中等置信度 | 122 | ? | ? | ? | 需要详细复核 |

---

## 🚨 遇到问题时

### 问题1：不确定是否是肺转移

**解决方法：**
1. 查看论文摘要
2. 查看Supplementary材料
3. 标记为"Uncertain"，先跳过
4. 优先确认明确的样本

### 问题2：CSV太大，Excel打开慢

**解决方法：**
1. 使用交互式工具 `review_v2_results.py`
2. 或先用工具导出子集
3. 或使用Google Sheets（在线）

### 问题3：GEO页面信息不全

**解决方法：**
1. 查看关联论文
2. 查看SRA页面的描述
3. 如果仍不清楚，标记"Uncertain"

---

## ✅ 完成后

1. **统计最终结果**
   ```bash
   wc -l Final_Accepted_Samples.csv
   ```

2. **提取SRR列表**（如果CSV包含SRX）
   ```python
   # 如果需要从SRX获取SRR，使用原脚本的get_srr_from_srx功能
   ```

3. **准备下载**
   将接受的样本传递给下载脚本

4. **记录复核日志**
   - 复核日期
   - 接受/拒绝的数量
   - 主要排除的数据集
   - 特别注意事项

---

## 📞 需要帮助？

如果复核过程中遇到技术问题或需要进一步分析特定数据集，可以：

1. 使用交互式工具的detail命令查看详情
2. 导出特定子集单独分析
3. 查阅GEO和SRA的官方文档

---

**祝复核顺利！** 🎉


