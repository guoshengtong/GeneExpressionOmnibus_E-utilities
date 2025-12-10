"""
配置文件 - GEO 肺部转移瘤数据挖掘流水线
Configuration file for GEO Lung Metastasis Mining Pipeline
"""

# ==================== NCBI E-utilities 配置 ====================
# 重要：请填入您的邮箱地址，NCBI以此追踪使用情况
# IMPORTANT: Please provide your email address for NCBI tracking
ENTREZ_EMAIL = "xiaotontong@outlook.com"

# NCBI API 限制通常为每秒3次请求
# NCBI API rate limit is typically 3 requests per second
API_DELAY = 0.5  # seconds

# 如果您有NCBI API Key，可以提高请求速率到每秒10次
# If you have an NCBI API Key, you can increase to 10 requests per second
ENTREZ_API_KEY = None  # Set to your API key string if available

# ==================== GEO 数据缓存配置 ====================
# GEO元数据缓存目录
GEO_CACHE_DIR = "./GEO_Cache"

# ==================== 搜索配置 ====================
# 技术术语 (scRNA-seq, snRNA-seq, 空间转录组, ATAC-seq)
TECH_TERMS = '("scRNA-seq" OR "single cell RNA-seq" OR "snRNA-seq" OR "single nucleus RNA-seq" OR "spatial transcriptomics" OR "ATAC-seq")'

# 生物学背景术语
BIOLOGY_TERMS = '(("lung" OR "pulmonary") AND ("metastasis" OR "secondary tumor"))'

# 基础过滤条件 (人类 和 高通量测序数据类型)
BASE_FILTERS = '"Homo sapiens"[Organism] AND ("Expression profiling by high throughput sequencing"[DataSet Type] OR "Genome binding/occupancy profiling by high throughput sequencing"[DataSet Type])'

# 完整搜索查询
SEARCH_QUERY = f"({TECH_TERMS}) AND ({BIOLOGY_TERMS}) AND ({BASE_FILTERS})"

# 最大返回结果数
MAX_SEARCH_RESULTS = 5000

# ==================== 过滤规则配置 ====================
# 已知的癌症原发部位（用于识别转移瘤来源）
KNOWN_PRIMARY_SITES = [
    "breast", "colon", "colorectal", "crc", "melanoma", 
    "kidney", "renal", "pancreas", "liver", "sarcoma", 
    "ovarian", "prostate", "gastric", "stomach", "esophageal",
    "thyroid", "bladder", "cervical", "uterine", "testicular"
]

# 原发性肺癌术语（需要排除）
PRIMARY_LUNG_CANCER_TERMS = [
    "nsclc", "sclc", "lung adenocarcinoma", 
    "squamous cell lung carcinoma", "lung squamous cell carcinoma"
]

# ==================== 输出配置 ====================
# 结果CSV文件名
OUTPUT_CSV = "GEO_Lung_Metastasis_Mining_Results.csv"

# SRR列表文件名
SRR_LIST_FILE = "SRR_accession_list.txt"

# 日志级别 (DEBUG, INFO, WARNING, ERROR, CRITICAL)
LOG_LEVEL = "INFO"

# ==================== SRA Toolkit 配置 ====================
# SRA Toolkit prefetch 批处理大小
SRA_BATCH_SIZE = 5

# SRA 下载目录
SRA_DOWNLOAD_DIR = "./SRA_Data"

# ==================== 高级配置 ====================
# ESummary批处理大小
ESUMMARY_BATCH_SIZE = 500

# SRA ESearch批处理大小
SRA_SEARCH_BATCH_SIZE = 200

