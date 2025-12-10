"""
改进的配置文件 - GEO 肺部转移瘤数据挖掘流水线 V2
Enhanced Configuration - GEO Lung Metastasis Mining Pipeline V2

主要改进：
1. 多阶段搜索策略
2. 扩展的癌症类型关键词
3. 外部知识库支持
4. 置信度评分系统
"""

# ==================== NCBI E-utilities 配置 ====================
ENTREZ_EMAIL = "xiaotontong@outlook.com"
API_DELAY = 0.5
ENTREZ_API_KEY = None

# ==================== GEO 数据缓存配置 ====================
GEO_CACHE_DIR = "./GEO_Cache"

# ==================== 搜索策略配置 ====================

# 搜索模式：'strict'（严格）, 'loose'（宽松）, 'multi_stage'（多阶段，推荐）
SEARCH_MODE = 'multi_stage'

# ==================== 阶段1：严格搜索（原始逻辑，高精度）====================

# 技术术语
TECH_TERMS = '("scRNA-seq" OR "single cell RNA-seq" OR "snRNA-seq" OR "single nucleus RNA-seq" OR "spatial transcriptomics" OR "ATAC-seq")'

# 生物学术语 - 原始版本
BIOLOGY_TERMS_STRICT = '(("lung" OR "pulmonary") AND ("metastasis" OR "secondary tumor"))'

# 基础过滤条件
BASE_FILTERS = '"Homo sapiens"[Organism] AND ("Expression profiling by high throughput sequencing"[DataSet Type] OR "Genome binding/occupancy profiling by high throughput sequencing"[DataSet Type])'

# 严格搜索查询
SEARCH_QUERY_STRICT = f"({TECH_TERMS}) AND ({BIOLOGY_TERMS_STRICT}) AND ({BASE_FILTERS})"

# ==================== 阶段2：宽松搜索（高召回率）====================

# 宽松的生物学术语 - 使用 OR 逻辑
BIOLOGY_TERMS_LOOSE = '''(
    ("lung" AND "metastasis") OR
    ("lung" AND "metastatic") OR
    ("pulmonary" AND "metastasis") OR
    ("osteosarcoma" AND ("lung" OR "metastasis")) OR
    ("breast cancer" AND ("lung" OR "metastasis" OR "metastatic")) OR
    ("melanoma" AND ("lung" OR "metastasis" OR "metastatic")) OR
    ("sarcoma" AND ("lung" OR "metastasis"))
)'''

# 宽松搜索查询
SEARCH_QUERY_LOOSE = f"({TECH_TERMS}) AND ({BIOLOGY_TERMS_LOOSE}) AND ({BASE_FILTERS})"

# ==================== 阶段3：癌症特异性搜索 ====================

# 针对特定癌症类型的搜索
CANCER_SPECIFIC_SEARCHES = {
    'Osteosarcoma': f'''
        ("osteosarcoma" OR "bone sarcoma") AND
        ("lung" OR "pulmonary" OR "metastasis" OR "metastatic") AND
        ({TECH_TERMS}) AND
        "Homo sapiens"[Organism]
    ''',
    
    'Breast_Cancer': f'''
        ("breast cancer" OR "breast carcinoma" OR "breast adenocarcinoma") AND
        ("lung metastasis" OR "lung metastatic" OR "pulmonary metastasis") AND
        ({TECH_TERMS}) AND
        "Homo sapiens"[Organism]
    ''',
    
    'Melanoma': f'''
        "melanoma" AND
        ("lung metastasis" OR "lung metastatic" OR "pulmonary metastasis") AND
        ({TECH_TERMS}) AND
        "Homo sapiens"[Organism]
    ''',
    
    'Colorectal': f'''
        ("colorectal cancer" OR "colon cancer" OR "CRC") AND
        ("lung metastasis" OR "lung metastatic" OR "pulmonary metastasis") AND
        ({TECH_TERMS}) AND
        "Homo sapiens"[Organism]
    ''',
    
    'Renal': f'''
        ("renal cell carcinoma" OR "kidney cancer" OR "RCC") AND
        ("lung metastasis" OR "lung metastatic" OR "pulmonary metastasis") AND
        ({TECH_TERMS}) AND
        "Homo sapiens"[Organism]
    '''
}

# ==================== 过滤器配置 ====================

# 过滤模式：'standard'（标准）, 'enhanced'（增强，推荐）
FILTER_MODE = 'enhanced'

# 置信度阈值
CONFIDENCE_THRESHOLDS = {
    'high': 0.8,        # 高置信度，自动接受
    'medium': 0.5,      # 中等置信度，需要复核
    'low': 0.3          # 低置信度，自动拒绝
}

# ==================== 已知癌症类型配置 ====================

# 癌症类型及其特征
CANCER_TYPE_PATTERNS = {
    'osteosarcoma': {
        'keywords': ['osteosarcoma', 'bone sarcoma', 'OS tumor', 'osteogenic sarcoma'],
        'abbreviations': ['OS'],
        'primary_site': 'bone',
        'common_metastasis_sites': ['lung', 'bone'],
        'weight': 1.0
    },
    'breast': {
        'keywords': ['breast cancer', 'breast carcinoma', 'breast adenocarcinoma', 'mammary carcinoma'],
        'abbreviations': ['BC', 'BRCA'],
        'primary_site': 'breast',
        'common_metastasis_sites': ['lung', 'bone', 'liver', 'brain'],
        'weight': 1.0
    },
    'melanoma': {
        'keywords': ['melanoma', 'malignant melanoma', 'cutaneous melanoma'],
        'abbreviations': ['MEL'],
        'primary_site': 'skin',
        'common_metastasis_sites': ['lung', 'liver', 'brain', 'bone'],
        'weight': 1.0
    },
    'colorectal': {
        'keywords': ['colorectal cancer', 'colon cancer', 'rectal cancer', 'CRC'],
        'abbreviations': ['CRC'],
        'primary_site': 'colon',
        'common_metastasis_sites': ['liver', 'lung', 'peritoneum'],
        'weight': 1.0
    },
    'renal': {
        'keywords': ['renal cell carcinoma', 'kidney cancer', 'RCC'],
        'abbreviations': ['RCC'],
        'primary_site': 'kidney',
        'common_metastasis_sites': ['lung', 'bone', 'liver', 'brain'],
        'weight': 1.0
    }
}

# ==================== 外部知识库配置 ====================

# 是否启用外部知识库
ENABLE_KNOWLEDGE_BASE = True

# 已知数据集（手动维护）
KNOWN_DATASETS = {
    'GSE234187': {
        'status': 'confirmed',
        'has_lung_metastasis': True,
        'lung_metastasis_samples': ['GSM7453693'],
        'cancer_type': 'osteosarcoma',
        'confidence': 1.0,
        'source': 'Sample title annotation',
        'notes': 'Clear annotation: OS tissue of lung metastatic patient',
        'verified_date': '2025-11-29'
    },
    
    'GSE152048': {
        'status': 'needs_annotation',
        'has_lung_metastasis': True,
        'cancer_type': 'osteosarcoma',
        'lung_metastasis_patients': None,  # 需要从论文补充材料获取
        'total_samples': 76,
        'paper': {
            'title': 'Single-cell RNA landscape of intratumoral heterogeneity',
            'journal': 'Nature Communications',
            'year': 2020,
            'doi': '10.1038/s41467-020-20059-6'
        },
        'notes': 'Patient BC codes need to be mapped from supplementary materials',
        'action_required': 'Download paper supplementary table to identify lung metastasis patients'
    },
    
    'GSE270231': {
        'status': 'found_in_search',
        'has_lung_metastasis': True,
        'cancer_type': 'osteosarcoma',
        'confidence': 0.9,
        'notes': 'Already found in original search results (7 samples)'
    }
}

# ==================== 过滤规则配置 ====================

# 已知的癌症原发部位（用于识别转移瘤来源）
KNOWN_PRIMARY_SITES = [
    "breast", "colon", "colorectal", "crc", "melanoma", 
    "kidney", "renal", "pancreas", "liver", "sarcoma", "osteosarcoma",
    "ovarian", "prostate", "gastric", "stomach", "esophageal",
    "thyroid", "bladder", "cervical", "uterine", "testicular"
]

# 原发性肺癌术语（需要排除）
PRIMARY_LUNG_CANCER_TERMS = [
    "nsclc", "sclc", "lung adenocarcinoma", 
    "squamous cell lung carcinoma", "lung squamous cell carcinoma",
    "primary lung cancer", "primary lung tumor"
]

# ==================== 输出配置 ====================

OUTPUT_CSV = "GEO_Lung_Metastasis_Mining_Results_V2.csv"
SRR_LIST_FILE = "SRR_accession_list_v2.txt"
LOG_LEVEL = "INFO"

# 置信度分组输出
OUTPUT_BY_CONFIDENCE = {
    'high': 'Results_High_Confidence.csv',
    'medium': 'Results_Needs_Review.csv',
    'low': 'Results_Low_Confidence.csv'
}

# ==================== 高级配置 ====================

MAX_SEARCH_RESULTS = 5000
ESUMMARY_BATCH_SIZE = 500
SRA_SEARCH_BATCH_SIZE = 200
SRA_BATCH_SIZE = 5
SRA_DOWNLOAD_DIR = "./SRA_Data"

# 详细日志
DETAILED_LOGGING = True
LOG_FILTER_DECISIONS = True  # 记录每个样本的过滤决策
LOG_SEARCH_QUERIES = True    # 记录搜索查询

# 性能优化
CACHE_SEARCH_RESULTS = True
PARALLEL_GSE_ANALYSIS = False  # 未来功能：并行分析


