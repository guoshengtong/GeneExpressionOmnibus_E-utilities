"""
配置文件 - GEO 小鼠骨髓单细胞B细胞发育和TAOK基因数据挖掘
Configuration - GEO Mouse Bone Marrow Single-Cell B Cell Development and TAOK Gene Mining

搜索目标：
1. 小鼠骨髓单细胞RNA测序数据
2. B细胞发育相关
3. TAOK基因相关（可选，建议先找数据再分析）
"""

# ==================== NCBI E-utilities 配置 ====================
ENTREZ_EMAIL = "xiaotontong@outlook.com"
API_DELAY = 0.5
ENTREZ_API_KEY = None

# ==================== GEO 数据缓存配置 ====================
GEO_CACHE_DIR = "./GEO_Cache"

# ==================== 搜索策略配置 ====================

# 搜索模式：'multi_stage'（多阶段，推荐）
SEARCH_MODE = 'multi_stage'

# ==================== 技术术语 ====================
TECH_TERMS = '("scRNA-seq" OR "single cell RNA-seq" OR "single-cell RNA-seq" OR "snRNA-seq" OR "single nucleus RNA-seq")'

# ==================== 阶段1：基础搜索（小鼠骨髓单细胞）====================

# 基础搜索查询 - 小鼠骨髓单细胞
SEARCH_QUERY_BASE = f'''
    ({TECH_TERMS}) AND
    ("mouse" OR "Mus musculus") AND
    ("bone marrow" OR "bone-marrow" OR "BM") AND
    "Mus musculus"[Organism] AND
    ("Expression profiling by high throughput sequencing"[DataSet Type])
'''

# ==================== 阶段2：B细胞相关搜索 ====================

# B细胞相关搜索查询
SEARCH_QUERY_B_CELL = f'''
    ({TECH_TERMS}) AND
    ("mouse" OR "Mus musculus") AND
    ("bone marrow" OR "bone-marrow" OR "BM") AND
    ("B cell" OR "B-cell" OR "B lymphocyte" OR "B-cell development" OR "B cell development" OR "B lymphopoiesis") AND
    "Mus musculus"[Organism] AND
    ("Expression profiling by high throughput sequencing"[DataSet Type])
'''

# ==================== 阶段3：B细胞发育相关搜索 ====================

# B细胞发育相关搜索查询
SEARCH_QUERY_B_CELL_DEVELOPMENT = f'''
    ({TECH_TERMS}) AND
    ("mouse" OR "Mus musculus") AND
    ("bone marrow" OR "bone-marrow" OR "BM") AND
    ("B cell development" OR "B-cell development" OR "B lymphopoiesis" OR "B-cell differentiation" OR "B cell differentiation") AND
    "Mus musculus"[Organism] AND
    ("Expression profiling by high throughput sequencing"[DataSet Type])
'''

# ==================== 阶段4：TAOK基因相关搜索（可选）====================

# TAOK基因相关搜索查询（注意：GEO中可能不会在标题中直接提到基因名）
SEARCH_QUERY_TAOK = f'''
    ({TECH_TERMS}) AND
    ("mouse" OR "Mus musculus") AND
    ("B cell" OR "B-cell" OR "B lymphocyte" OR "immune" OR "bone marrow") AND
    ("TAOK" OR "Taok" OR "TAO kinase" OR "thousand and one amino acid kinase") AND
    "Mus musculus"[Organism] AND
    ("Expression profiling by high throughput sequencing"[DataSet Type])
'''

# ==================== 多阶段搜索配置 ====================

# 定义所有搜索阶段
SEARCH_STAGES = {
    'base': {
        'name': '基础搜索（小鼠骨髓单细胞）',
        'query': SEARCH_QUERY_BASE
    },
    'b_cell': {
        'name': 'B细胞相关搜索',
        'query': SEARCH_QUERY_B_CELL
    },
    'b_cell_development': {
        'name': 'B细胞发育相关搜索',
        'query': SEARCH_QUERY_B_CELL_DEVELOPMENT
    },
    'taok': {
        'name': 'TAOK基因相关搜索（可选）',
        'query': SEARCH_QUERY_TAOK
    }
}

# ==================== 过滤器配置 ====================

# 过滤模式：'enhanced'（增强，推荐）
FILTER_MODE = 'enhanced'

# 置信度阈值
CONFIDENCE_THRESHOLDS = {
    'high': 0.8,        # 高置信度，自动接受
    'medium': 0.5,      # 中等置信度，需要复核
    'low': 0.3          # 低置信度，自动拒绝
}

# ==================== 关键词配置 ====================

# 必须包含的关键词（小鼠相关）
REQUIRED_MOUSE_KEYWORDS = [
    'mouse', 'mus musculus', 'murine'
]

# 必须包含的关键词（骨髓相关）
REQUIRED_BONE_MARROW_KEYWORDS = [
    'bone marrow', 'bone-marrow', 'bm'
]

# 必须包含的关键词（单细胞相关）
REQUIRED_SINGLE_CELL_KEYWORDS = [
    'scRNA-seq', 'single cell', 'single-cell', 'scRNAseq'
]

# B细胞相关关键词（加分项）
B_CELL_KEYWORDS = [
    'b cell', 'b-cell', 'b lymphocyte', 'b lymphopoiesis',
    'b cell development', 'b-cell development', 'b cell differentiation',
    'b-cell differentiation', 'b cell maturation', 'b-cell maturation',
    'pro-b', 'pre-b', 'immature b', 'mature b', 'plasma cell'
]

# B细胞发育阶段关键词
B_CELL_DEVELOPMENT_STAGES = [
    'pro-b', 'pre-b', 'immature b', 'mature b', 'plasma cell',
    'b lymphopoiesis', 'b cell development', 'b-cell development'
]

# TAOK基因相关关键词（可选）
TAOK_KEYWORDS = [
    'taok', 'tao kinase', 'thousand and one amino acid kinase',
    'taok1', 'taok2', 'taok3'
]

# 需要排除的关键词
EXCLUDE_KEYWORDS = [
    'human', 'homo sapiens',  # 排除人类
    'cell line', 'cell culture',  # 排除细胞系（除非是原代细胞）
    'tumor', 'cancer', 'carcinoma',  # 排除肿瘤（除非明确是B细胞发育研究）
    'in vitro', 'ivt'  # 排除体外培养（除非是原代细胞）
]

# ==================== 输出配置 ====================

OUTPUT_CSV = "GEO_Mouse_Bone_Marrow_Bcell_TAOK_Mining_Results.csv"
SRR_LIST_FILE = "SRR_accession_list_mouse_bcell.txt"
LOG_LEVEL = "INFO"

# 置信度分组输出
OUTPUT_BY_CONFIDENCE = {
    'high': 'Results_Mouse_Bcell_High_Confidence.csv',
    'medium': 'Results_Mouse_Bcell_Needs_Review.csv',
    'low': 'Results_Mouse_Bcell_Low_Confidence.csv'
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

