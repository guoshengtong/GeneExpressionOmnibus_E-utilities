#!/usr/bin/env python3
"""
运行 V2 改进版挖掘并对比结果
Run V2 Mining and Compare Results

这个脚本实现V2的核心改进功能，并与V1结果对比
"""

import pandas as pd
from Bio import Entrez
import GEOparse
import time
import re
import io
import os
import logging
from datetime import datetime
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import config_v2 as config

# 确保logs目录存在
logs_dir = Path(__file__).parent.parent / 'logs'
logs_dir.mkdir(exist_ok=True)

Entrez.email = config.ENTREZ_EMAIL

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(
            str(Path(__file__).parent.parent / 'logs' / f'geo_mining_v2_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
        ),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class ImprovedMiner:
    """改进的挖掘器"""
    
    def __init__(self):
        self.stats = {
            'stage1_strict': 0,
            'stage2_loose': 0,
            'stage3_cancer': 0,
            'total_gse': 0,
            'gsm_analyzed': 0,
            'gsm_relevant': 0
        }
    
    def multi_stage_search(self):
        """多阶段搜索"""
        all_gse = set()
        
        logger.info("="*70)
        logger.info("开始多阶段搜索")
        logger.info("="*70)
        
        # 阶段1：严格搜索（V1原始）
        logger.info("\n[阶段1] 严格搜索（V1逻辑）...")
        gse_strict = self._search_stage('strict', config.SEARCH_QUERY_STRICT)
        all_gse.update(gse_strict)
        self.stats['stage1_strict'] = len(gse_strict)
        logger.info(f"  找到: {len(gse_strict)} 个数据集")
        
        # 阶段2：宽松搜索
        logger.info("\n[阶段2] 宽松搜索（OR逻辑）...")
        gse_loose = self._search_stage('loose', config.SEARCH_QUERY_LOOSE)
        new_loose = gse_loose - all_gse
        all_gse.update(gse_loose)
        self.stats['stage2_loose'] = len(new_loose)
        logger.info(f"  找到: {len(gse_loose)} 个数据集")
        logger.info(f"  新增: {len(new_loose)} 个数据集")
        
        # 阶段3：癌症特异性搜索
        logger.info("\n[阶段3] 癌症特异性搜索...")
        cancer_gse = set()
        for cancer_type, query in config.CANCER_SPECIFIC_SEARCHES.items():
            logger.info(f"  - {cancer_type}...")
            gse_cancer = self._search_stage(cancer_type, query)
            new_cancer = gse_cancer - all_gse
            cancer_gse.update(gse_cancer)
            all_gse.update(gse_cancer)
            logger.info(f"    找到: {len(gse_cancer)} 个, 新增: {len(new_cancer)} 个")
        
        self.stats['stage3_cancer'] = len(cancer_gse - gse_loose - gse_strict)
        self.stats['total_gse'] = len(all_gse)
        
        logger.info(f"\n总计找到: {len(all_gse)} 个唯一数据集")
        
        return list(all_gse), {
            'strict': list(gse_strict),
            'loose': list(new_loose),
            'cancer': list(cancer_gse - gse_loose - gse_strict)
        }
    
    def _search_stage(self, stage_name, query):
        """执行单个搜索阶段"""
        try:
            handle = Entrez.esearch(db='gds', term=query, retmax=config.MAX_SEARCH_RESULTS, usehistory='y')
            search_results = Entrez.read(handle)
            handle.close()
            
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]
            count = int(search_results["Count"])
            
            # 获取 GSE 编号
            gse_list = []
            batch_size = config.ESUMMARY_BATCH_SIZE
            
            for start in range(0, count, batch_size):
                time.sleep(config.API_DELAY)
                try:
                    handle = Entrez.esummary(
                        db='gds',
                        webenv=webenv,
                        query_key=query_key,
                        retstart=start,
                        retmax=batch_size
                    )
                    summaries = Entrez.read(handle)
                    handle.close()
                    
                    # 处理不同格式
                    docs = []
                    if 'DocumentSummarySet' in summaries and 'DocumentSummary' in summaries['DocumentSummarySet']:
                        docs = summaries['DocumentSummarySet']['DocumentSummary']
                    elif isinstance(summaries, list):
                        docs = summaries
                    
                    for doc in docs:
                        accession = doc.get('Accession')
                        if accession and accession.startswith("GSE"):
                            gse_list.append(accession)
                
                except Exception as e:
                    logger.error(f"Error in batch {start}: {e}")
            
            return set(gse_list)
        
        except Exception as e:
            logger.error(f"Search stage {stage_name} failed: {e}")
            return set()
    
    def enhanced_filter(self, gsm_metadata, gse_id):
        """增强的过滤器"""
        # 1. 检查知识库
        if config.ENABLE_KNOWLEDGE_BASE and gse_id in config.KNOWN_DATASETS:
            kb_result = self._check_knowledge_base(gse_id, gsm_metadata)
            if kb_result:
                return True, 1.0, kb_result
        
        # 2. 标准检查
        text_to_analyze = (
            gsm_metadata.get('title', [''])[0] + " " +
            gsm_metadata.get('source_name_ch1', [''])[0] + " " +
            " ".join(gsm_metadata.get('characteristics_ch1', [])) + " " +
            " ".join(gsm_metadata.get('description', []))
        ).lower()
        
        # 检查人类
        organism = gsm_metadata.get('organism_ch1', [''])[0].lower()
        if 'homo sapiens' not in organism and 'human' not in organism:
            return False, 0.0, "Not human sample"
        
        # 检查肺部和转移
        has_lung = "lung" in text_to_analyze or "pulmonary" in text_to_analyze
        has_metastasis = "metastasis" in text_to_analyze or "metastatic" in text_to_analyze
        
        # 如果有明确证据
        if has_lung and has_metastasis:
            # 检查是否排除原发肺癌
            if re.search(r'(primary site:\s*lung|origin:\s*lung|primary lung tumor)', text_to_analyze):
                return False, 0.0, "Primary lung cancer"
            
            # 检查已知原发部位
            for site in config.KNOWN_PRIMARY_SITES:
                if site in text_to_analyze:
                    return True, 0.9, f"Detected {site} cancer with lung metastasis"
            
            return True, 0.7, "Lung metastasis (needs review)"
        
        # 3. 癌症类型推断
        detected_cancer = None
        for cancer_type, pattern in config.CANCER_TYPE_PATTERNS.items():
            for keyword in pattern['keywords']:
                if keyword.lower() in text_to_analyze:
                    detected_cancer = cancer_type
                    break
            if detected_cancer:
                break
        
        if detected_cancer:
            # 骨肉瘤特殊处理
            if detected_cancer in ['osteosarcoma']:
                if has_metastasis or has_lung:
                    return True, 0.6, f"{detected_cancer} with possible lung involvement (review required)"
        
        return False, 0.0, "No lung metastasis evidence"
    
    def _check_knowledge_base(self, gse_id, gsm_metadata):
        """检查知识库"""
        kb_entry = config.KNOWN_DATASETS.get(gse_id)
        if not kb_entry:
            return None
        
        gsm_id = gsm_metadata.get('geo_accession', [''])[0]
        
        # 检查样本ID
        if gsm_id in kb_entry.get('lung_metastasis_samples', []):
            return f"Knowledge base confirmed: {kb_entry.get('notes', 'Confirmed lung metastasis')}"
        
        # 检查患者编号
        if kb_entry.get('lung_metastasis_patients'):
            title = gsm_metadata.get('title', [''])[0]
            # 简单提取患者ID（如 BC2, BC3 等）
            import re
            match = re.match(r'^([A-Z]+\d+)', title)
            if match:
                patient_id = match.group(1)
                if patient_id in kb_entry['lung_metastasis_patients']:
                    return f"Knowledge base confirmed: Patient {patient_id} is lung metastasis"
        
        return None
    
    def analyze_gse(self, gse_id):
        """分析单个GSE"""
        logger.info(f"Parsing {gse_id}...")
        
        try:
            gse = GEOparse.get_GEO(geo=gse_id, destdir=config.GEO_CACHE_DIR, silent=True)
        except Exception as e:
            logger.error(f"Error processing {gse_id}: {e}")
            return []
        
        relevant_samples = []
        
        for gsm_name, gsm in gse.gsms.items():
            self.stats['gsm_analyzed'] += 1
            
            is_relevant, confidence, reason = self.enhanced_filter(gsm.metadata, gse_id)
            
            if is_relevant:
                # 查找 SRA 链接
                srx_id = None
                for relation in gsm.metadata.get('relation', []):
                    if 'SRA:' in relation or 'sra.cgi' in relation:
                        match = re.search(r'(SRX\d+)', relation)
                        if match:
                            srx_id = match.group(1)
                            break
                
                if srx_id:
                    self.stats['gsm_relevant'] += 1
                    relevant_samples.append({
                        "GSE": gse_id,
                        "GSM": gsm_name,
                        "SRX": srx_id,
                        "Library_Strategy": gsm.metadata.get('library_strategy', [''])[0],
                        "Title": gsm.metadata.get('title', [''])[0],
                        "Characteristics": "; ".join(gsm.metadata.get('characteristics_ch1', [])),
                        "Source_Name": gsm.metadata.get('source_name_ch1', [''])[0],
                        "Filter_Reason": reason,
                        "Confidence": confidence,
                        "Found_In_V1": False  # 稍后更新
                    })
        
        return relevant_samples
    
    def run(self):
        """运行完整挖掘"""
        logger.info("="*70)
        logger.info("GEO 肺部转移瘤挖掘 V2")
        logger.info("="*70)
        
        # 多阶段搜索
        all_gse, stage_details = self.multi_stage_search()
        
        # 分析所有数据集
        logger.info("\n" + "="*70)
        logger.info(f"开始分析 {len(all_gse)} 个数据集")
        logger.info("="*70)
        
        all_samples = []
        for idx, gse_id in enumerate(all_gse, 1):
            logger.info(f"\nProgress: {idx}/{len(all_gse)} - {gse_id}")
            samples = self.analyze_gse(gse_id)
            if samples:
                logger.info(f"  -> Found {len(samples)} relevant samples")
                all_samples.extend(samples)
            time.sleep(config.API_DELAY)
        
        return pd.DataFrame(all_samples), stage_details


def compare_with_v1(v2_results):
    """对比V1和V2结果"""
    logger.info("\n" + "="*70)
    logger.info("对比 V1 和 V2 结果")
    logger.info("="*70)
    
    # 读取V1结果
    try:
        v1_df = pd.read_csv('GEO_Lung_Metastasis_Mining_Results.csv')
        v1_gsm = set(v1_df['GSM'].values)
        v1_gse = set(v1_df['GSE'].values)
        logger.info(f"\nV1 结果: {len(v1_gsm)} 个样本, {len(v1_gse)} 个数据集")
    except FileNotFoundError:
        logger.warning("未找到V1结果文件")
        v1_gsm = set()
        v1_gse = set()
    
    # V2结果
    v2_gsm = set(v2_results['GSM'].values)
    v2_gse = set(v2_results['GSE'].values)
    logger.info(f"V2 结果: {len(v2_gsm)} 个样本, {len(v2_gse)} 个数据集")
    
    # 标记哪些在V1中找到
    v2_results['Found_In_V1'] = v2_results['GSM'].isin(v1_gsm)
    
    # 新增的样本
    new_gsm = v2_gsm - v1_gsm
    new_gse = v2_gse - v1_gse
    
    logger.info(f"\n新增样本: {len(new_gsm)} 个")
    logger.info(f"新增数据集: {len(new_gse)} 个")
    
    if new_gse:
        logger.info(f"\n新增的数据集:")
        for gse in sorted(new_gse):
            count = len(v2_results[v2_results['GSE'] == gse])
            logger.info(f"  {gse}: {count} 个样本")
    
    # 生成新增样本报告
    new_samples = v2_results[~v2_results['Found_In_V1']]
    
    if len(new_samples) > 0:
        # 按置信度分组
        high_conf = new_samples[new_samples['Confidence'] >= 0.8]
        med_conf = new_samples[(new_samples['Confidence'] >= 0.5) & (new_samples['Confidence'] < 0.8)]
        low_conf = new_samples[new_samples['Confidence'] < 0.5]
        
        logger.info(f"\n新增样本置信度分布:")
        logger.info(f"  高置信度 (≥0.8): {len(high_conf)} 个")
        logger.info(f"  中等置信度 (0.5-0.8): {len(med_conf)} 个")
        logger.info(f"  低置信度 (<0.5): {len(low_conf)} 个")
        
        return new_samples, {
            'high': high_conf,
            'medium': med_conf,
            'low': low_conf
        }
    
    return new_samples, None


def main():
    """主函数"""
    print("\n" + "="*70)
    print("运行 V2 改进版数据挖掘")
    print("="*70)
    print()
    
    # 创建挖掘器
    miner = ImprovedMiner()
    
    # 运行挖掘
    v2_results, stage_details = miner.run()
    
    # 保存V2结果
    logger.info("\n保存结果...")
    v2_results.to_csv('GEO_Lung_Metastasis_Mining_Results_V2.csv', index=False)
    logger.info(f"V2完整结果已保存: GEO_Lung_Metastasis_Mining_Results_V2.csv")
    
    # 对比V1和V2
    new_samples, conf_groups = compare_with_v1(v2_results)
    
    if len(new_samples) > 0:
        # 保存新增样本
        new_samples.to_csv('V2_New_Samples.csv', index=False)
        logger.info(f"\n新增样本已保存: V2_New_Samples.csv")
        
        # 按置信度保存
        if conf_groups:
            if len(conf_groups['high']) > 0:
                conf_groups['high'].to_csv('V2_New_Samples_High_Confidence.csv', index=False)
                logger.info(f"高置信度新增样本: V2_New_Samples_High_Confidence.csv")
            
            if len(conf_groups['medium']) > 0:
                conf_groups['medium'].to_csv('V2_New_Samples_Needs_Review.csv', index=False)
                logger.info(f"需要复核的新增样本: V2_New_Samples_Needs_Review.csv")
    
    # 打印统计
    logger.info("\n" + "="*70)
    logger.info("挖掘统计")
    logger.info("="*70)
    logger.info(f"搜索阶段1（严格）: {miner.stats['stage1_strict']} 个数据集")
    logger.info(f"搜索阶段2（宽松）: 新增 {miner.stats['stage2_loose']} 个数据集")
    logger.info(f"搜索阶段3（癌症）: 新增 {miner.stats['stage3_cancer']} 个数据集")
    logger.info(f"总数据集: {miner.stats['total_gse']} 个")
    logger.info(f"分析样本数: {miner.stats['gsm_analyzed']} 个")
    logger.info(f"符合条件样本: {miner.stats['gsm_relevant']} 个")
    
    logger.info("\n完成！")
    
    return v2_results, new_samples


if __name__ == "__main__":
    try:
        v2_results, new_samples = main()
        print("\n" + "="*70)
        print("✅ V2挖掘完成！")
        print("="*70)
        print("\n生成的文件:")
        print("  - GEO_Lung_Metastasis_Mining_Results_V2.csv  (V2完整结果)")
        print("  - V2_New_Samples.csv                         (V2新增样本)")
        print("  - V2_New_Samples_High_Confidence.csv         (高置信度)")
        print("  - V2_New_Samples_Needs_Review.csv            (需要复核)")
        print()
    except KeyboardInterrupt:
        print("\n\n用户中断")
    except Exception as e:
        logger.error(f"错误: {e}")
        import traceback
        traceback.print_exc()


