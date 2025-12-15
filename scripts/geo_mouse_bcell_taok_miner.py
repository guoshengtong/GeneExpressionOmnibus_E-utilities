#!/usr/bin/env python3
"""
GEO 小鼠骨髓单细胞B细胞发育和TAOK基因数据挖掘流水线
GEO Mouse Bone Marrow Single-Cell B Cell Development and TAOK Gene Mining Pipeline

主要功能：
1. 多阶段搜索策略（基础+B细胞+B细胞发育+TAOK）
2. 智能过滤器（置信度评分）
3. 详细的可追溯性日志
4. B细胞发育阶段识别

作者：AI Assistant
版本：1.0
日期：2025-12-10
"""

import pandas as pd
from Bio import Entrez
import GEOparse
import time
import re
import os
import logging
from datetime import datetime
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
# 添加config目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent / 'config'))
import config_mouse_bcell_taok as config

# 确保logs目录存在
logs_dir = Path(__file__).parent.parent / 'logs'
logs_dir.mkdir(exist_ok=True)

Entrez.email = config.ENTREZ_EMAIL
if config.ENTREZ_API_KEY:
    Entrez.api_key = config.ENTREZ_API_KEY

# 配置日志
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(
            str(Path(__file__).parent.parent / 'logs' / f'mouse_bcell_taok_mining_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
        ),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class MouseBcellTAOKMiner:
    """小鼠骨髓B细胞发育和TAOK基因挖掘器"""
    
    def __init__(self):
        self.stats = {
            'stage_base': 0,
            'stage_b_cell': 0,
            'stage_b_cell_development': 0,
            'stage_taok': 0,
            'total_gse': 0,
            'gsm_analyzed': 0,
            'gsm_relevant': 0
        }
    
    def multi_stage_search(self):
        """多阶段搜索"""
        all_gse = set()
        stage_details = {}
        
        logger.info("="*70)
        logger.info("开始多阶段搜索 - 小鼠骨髓单细胞B细胞发育和TAOK基因")
        logger.info("="*70)
        
        # 遍历所有搜索阶段
        for stage_key, stage_config in config.SEARCH_STAGES.items():
            logger.info(f"\n[{stage_key.upper()}] {stage_config['name']}...")
            gse_stage = self._search_stage(stage_key, stage_config['query'])
            new_gse = gse_stage - all_gse
            all_gse.update(gse_stage)
            self.stats[f'stage_{stage_key}'] = len(new_gse)
            stage_details[stage_key] = list(gse_stage)
            logger.info(f"  找到: {len(gse_stage)} 个数据集")
            logger.info(f"  新增: {len(new_gse)} 个数据集")
        
        self.stats['total_gse'] = len(all_gse)
        
        logger.info(f"\n总计找到: {len(all_gse)} 个唯一数据集")
        
        return list(all_gse), stage_details
    
    def _search_stage(self, stage_name, query):
        """执行单个搜索阶段"""
        try:
            if config.LOG_SEARCH_QUERIES:
                logger.debug(f"搜索查询: {query}")
            
            handle = Entrez.esearch(
                db='gds',
                term=query,
                retmax=config.MAX_SEARCH_RESULTS,
                usehistory='y'
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]
            count = int(search_results["Count"])
            
            logger.info(f"  搜索结果: {count} 条记录")
            
            if count == 0:
                return set()
            
            # 获取 GSE 编号
            gse_list = []
            batch_size = config.ESUMMARY_BATCH_SIZE
            
            for start in range(0, min(count, config.MAX_SEARCH_RESULTS), batch_size):
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
        """增强的过滤器 - 针对小鼠骨髓B细胞发育"""
        # 准备分析文本
        text_to_analyze = (
            gsm_metadata.get('title', [''])[0] + " " +
            gsm_metadata.get('source_name_ch1', [''])[0] + " " +
            " ".join(gsm_metadata.get('characteristics_ch1', [])) + " " +
            " ".join(gsm_metadata.get('description', []))
        ).lower()
        
        # 检查小鼠
        organism = gsm_metadata.get('organism_ch1', [''])[0].lower()
        if not any(keyword in organism for keyword in config.REQUIRED_MOUSE_KEYWORDS):
            return False, 0.0, "Not mouse sample"
        
        # 检查是否排除
        for exclude_keyword in config.EXCLUDE_KEYWORDS:
            if exclude_keyword in text_to_analyze:
                # 特殊处理：如果是B细胞发育研究中的肿瘤，可能是相关的
                if exclude_keyword in ['tumor', 'cancer', 'carcinoma']:
                    # 检查是否同时包含B细胞发育关键词
                    if any(bcell_keyword in text_to_analyze for bcell_keyword in config.B_CELL_KEYWORDS):
                        # 可能是B细胞淋巴瘤或B细胞发育异常研究，保留但降低置信度
                        pass
                    else:
                        return False, 0.0, f"Excluded keyword: {exclude_keyword}"
                else:
                    return False, 0.0, f"Excluded keyword: {exclude_keyword}"
        
        # 检查骨髓
        has_bone_marrow = any(keyword in text_to_analyze for keyword in config.REQUIRED_BONE_MARROW_KEYWORDS)
        if not has_bone_marrow:
            return False, 0.0, "No bone marrow keyword found"
        
        # 检查单细胞
        has_single_cell = any(keyword in text_to_analyze for keyword in config.REQUIRED_SINGLE_CELL_KEYWORDS)
        if not has_single_cell:
            return False, 0.0, "No single-cell keyword found"
        
        # 计算置信度
        confidence = 0.5  # 基础置信度（满足基本条件）
        reasons = []
        
        # B细胞相关（加分）
        b_cell_matches = [keyword for keyword in config.B_CELL_KEYWORDS if keyword in text_to_analyze]
        if b_cell_matches:
            confidence += 0.2
            reasons.append(f"B cell keywords: {', '.join(b_cell_matches[:3])}")
        
        # B细胞发育阶段（加分）
        development_matches = [keyword for keyword in config.B_CELL_DEVELOPMENT_STAGES if keyword in text_to_analyze]
        if development_matches:
            confidence += 0.2
            reasons.append(f"B cell development stages: {', '.join(development_matches[:3])}")
        
        # TAOK基因相关（加分，但可能不在标题中）
        taok_matches = [keyword for keyword in config.TAOK_KEYWORDS if keyword in text_to_analyze]
        if taok_matches:
            confidence += 0.3
            reasons.append(f"TAOK keywords: {', '.join(taok_matches)}")
        
        # 如果同时有B细胞和发育阶段，提高置信度
        if b_cell_matches and development_matches:
            confidence += 0.1
        
        # 限制置信度在0-1之间
        confidence = min(1.0, confidence)
        
        reason_str = "; ".join(reasons) if reasons else "Mouse bone marrow single-cell sample"
        
        # 根据置信度阈值决定是否接受
        if confidence >= config.CONFIDENCE_THRESHOLDS['low']:
            return True, confidence, reason_str
        else:
            return False, confidence, f"Low confidence: {reason_str}"
    
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
                srr_list = []
                
                for relation in gsm.metadata.get('relation', []):
                    if 'SRA:' in relation or 'sra.cgi' in relation:
                        match = re.search(r'(SRX\d+)', relation)
                        if match:
                            srx_id = match.group(1)
                            break
                
                # 尝试从relation中提取SRR
                for relation in gsm.metadata.get('relation', []):
                    srr_matches = re.findall(r'(SRR\d+)', relation)
                    srr_list.extend(srr_matches)
                
                # 如果没有找到SRR，尝试从supplementary_file中查找
                if not srr_list:
                    for supp_file in gsm.metadata.get('supplementary_file', []):
                        srr_matches = re.findall(r'(SRR\d+)', supp_file)
                        srr_list.extend(srr_matches)
                
                self.stats['gsm_relevant'] += 1
                
                # 提取B细胞发育阶段信息
                text_to_analyze = (
                    gsm.metadata.get('title', [''])[0] + " " +
                    " ".join(gsm.metadata.get('characteristics_ch1', []))
                ).lower()
                
                b_cell_stages = [stage for stage in config.B_CELL_DEVELOPMENT_STAGES if stage in text_to_analyze]
                has_taok = any(keyword in text_to_analyze for keyword in config.TAOK_KEYWORDS)
                
                relevant_samples.append({
                    "GSE": gse_id,
                    "GSM": gsm_name,
                    "SRX": srx_id if srx_id else "",
                    "SRR_List": "; ".join(srr_list) if srr_list else "",
                    "SRR_Count": len(srr_list),
                    "Library_Strategy": gsm.metadata.get('library_strategy', [''])[0],
                    "Title": gsm.metadata.get('title', [''])[0],
                    "Characteristics": "; ".join(gsm.metadata.get('characteristics_ch1', [])),
                    "Source_Name": gsm.metadata.get('source_name_ch1', [''])[0],
                    "Organism": gsm.metadata.get('organism_ch1', [''])[0],
                    "B_Cell_Stages": "; ".join(b_cell_stages) if b_cell_stages else "",
                    "Has_TAOK": "Yes" if has_taok else "No (may need expression analysis)",
                    "Filter_Reason": reason,
                    "Confidence": confidence
                })
        
        return relevant_samples
    
    def run(self):
        """运行完整挖掘"""
        logger.info("="*70)
        logger.info("GEO 小鼠骨髓单细胞B细胞发育和TAOK基因挖掘")
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


def save_results_by_confidence(df):
    """按置信度分组保存结果"""
    if len(df) == 0:
        logger.warning("No samples to save")
        return
    
    results_dir = Path(__file__).parent.parent / 'results'
    results_dir.mkdir(exist_ok=True)
    
    # 保存完整结果
    output_path = results_dir / config.OUTPUT_CSV
    df.to_csv(output_path, index=False)
    logger.info(f"完整结果已保存: {output_path}")
    
    # 按置信度分组
    high_conf = df[df['Confidence'] >= config.CONFIDENCE_THRESHOLDS['high']]
    med_conf = df[(df['Confidence'] >= config.CONFIDENCE_THRESHOLDS['medium']) & 
                  (df['Confidence'] < config.CONFIDENCE_THRESHOLDS['high'])]
    low_conf = df[df['Confidence'] < config.CONFIDENCE_THRESHOLDS['medium']]
    
    if len(high_conf) > 0:
        high_path = results_dir / config.OUTPUT_BY_CONFIDENCE['high']
        high_conf.to_csv(high_path, index=False)
        logger.info(f"高置信度结果: {high_path} ({len(high_conf)} 个样本)")
    
    if len(med_conf) > 0:
        med_path = results_dir / config.OUTPUT_BY_CONFIDENCE['medium']
        med_conf.to_csv(med_path, index=False)
        logger.info(f"需要复核的结果: {med_path} ({len(med_conf)} 个样本)")
    
    if len(low_conf) > 0:
        low_path = results_dir / config.OUTPUT_BY_CONFIDENCE['low']
        low_conf.to_csv(low_path, index=False)
        logger.info(f"低置信度结果: {low_path} ({len(low_conf)} 个样本)")
    
    # 保存SRR列表
    data_dir = Path(__file__).parent.parent / 'data'
    data_dir.mkdir(exist_ok=True)
    
    all_srr = []
    for srr_list_str in df['SRR_List'].values:
        if srr_list_str and srr_list_str.strip():
            srr_items = [srr.strip() for srr in srr_list_str.split(';') if srr.strip()]
            all_srr.extend(srr_items)
    
    if all_srr:
        srr_file = data_dir / config.SRR_LIST_FILE
        with open(srr_file, 'w') as f:
            for srr in sorted(set(all_srr)):
                f.write(f"{srr}\n")
        logger.info(f"SRR列表已保存: {srr_file} ({len(set(all_srr))} 个唯一SRR)")


def main():
    """主函数"""
    print("\n" + "="*70)
    print("运行 小鼠骨髓单细胞B细胞发育和TAOK基因数据挖掘")
    print("="*70)
    print()
    
    # 创建挖掘器
    miner = MouseBcellTAOKMiner()
    
    # 运行挖掘
    results_df, stage_details = miner.run()
    
    # 保存结果
    logger.info("\n" + "="*70)
    logger.info("保存结果...")
    logger.info("="*70)
    
    save_results_by_confidence(results_df)
    
    # 打印统计
    logger.info("\n" + "="*70)
    logger.info("挖掘统计")
    logger.info("="*70)
    logger.info(f"搜索阶段1（基础）: {miner.stats['stage_base']} 个数据集")
    logger.info(f"搜索阶段2（B细胞）: {miner.stats['stage_b_cell']} 个数据集")
    logger.info(f"搜索阶段3（B细胞发育）: {miner.stats['stage_b_cell_development']} 个数据集")
    logger.info(f"搜索阶段4（TAOK）: {miner.stats['stage_taok']} 个数据集")
    logger.info(f"总数据集: {miner.stats['total_gse']} 个")
    logger.info(f"分析样本数: {miner.stats['gsm_analyzed']} 个")
    logger.info(f"符合条件样本: {miner.stats['gsm_relevant']} 个")
    
    if len(results_df) > 0:
        logger.info(f"\n置信度分布:")
        logger.info(f"  高置信度 (≥{config.CONFIDENCE_THRESHOLDS['high']}): {len(results_df[results_df['Confidence'] >= config.CONFIDENCE_THRESHOLDS['high']])} 个")
        logger.info(f"  中等置信度 ({config.CONFIDENCE_THRESHOLDS['medium']}-{config.CONFIDENCE_THRESHOLDS['high']}): {len(results_df[(results_df['Confidence'] >= config.CONFIDENCE_THRESHOLDS['medium']) & (results_df['Confidence'] < config.CONFIDENCE_THRESHOLDS['high'])])} 个")
        logger.info(f"  低置信度 (<{config.CONFIDENCE_THRESHOLDS['medium']}): {len(results_df[results_df['Confidence'] < config.CONFIDENCE_THRESHOLDS['medium']])} 个")
        
        logger.info(f"\nB细胞发育阶段分布:")
        stages_df = results_df[results_df['B_Cell_Stages'] != '']
        logger.info(f"  包含B细胞发育阶段信息: {len(stages_df)} 个样本")
        
        logger.info(f"\nTAOK基因信息:")
        taok_df = results_df[results_df['Has_TAOK'] == 'Yes']
        logger.info(f"  元数据中提到TAOK: {len(taok_df)} 个样本")
        logger.info(f"  需要表达分析: {len(results_df[results_df['Has_TAOK'] == 'No (may need expression analysis)'])} 个样本")
    
    logger.info("\n完成！")
    
    return results_df, stage_details


if __name__ == "__main__":
    try:
        results_df, stage_details = main()
        print("\n" + "="*70)
        print("✅ 挖掘完成！")
        print("="*70)
        print("\n生成的文件:")
        print(f"  - results/{config.OUTPUT_CSV}  (完整结果)")
        print(f"  - results/{config.OUTPUT_BY_CONFIDENCE['high']}  (高置信度)")
        print(f"  - results/{config.OUTPUT_BY_CONFIDENCE['medium']}  (需要复核)")
        print(f"  - data/{config.SRR_LIST_FILE}  (SRR列表)")
        print()
    except KeyboardInterrupt:
        print("\n\n用户中断")
    except Exception as e:
        logger.error(f"错误: {e}")
        import traceback
        traceback.print_exc()

