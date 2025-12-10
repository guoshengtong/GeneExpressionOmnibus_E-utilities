#!/usr/bin/env python3
"""
GEO 肺部转移瘤数据挖掘流水线
GEO Lung Metastasis Data Mining Pipeline

该脚本通过NCBI E-utilities API自动化搜索、解析和过滤GEO数据库，
识别其他原发部位转移到肺部的肿瘤样本，并关联到SRA原始数据。

This script automates the search, parsing, and filtering of GEO database
via NCBI E-utilities API to identify tumor samples that metastasized to
the lung from other primary sites, and links them to SRA raw data.
"""

import pandas as pd
from Bio import Entrez
import GEOparse
import time
import re
import io
import os
import logging
from typing import List, Dict, Tuple
from datetime import datetime

# 导入配置
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import config

# ------------------------------------------------------------------
# 初始化和配置
# ------------------------------------------------------------------

class GEOLungMetastasisMiner:
    """GEO肺部转移瘤数据挖掘器"""
    
    def __init__(self):
        """初始化挖掘器"""
        # 配置 Entrez
        Entrez.email = config.ENTREZ_EMAIL
        if config.ENTREZ_API_KEY:
            Entrez.api_key = config.ENTREZ_API_KEY
        
        # 创建缓存目录
        if not os.path.exists(config.GEO_CACHE_DIR):
            os.makedirs(config.GEO_CACHE_DIR)
        
        # 确保logs目录存在
        logs_dir = Path(__file__).parent.parent / 'logs'
        logs_dir.mkdir(exist_ok=True)
        
        # 配置日志
        logging.basicConfig(
            level=getattr(logging, config.LOG_LEVEL),
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(
                    str(logs_dir / f'geo_mining_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
                ),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
        # 统计信息
        self.stats = {
            'gse_found': 0,
            'gse_analyzed': 0,
            'gsm_analyzed': 0,
            'gsm_relevant': 0,
            'srx_found': 0,
            'srr_found': 0
        }
    
    # ------------------------------------------------------------------
    # 阶段 1: 广泛搜索 (GEO)
    # ------------------------------------------------------------------
    
    def search_geo(self, query: str = None) -> List[str]:
        """
        使用 ESearch 搜索 GEO DataSets 数据库，并通过 ESummary 获取 GSE Accessions。
        
        Args:
            query: 搜索查询字符串，默认使用config中的SEARCH_QUERY
            
        Returns:
            List[str]: GSE编号列表
        """
        if query is None:
            query = config.SEARCH_QUERY
        
        self.logger.info(f"Searching GEO with query: {query}")
        
        try:
            # db='gds' 指向 GEO DataSets
            # 使用 usehistory='y' 来高效处理大量结果
            handle = Entrez.esearch(
                db='gds', 
                term=query, 
                retmax=config.MAX_SEARCH_RESULTS, 
                usehistory='y'
            )
            search_results = Entrez.read(handle)
            handle.close()
        except Exception as e:
            self.logger.error(f"Error during GEO search: {e}")
            return []
        
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        count = int(search_results["Count"])
        self.stats['gse_found'] = count
        
        self.logger.info(f"Found {count} potential datasets (GDS IDs). Fetching Accessions...")
        
        # 使用 ESummary 获取 GSE Accession Numbers，分批次进行
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
                
                # 处理不同的返回格式
                docs = []
                if 'DocumentSummarySet' in summaries and 'DocumentSummary' in summaries['DocumentSummarySet']:
                    docs = summaries['DocumentSummarySet']['DocumentSummary']
                elif isinstance(summaries, list):
                    docs = summaries
                
                # 遍历文档并提取 GSE 编号
                for summary in docs:
                    accession = summary.get('Accession')
                    if accession and accession.startswith("GSE"):
                        gse_list.append(accession)
            except Exception as e:
                self.logger.error(f"Error fetching summaries batch starting at {start}: {e}")
                time.sleep(5)  # 出错时等待更长时间
        
        unique_gse_list = list(set(gse_list))
        self.logger.info(f"Retrieved {len(unique_gse_list)} unique GSE Accessions.")
        return unique_gse_list
    
    # ------------------------------------------------------------------
    # 阶段 2: 深度解析与精准过滤 (样本级) - 核心步骤
    # ------------------------------------------------------------------
    
    def is_lung_metastasis_of_other_origin(self, gsm_metadata: Dict) -> Tuple[bool, str]:
        """
        核心过滤逻辑：判断样本是否为其他原位肿瘤转移到肺部的样本。
        这是一个启发式规则集，需要根据实际数据不断迭代优化。
        
        Args:
            gsm_metadata: GSM样本元数据字典
            
        Returns:
            Tuple[bool, str]: (是否符合条件, 判断理由)
        """
        # 合并所有关键字段进行分析 (Title, Source Name, Characteristics, Description)
        text_to_analyze = (
            gsm_metadata.get('title', [''])[0] + " " +
            gsm_metadata.get('source_name_ch1', [''])[0] + " " +
            " ".join(gsm_metadata.get('characteristics_ch1', [])) + " " +
            " ".join(gsm_metadata.get('description', []))
        ).lower()
        
        # 1. 必须是人类样本
        organism = gsm_metadata.get('organism_ch1', [''])[0].lower()
        if 'homo sapiens' not in organism and 'human' not in organism:
            return False, "Not human sample"
        
        # 2. 位置必须在肺部
        if "lung" not in text_to_analyze and "pulmonary" not in text_to_analyze:
            return False, "Not lung tissue"
        
        # 3. 状态必须是转移/继发
        has_metastasis_keyword = False
        if "metastasis" in text_to_analyze or "metastatic" in text_to_analyze:
            has_metastasis_keyword = True
        elif "secondary" in text_to_analyze:
            has_metastasis_keyword = True
        
        if not has_metastasis_keyword:
            return False, "No metastasis/secondary keywords"
        
        # 4. 排除细胞系 (Cell line) 或类器官 (Organoid)，除非是患者来源的异种移植瘤 (PDX)
        if ("cell line" in text_to_analyze or "organoid" in text_to_analyze):
            if "pdx" not in text_to_analyze and "xenograft" not in text_to_analyze and "patient-derived" not in text_to_analyze:
                return False, "Cell line or organoid (not PDX)"
        
        # 5. 排除原发性肺癌 (关键挑战)
        
        # 5a. 寻找明确指示原发灶在别处的证据 (强证据)
        # 例如: "metastasis from breast cancer", "origin: colon"
        known_origins_pattern = '|'.join(config.KNOWN_PRIMARY_SITES)
        known_origins_regex = rf'(from|origin:|primary site:|primary tumor:)\s*({known_origins_pattern})'
        match = re.search(known_origins_regex, text_to_analyze)
        if match:
            return True, f"Clear evidence of metastasis from {match.group(2)}"
        
        # 检查是否直接提到其他部位的癌症
        for site in config.KNOWN_PRIMARY_SITES:
            # 例如: "breast cancer metastasis to lung", "metastatic colon cancer"
            if re.search(rf'{site}\s+(cancer|carcinoma|tumor|adenocarcinoma).*metasta', text_to_analyze):
                return True, f"Detected {site} cancer with metastasis"
            if re.search(rf'metasta.*{site}\s+(cancer|carcinoma|tumor)', text_to_analyze):
                return True, f"Detected metastatic {site} cancer"
        
        # 5b. 如果没有明确证据，尝试排除常见的原发性肺癌术语
        
        # 如果明确说明原发灶在肺部，则排除
        if re.search(r'(primary site:\s*lung|origin:\s*lung|primary lung tumor|primary lung cancer)', text_to_analyze):
            return False, "Primary lung cancer (explicit)"
        
        # 如果检测到肺癌关键词，且没有其他原发灶证据，我们保守地排除
        for cancer in config.PRIMARY_LUNG_CANCER_TERMS:
            if cancer in text_to_analyze:
                # 但是如果同时有"metastasis from"或明确的原发部位，则可能是描述
                if not re.search(rf'metastasis from|primary:\s*(?!lung)', text_to_analyze):
                    return False, f"Likely primary lung cancer ({cancer})"
        
        # 6. 检查转移方向 (排除从肺部转移到其他部位的样本)
        # 例如描述为"转移到脑部的肺癌样本"
        if re.search(r'lung.*metastasis to (brain|liver|bone|adrenal)', text_to_analyze):
            return False, "Lung cancer metastasized to other organs"
        
        # 如果通过了所有检查，但没有找到明确的原发部位证据
        # 我们标记为"可能符合，需人工复核"
        return True, "Possible lung metastasis (manual review required)"
    
    def analyze_gse(self, gse_id: str) -> List[Dict]:
        """
        下载并解析GSE元数据，应用过滤逻辑。
        
        Args:
            gse_id: GSE编号
            
        Returns:
            List[Dict]: 符合条件的样本信息列表
        """
        self.logger.info(f"Parsing {gse_id}...")
        self.stats['gse_analyzed'] += 1
        
        try:
            # 使用 GEOparse 下载和解析 SOFT 文件
            gse = GEOparse.get_GEO(
                geo=gse_id, 
                destdir=config.GEO_CACHE_DIR, 
                silent=True
            )
        except Exception as e:
            self.logger.error(f"Error processing {gse_id}: {e}")
            return []
        
        relevant_samples = []
        
        for gsm_name, gsm in gse.gsms.items():
            self.stats['gsm_analyzed'] += 1
            
            is_relevant, reason = self.is_lung_metastasis_of_other_origin(gsm.metadata)
            
            if is_relevant:
                # 查找 SRA 链接 (SRX ID)
                srx_id = None
                # SRA链接通常在 'relation' 字段中
                for relation in gsm.metadata.get('relation', []):
                    if 'SRA:' in relation or 'sra.cgi' in relation:
                        # 提取 SRX ID (Experiment ID)，覆盖不同URL格式
                        match = re.search(r'term=(SRX\d+)', relation)
                        if match:
                            srx_id = match.group(1)
                            break
                        match = re.search(r'/(SRX\d+)/?$', relation)
                        if match:
                            srx_id = match.group(1)
                            break
                        # 有时候直接是SRX ID
                        match = re.search(r'(SRX\d+)', relation)
                        if match:
                            srx_id = match.group(1)
                            break
                
                # 只有找到SRX ID才能下载原始数据
                if srx_id:
                    self.stats['gsm_relevant'] += 1
                    self.stats['srx_found'] += 1
                    relevant_samples.append({
                        "GSE": gse_id,
                        "GSM": gsm_name,
                        "SRX": srx_id,
                        "Library_Strategy": gsm.metadata.get('library_strategy', [''])[0],
                        "Title": gsm.metadata.get('title', [''])[0],
                        "Characteristics": "; ".join(gsm.metadata.get('characteristics_ch1', [])),
                        "Filter_Reason": reason,
                        "Source_Name": gsm.metadata.get('source_name_ch1', [''])[0]
                    })
                else:
                    self.logger.warning(f"{gsm_name} matches criteria but has no SRX link")
        
        return relevant_samples
    
    # ------------------------------------------------------------------
    # 阶段 3: 关联原始数据 (SRA) - SRX 到 SRR
    # ------------------------------------------------------------------
    
    def get_srr_from_srx_list(self, srx_list: List[str]) -> Dict[str, List[str]]:
        """
        高效地将 SRX ID 列表转换为 SRR ID 列表。
        使用 EFetch 的 'runinfo' (CSV格式)，这是最快且最稳定的方法。
        
        Args:
            srx_list: SRX编号列表
            
        Returns:
            Dict[str, List[str]]: SRX到SRR列表的映射
        """
        if not srx_list:
            return {}
        
        self.logger.info(f"Fetching SRR IDs for {len(srx_list)} SRX experiments...")
        srx_to_srr = {}
        
        # 1. 使用 ESearch 找到 SRA 内部 ID (因为 EFetch 需要内部ID)
        # 我们需要分批处理以避免过长的 URL
        batch_size = config.SRA_SEARCH_BATCH_SIZE
        sra_ids = []
        
        for i in range(0, len(srx_list), batch_size):
            batch = srx_list[i:i+batch_size]
            time.sleep(config.API_DELAY)
            try:
                # 乘以系数以防一个SRX有多个SRR
                handle = Entrez.esearch(
                    db="sra", 
                    term=" OR ".join(batch), 
                    retmax=len(batch) * 20
                )
                record = Entrez.read(handle)
                sra_ids.extend(record["IdList"])
            except Exception as e:
                self.logger.error(f"Error during SRA ESearch for batch {i}: {e}")
                time.sleep(5)
        
        if not sra_ids:
            self.logger.warning("No SRA IDs found.")
            return {}
        
        self.logger.info(f"Found {len(sra_ids)} SRA internal IDs, fetching RunInfo...")
        
        # 2. 使用 EFetch 获取 RunInfo (CSV 格式)
        time.sleep(config.API_DELAY)
        try:
            handle = Entrez.efetch(
                db="sra", 
                id=",".join(sra_ids), 
                rettype="runinfo", 
                retmode="text"
            )
            data = handle.read()
            handle.close()
        except Exception as e:
            self.logger.error(f"Error during SRA EFetch: {e}")
            return {}
        
        # 使用 pandas 读取 CSV 格式的 runinfo (比手动解析更稳定)
        try:
            # 确保数据是字符串类型
            if isinstance(data, bytes):
                data = data.decode('utf-8')
            runinfo_df = pd.read_csv(io.StringIO(data))
        except pd.errors.EmptyDataError:
            self.logger.warning("Received empty RunInfo data from SRA.")
            return {}
        except Exception as e:
            self.logger.error(f"Error parsing RunInfo CSV: {e}")
            return {}
        
        # 创建 SRX 到 SRR 的映射
        for index, row in runinfo_df.iterrows():
            srx = row.get('Experiment')
            srr = row.get('Run')
            if pd.notna(srx) and pd.notna(srr):
                if srx not in srx_to_srr:
                    srx_to_srr[srx] = []
                srx_to_srr[srx].append(srr)
                self.stats['srr_found'] += 1
        
        self.logger.info(f"Mapped {len(srx_to_srr)} SRX to {self.stats['srr_found']} SRR runs")
        return srx_to_srr
    
    # ------------------------------------------------------------------
    # 主执行流程
    # ------------------------------------------------------------------
    
    def run(self) -> Tuple[pd.DataFrame, List[str]]:
        """
        执行完整的挖掘流程
        
        Returns:
            Tuple[pd.DataFrame, List[str]]: (结果DataFrame, SRR列表)
        """
        self.logger.info("="*70)
        self.logger.info("Starting GEO Lung Metastasis Mining Pipeline")
        self.logger.info("="*70)
        
        # 1. 执行搜索
        self.logger.info("\n[Stage 1] Searching GEO database...")
        gse_list = self.search_geo()
        
        if not gse_list:
            self.logger.info("Initial search returned no results.")
            return pd.DataFrame(), []
        
        # 2. 解析和过滤
        self.logger.info(f"\n[Stage 2] Deep metadata analysis of {len(gse_list)} GSE datasets...")
        all_relevant_samples = []
        
        # 遍历所有找到的 GSE
        for idx, gse_id in enumerate(gse_list, 1):
            self.logger.info(f"Progress: {idx}/{len(gse_list)} - {gse_id}")
            samples = self.analyze_gse(gse_id)
            if samples:
                self.logger.info(f"-> Found {len(samples)} relevant samples in {gse_id}.")
                all_relevant_samples.extend(samples)
            time.sleep(config.API_DELAY)  # 尊重API限制
        
        if not all_relevant_samples:
            self.logger.info("\nNo relevant samples found matching the strict criteria after metadata analysis.")
            return pd.DataFrame(), []
        
        results_df = pd.DataFrame(all_relevant_samples)
        self.logger.info(f"\n[Stage 2 Complete] Total relevant samples found: {len(results_df)}")
        
        # 3. 关联 SRA SRR ID
        self.logger.info(f"\n[Stage 3] Linking to SRA database...")
        srx_list = results_df['SRX'].unique().tolist()
        srx_to_srr_map = self.get_srr_from_srx_list(srx_list)
        
        # 将 SRR 映射回主 DataFrame
        results_df['SRR_List'] = results_df['SRX'].map(lambda x: srx_to_srr_map.get(x, []))
        results_df['SRR_Count'] = results_df['SRR_List'].apply(len)
        
        # 提取所有 SRR ID 用于下载
        all_srr_ids = [
            srr for sublist in srx_to_srr_map.values() 
            if sublist for srr in sublist
        ]
        
        # 打印统计信息
        self.print_statistics()
        
        return results_df, all_srr_ids
    
    def save_results(self, results_df: pd.DataFrame, all_srr_ids: List[str]):
        """
        保存结果到文件
        
        Args:
            results_df: 结果DataFrame
            all_srr_ids: 所有SRR ID列表
        """
        # 4. 保存结果
        self.logger.info("\n[Stage 4] Saving results...")
        
        # 保存详细结果CSV
        results_df.to_csv(config.OUTPUT_CSV, index=False)
        self.logger.info(f"✓ Results saved to {config.OUTPUT_CSV}")
        self.logger.info("  Please review this file manually to verify sample accuracy.")
        
        # 保存SRR列表
        with open(config.SRR_LIST_FILE, 'w') as f:
            for srr in all_srr_ids:
                f.write(f"{srr}\n")
        self.logger.info(f"✓ SRR accession list ({len(all_srr_ids)} runs) saved to {config.SRR_LIST_FILE}")
        
        self.logger.info("\n" + "="*70)
        self.logger.info("Pipeline completed successfully!")
        self.logger.info("="*70)
        self.logger.info("\n⚠️  IMPORTANT: Manual Review Required")
        self.logger.info("Due to the complexity of biomedical metadata, please manually")
        self.logger.info("review the results CSV before downloading large amounts of data.")
        self.logger.info("Verify that samples truly match your research criteria.")
    
    def print_statistics(self):
        """打印统计信息"""
        self.logger.info("\n" + "="*70)
        self.logger.info("Mining Statistics:")
        self.logger.info("="*70)
        self.logger.info(f"GSE datasets found in search:    {self.stats['gse_found']}")
        self.logger.info(f"GSE datasets analyzed:           {self.stats['gse_analyzed']}")
        self.logger.info(f"GSM samples analyzed:            {self.stats['gsm_analyzed']}")
        self.logger.info(f"GSM samples matching criteria:   {self.stats['gsm_relevant']}")
        self.logger.info(f"SRX experiments found:           {self.stats['srx_found']}")
        self.logger.info(f"SRR runs available for download: {self.stats['srr_found']}")
        self.logger.info("="*70)


def main():
    """主函数"""
    # 检查配置
    if config.ENTREZ_EMAIL == "your.email@example.com":
        print("="*70)
        print("⚠️  WARNING: Please configure your email in config.py")
        print("="*70)
        print("\nNCBI requires a valid email address for API usage tracking.")
        print("Please edit config.py and set ENTREZ_EMAIL to your email address.")
        print("\nExample: ENTREZ_EMAIL = 'your.name@institution.edu'")
        print("="*70)
        return
    
    # 创建挖掘器实例并运行
    miner = GEOLungMetastasisMiner()
    
    try:
        results_df, all_srr_ids = miner.run()
        
        if not results_df.empty:
            miner.save_results(results_df, all_srr_ids)
        else:
            miner.logger.info("\nNo results to save.")
    
    except KeyboardInterrupt:
        miner.logger.info("\n\nPipeline interrupted by user.")
    except Exception as e:
        miner.logger.error(f"\nUnexpected error: {e}", exc_info=True)


if __name__ == "__main__":
    main()

