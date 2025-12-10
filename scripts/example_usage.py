#!/usr/bin/env python3
"""
示例使用脚本 - GEO 肺部转移瘤数据挖掘流水线
Example Usage Script - GEO Lung Metastasis Mining Pipeline

这个脚本展示了如何使用挖掘器的高级功能。
This script demonstrates advanced usage of the miner.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from geo_lung_metastasis_miner import GEOLungMetastasisMiner
import config

def example_basic_usage():
    """
    示例 1: 基本使用
    Example 1: Basic Usage
    """
    print("\n=== 示例 1: 基本使用 | Example 1: Basic Usage ===\n")
    
    # 创建挖掘器实例
    miner = GEOLungMetastasisMiner()
    
    # 运行完整流程
    results_df, all_srr_ids = miner.run()
    
    # 保存结果
    if not results_df.empty:
        miner.save_results(results_df, all_srr_ids)
        
        # 显示一些统计信息
        print(f"\n找到 {len(results_df)} 个符合条件的样本")
        print(f"涉及 {results_df['GSE'].nunique()} 个不同的 GSE 数据集")
        print(f"共有 {len(all_srr_ids)} 个 SRR 运行可供下载")


def example_custom_query():
    """
    示例 2: 使用自定义查询
    Example 2: Using Custom Query
    """
    print("\n=== 示例 2: 自定义查询 | Example 2: Custom Query ===\n")
    
    # 创建挖掘器实例
    miner = GEOLungMetastasisMiner()
    
    # 自定义搜索查询 - 例如只搜索单细胞数据
    custom_query = '''
        ("scRNA-seq" OR "single cell RNA-seq") AND 
        ("lung" AND "metastasis") AND 
        "Homo sapiens"[Organism] AND 
        "Expression profiling by high throughput sequencing"[DataSet Type]
    '''
    
    # 使用自定义查询搜索
    gse_list = miner.search_geo(query=custom_query)
    
    print(f"自定义查询找到 {len(gse_list)} 个 GSE 数据集")
    
    # 继续分析这些数据集
    all_relevant_samples = []
    for gse_id in gse_list[:5]:  # 只处理前5个作为示例
        samples = miner.analyze_gse(gse_id)
        all_relevant_samples.extend(samples)
    
    print(f"前5个数据集中找到 {len(all_relevant_samples)} 个相关样本")


def example_analyze_specific_gse():
    """
    示例 3: 分析特定的 GSE 数据集
    Example 3: Analyze Specific GSE Dataset
    """
    print("\n=== 示例 3: 分析特定 GSE | Example 3: Analyze Specific GSE ===\n")
    
    # 创建挖掘器实例
    miner = GEOLungMetastasisMiner()
    
    # 指定要分析的 GSE ID
    target_gse = "GSE123456"  # 替换为实际的 GSE ID
    
    print(f"分析 {target_gse}...")
    
    # 分析单个 GSE
    samples = miner.analyze_gse(target_gse)
    
    if samples:
        print(f"\n找到 {len(samples)} 个相关样本:")
        for sample in samples:
            print(f"  - {sample['GSM']}: {sample['Title']}")
            print(f"    过滤理由: {sample['Filter_Reason']}")
    else:
        print(f"在 {target_gse} 中未找到符合条件的样本")


def example_filter_results():
    """
    示例 4: 过滤和筛选结果
    Example 4: Filter and Screen Results
    """
    print("\n=== 示例 4: 过滤结果 | Example 4: Filter Results ===\n")
    
    import pandas as pd
    
    # 假设已经运行过挖掘器并生成了结果文件
    try:
        results_df = pd.read_csv(config.OUTPUT_CSV)
        
        print(f"总样本数: {len(results_df)}")
        
        # 只选择 scRNA-seq 数据
        scrna_df = results_df[results_df['Library_Strategy'].str.contains('RNA', case=False, na=False)]
        print(f"scRNA-seq 样本数: {len(scrna_df)}")
        
        # 只选择有明确原发部位证据的样本
        high_confidence = results_df[
            results_df['Filter_Reason'].str.contains('Clear evidence|Detected', case=False, na=False)
        ]
        print(f"高置信度样本数: {len(high_confidence)}")
        
        # 按 GSE 分组统计
        print("\n每个 GSE 的样本数:")
        print(results_df['GSE'].value_counts().head(10))
        
        # 保存筛选后的结果
        high_confidence.to_csv("high_confidence_samples.csv", index=False)
        print("\n高置信度样本已保存到 high_confidence_samples.csv")
        
    except FileNotFoundError:
        print(f"未找到结果文件: {config.OUTPUT_CSV}")
        print("请先运行主程序生成结果")


def example_get_srr_info():
    """
    示例 5: 获取 SRR 信息
    Example 5: Get SRR Information
    """
    print("\n=== 示例 5: 获取 SRR 信息 | Example 5: Get SRR Info ===\n")
    
    # 创建挖掘器实例
    miner = GEOLungMetastasisMiner()
    
    # 指定一些 SRX ID
    srx_list = ["SRX1234567", "SRX7654321"]  # 替换为实际的 SRX ID
    
    # 获取对应的 SRR ID
    srx_to_srr = miner.get_srr_from_srx_list(srx_list)
    
    print("SRX 到 SRR 的映射:")
    for srx, srr_list in srx_to_srr.items():
        print(f"  {srx} -> {', '.join(srr_list)}")


def main():
    """
    主函数 - 运行所有示例
    Main Function - Run All Examples
    """
    print("\n" + "="*70)
    print("GEO 肺部转移瘤数据挖掘流水线 - 使用示例")
    print("GEO Lung Metastasis Mining Pipeline - Usage Examples")
    print("="*70)
    
    # 检查配置
    if config.ENTREZ_EMAIL == "your.email@example.com":
        print("\n⚠️  请先在 config.py 中配置您的邮箱地址")
        print("⚠️  Please configure your email in config.py first")
        return
    
    # 运行示例（取消注释你想运行的示例）
    
    # 示例 1: 基本使用
    # example_basic_usage()
    
    # 示例 2: 自定义查询
    # example_custom_query()
    
    # 示例 3: 分析特定 GSE
    # example_analyze_specific_gse()
    
    # 示例 4: 过滤结果
    # example_filter_results()
    
    # 示例 5: 获取 SRR 信息
    # example_get_srr_info()
    
    print("\n提示: 取消注释上面的示例函数来运行它们")
    print("Tip: Uncomment the example functions above to run them")
    print("\n或者直接运行主程序:")
    print("Or run the main program directly:")
    print("  python geo_lung_metastasis_miner.py")
    print()


if __name__ == "__main__":
    main()

