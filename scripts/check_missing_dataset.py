#!/usr/bin/env python3
"""
检查缺失的数据集 - 分析为什么 GSE152048 和 GSM7453693 没有被爬取
Check Missing Dataset - Analyze why GSE152048 and GSM7453693 were not retrieved
"""

from Bio import Entrez
import GEOparse
import time
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import config

Entrez.email = config.ENTREZ_EMAIL

def check_gse_info(gse_id):
    """检查 GSE 数据集的详细信息"""
    print(f"\n{'='*70}")
    print(f"检查 GSE 数据集: {gse_id}")
    print('='*70)
    
    # 1. 检查是否存在
    print("\n[1] 检查数据集是否存在...")
    try:
        handle = Entrez.esearch(db='gds', term=gse_id)
        result = Entrez.read(handle)
        handle.close()
        
        if result['Count'] == '0':
            print(f"  ❌ 在 GEO 数据库中找不到 {gse_id}")
            return False
        else:
            print(f"  ✓ 找到 {gse_id}，共 {result['Count']} 条记录")
            print(f"  GDS ID: {result['IdList']}")
    except Exception as e:
        print(f"  ❌ 查询出错: {e}")
        return False
    
    # 2. 获取详细信息
    print("\n[2] 获取数据集详细信息...")
    try:
        time.sleep(0.5)
        if result['IdList']:
            handle = Entrez.esummary(db='gds', id=result['IdList'][0])
            summary = Entrez.read(handle)
            handle.close()
            
            # 处理不同的返回格式
            if isinstance(summary, list) and len(summary) > 0:
                info = summary[0]
            else:
                info = summary
            
            print(f"  Accession: {info.get('Accession', 'N/A')}")
            print(f"  Title: {info.get('title', info.get('Title', 'N/A'))[:100]}")
            print(f"  Entry Type: {info.get('EntryType', 'N/A')}")
            print(f"  Platform: {info.get('GPL', 'N/A')}")
            print(f"  Organism: {info.get('taxon', info.get('Organism', 'N/A'))}")
            print(f"  Sample Count: {info.get('n_samples', 'N/A')}")
            
            # 检查测序类型
            gse_type = info.get('gdsType', info.get('EntryType', ''))
            print(f"  GDS Type: {gse_type}")
            
    except Exception as e:
        print(f"  ⚠ 获取详细信息时出错: {e}")
    
    # 3. 使用 GEOparse 获取完整元数据
    print("\n[3] 下载并解析 SOFT 文件...")
    try:
        from pathlib import Path
        cache_dir = Path(__file__).parent.parent / 'GEO_Cache'
        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(cache_dir), silent=True)
        print(f"  ✓ 成功解析 {gse_id}")
        print(f"  样本数量: {len(gse.gsms)}")
        
        # 显示前几个样本
        print("\n  前5个样本:")
        for i, (gsm_name, gsm) in enumerate(list(gse.gsms.items())[:5], 1):
            title = gsm.metadata.get('title', [''])[0]
            print(f"    {i}. {gsm_name}: {title[:70]}")
        
        # 检查是否包含目标 GSM
        if 'GSM7453693' in gse.gsms:
            print(f"\n  ✓ 找到目标样本 GSM7453693!")
            gsm = gse.gsms['GSM7453693']
            print(f"    Title: {gsm.metadata.get('title', [''])[0]}")
            print(f"    Organism: {gsm.metadata.get('organism_ch1', [''])[0]}")
            print(f"    Characteristics: {'; '.join(gsm.metadata.get('characteristics_ch1', []))[:100]}")
            print(f"    Library Strategy: {gsm.metadata.get('library_strategy', [''])[0]}")
        else:
            print(f"\n  ⚠ 未找到 GSM7453693")
            
        return True, gse
        
    except Exception as e:
        print(f"  ❌ 下载/解析出错: {e}")
        return False, None

def check_why_filtered_out(gse_id):
    """检查为什么数据集没有被搜索到"""
    print(f"\n{'='*70}")
    print(f"分析为什么 {gse_id} 没有被搜索到")
    print('='*70)
    
    # 测试不同的搜索查询
    test_queries = [
        # 原始查询
        ('完整搜索条件', config.SEARCH_QUERY),
        
        # 简化查询1 - 只要人类和肺和转移
        ('简化查询1: 人类+肺+转移', 
         f'{gse_id} AND "Homo sapiens"[Organism]'),
        
        # 简化查询2 - 只要 GSE 编号
        ('简化查询2: 仅GSE编号',
         gse_id),
         
        # 简化查询3 - 肺+转移+人类
        ('简化查询3: 肺+转移',
         '("lung" OR "pulmonary") AND ("metastasis" OR "metastatic") AND "Homo sapiens"[Organism]'),
         
        # 简化查询4 - 骨肉瘤
        ('简化查询4: 骨肉瘤',
         'osteosarcoma AND lung AND "Homo sapiens"[Organism]'),
    ]
    
    results = {}
    
    for name, query in test_queries:
        print(f"\n测试: {name}")
        print(f"查询: {query[:100]}...")
        
        try:
            time.sleep(0.5)
            handle = Entrez.esearch(db='gds', term=query, retmax=100)
            result = Entrez.read(handle)
            handle.close()
            
            count = int(result['Count'])
            found = False
            
            if count > 0:
                # 检查是否包含目标 GSE
                time.sleep(0.5)
                handle = Entrez.esummary(db='gds', id=','.join(result['IdList'][:min(100, len(result['IdList']))]))
                summaries = Entrez.read(handle)
                handle.close()
                
                # 处理返回格式
                docs = summaries if isinstance(summaries, list) else []
                
                for doc in docs:
                    if doc.get('Accession') == gse_id:
                        found = True
                        break
            
            results[name] = {'count': count, 'found': found}
            
            if found:
                print(f"  ✓ 找到 {gse_id}! (总计 {count} 个结果)")
            else:
                print(f"  ✗ 未找到 {gse_id} (总计 {count} 个结果)")
                
        except Exception as e:
            print(f"  ❌ 查询出错: {e}")
            results[name] = {'count': 0, 'found': False, 'error': str(e)}
    
    return results

def analyze_filtering(gse):
    """分析样本是否会被过滤器过滤"""
    print(f"\n{'='*70}")
    print("分析过滤器行为")
    print('='*70)
    
    from geo_lung_metastasis_miner import GEOLungMetastasisMiner
    
    miner = GEOLungMetastasisMiner()
    
    print("\n检查每个样本是否通过过滤器:")
    
    for gsm_name, gsm in list(gse.gsms.items())[:10]:  # 只检查前10个
        is_relevant, reason = miner.is_lung_metastasis_of_other_origin(gsm.metadata)
        
        status = "✓ 通过" if is_relevant else "✗ 未通过"
        print(f"\n  {gsm_name}: {status}")
        print(f"    Title: {gsm.metadata.get('title', [''])[0][:60]}")
        print(f"    Reason: {reason}")
        
        if gsm_name == 'GSM7453693':
            print(f"    >>> 这是目标样本!")
            print(f"    Characteristics: {'; '.join(gsm.metadata.get('characteristics_ch1', []))}")

def main():
    """主函数"""
    gse_id = 'GSE152048'
    
    print("="*70)
    print("分析缺失数据集: GSE152048 和 GSM7453693")
    print("="*70)
    
    # 1. 检查数据集信息
    exists, gse = check_gse_info(gse_id)
    
    if not exists:
        print("\n❌ 无法获取数据集信息，终止分析")
        return
    
    # 2. 检查为什么没被搜索到
    search_results = check_why_filtered_out(gse_id)
    
    # 3. 分析过滤器
    if gse:
        analyze_filtering(gse)
    
    # 4. 总结
    print("\n" + "="*70)
    print("分析总结")
    print("="*70)
    
    print("\n搜索测试结果:")
    for name, result in search_results.items():
        if result['found']:
            print(f"  ✓ {name}: 找到")
        else:
            print(f"  ✗ {name}: 未找到 (总计 {result['count']} 个结果)")
    
    print("\n可能的原因:")
    
    # 检查是否能通过任何查询找到
    found_in_any = any(r['found'] for r in search_results.values())
    
    if found_in_any:
        print("  原因1: ✗ 数据集存在且可以被某些查询找到")
        print("  原因2: ✓ 可能是搜索条件太严格，排除了这个数据集")
        print("  原因3: 可能需要调整 config.py 中的 SEARCH_QUERY")
    else:
        print("  原因1: ✓ 数据集可能不符合任何搜索条件")
        print("  原因2: ✓ 数据集的元数据中可能缺少关键词")
        print("  原因3: 数据集类型可能不是预期的类型")
    
    print("\n建议:")
    if not found_in_any:
        print("  1. 检查数据集的 'DataSet Type' 是否符合搜索条件")
        print("  2. 检查数据集标题和描述中是否包含搜索的关键词")
        print("  3. 考虑放宽搜索条件，或者直接指定 GSE 编号")
    else:
        print("  1. 修改 config.py 中的 SEARCH_QUERY，使用更宽松的条件")
        print("  2. 或者直接在搜索查询中添加 'OR GSE152048'")

if __name__ == "__main__":
    main()


