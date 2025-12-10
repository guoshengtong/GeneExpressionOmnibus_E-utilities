#!/usr/bin/env python3
"""
改进的搜索方法测试
Test improved search methods
"""

from Bio import Entrez
import time
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import config

Entrez.email = config.ENTREZ_EMAIL

def test_search_methods():
    """测试不同的搜索方法"""
    
    # 简化的查询 - 测试用
    simple_query = '("lung" AND "metastasis" AND "single cell") AND "Homo sapiens"[Organism]'
    
    print("="*70)
    print("测试不同的 NCBI 数据库搜索方法")
    print("="*70)
    print()
    
    # 方法 1: 搜索 gds 数据库
    print("方法 1: 搜索 'gds' 数据库")
    print("-"*70)
    try:
        handle = Entrez.esearch(db='gds', term=simple_query, retmax=10)
        results = Entrez.read(handle)
        handle.close()
        print(f"找到: {results['Count']} 条记录")
        print(f"ID 列表: {results['IdList'][:5]}")
        
        # 尝试获取详细信息
        if results['IdList']:
            time.sleep(0.5)
            handle = Entrez.esummary(db='gds', id=results['IdList'][0])
            summary = Entrez.read(handle)
            handle.close()
            print(f"\n第一条记录的详细信息:")
            if 'DocumentSummarySet' in summary:
                for doc in summary['DocumentSummarySet']['DocumentSummary']:
                    print(f"  Accession: {doc.get('Accession', 'N/A')}")
                    print(f"  Title: {doc.get('title', 'N/A')[:80]}")
                    print(f"  EntryType: {doc.get('EntryType', 'N/A')}")
            else:
                print("  格式:", type(summary))
                if isinstance(summary, list) and len(summary) > 0:
                    print(f"  Accession: {summary[0].get('Accession', 'N/A')}")
                    print(f"  Title: {summary[0].get('title', 'N/A')[:80]}")
    except Exception as e:
        print(f"错误: {e}")
    
    print("\n" + "="*70)
    
    # 方法 2: 直接搜索 GSE
    print("\n方法 2: 直接在标题中搜索 GSE 关键词")
    print("-"*70)
    try:
        # 添加 GSE 作为过滤条件
        gse_query = simple_query + ' AND GSE[Accession]'
        handle = Entrez.esearch(db='gds', term=gse_query, retmax=10)
        results = Entrez.read(handle)
        handle.close()
        print(f"找到: {results['Count']} 条记录")
        print(f"ID 列表: {results['IdList'][:5]}")
        
        if results['IdList']:
            time.sleep(0.5)
            handle = Entrez.esummary(db='gds', id=','.join(results['IdList'][:3]))
            summary = Entrez.read(handle)
            handle.close()
            
            print(f"\n前3条记录:")
            if 'DocumentSummarySet' in summary:
                for doc in summary['DocumentSummarySet']['DocumentSummary']:
                    acc = doc.get('Accession', 'N/A')
                    title = doc.get('title', 'N/A')[:80] if 'title' in doc else doc.get('Title', 'N/A')[:80]
                    print(f"  {acc}: {title}")
            elif isinstance(summary, list):
                for doc in summary[:3]:
                    acc = doc.get('Accession', 'N/A')
                    title = doc.get('title', 'N/A')[:80] if 'title' in doc else doc.get('Title', 'N/A')[:80]
                    print(f"  {acc}: {title}")
    except Exception as e:
        print(f"错误: {e}")
    
    print("\n" + "="*70)
    
    # 方法 3: 使用 ELink 从 GDS 到 GSE
    print("\n方法 3: 尝试使用更宽松的查询")
    print("-"*70)
    try:
        # 更宽松的查询
        loose_query = '(lung metastasis) AND (single cell OR scRNA-seq) AND "Homo sapiens"[Organism]'
        handle = Entrez.esearch(db='gds', term=loose_query, retmax=20)
        results = Entrez.read(handle)
        handle.close()
        print(f"找到: {results['Count']} 条记录")
        
        if results['IdList']:
            time.sleep(0.5)
            # 获取前5条的详细信息
            handle = Entrez.esummary(db='gds', id=','.join(results['IdList'][:5]))
            summary = Entrez.read(handle)
            handle.close()
            
            gse_list = []
            print(f"\n前5条记录的 Accession:")
            
            # 处理不同的返回格式
            docs = []
            if 'DocumentSummarySet' in summary and 'DocumentSummary' in summary['DocumentSummarySet']:
                docs = summary['DocumentSummarySet']['DocumentSummary']
            elif isinstance(summary, list):
                docs = summary
            
            for doc in docs:
                acc = doc.get('Accession', 'N/A')
                entry_type = doc.get('EntryType', 'N/A')
                print(f"  {acc} (类型: {entry_type})")
                if acc.startswith('GSE'):
                    gse_list.append(acc)
            
            print(f"\n提取到的 GSE: {gse_list}")
            
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_search_methods()


