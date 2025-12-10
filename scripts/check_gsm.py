#!/usr/bin/env python3
"""
检查 GSM 样本信息
"""

from Bio import Entrez
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import config

Entrez.email = config.ENTREZ_EMAIL

def check_gsm(gsm_id):
    """检查 GSM 样本属于哪个 GSE"""
    print(f"\n检查样本: {gsm_id}")
    print("="*70)
    
    try:
        # 搜索 GSM
        handle = Entrez.esearch(db='gds', term=gsm_id)
        result = Entrez.read(handle)
        handle.close()
        
        if result['Count'] == '0':
            print(f"❌ 找不到 {gsm_id}")
            return
        
        print(f"✓ 找到 {result['Count']} 条记录")
        
        # 获取详细信息
        import time
        time.sleep(0.5)
        handle = Entrez.esummary(db='gds', id=result['IdList'][0])
        summary = Entrez.read(handle)
        handle.close()
        
        info = summary[0] if isinstance(summary, list) else summary
        
        print(f"\n样本信息:")
        print(f"  Accession: {info.get('Accession', 'N/A')}")
        print(f"  Title: {info.get('title', info.get('Title', 'N/A'))}")
        
        # 从 Accession 中提取 GSE
        acc = info.get('Accession', '')
        if acc.startswith('GSM'):
            # 尝试使用 ELink 找到关联的 GSE
            time.sleep(0.5)
            handle = Entrez.elink(dbfrom='gds', db='gds', id=result['IdList'][0])
            link_result = Entrez.read(handle)
            handle.close()
            
            print(f"\n关联的数据集:")
            if link_result and len(link_result) > 0:
                for link in link_result[0].get('LinkSetDb', []):
                    for linked_id in link.get('Link', [])[:5]:
                        time.sleep(0.3)
                        handle2 = Entrez.esummary(db='gds', id=linked_id['Id'])
                        linked_summary = Entrez.read(handle2)
                        handle2.close()
                        
                        linked_info = linked_summary[0] if isinstance(linked_summary, list) else linked_summary
                        linked_acc = linked_info.get('Accession', '')
                        if linked_acc.startswith('GSE'):
                            print(f"  → {linked_acc}: {linked_info.get('title', linked_info.get('Title', 'N/A'))[:80]}")
        
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    check_gsm('GSM7453693')


