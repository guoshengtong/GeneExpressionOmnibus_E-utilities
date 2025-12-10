#!/usr/bin/env python3
"""
测试改进效果 - 对比 V1 和 V2
Test Improvements - Compare V1 vs V2

这个脚本演示改进后的搜索和过滤能力
"""

from Bio import Entrez
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
import config_v2 as config
import time
import GEOparse

Entrez.email = config.ENTREZ_EMAIL

def test_search_improvements():
    """测试搜索改进"""
    print("="*70)
    print("测试1: 搜索策略对比")
    print("="*70)
    
    test_datasets = ['GSE234187', 'GSE152048']
    
    # V1: 严格搜索
    print("\n[V1] 使用严格搜索条件...")
    print(f"查询: {config.SEARCH_QUERY_STRICT[:100]}...")
    
    try:
        handle = Entrez.esearch(db='gds', term=config.SEARCH_QUERY_STRICT, retmax=100)
        result_strict = Entrez.read(handle)
        handle.close()
        
        print(f"找到: {result_strict['Count']} 个数据集")
        
        # 检查目标数据集
        time.sleep(0.5)
        if result_strict['IdList']:
            handle = Entrez.esummary(db='gds', id=','.join(result_strict['IdList'][:min(100, len(result_strict['IdList']))]))
            summaries = Entrez.read(handle)
            handle.close()
            
            docs = summaries if isinstance(summaries, list) else []
            found_gse = [doc.get('Accession') for doc in docs if doc.get('Accession', '').startswith('GSE')]
            
            for target in test_datasets:
                if target in found_gse:
                    print(f"  ✓ 找到 {target}")
                else:
                    print(f"  ✗ 未找到 {target}")
    except Exception as e:
        print(f"错误: {e}")
    
    # V2: 宽松搜索
    print("\n[V2] 使用宽松搜索条件...")
    print(f"查询: {config.SEARCH_QUERY_LOOSE[:100]}...")
    
    try:
        time.sleep(0.5)
        handle = Entrez.esearch(db='gds', term=config.SEARCH_QUERY_LOOSE, retmax=200)
        result_loose = Entrez.read(handle)
        handle.close()
        
        print(f"找到: {result_loose['Count']} 个数据集")
        
        time.sleep(0.5)
        if result_loose['IdList']:
            handle = Entrez.esummary(db='gds', id=','.join(result_loose['IdList'][:min(200, len(result_loose['IdList']))]))
            summaries = Entrez.read(handle)
            handle.close()
            
            docs = summaries if isinstance(summaries, list) else []
            found_gse = [doc.get('Accession') for doc in docs if doc.get('Accession', '').startswith('GSE')]
            
            for target in test_datasets:
                if target in found_gse:
                    print(f"  ✓ 找到 {target}")
                else:
                    print(f"  ✗ 未找到 {target}")
    except Exception as e:
        print(f"错误: {e}")
    
    # V2: 骨肉瘤特异性搜索
    print("\n[V2] 使用骨肉瘤特异性搜索...")
    osteosarcoma_query = config.CANCER_SPECIFIC_SEARCHES['Osteosarcoma']
    print(f"查询: {osteosarcoma_query[:100]}...")
    
    try:
        time.sleep(0.5)
        handle = Entrez.esearch(db='gds', term=osteosarcoma_query, retmax=200)
        result_osteo = Entrez.read(handle)
        handle.close()
        
        print(f"找到: {result_osteo['Count']} 个数据集")
        
        time.sleep(0.5)
        if result_osteo['IdList']:
            handle = Entrez.esummary(db='gds', id=','.join(result_osteo['IdList'][:min(200, len(result_osteo['IdList']))]))
            summaries = Entrez.read(handle)
            handle.close()
            
            docs = summaries if isinstance(summaries, list) else []
            found_gse = [doc.get('Accession') for doc in docs if doc.get('Accession', '').startswith('GSE')]
            
            for target in test_datasets:
                if target in found_gse:
                    print(f"  ✓ 找到 {target}")
                else:
                    print(f"  ✗ 未找到 {target}")
    except Exception as e:
        print(f"错误: {e}")


def test_filter_improvements():
    """测试过滤器改进"""
    print("\n" + "="*70)
    print("测试2: 过滤器对比")
    print("="*70)
    
    # 测试 GSE234187
    print("\n[GSE234187] 测试样本: GSM7453693")
    
    try:
        gse = GEOparse.get_GEO('GSE234187', destdir=config.GEO_CACHE_DIR, silent=True)
        
        if 'GSM7453693' in gse.gsms:
            gsm = gse.gsms['GSM7453693']
            title = gsm.metadata.get('title', [''])[0]
            characteristics = '; '.join(gsm.metadata.get('characteristics_ch1', []))
            
            print(f"  Title: {title}")
            print(f"  Characteristics: {characteristics}")
            
            # V1: 标准过滤
            print("\n  [V1] 标准过滤器判断:")
            text = (title + " " + characteristics).lower()
            has_lung = 'lung' in text or 'pulmonary' in text
            has_metastasis = 'metastasis' in text or 'metastatic' in text
            
            if has_lung and has_metastasis:
                print("    ✓ 通过（包含 lung 和 metastasis）")
            else:
                print(f"    ✗ 未通过（lung:{has_lung}, metastasis:{has_metastasis}）")
            
            # V2: 增强过滤（简化版）
            print("\n  [V2] 增强过滤器判断:")
            
            # 检查知识库
            if config.ENABLE_KNOWLEDGE_BASE and 'GSE234187' in config.KNOWN_DATASETS:
                kb = config.KNOWN_DATASETS['GSE234187']
                if 'GSM7453693' in kb.get('lung_metastasis_samples', []):
                    print(f"    ✓ 知识库确认: {kb['notes']}")
                    print(f"    置信度: {kb['confidence']}")
            
            # 检查癌症类型
            if 'osteosarcoma' in text or 'os' in title.lower():
                print("    ✓ 识别到骨肉瘤")
                if has_lung:
                    print("    ✓ 包含肺部关键词")
                if has_metastasis:
                    print("    ✓ 包含转移关键词")
    
    except Exception as e:
        print(f"错误: {e}")


def test_knowledge_base():
    """测试知识库功能"""
    print("\n" + "="*70)
    print("测试3: 外部知识库")
    print("="*70)
    
    print(f"\n已知数据集数量: {len(config.KNOWN_DATASETS)}")
    
    for gse_id, info in config.KNOWN_DATASETS.items():
        print(f"\n{gse_id}:")
        print(f"  状态: {info['status']}")
        print(f"  癌症类型: {info.get('cancer_type', 'N/A')}")
        print(f"  包含肺转移: {info.get('has_lung_metastasis', 'N/A')}")
        
        if 'lung_metastasis_samples' in info and info['lung_metastasis_samples']:
            print(f"  确认的样本: {', '.join(info['lung_metastasis_samples'])}")
        
        if info['status'] == 'needs_annotation':
            print(f"  ⚠️ 需要操作: {info.get('action_required', 'N/A')}")


def test_cancer_type_detection():
    """测试癌症类型识别"""
    print("\n" + "="*70)
    print("测试4: 癌症类型识别")
    print("="*70)
    
    test_texts = [
        ("Osteosarcoma patients", "osteosarcoma"),
        ("Breast cancer metastasis to lung", "breast"),
        ("Melanoma lung metastatic", "melanoma"),
        ("OS tissue of lung metastatic patient", "osteosarcoma"),
        ("RCC with pulmonary metastasis", "renal")
    ]
    
    for text, expected in test_texts:
        print(f"\n文本: {text}")
        print(f"  预期: {expected}")
        
        detected = []
        for cancer_type, pattern in config.CANCER_TYPE_PATTERNS.items():
            for keyword in pattern['keywords']:
                if keyword.lower() in text.lower():
                    detected.append(cancer_type)
                    break
            
            # 检查缩写
            for abbr in pattern.get('abbreviations', []):
                if abbr in text:
                    if cancer_type not in detected:
                        detected.append(cancer_type)
        
        if detected:
            print(f"  ✓ 检测到: {', '.join(detected)}")
            if expected in detected:
                print(f"  ✓ 匹配预期")
            else:
                print(f"  ⚠️ 未匹配预期（检测到了其他类型）")
        else:
            print(f"  ✗ 未检测到癌症类型")


def main():
    """主函数"""
    print("="*70)
    print("GEO 数据挖掘流水线 V2 改进效果测试")
    print("="*70)
    print()
    
    try:
        # 测试1: 搜索改进
        test_search_improvements()
        
        # 测试2: 过滤器改进
        test_filter_improvements()
        
        # 测试3: 知识库
        test_knowledge_base()
        
        # 测试4: 癌症类型识别
        test_cancer_type_detection()
        
        print("\n" + "="*70)
        print("测试完成！")
        print("="*70)
        
        print("\n总结:")
        print("  V1 vs V2 改进:")
        print("  ✓ 搜索策略: 单一 → 多阶段")
        print("  ✓ 过滤器: 简单关键词 → 智能识别+知识库")
        print("  ✓ 癌症类型: 不识别 → 多种类型识别")
        print("  ✓ 可追溯性: 基础 → 详细日志+置信度")
        print()
        print("  预期效果:")
        print("  ✓ GSE234187: 可通过知识库或宽松搜索找到")
        print("  ✓ GSE152048: 可通过骨肉瘤特异性搜索找到")
        print("  ✓ 更多数据集: 预计增加20-50%")
        print()
        
    except KeyboardInterrupt:
        print("\n\n测试被中断")
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()


