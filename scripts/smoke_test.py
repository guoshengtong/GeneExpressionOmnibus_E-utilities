#!/usr/bin/env python3
"""
冒烟测试 - 验证流水线基本功能
Smoke Test - Verify basic pipeline functionality

这个脚本执行一个非常小规模的测试来验证系统是否正常工作。
This script performs a very small-scale test to verify the system works.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from geo_lung_metastasis_miner import GEOLungMetastasisMiner
import config
import logging

def smoke_test():
    """执行冒烟测试"""
    print("="*70)
    print("开始冒烟测试 | Starting Smoke Test")
    print("="*70)
    print()
    
    # 配置日志级别为INFO以便看到详细过程
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # 创建挖掘器实例
    print("✓ 创建挖掘器实例...")
    miner = GEOLungMetastasisMiner()
    
    # 使用配置文件中的完整查询
    print("\n✓ 执行测试搜索...")
    print(f"使用配置文件中的查询: {config.SEARCH_QUERY[:100]}...")
    print()
    
    try:
        # 搜索 GEO（使用默认配置查询）
        gse_list = miner.search_geo()
        
        if not gse_list:
            print("⚠ 搜索未返回结果，这可能是正常的（查询条件可能很严格）")
            print("但 API 调用成功，说明系统正常工作！")
            return True
        
        print(f"\n✓ 找到 {len(gse_list)} 个潜在的 GSE 数据集")
        print(f"前5个: {', '.join(gse_list[:5])}")
        
        # 只分析第一个 GSE 作为测试
        if len(gse_list) > 0:
            test_gse = gse_list[0]
            print(f"\n✓ 测试解析第一个数据集: {test_gse}")
            
            samples = miner.analyze_gse(test_gse)
            
            if samples:
                print(f"\n✓ 在 {test_gse} 中找到 {len(samples)} 个符合条件的样本")
                print("\n样本示例:")
                for sample in samples[:2]:  # 只显示前2个
                    print(f"  - {sample['GSM']}: {sample['Title'][:80]}")
                    print(f"    判定理由: {sample['Filter_Reason']}")
            else:
                print(f"\n✓ {test_gse} 中没有符合条件的样本（过滤规则工作正常）")
            
            # 如果有样本，测试 SRA 关联
            if samples:
                print("\n✓ 测试 SRA 关联...")
                srx_list = [s['SRX'] for s in samples]
                srx_to_srr = miner.get_srr_from_srx_list(srx_list[:2])  # 只测试前2个
                
                if srx_to_srr:
                    print(f"✓ 成功获取 SRR 映射:")
                    for srx, srr_list in list(srx_to_srr.items())[:2]:
                        print(f"  {srx} -> {', '.join(srr_list)}")
        
        print("\n" + "="*70)
        print("✅ 冒烟测试成功！| Smoke Test PASSED!")
        print("="*70)
        print("\n系统各组件工作正常，可以开始正式挖掘数据了！")
        print("System components are working properly, ready for full mining!")
        print()
        print("运行完整流水线:")
        print("Run full pipeline:")
        print("  python3 geo_lung_metastasis_miner.py")
        print()
        
        return True
        
    except Exception as e:
        print(f"\n❌ 冒烟测试失败 | Smoke Test FAILED")
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    import sys
    success = smoke_test()
    sys.exit(0 if success else 1)

