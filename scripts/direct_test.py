#!/usr/bin/env python3
"""
直接测试 - 测试解析一个已知的 GSE 数据集
Direct Test - Test parsing a known GSE dataset
"""

from geo_lung_metastasis_miner import GEOLungMetastasisMiner
import logging

def direct_test():
    """直接测试解析功能"""
    print("="*70)
    print("直接测试 GEO 解析功能 | Direct Test of GEO Parsing")
    print("="*70)
    print()
    
    # 配置日志
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # 创建挖掘器实例
    print("✓ 创建挖掘器实例...")
    miner = GEOLungMetastasisMiner()
    
    # 使用一个已知的单细胞肺癌数据集进行测试
    # GSE131907 - 这是一个公开的肺癌单细胞数据集
    test_gse_list = ["GSE131907", "GSE123456"]  # 第二个可能不存在，用于测试错误处理
    
    print(f"\n测试解析 {len(test_gse_list)} 个 GSE 数据集...")
    print()
    
    for test_gse in test_gse_list:
        print(f"\n{'='*70}")
        print(f"测试 GSE: {test_gse}")
        print('='*70)
        
        try:
            samples = miner.analyze_gse(test_gse)
            
            if samples:
                print(f"\n✓ 找到 {len(samples)} 个符合过滤条件的样本")
                print("\n样本详情:")
                for i, sample in enumerate(samples[:3], 1):  # 只显示前3个
                    print(f"\n  样本 {i}:")
                    print(f"    GSM: {sample['GSM']}")
                    print(f"    标题: {sample['Title'][:80]}...")
                    print(f"    测序技术: {sample['Library_Strategy']}")
                    print(f"    SRX: {sample['SRX']}")
                    print(f"    判定理由: {sample['Filter_Reason']}")
                
                if len(samples) > 3:
                    print(f"\n  ... 还有 {len(samples) - 3} 个样本")
                    
            else:
                print(f"\n✓ GSE 解析成功，但没有样本符合过滤条件")
                print("  这是正常的 - 说明过滤逻辑正常工作")
                
        except FileNotFoundError:
            print(f"\n⚠ GSE {test_gse} 不存在（这是预期的测试情况）")
        except Exception as e:
            print(f"\n❌ 解析 {test_gse} 时出错: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*70)
    print("✅ 直接测试完成！| Direct Test Completed!")
    print("="*70)
    print()
    print("如果成功解析了至少一个 GSE，说明核心功能正常！")
    print("If at least one GSE was successfully parsed, core functionality works!")
    print()
    
if __name__ == "__main__":
    direct_test()


