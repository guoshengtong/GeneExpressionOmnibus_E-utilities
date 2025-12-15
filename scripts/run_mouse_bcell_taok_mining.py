#!/usr/bin/env python3
"""
运行小鼠骨髓单细胞B细胞发育和TAOK基因数据挖掘
Run Mouse Bone Marrow Single-Cell B Cell Development and TAOK Gene Mining

这个脚本执行针对小鼠骨髓单细胞B细胞发育和TAOK基因的数据挖掘任务
"""

import sys
from pathlib import Path

# 添加项目根目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent))

# 导入挖掘脚本
from scripts.geo_mouse_bcell_taok_miner import main

if __name__ == "__main__":
    print("\n" + "="*70)
    print("小鼠骨髓单细胞B细胞发育和TAOK基因数据挖掘")
    print("Mouse Bone Marrow Single-Cell B Cell Development and TAOK Gene Mining")
    print("="*70)
    print()
    
    try:
        results_df, stage_details = main()
        print("\n" + "="*70)
        print("✅ 挖掘任务完成！")
        print("="*70)
        print("\n📊 结果摘要:")
        if len(results_df) > 0:
            print(f"  - 找到 {len(results_df)} 个相关样本")
            print(f"  - 涉及 {results_df['GSE'].nunique()} 个数据集")
            print(f"  - 包含 {results_df['SRR_Count'].sum()} 个SRR运行")
        else:
            print("  - 未找到符合条件的样本")
        print()
    except KeyboardInterrupt:
        print("\n\n⚠️  用户中断操作")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 发生错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

