#!/usr/bin/env python3
"""
复核辅助工具 - 帮助快速查看和分析结果
Review Helper Tool - Quick view and analysis of results
"""

import pandas as pd
import webbrowser
import sys

def load_results():
    """加载结果文件"""
    try:
        df = pd.read_csv('GEO_Lung_Metastasis_Mining_Results.csv')
        return df
    except FileNotFoundError:
        print("❌ 找不到结果文件: GEO_Lung_Metastasis_Mining_Results.csv")
        sys.exit(1)

def show_summary(df):
    """显示结果摘要"""
    print("="*70)
    print("数据挖掘结果摘要 | Mining Results Summary")
    print("="*70)
    print()
    
    print(f"总样本数: {len(df)}")
    print(f"不同的 GSE 数量: {df['GSE'].nunique()}")
    print(f"总 SRR 运行数: {df['SRR_Count'].sum()}")
    print()
    
    print("按 GSE 分组统计:")
    print("-"*70)
    gse_counts = df.groupby('GSE').size().sort_values(ascending=False)
    for gse, count in gse_counts.items():
        print(f"  {gse}: {count} 个样本")
    print()
    
    print("按过滤理由分组:")
    print("-"*70)
    reason_counts = df['Filter_Reason'].value_counts()
    for reason, count in reason_counts.items():
        print(f"  {reason}: {count} 个样本")
    print()

def show_detailed_samples(df, gse_id=None):
    """显示样本详细信息"""
    if gse_id:
        samples = df[df['GSE'] == gse_id]
        print(f"\n{'='*70}")
        print(f"GSE 数据集详情: {gse_id}")
        print(f"{'='*70}\n")
    else:
        samples = df
        print(f"\n{'='*70}")
        print(f"所有样本详情")
        print(f"{'='*70}\n")
    
    for idx, row in samples.iterrows():
        print(f"样本 #{idx + 1}")
        print(f"  GSE: {row['GSE']}")
        print(f"  GSM: {row['GSM']}")
        print(f"  标题: {row['Title'][:80]}...")
        print(f"  特征: {row['Characteristics'][:80]}...")
        print(f"  来源: {row['Source_Name']}")
        print(f"  判定: {row['Filter_Reason']}")
        print(f"  SRR数: {row['SRR_Count']}")
        print()

def analyze_problematic_samples(df):
    """分析可能有问题的样本"""
    print("\n" + "="*70)
    print("⚠️  需要特别注意的样本 | Samples Requiring Special Attention")
    print("="*70 + "\n")
    
    # 1. GSE266330 - 肺癌骨转移（应该排除）
    lung_to_bone = df[df['GSE'] == 'GSE266330']
    if not lung_to_bone.empty:
        print("❌ 建议排除 - 肺癌转移到骨（不符合要求）:")
        print("-"*70)
        print(f"  GSE266330: {len(lung_to_bone)} 个样本")
        print(f"  特征: {lung_to_bone.iloc[0]['Characteristics']}")
        print(f"  来源: {lung_to_bone.iloc[0]['Source_Name']}")
        print(f"  原因: 这些是肺癌转移到骨头，不是其他癌症转移到肺")
        print()
    
    # 2. GSE179373 - 脑转移（不是肺转移）
    brain_met = df[df['GSE'] == 'GSE179373']
    if not brain_met.empty:
        print("⚠️  需要确认 - 脑转移样本（不是肺转移）:")
        print("-"*70)
        print(f"  GSE179373: {len(brain_met)} 个样本")
        print(f"  特征: {brain_met.iloc[0]['Characteristics'][:80]}...")
        print(f"  原因: 这些是脑转移样本，包含乳腺癌和黑色素瘤")
        print(f"  建议: 如果您只需要肺转移样本，应该排除这些")
        print()
    
    # 3. 缺少明确原发部位的样本
    possible_only = df[df['Filter_Reason'] == 'Possible lung metastasis (manual review required)']
    print(f"⚠️  需要人工复核的样本: {len(possible_only)} 个")
    print("-"*70)
    gse_groups = possible_only.groupby('GSE').size()
    for gse, count in gse_groups.items():
        print(f"  {gse}: {count} 个样本")
    print()

def open_geo_pages(df, gse_list=None):
    """在浏览器中打开 GEO 页面"""
    if gse_list is None:
        gse_list = df['GSE'].unique()
    
    print(f"\n准备打开 {len(gse_list)} 个 GEO 页面...")
    print("这将在您的默认浏览器中打开多个标签页。")
    response = input("\n继续? (y/n): ")
    
    if response.lower() != 'y':
        print("已取消")
        return
    
    for gse in gse_list:
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
        print(f"打开: {url}")
        webbrowser.open(url)
        import time
        time.sleep(1)  # 避免同时打开太多

def export_filtered_results(df, output_file='GEO_Results_Filtered.csv'):
    """导出筛选后的结果"""
    print("\n" + "="*70)
    print("导出筛选后的结果")
    print("="*70 + "\n")
    
    # 排除 GSE266330（肺癌骨转移）
    df_filtered = df[df['GSE'] != 'GSE266330'].copy()
    print(f"✓ 已排除 GSE266330 (肺癌骨转移): {len(df) - len(df_filtered)} 个样本")
    
    # 可选：排除 GSE179373（脑转移）
    response = input("\n是否也排除 GSE179373 (脑转移)? (y/n): ")
    if response.lower() == 'y':
        brain_count = len(df_filtered[df_filtered['GSE'] == 'GSE179373'])
        df_filtered = df_filtered[df_filtered['GSE'] != 'GSE179373']
        print(f"✓ 已排除 GSE179373 (脑转移): {brain_count} 个样本")
    
    # 保存
    df_filtered.to_csv(output_file, index=False)
    print(f"\n✓ 筛选后的结果已保存到: {output_file}")
    print(f"  原始样本数: {len(df)}")
    print(f"  筛选后样本数: {len(df_filtered)}")
    print(f"  排除的样本数: {len(df) - len(df_filtered)}")
    
    # 生成新的 SRR 列表
    import ast
    all_srr = []
    for srr_list_str in df_filtered['SRR_List']:
        try:
            srr_list = ast.literal_eval(srr_list_str)
            all_srr.extend(srr_list)
        except:
            pass
    
    srr_file = 'SRR_accession_list_filtered.txt'
    with open(srr_file, 'w') as f:
        for srr in all_srr:
            f.write(f"{srr}\n")
    
    print(f"✓ 新的 SRR 列表已保存到: {srr_file}")
    print(f"  SRR 数量: {len(all_srr)}")

def interactive_menu():
    """交互式菜单"""
    df = load_results()
    
    while True:
        print("\n" + "="*70)
        print("GEO 数据挖掘结果复核工具 | Review Helper Tool")
        print("="*70)
        print("\n选择操作:")
        print("  1. 显示结果摘要")
        print("  2. 显示所有样本详情")
        print("  3. 按 GSE 查看样本")
        print("  4. 分析需要注意的样本")
        print("  5. 在浏览器中打开 GEO 页面")
        print("  6. 导出筛选后的结果")
        print("  0. 退出")
        print()
        
        choice = input("请选择 (0-6): ").strip()
        
        if choice == '0':
            print("\n再见!")
            break
        elif choice == '1':
            show_summary(df)
        elif choice == '2':
            show_detailed_samples(df)
        elif choice == '3':
            print("\n可用的 GSE:")
            for i, gse in enumerate(df['GSE'].unique(), 1):
                print(f"  {i}. {gse}")
            gse_input = input("\n输入 GSE 编号: ").strip()
            if gse_input in df['GSE'].values:
                show_detailed_samples(df, gse_input)
            else:
                print(f"❌ 找不到 GSE: {gse_input}")
        elif choice == '4':
            analyze_problematic_samples(df)
        elif choice == '5':
            print("\n打开哪些 GSE 的页面?")
            print("  1. 所有 GSE")
            print("  2. 需要复核的 GSE (前5个)")
            print("  3. 特定 GSE")
            sub_choice = input("选择 (1-3): ").strip()
            
            if sub_choice == '1':
                open_geo_pages(df)
            elif sub_choice == '2':
                priority = ['GSE270231', 'GSE109281', 'GSE193594', 'GSE156405', 'GSE179373']
                available = [gse for gse in priority if gse in df['GSE'].values]
                open_geo_pages(df, available)
            elif sub_choice == '3':
                gse = input("输入 GSE 编号: ").strip()
                if gse in df['GSE'].values:
                    open_geo_pages(df, [gse])
                else:
                    print(f"❌ 找不到 GSE: {gse}")
        elif choice == '6':
            export_filtered_results(df)
        else:
            print("❌ 无效的选择")
        
        input("\n按 Enter 继续...")

if __name__ == "__main__":
    try:
        interactive_menu()
    except KeyboardInterrupt:
        print("\n\n用户中断。再见!")
        sys.exit(0)


