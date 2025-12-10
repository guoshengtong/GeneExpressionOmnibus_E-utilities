#!/usr/bin/env python3
"""
V2结果复核工具 - 交互式复核界面
Review V2 Results - Interactive review interface
"""

import pandas as pd
import webbrowser
import sys
import os

class V2ResultsReviewer:
    def __init__(self):
        self.df = self.load_results()
        self.current_filter = None
        self.filtered_df = self.df.copy()
        
    def load_results(self):
        """加载V2结果文件"""
        try:
            df = pd.read_csv('GEO_Lung_Metastasis_Mining_Results_V2.csv')
            print(f"✓ 成功加载 {len(df)} 个样本")
            return df
        except FileNotFoundError:
            print("❌ 找不到结果文件: GEO_Lung_Metastasis_Mining_Results_V2.csv")
            sys.exit(1)
    
    def show_summary(self):
        """显示结果摘要"""
        df = self.filtered_df
        print("\n" + "="*80)
        print("V2 数据挖掘结果摘要")
        print("="*80)
        print()
        
        # 基本统计
        print("【整体统计】")
        print(f"  总样本数: {len(df)}")
        print(f"  数据集数量: {df['GSE'].nunique()}")
        print(f"  V1已有样本: {df['Found_In_V1'].sum()}")
        print(f"  V2新增样本: {(~df['Found_In_V1']).sum()}")
        print()
        
        # 置信度分布
        print("【置信度分布】")
        high_conf = len(df[df['Confidence'] >= 0.8])
        med_conf = len(df[(df['Confidence'] >= 0.6) & (df['Confidence'] < 0.8)])
        low_conf = len(df[df['Confidence'] < 0.6])
        print(f"  高置信度 (≥0.8): {high_conf} 个 ({high_conf/len(df)*100:.1f}%)")
        print(f"  中等置信度 (0.6-0.8): {med_conf} 个 ({med_conf/len(df)*100:.1f}%)")
        print(f"  低置信度 (<0.6): {low_conf} 个 ({low_conf/len(df)*100:.1f}%)")
        print()
        
        # 按数据集统计
        print("【按数据集统计】(前10)")
        print("-"*80)
        gse_counts = df.groupby('GSE').agg({
            'GSM': 'count',
            'Confidence': 'max',
            'Found_In_V1': lambda x: 'V1' if x.any() else 'NEW'
        }).rename(columns={'GSM': 'Sample_Count', 'Confidence': 'Max_Conf', 'Found_In_V1': 'Status'})
        gse_counts = gse_counts.sort_values('Sample_Count', ascending=False)
        
        for idx, (gse, row) in enumerate(gse_counts.head(10).iterrows(), 1):
            status = row['Status']
            conf = row['Max_Conf']
            conf_label = "高" if conf >= 0.8 else "中" if conf >= 0.6 else "低"
            print(f"  {idx:2d}. {gse}: {int(row['Sample_Count']):2d}个样本 | 置信度:{conf:.1f}({conf_label}) | {status}")
        
        if len(gse_counts) > 10:
            print(f"  ... 还有 {len(gse_counts) - 10} 个数据集")
        print()
        
        # 判定理由统计
        print("【判定理由统计】(前8)")
        print("-"*80)
        reason_counts = df['Filter_Reason'].value_counts()
        for idx, (reason, count) in enumerate(reason_counts.head(8).items(), 1):
            pct = count/len(df)*100
            print(f"  {idx}. {reason}: {count} ({pct:.1f}%)")
        if len(reason_counts) > 8:
            print(f"  ... 还有 {len(reason_counts) - 8} 种判定理由")
        print()
    
    def filter_by_confidence(self, min_conf=0.0, max_conf=1.0):
        """按置信度过滤"""
        self.filtered_df = self.df[(self.df['Confidence'] >= min_conf) & 
                                    (self.df['Confidence'] <= max_conf)]
        self.current_filter = f"Confidence: {min_conf}-{max_conf}"
        print(f"✓ 已过滤: {len(self.filtered_df)} 个样本 (置信度 {min_conf}-{max_conf})")
    
    def filter_by_status(self, new_only=True):
        """按新增/已有过滤"""
        if new_only:
            self.filtered_df = self.df[~self.df['Found_In_V1']]
            self.current_filter = "仅新增样本"
            print(f"✓ 已过滤: {len(self.filtered_df)} 个新增样本")
        else:
            self.filtered_df = self.df[self.df['Found_In_V1']]
            self.current_filter = "仅V1已有样本"
            print(f"✓ 已过滤: {len(self.filtered_df)} 个V1已有样本")
    
    def filter_by_gse(self, gse_id):
        """按数据集过滤"""
        self.filtered_df = self.df[self.df['GSE'] == gse_id]
        self.current_filter = f"GSE: {gse_id}"
        print(f"✓ 已过滤: {len(self.filtered_df)} 个样本来自 {gse_id}")
    
    def reset_filter(self):
        """重置过滤"""
        self.filtered_df = self.df.copy()
        self.current_filter = None
        print(f"✓ 已重置: 显示所有 {len(self.filtered_df)} 个样本")
    
    def show_samples(self, limit=10, sort_by='Confidence'):
        """显示样本列表"""
        df = self.filtered_df.copy()
        
        # 排序
        if sort_by == 'Confidence':
            df = df.sort_values('Confidence', ascending=False)
        elif sort_by == 'GSE':
            df = df.sort_values('GSE')
        
        print("\n" + "="*80)
        if self.current_filter:
            print(f"样本列表 (过滤: {self.current_filter}, 排序: {sort_by})")
        else:
            print(f"样本列表 (排序: {sort_by})")
        print("="*80)
        print()
        
        for idx, (_, row) in enumerate(df.head(limit).iterrows(), 1):
            status = "V1" if row['Found_In_V1'] else "NEW"
            conf = row['Confidence']
            conf_emoji = "⭐" if conf >= 0.9 else "✓" if conf >= 0.7 else "?"
            
            print(f"[{idx:2d}] {conf_emoji} {row['GSE']} - {row['GSM']} | 置信度:{conf:.1f} | {status}")
            print(f"     标题: {row['Title'][:60]}...")
            print(f"     判定: {row['Filter_Reason']}")
            print(f"     类型: {row['Library_Strategy']}")
            print()
        
        if len(df) > limit:
            print(f"... 还有 {len(df) - limit} 个样本 (使用 show_all 查看全部)")
        print()
    
    def show_detail(self, index):
        """显示样本详情"""
        if index < 1 or index > len(self.filtered_df):
            print(f"❌ 无效的索引: {index} (范围: 1-{len(self.filtered_df)})")
            return
        
        row = self.filtered_df.iloc[index - 1]
        
        print("\n" + "="*80)
        print(f"样本详情 #{index}")
        print("="*80)
        print()
        print(f"GSE编号:      {row['GSE']}")
        print(f"GSM编号:      {row['GSM']}")
        print(f"SRX编号:      {row['SRX']}")
        print(f"测序策略:     {row['Library_Strategy']}")
        print()
        print(f"样本标题:     {row['Title']}")
        print(f"来源名称:     {row['Source_Name']}")
        print()
        print(f"样本特征:")
        for line in str(row['Characteristics']).split('; ')[:5]:
            print(f"  - {line}")
        print()
        print(f"置信度:       {row['Confidence']:.2f}")
        print(f"判定理由:     {row['Filter_Reason']}")
        print(f"V1状态:       {'在V1中' if row['Found_In_V1'] else 'V2新增'}")
        print()
        print("链接:")
        print(f"  GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={row['GSE']}")
        print(f"  样本: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={row['GSM']}")
        print(f"  SRA: https://www.ncbi.nlm.nih.gov/sra/?term={row['SRX']}")
        print()
    
    def open_in_browser(self, index, target='geo'):
        """在浏览器中打开"""
        if index < 1 or index > len(self.filtered_df):
            print(f"❌ 无效的索引: {index}")
            return
        
        row = self.filtered_df.iloc[index - 1]
        
        if target == 'geo':
            url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={row['GSE']}"
        elif target == 'gsm':
            url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={row['GSM']}"
        elif target == 'sra':
            url = f"https://www.ncbi.nlm.nih.gov/sra/?term={row['SRX']}"
        else:
            print(f"❌ 无效的目标: {target} (可选: geo, gsm, sra)")
            return
        
        print(f"✓ 正在打开: {url}")
        webbrowser.open(url)
    
    def export_filtered(self, filename=None):
        """导出当前过滤的结果"""
        if filename is None:
            filename = f"filtered_results_{len(self.filtered_df)}_samples.csv"
        
        self.filtered_df.to_csv(filename, index=False)
        print(f"✓ 已导出 {len(self.filtered_df)} 个样本到: {filename}")
    
    def quick_review_new_high_conf(self):
        """快速复核：新增的高置信度样本"""
        print("\n" + "="*80)
        print("快速复核：V2新增的高置信度样本 (≥0.8)")
        print("="*80)
        
        high_conf_new = self.df[(~self.df['Found_In_V1']) & (self.df['Confidence'] >= 0.8)]
        
        print(f"\n找到 {len(high_conf_new)} 个高置信度新增样本\n")
        
        # 按数据集分组
        for gse in sorted(high_conf_new['GSE'].unique()):
            samples = high_conf_new[high_conf_new['GSE'] == gse]
            max_conf = samples['Confidence'].max()
            print(f"\n【{gse}】- {len(samples)}个样本 (最高置信度: {max_conf:.1f})")
            print("-"*80)
            
            for idx, (_, row) in enumerate(samples.iterrows(), 1):
                print(f"  {idx}. {row['GSM']}: {row['Title'][:50]}")
                print(f"     置信度: {row['Confidence']:.1f} | 判定: {row['Filter_Reason']}")
                print(f"     类型: {row['Library_Strategy']} | 来源: {row['Source_Name'][:40]}")
        
        print("\n" + "="*80)
        return high_conf_new

def interactive_mode():
    """交互式复核模式"""
    reviewer = V2ResultsReviewer()
    
    print("\n" + "╔" + "="*78 + "╗")
    print("║" + " "*25 + "V2 结果复核工具" + " "*38 + "║")
    print("╚" + "="*78 + "╝")
    
    reviewer.show_summary()
    
    print("\n可用命令:")
    print("  summary              - 显示摘要")
    print("  list [n]             - 显示前n个样本 (默认10)")
    print("  detail <n>           - 显示第n个样本详情")
    print("  open <n> [geo|gsm|sra] - 在浏览器中打开")
    print()
    print("  filter_conf <min> [max] - 按置信度过滤")
    print("  filter_new           - 仅显示新增样本")
    print("  filter_v1            - 仅显示V1已有样本")
    print("  filter_gse <GSE>     - 按数据集过滤")
    print("  reset                - 重置过滤")
    print()
    print("  quick_review         - 快速复核高置信度新增样本")
    print("  export [filename]    - 导出当前过滤结果")
    print()
    print("  help                 - 显示帮助")
    print("  quit / exit          - 退出")
    print()
    
    while True:
        try:
            cmd = input("\n复核> ").strip().lower()
            
            if not cmd:
                continue
            
            parts = cmd.split()
            command = parts[0]
            
            if command in ['quit', 'exit', 'q']:
                print("再见！")
                break
            
            elif command == 'help':
                print("\n详细帮助请查看上方命令列表")
            
            elif command == 'summary':
                reviewer.show_summary()
            
            elif command == 'list':
                limit = int(parts[1]) if len(parts) > 1 else 10
                reviewer.show_samples(limit=limit)
            
            elif command == 'detail':
                if len(parts) < 2:
                    print("用法: detail <索引>")
                else:
                    reviewer.show_detail(int(parts[1]))
            
            elif command == 'open':
                if len(parts) < 2:
                    print("用法: open <索引> [geo|gsm|sra]")
                else:
                    target = parts[2] if len(parts) > 2 else 'geo'
                    reviewer.open_in_browser(int(parts[1]), target)
            
            elif command == 'filter_conf':
                if len(parts) < 2:
                    print("用法: filter_conf <min> [max]")
                else:
                    min_conf = float(parts[1])
                    max_conf = float(parts[2]) if len(parts) > 2 else 1.0
                    reviewer.filter_by_confidence(min_conf, max_conf)
                    reviewer.show_samples()
            
            elif command == 'filter_new':
                reviewer.filter_by_status(new_only=True)
                reviewer.show_samples()
            
            elif command == 'filter_v1':
                reviewer.filter_by_status(new_only=False)
                reviewer.show_samples()
            
            elif command == 'filter_gse':
                if len(parts) < 2:
                    print("用法: filter_gse <GSE编号>")
                else:
                    reviewer.filter_by_gse(parts[1].upper())
                    reviewer.show_samples()
            
            elif command == 'reset':
                reviewer.reset_filter()
                reviewer.show_samples()
            
            elif command == 'quick_review':
                reviewer.quick_review_new_high_conf()
            
            elif command == 'export':
                filename = parts[1] if len(parts) > 1 else None
                reviewer.export_filtered(filename)
            
            else:
                print(f"❌ 未知命令: {command}")
                print("输入 'help' 查看可用命令")
        
        except KeyboardInterrupt:
            print("\n\n再见！")
            break
        except Exception as e:
            print(f"❌ 错误: {e}")
            print("请检查命令格式")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # 快速模式
        reviewer = V2ResultsReviewer()
        if sys.argv[1] == 'summary':
            reviewer.show_summary()
        elif sys.argv[1] == 'quick':
            reviewer.quick_review_new_high_conf()
        elif sys.argv[1] == 'new_high':
            reviewer.filter_by_status(new_only=True)
            reviewer.filter_by_confidence(0.8, 1.0)
            reviewer.show_samples(limit=50)
    else:
        # 交互式模式
        interactive_mode()


