#!/usr/bin/env python3
"""
安装测试脚本 - GEO 肺部转移瘤数据挖掘流水线
Installation Test Script - GEO Lung Metastasis Mining Pipeline

运行此脚本以验证所有依赖是否正确安装。
Run this script to verify all dependencies are correctly installed.
"""

import sys
import importlib

def test_python_version():
    """测试 Python 版本"""
    print("测试 Python 版本 | Testing Python version...")
    version = sys.version_info
    if version.major >= 3 and version.minor >= 8:
        print(f"  ✓ Python {version.major}.{version.minor}.{version.micro}")
        return True
    else:
        print(f"  ✗ Python {version.major}.{version.minor}.{version.micro} (需要 3.8+)")
        return False

def test_module_import(module_name, package_name=None):
    """测试模块导入"""
    if package_name is None:
        package_name = module_name
    
    try:
        module = importlib.import_module(module_name)
        version = getattr(module, '__version__', 'unknown')
        print(f"  ✓ {package_name} ({version})")
        return True
    except ImportError as e:
        print(f"  ✗ {package_name} - 未安装 | Not installed")
        print(f"    错误 | Error: {e}")
        return False

def test_config():
    """测试配置文件"""
    print("\n测试配置 | Testing configuration...")
    try:
        import config
        
        # 检查邮箱是否配置
        if config.ENTREZ_EMAIL == "your.email@example.com":
            print("  ⚠ 邮箱地址未配置 | Email not configured")
            print("    请编辑 config.py 设置 ENTREZ_EMAIL")
            print("    Please edit config.py to set ENTREZ_EMAIL")
            return False
        else:
            print(f"  ✓ 邮箱已配置 | Email configured: {config.ENTREZ_EMAIL}")
        
        # 检查其他配置
        print(f"  ✓ GEO 缓存目录 | GEO cache dir: {config.GEO_CACHE_DIR}")
        print(f"  ✓ API 延迟 | API delay: {config.API_DELAY}s")
        
        return True
        
    except ImportError as e:
        print(f"  ✗ 无法导入配置文件 | Cannot import config: {e}")
        return False
    except Exception as e:
        print(f"  ✗ 配置文件错误 | Config error: {e}")
        return False

def test_main_script():
    """测试主脚本"""
    print("\n测试主脚本 | Testing main script...")
    try:
        from geo_lung_metastasis_miner import GEOLungMetastasisMiner
        print("  ✓ 主脚本导入成功 | Main script imported successfully")
        
        # 尝试创建实例（不运行）
        try:
            miner = GEOLungMetastasisMiner()
            print("  ✓ 挖掘器实例创建成功 | Miner instance created successfully")
            return True
        except Exception as e:
            print(f"  ⚠ 实例创建警告 | Instance creation warning: {e}")
            return True  # 仍然返回True，因为导入成功
            
    except ImportError as e:
        print(f"  ✗ 主脚本导入失败 | Main script import failed: {e}")
        return False
    except Exception as e:
        print(f"  ✗ 主脚本错误 | Main script error: {e}")
        return False

def test_sra_toolkit():
    """测试 SRA Toolkit（可选）"""
    print("\n测试 SRA Toolkit (可选) | Testing SRA Toolkit (optional)...")
    import subprocess
    
    commands = ['prefetch', 'fasterq-dump']
    all_found = True
    
    for cmd in commands:
        try:
            result = subprocess.run(
                [cmd, '--help'], 
                capture_output=True, 
                timeout=5
            )
            if result.returncode == 0 or result.returncode == 1:  # 有些命令 --help 返回 1
                print(f"  ✓ {cmd} 已安装 | installed")
            else:
                print(f"  ✗ {cmd} 未正确安装 | not properly installed")
                all_found = False
        except FileNotFoundError:
            print(f"  ⚠ {cmd} 未找到 | not found")
            all_found = False
        except subprocess.TimeoutExpired:
            print(f"  ⚠ {cmd} 超时 | timeout")
            all_found = False
        except Exception as e:
            print(f"  ⚠ {cmd} 测试失败 | test failed: {e}")
            all_found = False
    
    if not all_found:
        print("\n  提示：SRA Toolkit 仅在下载数据时需要")
        print("  Tip: SRA Toolkit is only required for downloading data")
        print("  安装指南 | Installation guide:")
        print("    MacOS: brew install sratoolkit")
        print("    Linux: https://github.com/ncbi/sra-tools/wiki")
    
    return all_found

def main():
    """主测试函数"""
    print("="*70)
    print("GEO 肺部转移瘤数据挖掘流水线 - 安装测试")
    print("GEO Lung Metastasis Mining Pipeline - Installation Test")
    print("="*70)
    print()
    
    results = {}
    
    # 测试 Python 版本
    results['python'] = test_python_version()
    
    # 测试必需的 Python 模块
    print("\n测试 Python 依赖 | Testing Python dependencies...")
    required_modules = [
        ('Bio', 'biopython'),
        ('GEOparse', 'GEOparse'),
        ('pandas', 'pandas')
    ]
    
    all_modules_ok = True
    for module_name, package_name in required_modules:
        if not test_module_import(module_name, package_name):
            all_modules_ok = False
    
    results['modules'] = all_modules_ok
    
    # 测试可选模块
    print("\n测试可选依赖 | Testing optional dependencies...")
    optional_modules = [
        ('rich', 'rich')
    ]
    
    for module_name, package_name in optional_modules:
        test_module_import(module_name, package_name)
    
    # 测试配置
    results['config'] = test_config()
    
    # 测试主脚本
    results['main_script'] = test_main_script()
    
    # 测试 SRA Toolkit
    results['sra_toolkit'] = test_sra_toolkit()
    
    # 总结
    print("\n" + "="*70)
    print("测试总结 | Test Summary")
    print("="*70)
    
    if results['python'] and results['modules'] and results['config'] and results['main_script']:
        print("\n✓ 所有必需组件测试通过！")
        print("✓ All required components passed!")
        print("\n您可以开始使用流水线：")
        print("You can start using the pipeline:")
        print("  python geo_lung_metastasis_miner.py")
        
        if not results['sra_toolkit']:
            print("\n⚠ 提示：如需下载原始数据，请安装 SRA Toolkit")
            print("⚠ Tip: Install SRA Toolkit if you need to download raw data")
        
        print()
        return True
    else:
        print("\n✗ 部分测试失败，请解决以下问题：")
        print("✗ Some tests failed, please resolve the following issues:")
        print()
        
        if not results['python']:
            print("  - 升级 Python 到 3.8 或更高版本")
            print("    Upgrade Python to 3.8 or higher")
        
        if not results['modules']:
            print("  - 安装缺失的 Python 包:")
            print("    Install missing Python packages:")
            print("    pip install -r requirements.txt")
        
        if not results['config']:
            print("  - 配置 config.py 文件")
            print("    Configure config.py file")
        
        if not results['main_script']:
            print("  - 检查主脚本是否损坏")
            print("    Check if main script is corrupted")
        
        print()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)


