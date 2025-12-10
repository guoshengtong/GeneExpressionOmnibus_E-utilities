#!/bin/bash

################################################################################
# SRA 数据批量下载脚本
# SRA Data Batch Download Script
#
# 使用说明 | Usage:
#   1. 确保已安装 SRA Toolkit
#   2. 运行挖掘流水线生成 SRR_accession_list.txt
#   3. 执行此脚本: ./download_sra_data.sh
#
# 功能 | Features:
#   - 批量下载 FASTQ 文件
#   - 断点续传支持
#   - 错误日志记录
#   - 磁盘空间检查
################################################################################

set -e  # 遇到错误时退出

# ==================== 配置区域 | Configuration ====================

# SRR 列表文件
SRR_LIST_FILE="SRR_accession_list.txt"

# 下载目录
DOWNLOAD_DIR="./SRA_Data"

# 日志文件
LOG_FILE="sra_download_$(date +%Y%m%d_%H%M%S).log"
ERROR_LOG_FILE="sra_download_errors_$(date +%Y%m%d_%H%M%S).log"

# 每次下载的并行数（建议不超过3）
PARALLEL_JOBS=2

# 是否压缩 FASTQ 文件 (0=否, 1=是)
COMPRESS_FASTQ=1

# 最小磁盘空间要求 (GB)
MIN_DISK_SPACE_GB=100

# ==================== 颜色输出 | Color Output ====================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ==================== 辅助函数 | Helper Functions ====================

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1" | tee -a "$LOG_FILE"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" | tee -a "$LOG_FILE"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" | tee -a "$LOG_FILE"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE" "$ERROR_LOG_FILE"
}

check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "命令 '$1' 未找到，请先安装 SRA Toolkit"
        log_error "Command '$1' not found, please install SRA Toolkit first"
        log_info "安装指南 | Installation guide: https://github.com/ncbi/sra-tools/wiki"
        exit 1
    fi
}

check_disk_space() {
    # 检查可用磁盘空间
    available_space_kb=$(df -k "$DOWNLOAD_DIR" | tail -1 | awk '{print $4}')
    available_space_gb=$((available_space_kb / 1024 / 1024))
    
    log_info "可用磁盘空间 | Available disk space: ${available_space_gb} GB"
    
    if [ "$available_space_gb" -lt "$MIN_DISK_SPACE_GB" ]; then
        log_warning "磁盘空间不足 ${MIN_DISK_SPACE_GB} GB"
        log_warning "Disk space is less than ${MIN_DISK_SPACE_GB} GB"
        read -p "是否继续? Continue? (y/n): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log_info "用户取消下载 | User cancelled download"
            exit 0
        fi
    fi
}

download_srr() {
    local srr=$1
    local srr_dir="${DOWNLOAD_DIR}/${srr}"
    
    log_info "处理 | Processing: $srr"
    
    # 检查是否已下载
    if [ -d "$srr_dir" ] && [ "$(ls -A $srr_dir/*.fastq* 2>/dev/null)" ]; then
        log_success "$srr 已存在，跳过 | Already exists, skipping"
        return 0
    fi
    
    # 创建临时目录
    mkdir -p "$srr_dir"
    
    # 步骤 1: prefetch - 下载 SRA 文件
    log_info "  [1/2] Downloading SRA file for $srr..."
    if prefetch "$srr" -O "$srr_dir" 2>> "$ERROR_LOG_FILE"; then
        log_success "  SRA file downloaded: $srr"
    else
        log_error "  Failed to download SRA file: $srr"
        return 1
    fi
    
    # 步骤 2: fasterq-dump - 转换为 FASTQ
    log_info "  [2/2] Converting to FASTQ for $srr..."
    if fasterq-dump "$srr_dir/$srr/$srr.sra" -O "$srr_dir" -e $PARALLEL_JOBS 2>> "$ERROR_LOG_FILE"; then
        log_success "  FASTQ files generated: $srr"
        
        # 删除原始 SRA 文件以节省空间
        rm -rf "$srr_dir/$srr"
        
        # 压缩 FASTQ 文件
        if [ "$COMPRESS_FASTQ" -eq 1 ]; then
            log_info "  Compressing FASTQ files for $srr..."
            gzip "$srr_dir"/*.fastq 2>> "$ERROR_LOG_FILE" || log_warning "  Compression failed for $srr"
        fi
        
        log_success "✓ Completed: $srr"
        return 0
    else
        log_error "  Failed to convert to FASTQ: $srr"
        return 1
    fi
}

# ==================== 主程序 | Main Program ====================

main() {
    echo ""
    echo "=========================================="
    echo "   SRA 数据批量下载工具"
    echo "   SRA Data Batch Downloader"
    echo "=========================================="
    echo ""
    
    # 检查必要的命令
    log_info "检查依赖 | Checking dependencies..."
    check_command "prefetch"
    check_command "fasterq-dump"
    check_command "gzip"
    
    # 检查 SRR 列表文件
    if [ ! -f "$SRR_LIST_FILE" ]; then
        log_error "找不到文件: $SRR_LIST_FILE"
        log_error "File not found: $SRR_LIST_FILE"
        log_info "请先运行挖掘流水线生成此文件"
        log_info "Please run the mining pipeline first to generate this file"
        exit 1
    fi
    
    # 读取 SRR 列表
    mapfile -t srr_list < "$SRR_LIST_FILE"
    total_count=${#srr_list[@]}
    
    if [ "$total_count" -eq 0 ]; then
        log_error "SRR 列表为空 | SRR list is empty"
        exit 1
    fi
    
    log_info "找到 ${total_count} 个 SRR 编号 | Found ${total_count} SRR accessions"
    
    # 创建下载目录
    mkdir -p "$DOWNLOAD_DIR"
    
    # 检查磁盘空间
    check_disk_space
    
    # 询问用户确认
    echo ""
    log_warning "即将下载 ${total_count} 个样本的原始数据"
    log_warning "About to download raw data for ${total_count} samples"
    log_warning "这可能需要大量时间和磁盘空间！"
    log_warning "This may take considerable time and disk space!"
    echo ""
    read -p "确认开始下载？ Confirm to start download? (y/n): " -n 1 -r
    echo ""
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_info "用户取消下载 | User cancelled download"
        exit 0
    fi
    
    # 开始下载
    log_info "开始下载 | Starting download..."
    log_info "日志文件 | Log file: $LOG_FILE"
    log_info "错误日志 | Error log: $ERROR_LOG_FILE"
    echo ""
    
    success_count=0
    fail_count=0
    
    # 遍历所有 SRR
    for i in "${!srr_list[@]}"; do
        srr="${srr_list[$i]}"
        
        # 跳过空行
        if [ -z "$srr" ]; then
            continue
        fi
        
        current=$((i + 1))
        log_info "========================================"
        log_info "进度 | Progress: [$current/$total_count]"
        log_info "========================================"
        
        if download_srr "$srr"; then
            ((success_count++))
        else
            ((fail_count++))
        fi
        
        echo ""
    done
    
    # 打印总结
    echo ""
    echo "=========================================="
    echo "          下载完成 | Download Complete"
    echo "=========================================="
    log_success "成功 | Success: $success_count"
    if [ "$fail_count" -gt 0 ]; then
        log_error "失败 | Failed: $fail_count"
        log_info "查看错误日志 | Check error log: $ERROR_LOG_FILE"
    fi
    log_info "数据目录 | Data directory: $DOWNLOAD_DIR"
    log_info "完整日志 | Full log: $LOG_FILE"
    echo "=========================================="
    echo ""
    
    # 如果有失败的，提供重试建议
    if [ "$fail_count" -gt 0 ]; then
        log_info "提示：可以重新运行此脚本，已下载的文件会自动跳过"
        log_info "Tip: You can re-run this script, already downloaded files will be skipped"
    fi
}

# ==================== 信号处理 | Signal Handling ====================

trap 'echo ""; log_warning "下载被中断 | Download interrupted"; exit 130' INT TERM

# ==================== 执行主程序 | Execute Main ====================

main "$@"

