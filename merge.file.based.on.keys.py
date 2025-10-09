#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
from collections import defaultdict
import os
import re

def parse_column_spec(spec):
    """解析列规范字符串（如'1,3'或'1-3'），返回从0开始的列索引列表"""
    indices = []
    parts = spec.split(',')
    for part in parts:
        if '-' in part:
            start, end = part.split('-')
            start = int(start.strip())
            end = int(end.strip())
            indices.extend(range(start-1, end))  # 从1开始的索引转为0开始
        else:
            indices.append(int(part.strip())-1)  # 从1开始的索引转为0开始
    return sorted(set(indices))  # 去重并排序

def check_columns_in_file(file_obj, column_indices, file_name, file_type, delimiter):
    """检查文件是否包含所有指定的列索引"""
    file_obj.seek(0)
    found_valid_line = False
    
    for line in file_obj:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
            
        # 分割行 - 根据分隔符类型
        if delimiter == "whitespace":
            parts = re.split(r'\s+', stripped)
        else:
            parts = stripped.split(delimiter)
        
        # 检查所有列索引是否都存在
        if all(idx < len(parts) for idx in column_indices):
            found_valid_line = True
            break
    
    file_obj.seek(0)
    return found_valid_line

def process_reference_file(file_obj, key_indices, delimiter):
    """处理参考文件，返回字典和错误行（支持多列键值）"""
    ref_dict = defaultdict(list)
    blank_lines = []
    comment_lines = []
    error_lines = []
    unmatched_lines = []
    total_lines = 0
    valid_lines = 0
    collected_keys = []  # 收集前三个不同的键
    seen_keys = set()
    
    for line in file_obj:
        total_lines += 1
        stripped = line.strip()
        
        # 处理空行
        if not stripped:
            blank_lines.append(line.rstrip())
            continue
            
        # 处理注释行
        if stripped.startswith('#'):
            comment_lines.append(stripped)
            continue
            
        # 分割行
        if delimiter == "whitespace":
            parts = re.split(r'\s+', stripped)
        else:
            parts = stripped.split(delimiter)
        
        # 检查所有键列是否存在
        if any(idx >= len(parts) for idx in key_indices):
            error_lines.append(stripped)
            continue
            
        # 构建复合键
        key = tuple(parts[idx] for idx in key_indices)
        ref_dict[key].append(stripped)
        valid_lines += 1
        
        # 收集前三个不同的键
        if key not in seen_keys and len(collected_keys) < 3:
            seen_keys.add(key)
            collected_keys.append(key)

    return ref_dict, blank_lines, comment_lines, error_lines, unmatched_lines, total_lines, valid_lines, collected_keys

def process_query_file(file_obj, key_indices, delimiter, ref_dict):
    """处理查询文件，返回组合行和错误行（支持多列键值）"""
    combine_lines = []
    blank_lines = []
    comment_lines = []
    unmatched_lines = []
    error_lines = []
    total_lines = 0
    valid_lines = 0
    matched_lines = 0
    matched_ref_keys = set()  # 记录匹配到的参考文件键
    collected_keys = []  # 收集前三个不同的键
    seen_keys = set()
    
    for line in file_obj:
        total_lines += 1
        stripped = line.strip()
        
        # 处理空行
        if not stripped:
            blank_lines.append(line.rstrip())
            continue
            
        # 处理注释行
        if stripped.startswith('#'):
            comment_lines.append(stripped)
            continue
            
        # 分割行
        if delimiter == "whitespace":
            parts = re.split(r'\s+', stripped)
        else:
            parts = stripped.split(delimiter)
        
        # 检查所有键列是否存在
        if any(idx >= len(parts) for idx in key_indices):
            error_lines.append(stripped)
            continue
            
        # 构建复合键
        key = tuple(parts[idx] for idx in key_indices)
        valid_lines += 1
        
        # 收集前三个不同的键
        if key not in seen_keys and len(collected_keys) < 3:
            seen_keys.add(key)
            collected_keys.append(key)
        
        if key in ref_dict:
            matched_lines += 1
            # 记录匹配到的参考文件键
            matched_ref_keys.add(key)
            # 匹配成功：先输出参考文件内容，再输出查询文件内容
            for ref_line in ref_dict[key]:
                if delimiter == "whitespace":
                    combine_lines.append(f"{ref_line} {stripped}")
                else:
                    combine_lines.append(f"{ref_line}{delimiter}{stripped}")
        else:
            unmatched_lines.append(stripped)
    
    return combine_lines, blank_lines, comment_lines, unmatched_lines, error_lines, total_lines, valid_lines, matched_lines, collected_keys, matched_ref_keys

def write_output_files(prefix, combine_lines, rf_comment_lines, rf_unmatched_lines, rf_error_lines, qf_comment_lines, qf_unmatched_lines, qf_error_lines):
    """写入输出文件"""
    # 组合结果文件
    with open(f"{prefix}.combine.file.tsv", 'w') as file:
        for line in combine_lines:
            file.write(line + '\n')
    
    # 参考文件未匹配行文件
    with open(f"{prefix}.rf.unmatched.line.tsv", 'w') as file:
        for line in rf_comment_lines:
            file.write(line + '\n')
        for line in rf_unmatched_lines:
            file.write(line + '\n')
    
    # 参考文件错误行文件
    with open(f"{prefix}.rf.error.line.tsv", 'w') as file:
        for line in rf_error_lines:
            file.write(line + '\n')
    
    # 查询文件未匹配行文件
    with open(f"{prefix}.qf.unmatched.line.tsv", 'w') as file:
        for line in qf_comment_lines:
            file.write(line + '\n')
        for line in qf_unmatched_lines:
            file.write(line + '\n')
    
    # 查询文件错误行文件
    with open(f"{prefix}.qf.error.line.tsv", 'w') as file:
        for line in qf_error_lines:
            file.write(line + '\n')

def generate_statistics(args, ref_stats, qry_stats):
    """生成统计信息并输出到屏幕和日志"""
    # 解包统计信息
    ref_dict, ref_blanks, ref_comments, ref_errors, ref_unmatched, ref_total, ref_valid, ref_keys = ref_stats
    combine_lines, qry_blanks, qry_comments, qry_unmatched, qry_errors, qry_total, qry_valid, qry_matched, qry_keys, matched_ref_keys = qry_stats
    
    # 计算参考文件匹配情况
    ref_unique_keys = set(ref_dict.keys())
    
    ref_matched_keys = len(matched_ref_keys)
    ref_matched_lines = len(combine_lines)
    ref_unmatched_lines = len(ref_dict) - ref_matched_keys
    
    # 获取分隔符显示名称
    sep_display = "空格或制表符" if args.separator == "whitespace" else f"'{args.separator}'"
    
    # 创建统计信息文本
    stats = [
        "=" * 50,
        "文件处理统计信息",
        "=" * 50,
        f"参考文件: {args.ref_file.name}",
        f"  总行数: {ref_total}",
        f"  有效行数: {ref_valid} (去除空行和注释行)",
        f"  空白行: {len(ref_blanks)} (仅记录，不写入文件)",
        f"  注释行: {len(ref_comments)} (已写入未匹配文件)",
        f"  错误行数: {len(ref_errors)} (键值列不存在)",
        f"  唯一键值数: {len(ref_unique_keys)}",
        f"  匹配键值数: {ref_matched_keys}",
        f"  匹配行数: {ref_matched_lines}",
        f"  未匹配行数: {ref_unmatched_lines}",
        "",
        f"查询文件: {args.query_file.name}",
        f"  总行数: {qry_total}",
        f"  有效行数: {qry_valid} (去除空行和注释行)",
        f"  空白行: {len(qry_blanks)} (仅记录，不写入文件)",
        f"  注释行: {len(qry_comments)} (已写入未匹配文件)",
        f"  错误行数: {len(qry_errors)} (键值列不存在)",
        f"  匹配行数: {qry_matched}",
        f"  未匹配行数: {len(qry_unmatched)}",
        "",
        "输出文件:",
        f"  匹配结果: {args.prefix}.combine.file.tsv ({len(combine_lines)} 行)",
        f"  参考文件未匹配行: {args.prefix}.rf.unmatched.line.tsv ({len(ref_comments) + len(ref_unmatched)} 行)",
        f"  参考文件错误行: {args.prefix}.rf.error.line.tsv ({len(ref_errors)} 行)",
        f"  查询文件未匹配行: {args.prefix}.qf.unmatched.line.tsv ({len(qry_comments) + len(qry_unmatched)} 行)",
        f"  查询文件错误行: {args.prefix}.qf.error.line.tsv ({len(qry_errors)} 行)",
        f"  使用的分隔符: {sep_display}",
        f"  键值列: 参考文件={args.ref_column}, 查询文件={args.query_column}",
        f"  合并顺序: 参考文件内容 + 查询文件内容",
        "=" * 50
    ]
    
    # 输出到屏幕
    for line in stats:
        print(line)
    
    # 检查键值是否完全不同
    if ref_matched_keys == 0 and len(ref_dict) > 0 and qry_valid > 0:
        error_msg = "错误：给定键值完全不同，无法做匹配识别，无法完成文件合并"
        print(f"\n{error_msg}", file=sys.stderr)
        print(f"参考文件键值示例: {ref_keys}", file=sys.stderr)
        print(f"查询文件键值示例: {qry_keys}", file=sys.stderr)
        
        with open(f"{args.prefix}.log", 'a') as log_file:
            log_file.write(f"\n{error_msg}\n")
            log_file.write(f"参考文件键值示例: {ref_keys}\n")
            log_file.write(f"查询文件键值示例: {qry_keys}\n")
        
        sys.exit(1)
    
    print(f"\n相关提示:",
          f"1. {args.prefix}.combine.file.tsv 只包含匹配成功的行",
          f"2. {args.prefix}.rf.unmatched.line.tsv 包含参考文件中的注释行和未匹配行",
          f"3. {args.prefix}.qf.unmatched.line.tsv 包含查询文件中的注释行和未匹配行",
          f"4. {args.prefix}.rf.error.line.tsv 包含参考文件中键值列不存在的行",
          f"5. {args.prefix}.qf.error.line.tsv 包含查询文件中键值列不存在的行",
          f"6. 空白行仅记录数量，不写入任何文件",
          sep="\n")
    
    # 输出到日志文件
    with open(f"{args.prefix}.log", 'w') as log_file:
        for line in stats:
            log_file.write(line + '\n')
    
    return stats

def runcominbefile(args):
    """主逻辑处理函数"""
    # 解析列规范
    ref_key_indices = parse_column_spec(args.ref_column)
    qry_key_indices = parse_column_spec(args.query_column)
    
    # 检查键值列是否超出文件范围
    if not check_columns_in_file(args.ref_file, ref_key_indices, args.ref_file.name, "参考文件", args.separator):
        print(f"错误：参考文件 '{args.ref_file.name}' 中不存在指定的列 {args.ref_column}", file=sys.stderr)
        sep_display = "空格或制表符" if args.separator == "whitespace" else f"'{args.separator}'"
        print(f"请确认文件格式和列分隔符（当前识别的分隔符为 {sep_display}）", file=sys.stderr)
        sys.exit(1)
    
    if not check_columns_in_file(args.query_file, qry_key_indices, args.query_file.name, "查询文件", args.separator):
        print(f"错误：查询文件 '{args.query_file.name}' 中不存在指定的列 {args.query_column}", file=sys.stderr)
        sep_display = "空格或制表符" if args.separator == "whitespace" else f"'{args.separator}'"
        print(f"请确认文件格式和列分隔符（当前识别的分隔符: {sep_display}）", file=sys.stderr)
        sys.exit(1)
    
    # 重置文件指针
    args.ref_file.seek(0)
    args.query_file.seek(0)
    
    # 处理参考文件
    ref_stats = process_reference_file(args.ref_file, ref_key_indices, args.separator)
    ref_dict = ref_stats[0]
    
    # 处理查询文件
    qry_stats = process_query_file(args.query_file, qry_key_indices, args.separator, ref_dict)
    
    # 提取参考文件中未匹配的行
    matched_ref_keys = qry_stats[-1]  # 最后一个元素是匹配的参考文件键
    ref_unmatched_lines = []
    
    for key, lines in ref_dict.items():
        if key not in matched_ref_keys:
            ref_unmatched_lines.extend(lines)
    
    # 更新参考文件统计中的未匹配行
    ref_stats = list(ref_stats)
    ref_stats[4] = ref_unmatched_lines  # 索引4是未匹配行
    
    # 写入输出文件
    write_output_files(
        args.prefix, 
        qry_stats[0],  # combine_lines
        ref_stats[2],  # 参考文件注释行
        ref_stats[4],  # 参考文件未匹配行
        ref_stats[3],  # 参考文件错误行
        qry_stats[2],  # 查询文件注释行
        qry_stats[3],  # 查询文件未匹配行
        qry_stats[4]   # 查询文件错误行
    )
    
    # 生成统计信息
    generate_statistics(args, tuple(ref_stats), qry_stats)

def main():
    """主函数，解析命令行参数并调用处理函数"""
    parser = argparse.ArgumentParser(
        description='文件键值匹配工具：基于指定列匹配两个文件（支持多列键值）',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-rf', '--ref_file', required=True, type=argparse.FileType('r'),
                        help='参考文件路径')
    parser.add_argument('-rc', '--ref_column', type=str, default="1",
                        help='参考文件键列索引（从1开始），支持多列（如：1,3或1-3或1,4-5）')
    parser.add_argument('-qf', '--query_file', required=True, type=argparse.FileType('r'),
                        help='查询文件路径')
    parser.add_argument('-qc', '--query_column', type=str, default="1",
                        help='查询文件键列索引（从1开始），支持多列（如：1,3或1-3或1,4-5）')
    parser.add_argument('-sp', '--separator', default="\t", 
                        help='列分隔符（默认: 制表符, 或指定特定分隔符如","，或"whitespace"表示空格/制表符）')
    parser.add_argument('-pf', '--prefix', required=True, 
                        help='输出文件前缀')
    
    args = parser.parse_args()
    
    # 处理分隔符参数的特殊值
    if args.separator.lower() in ["tab", "t"]:
        args.separator = "\t"
    elif args.separator.lower() in ["space", "whitespace", "ws"]:
        args.separator = "whitespace"
    
    runcominbefile(args)

if __name__ == "__main__":
    main()
