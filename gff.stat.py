import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import numpy as np

def parse_gff3(gff3_file):
    gene_dict = {}
    mrna_to_gene = {}
    mrna_ids = set()  # 存储所有mRNA ID
    
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            seqid = fields[0]
            feature_type = fields[2].lower()
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]
            
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val
            
            if feature_type == 'gene':
                gene_id = attr_dict.get('ID')
                if gene_id:
                    gene_dict[gene_id] = {
                        'seqid': seqid,
                        'start': start,
                        'end': end,
                        'mRNAs': {}
                    }
                    
            elif feature_type == 'mrna':
                parent = attr_dict.get('Parent')
                if not parent:
                    continue
                if ',' in parent:
                    parent = parent.split(',')[0]
                if parent in gene_dict:
                    mrna_id = attr_dict.get('ID')
                    if mrna_id:
                        gene_dict[parent]['mRNAs'][mrna_id] = {
                            'start': start,
                            'end': end,
                            'exons': [],
                            'cds': [],
                            'three_prime_utr': [],
                            'five_prime_utr': [],
                            'exon_details': [],
                            'cds_details': [],
                            'introns': []
                        }
                        mrna_to_gene[mrna_id] = parent
                        mrna_ids.add(mrna_id)  # 添加到mRNA ID集合
            
            elif feature_type in ['exon', 'cds']:
                parent = attr_dict.get('Parent')
                if not parent:
                    continue
                if ',' in parent:
                    parents = parent.split(',')
                else:
                    parents = [parent]
                
                for p in parents:
                    if p in mrna_to_gene:
                        gene_id = mrna_to_gene[p]
                        mrna_info = gene_dict[gene_id]['mRNAs'].get(p)
                        if mrna_info:
                            feature_id = attr_dict.get('ID')
                            if feature_type == 'exon':
                                mrna_info['exons'].append((start, end))
                                mrna_info['exon_details'].append({'start': start, 'end': end, 'id': feature_id})
                            elif feature_type == 'cds':
                                mrna_info['cds'].append((start, end))
                                mrna_info['cds_details'].append({'start': start, 'end': end, 'id': feature_id})
            
            elif feature_type == 'five_prime_utr':
                parent = attr_dict.get('Parent')
                if not parent:
                    continue
                if ',' in parent:
                    parents = parent.split(',')
                else:
                    parents = [parent]
                
                for p in parents:
                    if p in mrna_to_gene:
                        gene_id = mrna_to_gene[p]
                        mrna_info = gene_dict[gene_id]['mRNAs'].get(p)
                        if mrna_info:
                            mrna_info['five_prime_utr'].append((start, end))
            
            elif feature_type == 'three_prime_utr':
                parent = attr_dict.get('Parent')
                if not parent:
                    continue
                if ',' in parent:
                    parents = parent.split(',')
                else:
                    parents = [parent]
                
                for p in parents:
                    if p in mrna_to_gene:
                        gene_id = mrna_to_gene[p]
                        mrna_info = gene_dict[gene_id]['mRNAs'].get(p)
                        if mrna_info:
                            mrna_info['three_prime_utr'].append((start, end))
            
            # 增强intron特征处理
            elif feature_type == 'intron':
                intron_id = attr_dict.get('ID')
                if not intron_id:
                    continue
                
                # 尝试从ID推断mRNA ID
                mrna_id_candidate = None
                if '_intron' in intron_id:
                    # 从intron ID中提取mRNA ID部分
                    mrna_id_candidate = intron_id.rsplit('_intron', 1)[0]
                
                # 检查候选mRNA ID是否有效
                if mrna_id_candidate and mrna_id_candidate in mrna_ids:
                    gene_id = mrna_to_gene.get(mrna_id_candidate)
                    if gene_id:
                        mrna_info = gene_dict[gene_id]['mRNAs'].get(mrna_id_candidate)
                        if mrna_info:
                            mrna_info['introns'].append({'start': start, 'end': end, 'id': intron_id})
    
    return gene_dict

def calculate_feature_length(feature_list):
    return sum(end - start + 1 for start, end in feature_list)

def process_gff3(gff3_file, prefix):
    gene_dict = parse_gff3(gff3_file)
    
    detail_file = f"{prefix}.gene.information.stat.tsv"
    summary_file = f"{prefix}.summary.information.stat.tsv"
    feature_file = f"{prefix}.intron.exon.cds.stat.tsv"
    
    # Data for summary stats
    gene_lengths = []
    mrna_lengths = []
    exon_counts = []
    exon_total_lengths = []
    intron_total_lengths = []
    
    # Data for boxplots and histograms
    gene_lengths_bp = []
    mrna_lengths_bp = []
    exon_total_per_mrna = []
    intron_total_per_mrna = []
    
    # Data for histograms
    all_gene_lengths = []
    all_intron_lengths = []
    all_exon_lengths = []
    all_mrna_lengths = []
    
    processed_genes = set()
    
    with open(detail_file, 'w') as f_detail, open(feature_file, 'w') as f_feature:
        # 写入特征文件表头
        feature_headers = [
            "chr_id", "start", "end", "feature_id", "feature_type",
            "gene_id", "mrna_id", "exon_count", "cds_count", "intron_count", "length"
        ]
        f_feature.write("\t".join(feature_headers) + "\n")
        
        # 原有detail文件表头
        detail_headers = [
            "seqid", "gene_start", "gene_end", "gene_id", "gene_length",
            "mrna_id", "mrna_length", "exon_count", "intron_count", "cds_count",
            "three_utr_count", "five_utr_count", "exon_total_length", "intron_total_length",
            "cds_total_length", "three_utr_total_length", "five_utr_total_length",
            "exon_avg_length", "intron_avg_length", "cds_avg_length", "three_utr_avg_length",
            "five_utr_avg_length", "intron_per_mrna", "intron_per_gene"
        ]
        f_detail.write("\t".join(detail_headers) + "\n")
        
        for gene_id, gene_info in gene_dict.items():
            seqid = gene_info['seqid']
            gene_start = gene_info['start']
            gene_end = gene_info['end']
            gene_length = gene_end - gene_start + 1
            
            if gene_id not in processed_genes:
                gene_lengths.append(gene_length)
                gene_lengths_bp.append(gene_length)
                all_gene_lengths.append(gene_length)
                processed_genes.add(gene_id)
            
            for mrna_id, mrna_info in gene_info['mRNAs'].items():
                mrna_start = mrna_info['start']
                mrna_end = mrna_info['end']
                mrna_length = mrna_end - mrna_start + 1
                mrna_lengths.append(mrna_length)
                mrna_lengths_bp.append(mrna_length)
                all_mrna_lengths.append(mrna_length)
                
                # Feature counts
                exon_count = len(mrna_info['exons'])
                cds_count = len(mrna_info['cds'])
                three_utr_count = len(mrna_info['three_prime_utr'])
                five_utr_count = len(mrna_info['five_prime_utr'])
                intron_count = len(mrna_info['introns'])  # 直接从introns列表获取
                
                # Total lengths
                exon_total = calculate_feature_length(mrna_info['exons'])
                cds_total = calculate_feature_length(mrna_info['cds'])
                three_utr_total = calculate_feature_length(mrna_info['three_prime_utr'])
                five_utr_total = calculate_feature_length(mrna_info['five_prime_utr'])
                intron_total = calculate_feature_length([(i['start'], i['end']) for i in mrna_info['introns']])
                
                # 如果没有直接解析到intron，则通过exon计算
                if intron_total == 0 and exon_count > 1:
                    intron_total = mrna_length - exon_total
                
                exon_total_per_mrna.append(exon_total)
                intron_total_per_mrna.append(intron_total)
                exon_counts.append(exon_count)
                exon_total_lengths.append(exon_total)
                intron_total_lengths.append(intron_total)
                
                # Collect exon and intron lengths for histograms
                for exon in mrna_info['exons']:
                    all_exon_lengths.append(exon[1] - exon[0] + 1)
                
                # Calculate intron lengths
                sorted_exons = sorted(mrna_info['exons'], key=lambda x: x[0])
                for i in range(len(sorted_exons) - 1):
                    intron_length = sorted_exons[i+1][0] - sorted_exons[i][1] - 1
                    if intron_length > 0:
                        all_intron_lengths.append(intron_length)
                
                # Average lengths
                exon_avg = exon_total / exon_count if exon_count > 0 else 0
                intron_avg = intron_total / intron_count if intron_count > 0 else 0
                cds_avg = cds_total / cds_count if cds_count > 0 else 0
                three_utr_avg = three_utr_total / three_utr_count if three_utr_count > 0 else 0
                five_utr_avg = five_utr_total / five_utr_count if five_utr_count > 0 else 0
                
                # Intron ratios
                intron_per_mrna = intron_total / mrna_length if mrna_length > 0 else 0
                intron_per_gene = intron_total / gene_length if gene_length > 0 else 0
                
                # Write detail row
                row = [
                    seqid, str(gene_start), str(gene_end), gene_id, str(gene_length),
                    mrna_id, str(mrna_length), str(exon_count), str(intron_count), str(cds_count),
                    str(three_utr_count), str(five_utr_count), str(exon_total), str(intron_total),
                    str(cds_total), str(three_utr_total), str(five_utr_total),
                    f"{exon_avg:.2f}", f"{intron_avg:.2f}", f"{cds_avg:.2f}", 
                    f"{three_utr_avg:.2f}", f"{five_utr_avg:.2f}",
                    f"{intron_per_mrna:.4f}", f"{intron_per_gene:.4f}"
                ]
                f_detail.write("\t".join(row) + "\n")
                
                # 写入特征统计信息
                # 处理exon特征
                for exon in mrna_info['exon_details']:
                    if exon['id']:
                        length = exon['end'] - exon['start'] + 1
                        feature_row = [
                            seqid,
                            str(exon['start']),
                            str(exon['end']),
                            exon['id'],
                            "exon",
                            gene_id,
                            mrna_id,
                            str(exon_count),
                            str(cds_count),
                            str(intron_count),
                            str(length)
                        ]
                        f_feature.write("\t".join(feature_row) + "\n")
                
                # 处理cds特征
                for cds in mrna_info['cds_details']:
                    if cds['id']:
                        length = cds['end'] - cds['start'] + 1
                        feature_row = [
                            seqid,
                            str(cds['start']),
                            str(cds['end']),
                            cds['id'],
                            "cds",
                            gene_id,
                            mrna_id,
                            str(exon_count),
                            str(cds_count),
                            str(intron_count),
                            str(length)
                        ]
                        f_feature.write("\t".join(feature_row) + "\n")
                
                # 处理intron特征
                for intron in mrna_info['introns']:
                    if intron['id']:
                        length = intron['end'] - intron['start'] + 1
                        feature_row = [
                            seqid,
                            str(intron['start']),
                            str(intron['end']),
                            intron['id'],
                            "intron",
                            gene_id,
                            mrna_id,
                            str(exon_count),
                            str(cds_count),
                            str(intron_count),
                            str(length)
                        ]
                        f_feature.write("\t".join(feature_row) + "\n")
    
    # Calculate summary statistics
    num_genes = len(gene_lengths)
    num_mrnas = len(mrna_lengths)
    avg_gene_len = sum(gene_lengths) / num_genes if num_genes > 0 else 0
    avg_mrna_len = sum(mrna_lengths) / num_mrnas if num_mrnas > 0 else 0
    avg_exon_count = sum(exon_counts) / num_mrnas if num_mrnas > 0 else 0
    avg_intron_count = sum(intron_total_lengths) / num_mrnas if num_mrnas > 0 else 0
    
    total_exon_len = sum(exon_total_lengths)
    total_exon_count = sum(exon_counts)
    avg_exon_len = total_exon_len / total_exon_count if total_exon_count > 0 else 0
    
    total_intron_len = sum(intron_total_lengths)
    total_intron_count = total_exon_count - num_mrnas
    avg_intron_len = total_intron_len / total_intron_count if total_intron_count > 0 else 0
    
    # Write summary file
    with open(summary_file, 'w') as f_summary:
        headers = ["sample", "num_genes", "avg_gene_length", "avg_mrna_length", 
                   "avg_exon_length", "avg_exon_count", "avg_intron_length", "avg_intron_count"]
        f_summary.write("\t".join(headers) + "\n")
        row = [
            prefix, str(num_genes), f"{avg_gene_len:.2f}", f"{avg_mrna_len:.2f}",
            f"{avg_exon_len:.2f}", f"{avg_exon_count:.2f}", f"{avg_intron_len:.2f}", f"{avg_intron_count:.2f}"
        ]
        f_summary.write("\t".join(row) + "\n")
    

def main():
    parser = argparse.ArgumentParser(description='Process GFF3 files and generate statistics.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', '--config', help='Configuration file with sample and GFF3 paths')
    group.add_argument('-g', '--gff3', help='Single GFF3 file to process')
    parser.add_argument('-p', '--prefix', help='Output prefix (used with -g)')
    
    args = parser.parse_args()
    
    if args.config:
        with open(args.config, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                sample = parts[0]
                gff3_path = parts[1]
                process_gff3(gff3_path, sample)
                
    elif args.gff3:
        if not args.prefix:
            parser.error("Prefix (-p) is required when using -g")
        process_gff3(args.gff3, args.prefix)

if __name__ == "__main__":
    main()
