import argparse
import sys

def parse_attributes(attr_str):
    attributes = {}
    for part in attr_str.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            attributes[key] = value
    return attributes

def format_attributes(attrs):
    return ';'.join(f"{key}={value}" for key, value in attrs.items())

def process_gene_module(module_lines):
    output_lines = []
    current_mrna = None
    current_strand = None
    exons = []
    cdss = []
    
    for fields in module_lines:
        # 如果是注释行（字符串）
        if isinstance(fields, str):
            output_lines.append(fields)
            continue
            
        # 检查字段数量
        if len(fields) < 9:
            output_lines.append('\t'.join(fields))
            continue
            
        seqid, source, feature_type, start, end, score, strand, phase, attributes_str = fields[:9]
        attrs = parse_attributes(attributes_str)
        
        # 处理基因和mRNA
        if feature_type == 'gene' or feature_type == 'mRNA':
            # 处理前一个mRNA的特征
            if current_mrna:
                process_mrna_features(output_lines, current_mrna, current_strand, exons, cdss)
                exons = []
                cdss = []
            
            output_lines.append('\t'.join(fields))
            
            if feature_type == 'mRNA' and 'ID' in attrs:
                current_mrna = attrs['ID']
                current_strand = strand
            else:
                current_mrna = None
                current_strand = None
            continue
        
        # 收集外显子和CDS
        if current_mrna:
            if feature_type == 'exon':
                exons.append(fields)
            elif feature_type == 'CDS':
                cdss.append(fields)
            else:
                output_lines.append('\t'.join(fields))
        else:
            output_lines.append('\t'.join(fields))
    
    # 处理最后一个mRNA
    if current_mrna:
        process_mrna_features(output_lines, current_mrna, current_strand, exons, cdss)
    
    return output_lines

def process_mrna_features(output_lines, mrna_id, strand, exons, cdss):
    if not exons:
        return
    
    # 按起始位置排序
    exons_sorted = sorted(exons, key=lambda x: int(x[3]))
    cdss_sorted = sorted(cdss, key=lambda x: int(x[3])) if cdss else []
    
    # 根据链方向调整顺序
    if strand == '-':
        exons_sorted = exons_sorted[::-1]
        if cdss_sorted:
            cdss_sorted = cdss_sorted[::-1]
    
    # 处理外显子
    for i, exon in enumerate(exons_sorted, 1):
        attrs = parse_attributes(exon[8])
        attrs['ID'] = f"{mrna_id}_exon{i}"
        exon[8] = format_attributes(attrs)
        output_lines.append('\t'.join(exon))
    
    # 添加内含子（需要至少2个外显子）
    if len(exons_sorted) > 1:
        # 使用基因组顺序计算内含子
        genomic_exons = sorted(exons, key=lambda x: int(x[3]))
        for i in range(len(genomic_exons) - 1):
            prev_end = int(genomic_exons[i][4])
            next_start = int(genomic_exons[i+1][3])
            
            # 确保内含子坐标有效
            if prev_end + 1 <= next_start - 1:
                intron_start = str(prev_end + 1)
                intron_end = str(next_start - 1)
                
                # 创建内含子记录
                intron_fields = [
                    genomic_exons[0][0],  # seqid
                    genomic_exons[0][1],  # source
                    'intron',             # type
                    intron_start,
                    intron_end,
                    '.',                 # score
                    strand,               # strand
                    '.',                 # phase
                    f'ID={mrna_id}_intron{i+1}'
                ]
                output_lines.append('\t'.join(intron_fields))
    
    # 处理CDS
    for i, cds in enumerate(cdss_sorted, 1):
        attrs = parse_attributes(cds[8])
        attrs['ID'] = f"{mrna_id}_cds{i}"
        cds[8] = format_attributes(attrs)
        output_lines.append('\t'.join(cds))

def main():
    parser = argparse.ArgumentParser(description='Process GFF file to add introns and rename features.')
    parser.add_argument('-i', '--input', required=True, help='Input GFF file')
    parser.add_argument('-o', '--output', required=True, help='Output GFF file')
    args = parser.parse_args()
    
    # 读取输入文件
    with open(args.input, 'r') as f:
        input_lines = f.readlines()
    
    all_output_lines = []
    current_module = []
    
    # 按基因模块处理
    for line in input_lines:
        stripped = line.strip()
        
        # 处理注释行
        if line.startswith('#'):
            if current_module:
                processed = process_gene_module(current_module)
                all_output_lines.extend(processed)
                current_module = []
            all_output_lines.append(line)
            continue
            
        # 处理空行
        if not stripped:
            if current_module:
                processed = process_gene_module(current_module)
                all_output_lines.extend(processed)
                current_module = []
            all_output_lines.append(line)
            continue
            
        fields = stripped.split('\t')
        
        # 发现新基因时处理当前模块
        if len(fields) >= 3 and fields[2] == 'gene':
            if current_module:
                processed = process_gene_module(current_module)
                all_output_lines.extend(processed)
            current_module = [fields]
        else:
            current_module.append(fields)
    
    # 处理最后一个模块
    if current_module:
        processed = process_gene_module(current_module)
        all_output_lines.extend(processed)
    
    # 写入输出文件
    with open(args.output, 'w') as f:
        for line in all_output_lines:
            f.write(line + '\n' if not line.endswith('\n') else line)

if __name__ == "__main__":
    main()
