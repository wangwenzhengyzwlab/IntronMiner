import pandas as pd
from plotnine import (
    ggplot, aes, geom_point, geom_hline, facet_wrap, scale_color_manual,
    scale_y_continuous, labs, theme_classic, theme, element_blank,
    element_text, element_rect
)
from mizani.formatters import number_format

# -------- 1) 读取文件 --------
infile = "/share/org/YZWL/yzwl_guoxp/Documents/Jeffrey/CML304.chr.tsv"
outfile = "/share/org/YZWL/yzwl_guoxp/Documents/Jeffrey/CML304_Intron_Diff_ByChr_Horizontal_PosNeg.pdf"

colnames_in = [
    "seqid", "gene_start", "gene_end", "mRNA_id", "type", "gene_id",
    "mRNA_id_dup", "exon_number", "length.bp", "ref_chr_id",
    "ref_gene_start", "ref_gene_end", "ref_exon_number", "ref_length.bp", "dif.length.bp"
]


df_raw = pd.read_csv(infile, sep="\t", header=None, names=colnames_in, 
                     dtype=str, quoting=3)  # quoting=3 -> quote="" in R

# -------- 2) 数据处理 --------
# 转换数值
df_raw['Diff'] = pd.to_numeric(df_raw['dif.length.bp'], errors='coerce')

# 删除缺失值，得到干净的数据
df = df_raw.dropna(subset=['Diff'])

# 添加方向标签
df['Direction'] = df['Diff'].apply(
    lambda x: 'Positive' if x > 0 else ('Negative' if x < 0 else 'Zero')
)

# 只保留非零的
df = df[df['Direction'] != 'Zero']


# -------- 3) 仅保留 Chr01–Chr10，按顺序排列 --------
chr_levels = [f'Chr{str(i).zfill(2)}' for i in range(1, 11)]
df = df[df['seqid'].isin(chr_levels)]
df['seqid'] = pd.Categorical(df['seqid'], categories=chr_levels, ordered=True)
df['Diff'] = df['Diff'].astype(float)

# 删除缺失值
df = df.dropna(subset=['Diff'])
# 处理 Diff 列
df['Diff'] = pd.to_numeric(df['Diff'], errors='coerce')  # 转为数字，无法转的变 NaN
df = df.dropna(subset=['Diff'])                         # 删除 NaN
df['Diff'] = df['Diff'].astype(float)                  # 强制为标准 float64

p = (ggplot(df, aes(x='gene_start', y='Diff', color='Direction'))
     + geom_hline(yintercept=0, linetype='dashed', size=0.4, alpha=0.6)
     + geom_point(size=1.5, alpha=0.7)
     + facet_wrap('~seqid', ncol=10, scales='free_x')
     + scale_color_manual({'Positive': '#CB5979', 'Negative': '#5494BE'})
     + labs(y='Length Difference (CML304 - B73)')
     + theme_classic(base_size=14)
     + theme(
         axis_title_x=element_blank(),
         axis_text_x=element_blank(),
         axis_ticks_major_x=element_blank(),
         axis_title_y=element_text(face='bold', size=14),
         axis_text_y=element_text(color='black', size=11),
         strip_background=element_rect(fill='white', color='black', size=0.6),
         strip_text=element_text(face='bold', size=12),
         #strip_placement='outside',
         panel_border=element_rect(fill=None, color='black', size=0.6),
         panel_grid_major=element_blank(),
         panel_grid_minor=element_blank(),
         legend_position='none'
     ))

p.save(outfile, width=16, height=5, dpi=600)
