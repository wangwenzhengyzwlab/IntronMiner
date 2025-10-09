#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from plotnine import (
    ggplot, aes, geom_point, geom_hline, facet_wrap, scale_color_manual,
    scale_y_continuous, labs, theme_classic, theme, element_blank,
    element_text, element_rect
)
from mizani.formatters import number_format

ALIASES = {
    "seq": ["seqid", "Chromosome", "chromosome", "chr", "chr_id", "scaffold", "contig"],
    "diff": ["dif.length.bp", "Diff", "diff", "dif", "dif_length_bp", "delta_len", "delta"],
    "start": ["gene_start", "start", "start_bp", "begin", "pos_start"]
}

CHR_LEVELS = [f"Chr{str(i).zfill(2)}" for i in range(1, 11)]

def choose_col(cols, prefer):
    cols_lower = [c.lower() for c in cols]
    for cand in ALIASES[prefer]:
        if cand.lower() in cols_lower:
            return cols[cols_lower.index(cand.lower())]
    return None

def normalize_chr(x):
    # 将 1/01/chr1/Chr1 -> Chr01；其余值保持原样
    v = str(x)
    v2 = v
    if v2.lower().startswith("chr"):
        v2 = v2[3:]
    if v2.isdigit():
        n = int(v2)
        if 1 <= n <= 10:
            return f"Chr{n:02d}"
    # 已经是 Chr01..Chr10
    if v.startswith("Chr") and len(v) == 5 and v[3:].isdigit():
        return v
    return v

def parse_args():
    ap = argparse.ArgumentParser(description="Plot intron length differences by chromosome (robust column detection).")
    ap.add_argument("-i", "--input", required=True, help="Input TSV (the *.chr.tsv)")
    ap.add_argument("-o", "--output", required=True, help="Output PDF path")
    ap.add_argument("--seq-col", help="Name of sequence/chr column (optional)")
    ap.add_argument("--diff-col", help="Name of difference column (optional)")
    ap.add_argument("--start-col", help="Name of gene start column (optional)")
    ap.add_argument("--ylim_pos", type=float, default=40000, help="Upper limit for positive differences")
    ap.add_argument("--ylim_neg_step", type=int, default=20000, help="Step for flooring negative limit")
    ap.add_argument("--facet_all", action="store_true",
                    help="Facet by ALL unique seq values instead of forcing Chr01..Chr10")
    return ap.parse_args()

def main():
    args = parse_args()

    # 读入：尝试 header=0；若第一行不像表头，改用 header=None 并赋标准名
    try:
        df = pd.read_csv(args.input, sep="\t", header=0, dtype=str, quoting=3)
    except Exception:
        df = pd.read_csv(args.input, sep="\t", header=None, dtype=str, quoting=3)

    # 列名选择
    seq_col = args.seq_col or choose_col(list(df.columns), "seq")
    diff_col = args.diff_col or choose_col(list(df.columns), "diff")
    start_col = args.start_col or choose_col(list(df.columns), "start")

    missing = [name for name, col in [("seq", seq_col), ("diff", diff_col), ("start", start_col)] if col is None]
    if missing:
        raise SystemExit(
            f"[ERROR] Cannot find columns for: {', '.join(missing)}.\n"
            f"  Available columns: {list(df.columns)}\n"
            f"  Try specifying --seq-col/--diff-col/--start-col explicitly."
        )

    # 转数值 & 规范化
    df["Diff"] = pd.to_numeric(df[diff_col], errors="coerce")
    df["gene_start_plot"] = pd.to_numeric(df[start_col], errors="coerce")
    df["seq_plot"] = df[seq_col].map(normalize_chr)

    # 丢掉 Diff 为 NA 的
    df = df.dropna(subset=["Diff"]).copy()

    # 方向
    df["Direction"] = df["Diff"].apply(lambda x: "Positive" if x > 0 else ("Negative" if x < 0 else "Zero"))
    df = df[df["Direction"] != "Zero"].copy()

    # 染色体集合
    if args.facet_all:
        # 全部分面（顺序按自然排序）
        cats = sorted(df["seq_plot"].dropna().unique().tolist())
    else:
        cats = CHR_LEVELS

    df = df[df["seq_plot"].isin(cats)].copy()
    if df.empty:
        # 打印可用值，帮助定位
        uniq = df.assign(_tmp="x")
        # 如果过滤后为空，uniq 也会空；改用原数据的唯一值提示
        uniq_vals = sorted(pd.Series([normalize_chr(v) for v in (pd.unique(df.assign(seq=df["seq_plot"])["seq_plot"]) if not df.empty else pd.unique(df.assign(seq=df[seq_col].map(normalize_chr))["seq"]))]).dropna().astype(str).unique())
        raise SystemExit(
            "[ERROR] No rows left after filtering by chromosomes.\n"
            f"  Expected: {cats}\n"
            f"  Found (examples from file): {uniq_vals[:20]}"
        )

    df["seq_plot"] = pd.Categorical(df["seq_plot"], categories=cats, ordered=True)

    # 负向轴下限
    if (df["Diff"] < 0).any():
        min_neg = df.loc[df["Diff"] < 0, "Diff"].min()
        floor_neg = int((min_neg // args.ylim_neg_step) * args.ylim_neg_step)
    else:
        floor_neg = -args.ylim_neg_step

    p = (ggplot(df, aes(x="gene_start_plot", y="Diff", color="Direction"))
         + geom_hline(yintercept=0, linetype="dashed", size=0.4, alpha=0.6)
         + geom_point(size=1.0, alpha=0.6)
         + facet_wrap("~seq_plot", ncol=min(10, len(cats)), scales="free_x")
         + scale_color_manual({"Positive": "#CB5979", "Negative": "#5494BE"})
         + scale_y_continuous(limits=(floor_neg, args.ylim_pos), labels=number_format(accuracy=1))
         + labs(y="Length Difference (Target - B73)")
         + theme_classic(base_size=13)
         + theme(
             axis_title_x=element_blank(),
             axis_text_x=element_blank(),
             axis_ticks_major_x=element_blank(),
             axis_title_y=element_text(face='bold', size=13),
             axis_text_y=element_text(color='black', size=10),
             strip_background=element_rect(fill='white', color='black', size=0.6),
             strip_text=element_text(face='bold', size=11),
             panel_border=element_rect(fill=None, color='black', size=0.6),
             panel_grid_major=element_blank(),
             panel_grid_minor=element_blank(),
             legend_position='none'
         ))

    p.save(args.output, width=16, height=5, dpi=600)
    print(f"[OK] Plot saved -> {args.output}")

if __name__ == "__main__":
    main()

