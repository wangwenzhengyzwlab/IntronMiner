
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
intron_pipeline.py  (relative-path edition)

将所有子脚本以“同目录相对路径”调用：
- change.gff3.add.intron(2).py
- gff.stat(2).py
- merge.file.based.on.keys(2).py
- plot_introns_v2.py

放置方式：把以上 5 个文件放到同一个目录下，直接运行本脚本即可。
"""

import argparse
import os
import sys
import shutil
import subprocess
from pathlib import Path

# -------- Helpers --------
def run_sh(cmd, cwd=None):
    print(f"[run/sh] {cmd}")
    result = subprocess.run(cmd, cwd=cwd, shell=True)
    if result.returncode != 0:
        sys.exit(result.returncode)

def run_cmd(args_list, cwd=None):
    print(f"[run/cmd] {' '.join(map(str, args_list))}")
    result = subprocess.run(args_list, cwd=cwd)
    if result.returncode != 0:
        sys.exit(result.returncode)

def which(bin_name):
    return shutil.which(bin_name)

def parse_args():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--sample", required=True, help="样本名称（文件前缀）")
    ap.add_argument("--workdir", required=True, help="工作目录（将创建 {sample}.liftoff 子目录）")
    # Liftoff 相关
    ap.add_argument("--run-liftoff", action="store_true", help="是否执行 liftoff（需要 liftoff & minimap2）")
    ap.add_argument("--liftoff-bin", default="liftoff", help="liftoff 可执行程序名或路径")
    ap.add_argument("--minimap2-bin", default="minimap2", help="minimap2 可执行程序名或路径")
    ap.add_argument("--target-fasta", help="目标物种的基因组 FASTA（liftoff 目标）")
    ap.add_argument("--ref-fasta", help="B73 参考基因组 FASTA （liftoff 参考）")
    ap.add_argument("--ref-gff", help="B73 参考注释（含 intron 的 gff3），liftoff -g 输入")
    ap.add_argument("--liftoff-mapped-gff", help="liftoff 输出的 *.mapped.gff3_polished 或等价文件（作为 change.gff3.add.intron 的输入）")
    # 合并与作图
    ap.add_argument("--ref-feature-tsv", dest="ref_feature_tsv", required=True, help="B73 端 intron/exon/cds 特征统计 TSV（合并参考文件）")
    ap.add_argument("--skip-plot", action="store_true", help="仅生成 TSV，不绘图")
    ap.add_argument("--threads", type=int, default=8, help="liftoff 线程数（仅在 --run-liftoff 生效）")
    return ap.parse_args()

def main():
    args = parse_args()

    script_dir = Path(__file__).parent.resolve()

    # 工作路径
    workdir = Path(args.workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    liftoff_dir = workdir / f"{args.sample}.liftoff"
    liftoff_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(liftoff_dir)

    # 1) Liftoff（可选）
    mapped_polished = None
    if args.run_liftoff:
        if not (args.target_fasta and args.ref_fasta and args.ref_gff):
            sys.exit("运行 liftoff 需要 --target-fasta/--ref-fasta/--ref-gff")
        if not which(args.liftoff_bin):
            sys.exit(f"找不到 liftoff 可执行文件：{args.liftoff_bin}")
        if not which(args.minimap2_bin):
            sys.exit(f"找不到 minimap2 可执行文件：{args.minimap2_bin}")

        mapped_base = f"{args.sample}.liftoff.B73.mapped.gff3"
        mapped = liftoff_dir / mapped_base
        unmapped = liftoff_dir / f"{args.sample}.liftoff.B73.unmapped.gff3"

        cmd = (
            f"{args.liftoff_bin} -copies -p {args.threads} "
            f"-g {args.ref_gff} -m {args.minimap2_bin} -polish -cds "
            f"-o {mapped} -u {unmapped} {args.target_fasta} {args.ref_fasta}"
        )
        run_sh(cmd)

        mapped_polished = str(mapped) + "_polished"
        if not Path(mapped_polished).exists():
            mapped_polished = str(mapped)
    else:
        if not args.liftoff_mapped_gff:
            sys.exit("未运行 liftoff 时，必须提供 --liftoff-mapped-gff")
        mapped_polished = args.liftoff_mapped_gff
        if not Path(mapped_polished).exists():
            sys.exit(f"找不到 liftoff 映射注释：{mapped_polished}")

    # 2) 添加 intron
    gff_with_intron = f"{args.sample}.liftoff.B73.mapped.gff3_polished.gff3"
    run_cmd(["python", str(script_dir / "change.gff3.add.intron(2).py"),
             "-i", mapped_polished, "-o", gff_with_intron])

    # 3) 统计
    run_cmd(["python", str(script_dir / "gff.stat(2).py"),
             "-g", gff_with_intron, "-p", f"{args.sample}.liftoff"])

    feature_stat = f"{args.sample}.liftoff.intron.exon.cds.stat.tsv"
    if not Path(feature_stat).exists():
        sys.exit(f"未找到特征统计文件：{feature_stat}")

    # 4) 重写 4/6/7 列
    change_tsv = f"{args.sample}.liftoff.intron.exon.cds.stat.change.tsv"
    with open(feature_stat, "r") as fin, open(change_tsv, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        fout.write("\t".join(header) + "\n")
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            def get_or_blank(i):
                return cols[i] if i < len(cols) else ""

            f4 = get_or_blank(3)
            parts4 = f4.split("_")
            f4_new = f"{parts4[0]}_{parts4[1]}_{parts4[-1]}" if len(parts4) >= 3 else f4

            f6 = get_or_blank(5)
            f6_new = f6.split("_")[0] if f6 else f6

            f7 = get_or_blank(6)
            p7 = f7.split("_")
            f7_new = f"{p7[0]}_{p7[1]}" if len(p7) >= 2 else f7

            cols[3] = f4_new
            if len(cols) >= 6:
                cols[5] = f6_new
            if len(cols) >= 7:
                cols[6] = f7_new
            fout.write("\t".join(cols) + "\n")

    # 5) 合并（以第4列为 key）
    run_cmd([
        "python", str(script_dir / "merge.file.based.on.keys(2).py"),
        "-rf", change_tsv, "-rc", "4",
        "-qf", args.ref_feature_tsv, "-qc", "4",
        "-pf", f"{args.sample}.liftoff.B73"
    ])

    combined_tsv = f"{args.sample}.liftoff.B73.combine.file.tsv"
    if not Path(combined_tsv).exists():
        sys.exit("未找到合并输出 .combine.file.tsv")

    # 6) 筛选等位内含子并计算长度差
    equal_intron = f"{args.sample}.liftoff.B73.combine.equal.intron.dif.tsv"
    with open(combined_tsv, "r") as fin, open(equal_intron, "w") as fout:
        for raw in fin:
            row = raw.rstrip("\n").split("\t")
            if len(row) < 14:
                continue
            if row[4] != "intron":
                continue
            picked = row[0:7] + row[9:14] + row[20:]
            if len(picked) < 10:
                continue
            if picked[7] == picked[-2]:
                try:
                    length = int(picked[-1]) - int(picked[8])
                except ValueError:
                    continue
                fout.write("\t".join(picked + [str(length)]) + "\n")

    # 7) 规范化到作图输入并绘图
    in_plot = f"{args.sample}.chr.tsv"
    with open(in_plot, "w") as fout, open(equal_intron, "r") as fin:
        fout.write("\t".join([
            "seqid","gene_start","gene_end","mRNA_id","type","gene_id",
            "mRNA_id_dup","exon_number","length.bp",
            "ref_chr_id","ref_gene_start","ref_gene_end","ref_exon_number","ref_length.bp",
            "dif.length.bp"
        ]) + "\n")
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 15:
                cols = cols + [""]*(15-len(cols))
            row = [
                cols[0], cols[1], cols[2],
                cols[3], cols[4], cols[5],
                cols[6], cols[7],
                cols[9] if len(cols)>9 else "",
                cols[10] if len(cols)>10 else "",
                cols[11] if len(cols)>11 else "",
                cols[12] if len(cols)>12 else "",
                cols[-3] if len(cols)>3 else "",
                cols[-2] if len(cols)>2 else "",
                cols[-1],
            ]
            fout.write("\t".join(row) + "\n")

    if not args.skip_plot:
        pdf = f"{args.sample}_Intron_Diff_ByChr_Horizontal_PosNeg.pdf"
        run_cmd(["python", str(script_dir / "plot_introns_v2.py"),
                 "-i", in_plot, "-o", pdf])
        print(f"[OK] Plot saved -> {pdf}")
    print("[DONE]")

if __name__ == "__main__":
    main()
