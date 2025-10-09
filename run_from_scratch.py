#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, subprocess, sys
from pathlib import Path

here = Path(__file__).parent.resolve()

def run(cmd, cwd=None):
    print("[run]", cmd)
    ret = subprocess.run(cmd, shell=True, cwd=cwd)
    if ret.returncode != 0:
        sys.exit(ret.returncode)

def main():
    ap = argparse.ArgumentParser(description="Run intron pipeline from scratch: build B73 stats + liftoff + plot")
    ap.add_argument("--workdir", required=True, help="工作目录（会在里面生成中间与结果文件）")
    ap.add_argument("--sample",   required=True, help="目标样本名（文件前缀）")
    # Reference inputs
    ap.add_argument("--ref-fasta", required=True, help="B73 参考基因组 FASTA")
    ap.add_argument("--ref-gff",   required=True, help="B73 参考 GFF3（原始，尚未补 intron）")
    # Target input
    ap.add_argument("--target-fasta", required=True, help="目标物种基因组 FASTA")
    # Tools
    ap.add_argument("--liftoff-bin", default="liftoff", help="liftoff 路径")
    ap.add_argument("--minimap2-bin", default="minimap2", help="minimap2 路径")
    ap.add_argument("--threads", type=int, default=8, help="liftoff 线程")
    ap.add_argument("--skip-plot", action="store_true", help="仅生成 TSV，不出图")
    args = ap.parse_args()

    work = Path(args.workdir).expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)

    # === 1) 为 B73 注释补 intron 并统计，得到 B73.intron.exon.cds.stat.tsv ===
    b73_dir = work / "B73.ref"
    b73_dir.mkdir(exist_ok=True, parents=True)
    b73_with_intron = b73_dir / "B73.with_intron.gff3"
    run(f"python {here/'change.gff3.add.intron.py'} -i {args.ref_gff} -o {b73_with_intron}")
    # 统计会生成三份，目标表名为 *.intron.exon.cds.stat.tsv
    run(f"python {here/'gff.stat.py'} -g {b73_with_intron} -p {b73_dir/'B73'}")
    b73_feature_tsv = b73_dir / "B73.intron.exon.cds.stat.tsv"

    # === 2) 运行主流水线（含 liftoff + 合并 + 差值 + 作图）===
    # intron_pipeline.py 已支持 --run-liftoff / --ref-feature-tsv / --target-fasta 等参数
    cmd = (
        f"python {here/'intron_pipeline.py'} "
        f"--sample {args.sample} --workdir {work} "
        f"--run-liftoff --liftoff-bin {args.liftoff_bin} --minimap2-bin {args.minimap2_bin} "
        f"--target-fasta {args.target_fasta} --ref-fasta {args.ref_fasta} --ref-gff {b73_with_intron} "
        f"--ref-feature-tsv {b73_feature_tsv} " + ("--skip-plot" if args.skip_plot else "")
    )
    run(cmd)

    print("\n[ALL DONE]")
    print(f"- B73 feature TSV : {b73_feature_tsv}")
    out_pdf = (Path(args.workdir)/f"{args.sample}.liftoff"/f"{args.sample}_Intron_Diff_ByChr_Horizontal_PosNeg.pdf")
    print(f"- Plot (if enabled): {out_pdf}")

if __name__ == "__main__":
    main()

