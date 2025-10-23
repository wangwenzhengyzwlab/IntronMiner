#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fusion_collect.py

主流工具（Arriba/STAR-Fusion/FusionCatcher/JAFFA）输出整合

useage：
  python fusion_collect.py \
    --inputs arriba.tsv star-fusion.fusion_predictions.abridged.tsv fusioncatcher.txt jaffa.csv \
    -o sample.fusion.external.tsv

  合并去冗余
  python fusion_collect_external.py \
    --merge sample.fusion.external.tsv -o sample.fusion.external.merged.tsv
"""

import argparse, csv, re, sys
from collections import defaultdict, namedtuple

Row = namedtuple("Row", [
    "geneA","geneB","chrA","posA","strandA","chrB","posB","strandB",
    "support_split","support_spanning","tool","extra"
])

def try_int(x, default=0):
    try:
        return int(x)
    except:
        try:
            return int(float(x))
        except:
            return default

def parse_arriba(path):
    rows=[]
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for d in r:
            g1 = d.get("gene1", d.get("Gene1",""))
            g2 = d.get("gene2", d.get("Gene2",""))
            b1 = d.get("breakpoint1","")
            b2 = d.get("breakpoint2","")
            chrA,posA = (b1.split(":")+[".","0"])[:2]
            chrB,posB = (b2.split(":")+[".","0"])[:2]
            sp  = try_int(d.get("split_reads", d.get("splitreads",0)))
            spn = try_int(d.get("discordant_mates", d.get("spanning_pairs",0)))
            sA = d.get("strand1","."); sB = d.get("strand2",".")
            rows.append(Row(g1,g2,chrA,int(posA),sA,chrB,int(posB),sB,sp,spn,"Arriba",""))
    return rows

def parse_star_fusion(path):
    rows=[]
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for d in r:
            # Arriba列名：#FusionName, LeftGene, RightGene, LeftBreakpoint, RightBreakpoint, JunctionReadCount, SpanningFragCount
            g1 = d.get("LeftGene","").split("^")[0]
            g2 = d.get("RightGene","").split("^")[0]
            b1 = d.get("LeftBreakpoint","")
            b2 = d.get("RightBreakpoint","")
            chrA,posA = (b1.split(":")+[".","0"])[:2]
            chrB,posB = (b2.split(":")+[".","0"])[:2]
            sp  = try_int(d.get("JunctionReadCount",0))
            spn = try_int(d.get("SpanningFragCount",0))
            rows.append(Row(g1,g2,chrA,int(posA),".",chrB,int(posB),".",sp,spn,"STAR-Fusion",""))
    return rows

def parse_fusioncatcher(path):
    rows=[]
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): 
                continue
            parts = re.split(r"\t|,", line.strip())
            if len(parts) < 5:
                continue
            # 粗配
            g1,g2 = parts[0],parts[1]
            b1 = next((p for p in parts if re.match(r"^\w+:\d+$", p)), None)
            b2 = None
            if b1:
                for p in parts:
                    if p!=b1 and re.match(r"^\w+:\d+$", p):
                        b2 = p; break
            if not (b1 and b2): 
                continue
            chrA,posA = b1.split(":")
            chrB,posB = b2.split(":")
            rows.append(Row(g1,g2,chrA,int(posA),".",chrB,int(posB),".",0,0,"FusionCatcher",""))
    return rows

def parse_jaffa(path):
    rows=[]
    with open(path) as f:
        head = f.readline()
        delim = "," if (head.count(",")>head.count("\t")) else "\t"
        f.seek(0)
        r = csv.DictReader(f, delimiter=delim)
        for d in r:
            # 常见列名：gene1, gene2, chr1, pos1, chr2, pos2, spanning, split
            g1 = d.get("gene1","")
            g2 = d.get("gene2","")
            chrA = d.get("chr1", d.get("chrom1",".")); chrB = d.get("chr2", d.get("chrom2","."))
            posA = try_int(d.get("pos1", d.get("break1",0)))
            posB = try_int(d.get("pos2", d.get("break2",0)))
            spn  = try_int(d.get("spanning",0)); sp = try_int(d.get("split",0))
            rows.append(Row(g1,g2,chrA,posA,".",chrB,posB,".",sp,spn,"JAFFA",""))
    return rows

def write_rows(rows, out_tsv):
    with open(out_tsv,"w") as w:
        w.write("\t".join(Row._fields)+"\n")
        for r in rows:
            w.write("\t".join(map(str,r))+"\n")

def merge_rows(in_tsv, out_tsv, cluster_win=10):
    """
    combine&vote
    """
    from collections import defaultdict
    buckets = defaultdict(lambda: {"sp":0,"spn":0,"tools":set(),"members":0})
    for line in open(in_tsv):
        if line.startswith("geneA"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 12: continue
        gA,gB,chrA,posA,sA,chrB,posB,sB,sp,spn,tool,extra = c[:12]
        posA,posB = int(posA), int(posB)
        key = (gA,gB,chrA,round(posA/cluster_win),chrB,round(posB/cluster_win),sA,sB)
        buckets[key]["sp"] = max(buckets[key]["sp"], try_int(sp))
        buckets[key]["spn"] = max(buckets[key]["spn"], try_int(spn))
        buckets[key]["tools"].update(tool.split(","))
        buckets[key]["members"] += 1
    with open(out_tsv,"w") as w:
        w.write("\t".join(["geneA","geneB","chrA","posA_bin","chrB","posB_bin","strandA","strandB","max_split","max_spanning","tools_support","n_items"])+"\n")
        for key, val in buckets.items():
            gA,gB,chrA,bA,chrB,bB,sA,sB = key
            w.write("\t".join(map(str,[
                gA,gB,chrA,bA,chrB,bB,sA,sB,val["sp"],val["spn"],",".join(sorted(val["tools"])),val["members"]
            ]))+"\n")

def detect_tool_type(path):
    p = path.lower()
    if p.endswith(".abridged.tsv") or "star-fusion" in p:
        return "starfusion"
    if p.endswith(".tsv"): # 多数 Arriba
        
        with open(path) as f:
            head = f.readline().lower()
            if "gene1" in head and "breakpoint1" in head:
                return "arriba"
    if "fusioncatcher" in p or p.endswith(".txt") or p.endswith(".summary"):
        return "fusioncatcher"
    if p.endswith(".csv") or "jaffa" in p:
        return "jaffa"
    return "unknown"

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--inputs", nargs="+", help="外部工具的结果文件（Arriba/STAR-Fusion/FusionCatcher/JAFFA）")
    ap.add_argument("--merge", help="上一步统一表（*.fusion.external.tsv），做聚类合并")
    ap.add_argument("-o","--out", required=True, help="输出文件")
    ap.add_argument("--cluster-win", type=int, default=10)
    args = ap.parse_args()

    if args.inputs and not args.merge:
        rows=[]
        for p in args.inputs:
            tp = detect_tool_type(p)
            if tp=="arriba": rows += parse_arriba(p)
            elif tp=="starfusion": rows += parse_star_fusion(p)
            elif tp=="fusioncatcher": rows += parse_fusioncatcher(p)
            elif tp=="jaffa": rows += parse_jaffa(p)
            else:
                print(f"[WARN] 无法识别工具类型，跳过：{p}", file=sys.stderr)
        write_rows(rows, args.out)
    elif args.merge and not args.inputs:
        merge_rows(args.merge, args.out, args.cluster_win)
    else:
        ap.print_help(); sys.exit(1)

if __name__ == "__main__":
    main()

