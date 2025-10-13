#  IntronDiff-Pipeline
**A one-stop pipeline for comparing and visualizing intron length differences between genomes**

---

##  Overview
`IntronDiff-Pipeline` is an automated workflow designed to **compare intron length differences across maize lines (or other species)**.  
It integrates:
1. **Liftoff-based gene annotation mapping**
2. **Automated intron completion & GFF normalization**
3. **Intron/Exon/CDS statistics**
4. **File merging and delta computation by gene features**
5. **Visualization of intron length differences per chromosome**

---

##  Components
| File | Description |
|:--|:--|
| `change.gff3.add.intron.py` | Fill missing introns and normalize exon/CDS IDs in GFF3. |
| `gff.stat.py` | Summarize gene/transcript/exon/CDS/intron positions and lengths. |
| `merge.file.based.on.keys.py` | Merge reference and target statistics by key columns. |
| `plot_introns_v2.py` | Robust plotting script by chromosome facets. |
| `intron_pipeline.py` | Main pipeline controller, integrates the full process. |
| `run_from_scratch.py` *(optional)* | One-click workflow from annotation to plot. |

---

##  Environment
Python ≥ 3.8  
```bash
pip install pandas plotnine mizani
```
External dependencies: Liftoff, minimap2 (must be in `$PATH`).

---

##  Input Files
| File Type | Description |
|:--|:--|
| Reference genome FASTA (`--ref-fasta`) | e.g., `B73.fa` |
| Reference annotation GFF3 (`--ref-gff`) | Raw annotation; introns are auto-completed. |
| Target genome FASTA (`--target-fasta`) | Target genome sequence (e.g., `Mo17.fa`). |

> If Liftoff results already exist, use `--liftoff-mapped-gff` to skip mapping.

---

##  Usage

### Option 1: From Liftoff to final visualization
```bash
python intron_pipeline.py --sample Mo17 --workdir ./work --run-liftoff   --liftoff-bin liftoff --minimap2-bin minimap2 --threads 16   --ref-fasta /path/to/B73.fa --ref-gff /path/to/B73.gff3   --target-fasta /path/to/Mo17.fa --ref-feature-tsv /path/to/B73.intron.exon.cds.stat.tsv   --plot-script /path/to/plot_introns_v2.py --ylim-pos 40000 --ylim-neg-step 20000
```

### Option 2: Using existing Liftoff results
```bash
python intron_pipeline.py --sample Mo17 --workdir ./work   --liftoff-mapped-gff ./work/Mo17.liftoff.B73.mapped.gff3_polished.gff3   --ref-feature-tsv ./work/B73.ref/B73.intron.exon.cds.stat.tsv   --plot-script ./plot_introns_v2.py
```

### Option 3: Plot only
```bash
python plot_introns_v2.py -i ./work/Mo17.liftoff/Mo17.chr.tsv   -o ./work/Mo17.liftoff/Mo17_Intron_Diff_ByChr_Horizontal_PosNeg.pdf   --ylim_pos 40000 --ylim_neg_step 20000
```

---

##  Output Files
| File | Description |
|:--|:--|
| `<sample>.liftoff.B73.mapped.gff3_polished.gff3` | Liftoff-mapped annotation |
| `<sample>.liftoff.intron.exon.cds.stat.tsv` | Target genome feature summary |
| `<sample>.liftoff.B73.combine.file.tsv` | Combined reference–target table |
| `<sample>.chr.tsv` | Input table for plotting |
| `<sample>_Intron_Diff_ByChr_Horizontal_PosNeg.pdf` | Visualization result |

---

## 🎨 Plot Parameters
| Parameter | Description |
|:--|:--|
| `-i, --input` | Input table (`*.chr.tsv`) |
| `-o, --output` | Output PDF |
| `--ylim_pos` | Positive Y-axis limit |
| `--ylim_neg_step` | Negative Y-axis scaling step |
| `--facet_all` | Include all chromosomes |

---

##  Plot Features
- One facet per chromosome  
- Red: target > B73; Blue: target < B73  
- X-axis: gene start position  
- Y-axis: intron length difference (Target - B73)

---

##  Troubleshooting
| Issue | Cause | Solution |
|:--|:--|:--|
| No PDF output | Outdated script | Update `plot_introns_v2.py`. |
| facet_wrap error | Chromosome format mismatch | Use `--facet-all`. |
| `.chr.tsv` empty | Key mismatch | Ensure introns are completed on both sides. |
| ModuleNotFoundError | Missing dependency | `pip install plotnine mizani`. |

---

## 📂 Suggested Folder Structure
```
intron_diff_pipeline/
├── change.gff3.add.intron.py
├── gff.stat.py
├── merge.file.based.on.keys.py
├── intron_pipeline.py
├── plot_introns_v2.py
└── README.md
```

---

##  Quick Test
```bash
python intron_pipeline.py --sample Mo17 --workdir ./work --run-liftoff   --ref-fasta B73.fa --ref-gff B73.gff3 --target-fasta Mo17.fa   --ref-feature-tsv ./B73.ref/B73.intron.exon.cds.stat.tsv --plot-script ./plot_introns_v2.py
```
Output: `./work/Mo17.liftoff/Mo17_Intron_Diff_ByChr_Horizontal_PosNeg.pdf`
#  IntronMiner
**从基因组比对到 Intron 差异可视化的一站式流程**

---

##  一、功能概述
`IntronDiff-Pipeline` 是一套用于 **比较不同玉米材料（或其他物种）间 intron 长度差异** 的自动化流程。  
它整合了：
1. **Liftoff 基因注释映射**  
2. **Intron 自动补全与 GFF 规范化**  
3. **Intron/Exon/CDS 统计**  
4. **按基因特征键值的文件合并与差值计算**  
5. **差异结果可视化（按染色体分面）**

适用于任意具备参考注释（如 B73）与目标基因组（如 Mo17、CML304 等）的对比分析。

---

##  二、软件组成
| 文件名 | 功能说明 |
|:--|:--|
| **`change.gff3.add.intron.py`** | 为 GFF3 文件补齐 intron，规范 exon/CDS 的 ID（支持多转录本）。 |
| **`gff.stat.py`** | 从 GFF3 统计基因、转录本、exon、CDS、intron 的起止、长度信息，输出 3 个表。 |
| **`merge.file.based.on.keys.py`** | 按指定列（默认第 4 列）合并目标端与参考端的统计结果。 |
| **`plot_introns_v2.py`** | 绘图脚本：按染色体分面绘制 Intron 长度差分布图。 |
| **`intron_pipeline.py`** | 主控脚本：整合全流程，在生成 `<sample>.chr.tsv` 后自动调用 `plot_introns_v2.py` 出图。 |
| **（可选）run_from_scratch.py** | 驱动脚本：可从 B73 注释开始直至绘图。 |

---

##  三、环境依赖

### 1. Python 环境
推荐 Python ≥ 3.8  
必要包：
```bash
pip install pandas plotnine mizani
```

### 2. 外部依赖
- Liftoff  
- minimap2

确保可执行文件在 `$PATH` 或用参数指定。

---

##  四、输入文件

| 文件类型 | 说明 |
|:--|:--|
| **参考基因组 FASTA** (`--ref-fasta`) | 如 `B73.fa` |
| **参考注释 GFF3** (`--ref-gff`) | 原始注释；脚本会自动补 intron。 |
| **目标基因组 FASTA** (`--target-fasta`) | 目标物种序列（如 Mo17.fa） |

> 若已有 liftoff 的映射结果，可直接指定 `--liftoff-mapped-gff` 跳过 liftoff。

---

##  五、使用流程

### 方式 1 ：一条命令从 Liftoff → 出图
```bash
python intron_pipeline.py   --sample Mo17   --workdir ./work   --run-liftoff   --liftoff-bin liftoff   --minimap2-bin minimap2   --threads 16   --ref-fasta /path/to/B73.fa   --ref-gff /path/to/B73.gff3   --target-fasta /path/to/Mo17.fa   --ref-feature-tsv /path/to/B73.intron.exon.cds.stat.tsv   --plot-script /path/to/plot_introns_v2.py   --ylim-pos 40000 --ylim-neg-step 20000
```

**输出路径（示例）**
```
work/
└── Mo17.liftoff/
    ├── Mo17.liftoff.B73.mapped.gff3_polished.gff3
    ├── Mo17.liftoff.intron.exon.cds.stat.tsv
    ├── Mo17.liftoff.B73.combine.file.tsv
    ├── Mo17.chr.tsv
    └── Mo17_Intron_Diff_ByChr_Horizontal_PosNeg.pdf
```

---

### 方式 2 ：已有 liftoff 结果，仅出图
```bash
python intron_pipeline.py   --sample Mo17   --workdir ./work   --liftoff-mapped-gff ./work/Mo17.liftoff.B73.mapped.gff3_polished.gff3   --ref-feature-tsv ./work/B73.ref/B73.intron.exon.cds.stat.tsv   --plot-script ./plot_introns_v2.py
```

---

### 方式 3 ：仅绘图
```bash
python plot_introns_v2.py   -i ./work/Mo17.liftoff/Mo17.chr.tsv   -o ./work/Mo17.liftoff/Mo17_Intron_Diff_ByChr_Horizontal_PosNeg.pdf   --ylim_pos 40000 --ylim_neg_step 20000
```

---

##  六、输出文件说明

| 文件名 | 含义 |
|:--|:--|
| `<sample>.liftoff.B73.mapped.gff3_polished.gff3` | Liftoff 映射结果 |
| `<sample>.liftoff.intron.exon.cds.stat.tsv` | 目标端特征统计表 |
| `<sample>.liftoff.B73.combine.file.tsv` | 合并匹配结果 |
| `<sample>.chr.tsv` | 绘图输入表 |
| `<sample>_Intron_Diff_ByChr_Horizontal_PosNeg.pdf` | 可视化结果 |

---

##  七、绘图参数（plot_introns_v2.py）

| 参数 | 说明 |
|:--|:--|
| `-i, --input` | 输入表 (`*.chr.tsv`) |
| `-o, --output` | 输出 PDF |
| `--ylim_pos` | 正向 y 轴上限 |
| `--ylim_neg_step` | 负向 y 轴下限步长 |
| `--facet_all` | 分面包含所有染色体 |

---

##  八、图像特征
- 每个染色体为一个分面  
- 红色：目标 > B73；蓝色：目标 < B73  
- 横坐标：基因起始位置  
- 纵坐标：intron 长度差 (Target - B73)

---

##  九、常见问题

| 问题 | 原因 | 解决办法 |
|:--|:--|:--|
| 无 PDF 输出 | 绘图脚本未更新 | 替换为新版 plot_introns_v2.py |
| 报错 facet_wrap 无数据 | 染色体列格式不同 | 加 `--facet-all` |
| *.chr.tsv 为空 | 键值不匹配 | 确保双方都补齐 intron |
| ModuleNotFoundError: plotnine | 环境缺包 | `pip install plotnine mizani` |

---

##  十、推荐目录结构
```
IntronMinner/
├── change.gff3.add.intron.py
├── gff.stat.py
├── merge.file.based.on.keys.py
├── intron_pipeline.py
├── plot_introns_v2.py
└── README.md
```

---

##  十一、快速测试
```bash
python intron_pipeline.py   --sample Mo17   --workdir ./work   --run-liftoff   --ref-fasta B73.fa --ref-gff B73.gff3 --target-fasta Mo17.fa   --ref-feature-tsv ./B73.ref/B73.intron.exon.cds.stat.tsv   --plot-script ./plot_introns_v2.py
```

输出文件：`./work/Mo17.liftoff/Mo17_Intron_Diff_ByChr_Horizontal_PosNeg.pdf`
