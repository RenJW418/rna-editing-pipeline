# RNA 编辑率分析流程

**语言：中文 | [English](README_EN.md)**

本仓库包含一个用于分析单个双端测序样本 RNA 编辑情况的独立 Bash 流程，以及一个用于统计相邻靶位点双核苷酸组成的 R 脚本。Bash 流程使用一对双端 FASTQ 文件和一条参考序列，并提取指定目标位点的 REDItools2 分析结果；R 脚本进一步从坐标排序后的 BAM 文件中统计同一条 read 上的双核苷酸类型。

## 分析步骤

流程依次执行：

1. 使用 FastQC 对原始 reads 进行质量检查
2. 使用 fastp 进行 reads 过滤
3. 使用 FastQC 对过滤后的 reads 再次进行质量检查
4. 使用 Bowtie2 将双端 reads 比对到参考序列
5. 使用 SAMtools 对 BAM 文件进行排序和索引
6. 使用 REDItools2 统计各位置的碱基组成
7. 从完整结果中提取指定目标位点

## 软件要求

- Bash 4 或更高版本
- FastQC
- fastp
- Bowtie2
- SAMtools
- Python
- REDItools2
- R 4.0 或更高版本（运行双核苷酸分析时需要）
- Bioconductor Rsamtools（运行双核苷酸分析时需要）

除 REDItools2 主脚本外，其余命令需要能够通过系统的 `PATH` 环境变量直接调用。

## 输入文件

该流程需要：

- 一份 R1 FASTQ 文件
- 一份与其配对的 R2 FASTQ 文件
- 一个只包含一条参考序列的 FASTA 文件
- 一个目标位点

目标位点使用从 1 开始的坐标，并按照 FASTA 中参考序列的完整长度计算。

## 使用方法

首先赋予脚本执行权限：

```bash
chmod +x pipeline.sh
```

运行流程：

```bash
REDITOOLS_SCRIPT=/path/to/REDItools2/src/cineca/reditools.py \
  ./pipeline.sh \
  sample_R1.fastq.gz \
  sample_R2.fastq.gz \
  reference.fa \
  132 \
  results \
  sample
```

参数说明：

| 顺序 | 参数 | 说明 |
| --- | --- | --- |
| 1 | `R1.fastq.gz` | R1 FASTQ 文件 |
| 2 | `R2.fastq.gz` | R2 FASTQ 文件 |
| 3 | `reference.fa` | 只包含一条序列的参考 FASTA |
| 4 | `target_position` | 参考序列上的目标位点，从 1 开始计数 |
| 5 | `output_dir` | 输出目录 |
| 6 | `sample_name` | 可选的样本名称；省略时根据 R1 文件名自动生成 |

脚本也可以通过 `bash` 直接运行：

```bash
bash pipeline.sh R1.fq.gz R2.fq.gz reference.fa 132 output sample01
```

## 可选配置

可以通过环境变量调整部分参数：

| 环境变量 | 默认值 | 说明 |
| --- | --- | --- |
| `THREADS` | `6` | 分析使用的线程数 |
| `FASTP_QUALITY` | `30` | fastp 的合格碱基质量阈值 |
| `REDITOOLS_SPLITS` | `10` | REDItools2 拆分 BAM 的数量 |
| `PYTHON_BIN` | `python` | 运行 REDItools2 使用的 Python 命令 |
| `REDITOOLS_SCRIPT` | `/root/software/reditools2.0/src/cineca/reditools.py` | REDItools2 主脚本路径 |

自定义参数示例：

```bash
THREADS=12 \
FASTP_QUALITY=25 \
REDITOOLS_SPLITS=12 \
PYTHON_BIN=python3 \
REDITOOLS_SCRIPT=/opt/REDItools2/src/cineca/reditools.py \
  ./pipeline.sh R1.fq.gz R2.fq.gz reference.fa 132 output sample01
```

## Read 水平双核苷酸分析

`dinucleotide_analysis.R` 对一个目录中的坐标排序 BAM 文件进行逐 read 分析。脚本仅保留同时覆盖预设靶位点处两个相邻参考位置的 reads，按照 CIGAR 字符串将参考坐标映射到 read 序列，并统计 AA、AG、GG 和 GA 等双核苷酸类型。该映射方式能够正确处理软剪切、插入和缺失，避免直接根据比对起点和 read 长度截取序列造成坐标偏移。

安装 Rsamtools：

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Rsamtools")
```

基本用法：

```bash
Rscript dinucleotide_analysis.R \
  --input-dir results/04-bowtie2 \
  --target-start 112 \
  --output dinucleotide_summary.csv
```

`--target-start 112` 默认分析相邻的第 112 和 113 位。所有位置均为参考序列上从 1 开始的坐标。

带自定义参数的示例：

```bash
Rscript dinucleotide_analysis.R \
  --input-dir results/04-bowtie2 \
  --target-start 112 \
  --target-end 113 \
  --categories AA,AG,GG,GA \
  --reference ref1 \
  --min-mapq 20 \
  --yield-size 200000 \
  --exclude-duplicates true \
  --output dinucleotide_summary.csv
```

可用参数：

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--input-dir` | 必填 | BAM 文件所在目录 |
| `--target-start` | 必填 | 相邻靶位点的第一个参考坐标 |
| `--target-end` | `target-start + 1` | 靶位点窗口的最后一个参考坐标 |
| `--output` | `<input-dir>/dinucleotide_summary.csv` | 输出 CSV 文件 |
| `--bam-pattern` | `\.bam$` | 用于筛选 BAM 文件名的正则表达式 |
| `--categories` | `AA,AG,GG,GA` | 纳入统计的序列类型，长度必须与靶位点窗口一致 |
| `--reference` | 不限制 | 仅统计指定参考序列上的比对 |
| `--min-mapq` | `0` | 最低比对质量 |
| `--yield-size` | `100000` | 每批读取的比对记录数，用于控制内存占用 |
| `--include-secondary` | `false` | 是否纳入 secondary alignments |
| `--include-supplementary` | `false` | 是否纳入 supplementary alignments |
| `--include-qcfail` | `false` | 是否纳入未通过测序平台质控的 reads |
| `--exclude-duplicates` | `false` | 是否排除被标记为 duplicate 的 reads |

输出表中每个 BAM 文件对应一行，包含：

- `total_reads`：通过参考序列、MAPQ 和 alignment flag 条件的 reads 数量；
- `spanning_reads`：能够在两个靶位点提取完整序列的 reads 数量；
- `filtered_reads`：属于指定双核苷酸类型的 reads 总数；
- `<类型>_count`：各类型的 reads 数量；
- `<类型>_percent`：该类型占 `filtered_reads` 的百分比。

当 `filtered_reads` 为 0 时，各类型百分比留空。使用 `Rscript dinucleotide_analysis.R --help` 可查看完整帮助。

## 输出目录

```text
output/
├── 01-fastqc/        原始 reads 的 FastQC 报告
├── 02-fastp/         过滤后的 reads 及 fastp HTML、JSON 报告
├── 03-fastqc/        过滤后 reads 的 FastQC 报告
├── 04-bowtie2/       排序并建立索引的 BAM 文件
├── 05-reditools2/    REDItools2 完整结果和目标位点结果
├── logs/             各分析步骤的运行日志
└── ref/              参考序列副本和 Bowtie2 索引
```

主要结果文件：

- `05-reditools2/<样本名>`：REDItools2 生成的完整位点结果表。
- `05-reditools2/<样本名>_target.tsv`：表头以及指定参考序列和目标位点对应的结果。
- `04-bowtie2/<样本名>.sorted.bam`：排序后的比对结果。
- `04-bowtie2/<样本名>.sorted.bam.bai`：BAM 索引文件。

如果 REDItools2 没有报告指定目标位点，目标位点文件将只包含表头，同时流程会输出警告信息。

## 目标位点结果说明

REDItools2 输出中通常包含参考序列名称、位置、参考碱基、覆盖深度、平均质量、各类碱基计数及观察到的碱基变化。具体列名和计算方式可能随 REDItools2 版本变化，应以实际输出表头为准。

目标位点的某种碱基变化比例可根据该替代碱基的 reads 数量和位点总覆盖量计算：

```text
编辑率 = 替代碱基 reads 数量 / 目标位点总覆盖量 × 100%
```

## 方法概述

首先使用 FastQC 评估原始测序数据质量，再通过 fastp 对双端 reads 进行质量过滤。过滤后的 reads 使用 Bowtie2 比对到输入的参考序列，并通过 SAMtools 对比对结果进行坐标排序和索引。随后使用 REDItools2 统计参考序列各位置的碱基组成，最后根据参考序列名称和目标坐标提取目标位点结果。

## 参考文献

- Andrews S. (2010). [FastQC: A Quality Control Tool for High Throughput Sequence Data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- Chen S, Zhou Y, Chen Y, Gu J. (2018). [fastp: an ultra-fast all-in-one FASTQ preprocessor](https://doi.org/10.1093/bioinformatics/bty560). *Bioinformatics* 34: i884–i890.
- Langmead B, Salzberg SL. (2012). [Fast gapped-read alignment with Bowtie 2](https://doi.org/10.1038/nmeth.1923). *Nature Methods* 9: 357–359.
- Li H, Handsaker B, Wysoker A, et al. (2009). [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352). *Bioinformatics* 25: 2078–2079.
- Picardi E, Pesole G. (2013). [REDItools: high-throughput RNA editing detection made easy](https://doi.org/10.1093/bioinformatics/btt287). *Bioinformatics* 29: 1813–1814.

## 许可证

本项目使用 [MIT License](LICENSE)。
