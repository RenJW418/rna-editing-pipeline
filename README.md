# RNA 编辑率分析流程

**语言：中文 | [English](README_EN.md)**

这是一个用于分析单个双端测序样本 RNA 编辑情况的独立 Bash 流程。每次运行使用一对双端 FASTQ 文件和一条参考序列，并提取指定目标位点的 REDItools2 分析结果。

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
