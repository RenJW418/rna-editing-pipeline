# RNA Editing Pipeline

**Language: [中文](README.md) | English**

This is a standalone Bash pipeline for analyzing RNA editing in one paired-end sequencing sample against one reference sequence. Each run accepts one pair of FASTQ files and one reference sequence, then extracts the REDItools2 result at a specified target position.

## Workflow

The pipeline performs the following steps:

1. Quality control of raw reads with FastQC
2. Read filtering with fastp
3. Quality control of filtered reads with FastQC
4. Paired-end alignment to the reference sequence with Bowtie2
5. BAM sorting and indexing with SAMtools
6. Per-position nucleotide counting with REDItools2
7. Extraction of the specified target position

## Requirements

- Bash 4 or newer
- FastQC
- fastp
- Bowtie2
- SAMtools
- Python
- REDItools2

All required commands except the REDItools2 entry script must be available through the system `PATH`.

## Input files

The pipeline requires:

- One R1 FASTQ file
- One paired R2 FASTQ file
- One FASTA file containing exactly one reference sequence
- One target position

The target position uses 1-based coordinates and is calculated on the complete reference sequence in the FASTA file.

## Usage

Make the script executable:

```bash
chmod +x pipeline.sh
```

Run the pipeline:

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

Arguments:

| Position | Argument | Description |
| --- | --- | --- |
| 1 | `R1.fastq.gz` | R1 FASTQ file |
| 2 | `R2.fastq.gz` | R2 FASTQ file |
| 3 | `reference.fa` | Reference FASTA containing one sequence |
| 4 | `target_position` | 1-based target position on the reference sequence |
| 5 | `output_dir` | Output directory |
| 6 | `sample_name` | Optional sample name; inferred from the R1 filename when omitted |

The script can also be started explicitly with Bash:

```bash
bash pipeline.sh R1.fq.gz R2.fq.gz reference.fa 132 output sample01
```

## Optional configuration

The following environment variables can be used to change pipeline settings:

| Environment variable | Default | Description |
| --- | --- | --- |
| `THREADS` | `6` | Number of processing threads |
| `FASTP_QUALITY` | `30` | fastp qualified-quality threshold |
| `REDITOOLS_SPLITS` | `10` | Number of BAM splits used by REDItools2 |
| `PYTHON_BIN` | `python` | Python executable used to run REDItools2 |
| `REDITOOLS_SCRIPT` | `/root/software/reditools2.0/src/cineca/reditools.py` | Path to the REDItools2 entry script |

Example with custom settings:

```bash
THREADS=12 \
FASTP_QUALITY=25 \
REDITOOLS_SPLITS=12 \
PYTHON_BIN=python3 \
REDITOOLS_SCRIPT=/opt/REDItools2/src/cineca/reditools.py \
  ./pipeline.sh R1.fq.gz R2.fq.gz reference.fa 132 output sample01
```

## Output directory

```text
output/
├── 01-fastqc/        FastQC reports for raw reads
├── 02-fastp/         filtered reads and fastp HTML/JSON reports
├── 03-fastqc/        FastQC reports for filtered reads
├── 04-bowtie2/       sorted and indexed BAM files
├── 05-reditools2/    complete REDItools2 output and target-site result
├── logs/             logs from each analysis step
└── ref/              reference copy and Bowtie2 index
```

Main result files:

- `05-reditools2/<sample>`: complete REDItools2 position table.
- `05-reditools2/<sample>_target.tsv`: the header and row matching the specified reference and target position.
- `04-bowtie2/<sample>.sorted.bam`: sorted alignment file.
- `04-bowtie2/<sample>.sorted.bam.bai`: BAM index file.

If REDItools2 does not report the specified target position, the target file contains only the header and the pipeline prints a warning.

## Target-site result

REDItools2 output usually includes the reference name, position, reference base, coverage, mean quality, nucleotide counts, and observed substitutions. Exact columns and calculations may differ between REDItools2 versions; refer to the header of the generated output.

The percentage of a nucleotide substitution at the target position can be calculated as:

```text
Editing rate = reads containing the alternative base / total target-site coverage × 100%
```

## Methods

Raw sequencing reads are first assessed with FastQC and filtered with fastp. The filtered paired-end reads are aligned to the supplied reference sequence using Bowtie2. SAMtools then coordinate-sorts and indexes the alignment. REDItools2 calculates the nucleotide composition at each reference position, after which the row matching the reference name and target coordinate is extracted.

## References

- Andrews S. (2010). [FastQC: A Quality Control Tool for High Throughput Sequence Data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- Chen S, Zhou Y, Chen Y, Gu J. (2018). [fastp: an ultra-fast all-in-one FASTQ preprocessor](https://doi.org/10.1093/bioinformatics/bty560). *Bioinformatics* 34: i884–i890.
- Langmead B, Salzberg SL. (2012). [Fast gapped-read alignment with Bowtie 2](https://doi.org/10.1038/nmeth.1923). *Nature Methods* 9: 357–359.
- Li H, Handsaker B, Wysoker A, et al. (2009). [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352). *Bioinformatics* 25: 2078–2079.
- Picardi E, Pesole G. (2013). [REDItools: high-throughput RNA editing detection made easy](https://doi.org/10.1093/bioinformatics/btt287). *Bioinformatics* 29: 1813–1814.

## License

This project is released under the [MIT License](LICENSE).
