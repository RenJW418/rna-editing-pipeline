# RNA Editing Pipeline

A standalone Bash pipeline for RNA editing analysis of one paired-end NGS sample against one reference sequence.

The workflow runs:

1. FastQC on raw reads
2. fastp read filtering
3. FastQC on filtered reads
4. Bowtie2 paired-end alignment
5. SAMtools sorting and indexing
6. REDItools2 nucleotide counting
7. Extraction of the requested target position

## Requirements

- Bash 4 or newer
- FastQC
- fastp
- Bowtie2
- SAMtools
- Python
- REDItools2

All required commands except the REDItools2 script must be available in `PATH`.

## Input

The script accepts one paired-end FASTQ sample, one single-record reference FASTA, and one 1-based target position. The reference FASTA must contain exactly one sequence.

## Usage

```bash
chmod +x pipeline.sh

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
| 1 | `R1.fastq.gz` | Forward-read FASTQ |
| 2 | `R2.fastq.gz` | Reverse-read FASTQ |
| 3 | `reference.fa` | Single-record reference FASTA |
| 4 | `target_position` | 1-based position on the complete reference sequence |
| 5 | `output_dir` | Output directory |
| 6 | `sample_name` | Optional sample name; inferred from R1 when omitted |

## Configuration

| Variable | Default | Description |
| --- | --- | --- |
| `THREADS` | `6` | Number of processing threads |
| `FASTP_QUALITY` | `30` | fastp qualified-quality threshold |
| `REDITOOLS_SPLITS` | `10` | Number of BAM splits used by REDItools2 |
| `PYTHON_BIN` | `python` | Python executable used for REDItools2 |
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

## Output

```text
output/
├── 01-fastqc/        FastQC reports for raw reads
├── 02-fastp/         filtered reads and fastp HTML/JSON reports
├── 03-fastqc/        FastQC reports for filtered reads
├── 04-bowtie2/       coordinate-sorted and indexed BAM files
├── 05-reditools2/    complete REDItools2 output and target-site table
├── logs/             logs from every analysis step
└── ref/              copied reference and Bowtie2 index
```

The main results are:

- `05-reditools2/<sample>`: complete REDItools2 position table.
- `05-reditools2/<sample>_target.tsv`: header plus the row matching the requested reference ID and position.
- `04-bowtie2/<sample>.sorted.bam`: sorted alignment file.

If REDItools2 does not report the requested position, the target file contains only the header and the pipeline prints a warning.

## Methods

Raw reads are assessed with FastQC and filtered with fastp. Filtered paired-end reads are aligned to the supplied reference sequence using Bowtie2. The alignment is coordinate-sorted and indexed with SAMtools. REDItools2 then produces per-position nucleotide counts, from which the requested target-position row is extracted.

## References

- Andrews S. (2010). [FastQC: A Quality Control Tool for High Throughput Sequence Data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
- Chen S, Zhou Y, Chen Y, Gu J. (2018). [fastp: an ultra-fast all-in-one FASTQ preprocessor](https://doi.org/10.1093/bioinformatics/bty560). *Bioinformatics* 34: i884–i890.
- Langmead B, Salzberg SL. (2012). [Fast gapped-read alignment with Bowtie 2](https://doi.org/10.1038/nmeth.1923). *Nature Methods* 9: 357–359.
- Li H, Handsaker B, Wysoker A, et al. (2009). [The Sequence Alignment/Map format and SAMtools](https://doi.org/10.1093/bioinformatics/btp352). *Bioinformatics* 25: 2078–2079.
- Picardi E, Pesole G. (2013). [REDItools: high-throughput RNA editing detection made easy](https://doi.org/10.1093/bioinformatics/btt287). *Bioinformatics* 29: 1813–1814.

## License

This project is released under the [MIT License](LICENSE).
