# PRIME-In_for_T_cell_engineering

This code implements our data preprocessing and analysis pipeline for PRIME-In_for_T_cell_engineering data. It takes a parameter manifest file (`.yaml`) with raw sequencing reads (FASTQ) info and produces a table of annotated gene knock-in types as output.

## Dependencies<a name="dependencies"></a>

To ensure successful execution of the pipeline, the following dependencies are required:

- **Python 3** — Programming language for executing the pipeline scripts  
- **Reference genome FASTA file**  ([Example](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz)) — Used for alignment and annotation  
- [`bwa`](http://bio-bwa.sourceforge.net/) — Alignment tool for mapping sequencing reads to the reference genome  
- [`trimmomatic`](https://github.com/timflutre/trimmomatic/) — Flexible read trimming tool for Illumina NGS data  
- [`samtools`](https://github.com/samtools/samtools/) — Tool for manipulating alignments in the SAM/BAM format  
- [`bedtools`](https://bedtools.readthedocs.io/en/latest/) — Genome arithmetic utility  
- [`blastn`](https://github.com/asadprodhan/blastn/blob/main/README.md) — Nucleotide-Nucleotide BLAST  

## Running the Knock-in Type Classifier

To run the full knock-in type classification pipeline, you must first create a manifest YAML file that describes all pipeline inputs, including a `data_info.xlsx` file which contains sample-specific metadata.

Once your files are ready, simply run:

```
python batch_reads_type_classifier.py --config ./config.yaml
```


An example `config.yaml` and `data_info.xlsx` can be found in the `src` directory.

---

## Required Columns in `data_info.xlsx`

The `data_info.xlsx` file should contain the following columns for each sample:

| Column Name              | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| `sample`                 | Sample name / identifier                                                    |
| `chr_seq`                | Chromosome or sequence name to align against (e.g., `chr1`, `chr7`)         |
| `cut_site`               | Genomic coordinate (1-based) of the CRISPR cut site                         |
| `strand`                 | Strand direction: `+` or `-`                                                |
| `primer_plus_5bp`        | Primer sequence including 5 bp upstream flanking sequence                   |
| `primer_to_cut_plus20bp` | Subsequence from primer to cut site, plus 20 bp downstream from the cut     |

Make sure all values are correctly formatted for each sample; otherwise, the pipeline may raise input parsing errors.

---

