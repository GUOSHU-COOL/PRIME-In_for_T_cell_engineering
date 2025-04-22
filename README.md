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


**Note on Path Configuration**             

If you see a warning when running the script saying that certain software cannot be found, it's likely due to the use of relative paths in the src directory.
To fix this, you can manually replace the relative paths in the script with absolute paths specific to your environment.




## Getting Set Up
Installation
It's recommended (but not essential) to set up a conda environment to manage dependencies
```
conda create -n PEM python=3.8
conda activate PEM        

git clone --recursive https://github.com/GUOSHU-COOL/PRIME-In_for_T_cell_engineering
cd PRIME-In_for_T_cell_engineering

pip install -r requirements.txt
```
## remove redundancy seq avoid multi alignment with transgene
```
cd src
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
python remove_specfic_seq_infa.py
```
## Create blastn db
To generate the BLASTN database, run:
```
python generate_index_add_transgene_fa.py
python generate_index_add_transgene_KI_fa.py     
```

---
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
After the Knock-in Type classify, we analysis sgRNA-caused Off-target type

## Running the sgRNA-caused Off-target analysis  

```python
python sgRNA_analysis.py
```
