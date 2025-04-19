#     PRIME-In_for_T_cell_engineering

This code implements our data preprocessing and analysis pipeline for PRIME-In_for_T_cell_engineering data. It takes parameter manifest file (.yaml) with raw sequencing reads (FASTQ) info and and produces a table of annotated gene knock in type and as output.

 
## Dependencies<a name="dependencies"></a>
* Python 3
* Reference genome fasta file ([Example](http://www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta))
* [`bwa`](<http://bio-bwa.sourceforge.net/>) alignment tool
* [`trimmomatic`](<https://github.com/timflutre/trimmomatic/>) a  flexible read trimming tool for Illumina NGS data.
* [`samtools`](<https://github.com/samtools/samtools/>) genome arithmetic utility 
* [`blastn`](<https://github.com/asadprodhan/blastn/blob/main/README.md>)Nucleotide-Nucleotide BLAST
