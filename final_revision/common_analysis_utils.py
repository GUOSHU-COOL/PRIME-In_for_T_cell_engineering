import os
import pandas as pd
from typing import Optional, Tuple
from Bio import SeqIO
import pysam


def flash_stitch(cm_dir: str, sample: str) -> None:
    """
    Use the FLASH tool to stitch paired-end sequencing files into single-end.

    Args:
        cm_dir (str): Main directory path.
        sample (str): Sample name.
        clean_data_dir (str): Subdirectory path of cleaned data.
    """
    output_dir = os.path.join(cm_dir, sample, 'flash/')
    os.makedirs(output_dir, exist_ok=True)

    fq_r1 = os.path.join(cm_dir, f"{sample}_R1.fq.gz")
    fq_r2 = os.path.join(cm_dir, f"{sample}_R2.fq.gz")

    cmd = f"/home/wangxin/miniconda3_1/envs/PEM-Q/bin/flash {fq_r1} {fq_r2} -o {sample} -t 16 -d {output_dir} > flash.log"
    print(cmd)
    os.system(cmd)

    f1 = os.path.join(output_dir, f"{sample}.extendedFrags.fastq")
    f2 = os.path.join(output_dir, f"{sample}.notCombined_1.fastq")
    merge_f = os.path.join(output_dir, f"{sample}.merge.fastq")

    cmd = f"cat {f1} {f2} > {merge_f}"
    print(cmd)
    os.system(cmd)

    print("Merging fastq files... Done!")


def create_db(fa: str, db: str) -> None:
    """
    Create a BLASTN nucleotide database.

    Args:
        fa (str): Input FASTA sequence file.
        db (str): Output database file prefix.
    """
    cmd = f"makeblastdb -dbtype nucl -in {fa} -input_type fasta -parse_seqids -out {db}"
    os.system(cmd)
    print('Create blastn ref database, done!')


def blastn(test_fa: str, db: str, blastn_result_txt: str) -> None:
    """
    Perform BLASTN comparison on test sequences and output results.

    Args:
        test_fa (str): Query sequence file in FASTA format.
        db (str): Established BLAST database name.
        blastn_result_txt (str): Output BLASTN results file path.
    """
    cmd = f"blastn -query {test_fa} -db {db} -out {blastn_result_txt} -evalue 0.05 -task blastn -num_threads 25 -min_raw_gapped_score 20 -outfmt 6"
    os.system(cmd)
    print('Blastn analysis, done!')


def get_sample_info_from_excel(excel_file: str, sample_name: str) -> Optional[pd.Series]:
    """
    Extract sample information from an Excel file.

    Args:
        excel_file (str): Path to the Excel file.
        sample_name (str): Sample name.

    Returns:
        Optional[pd.Series]: Return the row data if found, otherwise return None.
    """
    df = pd.read_excel(excel_file, dtype=str)
    sample_row = df[df['sample'] == sample_name]
    if not sample_row.empty:
        return sample_row.iloc[0]
    print(f"Sample {sample_name} not found.")
    return None


def extract_umi_and_ratio(record_id: str) -> Tuple[Optional[str], Optional[int]]:
    """
    Extract UMI and ratio information from FASTQ record ID.

    Args:
        record_id (str): Record ID in the FASTQ file.

    Returns:
        Tuple[Optional[str], Optional[int]]: Return UMI string and ratio as integer.
    """
    parts = record_id.split(':')
    if len(parts) >= 3:
        return parts[0], int(parts[-1])
    return None, None


def is_umi_similar(umi1: str, umi2: str) -> bool:
    """
    Check if two UMIs differ by only one base.

    Args:
        umi1 (str): First UMI.
        umi2 (str): Second UMI.

    Returns:
        bool: Whether the UMIs are similar (differ by only 1 base).
    """
    if len(umi1) != len(umi2):
        return False
    return sum(1 for a, b in zip(umi1, umi2) if a != b) == 1


def is_seq_similar(seq1: str, seq2: str, max_mismatch: int = 5) -> bool:
    """
    Check if two sequences are similar based on the number of mismatches in the shortest part.

    Args:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
        max_mismatch (int): Maximum allowed mismatch count.

    Returns:
        bool: Whether the sequences are similar.
    """
    min_len = min(len(seq1), len(seq2))
    return sum(1 for a, b in zip(seq1[:min_len], seq2[:min_len]) if a != b) <= max_mismatch


def filter_fastq_by_length(input_fastq: str, output_fastq: str, primer_to_cut_plus20bp: int) -> None:
    """
    Filter FASTQ records, collapse similar UMIs, and enforce length and N content limits.

    Args:
        input_fastq (str): Input FASTQ file path.
        output_fastq (str): Output filtered FASTQ file path.
        primer_to_cut_plus20bp (int): Minimum sequence length.
    """
    umi_dict = {}
    umi_seq_dict = {}

    for record in SeqIO.parse(input_fastq, "fastq"):
        umi, ratio = extract_umi_and_ratio(record.id)
        if umi and ratio:
            seq = str(record.seq)
            if umi in umi_dict:
                if ratio > umi_dict[umi]:
                    umi_dict[umi] = ratio
                    umi_seq_dict[umi] = seq
            else:
                umi_dict[umi] = ratio
                umi_seq_dict[umi] = seq

    umis_to_remove = set()
    umis = list(umi_dict.keys())

    for i in range(len(umis)):
        for j in range(i + 1, len(umis)):
            umi1, umi2 = umis[i], umis[j]
            if is_umi_similar(umi1, umi2):
                if umi_dict[umi2] > 5:
                    continue
                if is_seq_similar(umi_seq_dict[umi1], umi_seq_dict[umi2]):
                    if umi_dict[umi1] > umi_dict[umi2]:
                        umis_to_remove.add(umi2)
                    else:
                        umis_to_remove.add(umi1)

    for umi in umis_to_remove:
        umi_dict.pop(umi, None)
        umi_seq_dict.pop(umi, None)

    filtered_records = []
    for record in SeqIO.parse(input_fastq, "fastq"):
        umi, ratio = extract_umi_and_ratio(record.id)
        is_N_valid = record.seq.count('N') <= len(record.seq) * 0.05
        if umi and umi in umi_dict and primer_to_cut_plus20bp <= len(record.seq) <= 300 and is_N_valid:
            filtered_records.append(record)

    SeqIO.write(filtered_records, output_fastq, "fastq")
    print(f"Filtered and compressed {len(filtered_records)} records.")


def reverse_complement(sequence: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.

    Args:
        sequence (str): Input sequence.

    Returns:
        str: Reverse complement sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(sequence))


def reverse_cigar(cigar: str) -> str:
    """
    Reverse the operation order in a CIGAR string.

    Args:
        cigar (str): CIGAR string.

    Returns:
        str: Reversed CIGAR string.
    """
    new_cigar = []
    i = 0
    while i < len(cigar):
        num = ""
        while i < len(cigar) and cigar[i].isdigit():
            num += cigar[i]
            i += 1
        if i < len(cigar):
            new_cigar.insert(0, num + cigar[i])
            i += 1
    return ''.join(new_cigar)


def filter_and_reverse_sam(input_sam: str, output_sam: str) -> None:
    """
    Process SAM file: retain forward/reverse alignment results, reverse complement and reverse CIGAR for reverse sequences.

    Args:
        input_sam (str): Input SAM file.
        output_sam (str): Output processed SAM file.
    """
    with pysam.AlignmentFile(input_sam, "r") as sam_in, pysam.AlignmentFile(output_sam, "w", header=sam_in.header) as sam_out:
        for read in sam_in:
            if read.flag in [0, 16]:
                if read.flag == 16:
                    read.query_sequence = reverse_complement(read.query_sequence)
                    read.cigarstring = reverse_cigar(read.cigarstring)
                sam_out.write(read)
group_genome_mapping = {
    # AAVS1-P group
    "AAVS1-P1-T": "AAVS1-P-T",
    "AAVS1-P1-OT1": "AAVS1-P-T",
    "AAVS1-P1-OT2": "AAVS1-P-T",
    "AAVS1-P2-T": "AAVS1-P-T",
    "AAVS1-P2-OT1": "AAVS1-P-T",
    "AAVS1-P2-OT2": "AAVS1-P-T",

    # AAVS1-H group
    "AAVS1-H-T": "AAVS1-H-T",
    "AAVS1-H-OT1": "AAVS1-H-T",
    "AAVS1-H-OT2": "AAVS1-H-T",

    # CCR5-P group
    "CCR5-P1-T": "CCR5-P-T",
    "CCR5-P1-OT1": "CCR5-P-T",
    "CCR5-P1-OT2": "CCR5-P-T",
    "CCR5-P2-T": "CCR5-P-T",
    "CCR5-P2-OT1": "CCR5-P-T",
    "CCR5-P2-OT2": "CCR5-P-T",

    # CCR5-H group
    "CCR5-H-T": "CCR5-H-T",
    "CCR5-H-OT1": "CCR5-H-T",
    "CCR5-H-OT2": "CCR5-H-T",

    # KI type
    "AAVS1-P1-KI": "AAVS1-P-KI",
    "AAVS1-P2-KI": "AAVS1-P-KI",
    "AAVS1-H-KI": "AAVS1-H-KI",
    "CCR5-P1-KI": "CCR5-P-KI",
    "CCR5-P2-KI": "CCR5-P-KI",
    "CCR5-H-KI": "CCR5-H-KI",
}
