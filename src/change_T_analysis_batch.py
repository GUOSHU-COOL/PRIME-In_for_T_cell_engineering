# -*- coding: utf-8 -*-
#!usr/bin/python3


import sys
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import holoviews as hv
from holoviews import opts, dim
import os
import pysam
import pandas as pd
import re
from common_analysis_utils import (
    filter_and_reverse_sam,
    flash_stitch,
    blastn,
    filter_fastq_by_length,
    get_sample_info_from_excel,
    group_genome_mapping
)
def get_edit_count(bam, chr_seq, consolidated_fq_tri_filtered_file, output_dir, sample, db,strand):
    # Initialize the dictionary
    sm_main_chr_dic = {}  # Store the SM type alignment on the main chromosome
    ms_dic = {}  # Store the comparison of MS types
    sm_ki_dic = {}
    sm_translocation_dic = {}
    small_sm_deletion = {}
    medium_sm_deletion = {}
    large_sm_deletion = {}
    small_indels_dic = {}
    wt_dic = {}
    small_ms_deletion = {}
    medium_ms_deletion = {}
    large_ms_deletion = {}
    ms_ki_dic = {}
    ms_translocation_dic = {}
    other_dic = {}
    ms_other_dict = {}
    ms_translocation_type = {}
    sm_translocation_type = {}
    deletion_type_length = {}

    # A dictionary for storing sequences
    raw_fa_dic = {}
    print("consolidated_fq_tri_filtered_file",consolidated_fq_tri_filtered_file)
    # Read FASTQ
    with open(consolidated_fq_tri_filtered_file, "r") as raw_fq:
        for record in SeqIO.parse(raw_fq, "fastq"):
            if record.seq is None:
                print(f"empty序列: {record.id}")
            else:
                raw_fa_dic[record.id] = str(record.seq)  # record ID and seq

    # read BAM file
    bam_file = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    for read in bam_file:
        if read.qname in raw_fa_dic:
            cigar = str(read.cigarstring)
            letter = re.findall('\D', cigar)
            number = re.findall('\d+', cigar)
            number = list(map(int, number))  # Convert the numeric part to an integer

            # Filter CIGAR strings containing 'H'
            if 'H' in letter:
                print("H in CIGAR")
                continue


            condition2 = read.query_sequence.count('N') <= len(read.query_sequence) * 0.05  # The quantity of N does not exceed 5%

            if  condition2:
                # Determine the string type of the CIGAR
                if 'S' in letter and 'M' in letter:
                    if letter.count('S') == 1 and letter.count('M') == 1:
                        if letter.index('S') < letter.index('M'):
                            # Extract the sequence of part S
                            s_length = number[letter.index('S')]
                            s_sequence = read.query_sequence[:s_length]
                            match_start_index = read.blocks[0][0] # Obtain the starting point of the match sequence in the sam sequence # plus 1 would match start from 0
                            match_end_index = read.blocks[0][1]

                            # Obtain the alignment chromosome of part M
                            ref_name = read.reference_name  # Name of the comparison chromosome

                            # Determine the alignment chromosome of part M
                            if ref_name == "chr_transgene":
                                # if is transgene，categorize into KI
                                #sm_ki_dic[read.qname] = [read.query_sequence, match_start_index, match_end_index, cigar,ref_name]
                                sm_ki_dic[read.qname] = [read.query_sequence, cigar] # for unify
                            elif ref_name != chr_seq:
                                # if is other chrom，categorize into  Translocation
                                #sm_translocation_dic[read.qname] = [read.query_sequence, match_start_index, match_end_index, cigar,ref_name]
                                sm_translocation_dic[read.qname] = [read.query_sequence, cigar]
                                if strand =="+":
                                    sm_translocation_type[read.qname]={"ref_name" : ref_name,
                                                                       "chr_seq" : chr_seq,
                                                                        "reference_start":read.reference_start,
                                                                       "reference_end":read.reference_end,
                                                                       } # both strand + - read.reference_start<reference_end
                                else:
                                    sm_translocation_type[read.qname] = {"ref_name": ref_name,
                                                                        "chr_seq": chr_seq,
                                                                        "reference_start": read.reference_end,
                                                                         "reference_end": read.reference_end,
                                                                        }
                            else:
                                # If it is the main chromosome, it is categorize into SM
                                #sm_main_chr_dic[read.qname] = [read.query_sequence, cigar]
                                sm_main_chr_dic[read.qname] = [read.query_sequence, s_sequence, match_start_index,match_end_index, cigar]
                        else:
                            # deal MS type
                            s_length = number[letter.index('S')]
                            s_sequence = read.query_sequence[-s_length:]
                            match_start_index = read.blocks[0][0]
                            match_end_index = read.blocks[0][1]
                            ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
                            #ms_dic[read.qname] = [read.query_sequence,cigar]
                elif 'M' in letter and ('D' in letter or 'I' in letter):
                    # Strictly determine whether it is M... M-form
                    if len(letter) >= 3 and letter[0] == 'M' and letter[-1] == 'M':
                        # As long as it contains D or I in the middle, it is categorized into small_indels_dic
                        if any(mid in letter[1:-1] for mid in ['D', 'I', 'M']):
                            small_indels_dic[read.qname] = [read.query_sequence, cigar]
                            match = re.findall(r'(\d+)D', cigar)
                            deletion_length = sum(map(int, match))
                            # save to dict
                            deletion_type_length[read.qname] = {
                                "deletion_type": "small_indels_dic",  # small_indels_dic
                                "deletion_length": deletion_length
                            }
                        else:
                            other_dic[read.qname] = [read.query_sequence, cigar]
                    else:
                        other_dic[read.qname] = [read.query_sequence, cigar]
                elif letter == ['M']:
                    wt_dic[read.qname] = [read.query_sequence, cigar]
                else:
                    other_dic[read.qname] = [read.query_sequence, cigar]
    # 输出统计结果
    print('SM type align number:', len(sm_main_chr_dic))
    print('KI type align number:', len(sm_ki_dic))
    print('Translocation type align number:', len(sm_translocation_dic))
    print('Small Indels type align number:', len(small_indels_dic))
    print('WT type align number:', len(wt_dic))
    print('MS type align number:', len(ms_dic))
    print('Other type align number:', len(other_dic))

    # Write the sequence of the S part of the SM type to the FASTA file
    if len(sm_main_chr_dic)!=0:
        sm_main_chr_s_fa = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_primeradd5bp_noadapter_consolidated_tri_SM_S.fa')
        os.makedirs(os.path.join(output_dir, 'blastn'), exist_ok=True) # create first
        with open(sm_main_chr_s_fa, 'w') as f:
            for key in sm_main_chr_dic.keys():
                s_seq = sm_main_chr_dic[key][1]
                f.write(f'>{key}\n{s_seq}\n')
        sm_main_chr_s_fa_output= os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_noadapter_consolidated_tri_SM_S_fa_to_hg38fa_blastn_result.txt')
        blastn(sm_main_chr_s_fa, db, sm_main_chr_s_fa_output)
        sm_main_chr_best_alignment_dict=extract_best_alignment_qstart(sm_main_chr_s_fa_output)
        small_sm_deletion, medium_sm_deletion, large_sm_deletion=categorize_sm(sm_main_chr_best_alignment_dict, sm_main_chr_dic,deletion_type_length,strand)
        # using blastn to align（only to SM type in main chromosome）
    # Write the S part sequence of the MS type to the FASTA file
    os.makedirs(os.path.join(output_dir, 'blastn'), exist_ok=True)

    # Write the S part sequence of the MS type to the FASTA file
    if len(ms_dic)!=0:
        ms_s_fa = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_primeradd5bp_noadapter_consolidated_tri_MS_S.fa')
        with open(ms_s_fa, 'w') as f:
            for key in ms_dic.keys():
                s_seq = ms_dic[key][1]
                f.write(f'>{key}\n{s_seq}\n')

        # define blastn output file
        ms_s_fa_output = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_noadapter_consolidated_tri_MS_S_fa_to_hg38fa_blastn_result.txt')

        # Call blastn for comparison (only for MS types on the main chromosome)
        blastn(ms_s_fa, db, ms_s_fa_output)

        # Extract the best alignment information from the blastn comparison results
        ms_main_chr_best_alignment_dict, ms_other_dict = extract_best_alignment_qend(ms_s_fa_output,ms_dic)

        # Category deletion type (classify based on the comparison results)
        #ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
        small_ms_deletion, medium_ms_deletion, large_ms_deletion, ms_translocation_dic,ms_ki_dic,ms_translocation_type= (
            categorize_ms(ms_main_chr_best_alignment_dict, ms_dic,chr_seq,deletion_type_length,strand))

    return sm_ki_dic, sm_translocation_dic,small_sm_deletion, medium_sm_deletion, large_sm_deletion,small_indels_dic, wt_dic,small_ms_deletion, medium_ms_deletion, large_ms_deletion, ms_ki_dic, ms_translocation_dic, other_dic,ms_other_dict,ms_translocation_type,sm_translocation_type,deletion_type_length

def extract_best_alignment_qstart(blastn_result):
    """
    Extract the best alignment for each query sequence in the BLASTN alignment results (based on qstart == 0).

    parameters：
    - blastn_result: BLASTN alignment file path

    返回：
    - best_alignment_dict: Store the best alignment information for each query sequence
    """
    # read BLASTN alignment
    try:
        df = pd.read_csv(blastn_result, delimiter='\t', header=None)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        # groupby qseqid 
        grouped = df.groupby('qseqid')

        # A dictionary for storing the results
        best_alignment_dict = {}

        # Traverse each group and extract the ones where qstart == 0
        for name, group in grouped:
            group = group[group['qstart'] == 1]  # Only keep the rows where qstart == 1

            if not group.empty:
                # Align the rows with the maximum percentage in sequence
                max_row = group.loc[group['pident'].idxmax()]

                # Parse the required fields
                qseqid = max_row['qseqid']  # as key
                s_sseqid = max_row['sseqid']
                s_sstart = int(max_row['sstart'])
                s_send = int(max_row['send'])

                # Store in the dictionary
                best_alignment_dict[qseqid] = {
                    'sseqid': s_sseqid,
                    'sstart': s_sstart,
                    'send': s_send
                }

        return best_alignment_dict
    except Exception as e:
        print(f"read file  {blastn_result} have mistake: {e}")
        return {}

def extract_best_alignment_qend(blastn_result,ms_dic):
    """
    Extract the best alignment for each query sequence in the BLASTN alignment results (based on qend == s_length).
    If there is no matching of qend == s_length, put it in the 'other' dictionary.

    parameters：
    - blastn_result: BLASTN alignment file path
    - s_length: The sequence length used for comparison, filter the alignment of qend == s_length

    output：
    - best_alignment_dict: Store the best alignment information for each query sequence
    - other_dict: Store the query sequence that does not match qend == s_length
    """
    # read BLASTN alignment
    try:
        df = pd.read_csv(blastn_result, delimiter='\t', header=None)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        # groupby qseqid 
        grouped = df.groupby('qseqid')

        # A dict for storing the results
        best_alignment_dict = {}
        other_dict = {}

        # Traverse each group and extract the lines qend == s_length
        for name, group in grouped:
            group_id=group['qseqid'].iloc[0]
            s_length=len(ms_dic[group_id][1])
            origin_seq=ms_dic[group_id][0]
            cigar=ms_dic[group_id][4]
            qseqid = group['qseqid'].iloc[0]
            group = group[group['qend'] == s_length]  # all qseqid is same

            if not group.empty:
                # Align the rows with the maximum percentage in sequence
                max_row = group.loc[group['pident'].idxmax()]

                # Parse the required fields
                qseqid = max_row['qseqid']  # as key
                s_sseqid = max_row['sseqid']
                s_sstart = int(max_row['sstart'])
                s_send = int(max_row['send'])

                # Store in the dictionary
                best_alignment_dict[qseqid] = {
                    'sseqid': s_sseqid,
                    'sstart': s_sstart,
                    'send': s_send
                }
            else:
                # If qend == s_length is not matched, place it in the 'other' dictionary
                other_dict[qseqid] = [origin_seq,cigar]

        return best_alignment_dict, other_dict
    except Exception as e:
        print(f"read file  {blastn_result} have mistake: {e}")
        return {}, {}

def categorize_sm(sm_main_chr_best_alignment_dict, sm_main_chr_dic,deletion_type_length,strand):
    # sm_main_chr_best_alignment_dict blast seq
    small_deletion = {}
    medium_deletion = {}
    large_deletion = {}

    for key in sm_main_chr_best_alignment_dict:
        if key in sm_main_chr_dic:
            # Calculate the missing size
            match_start = sm_main_chr_dic[key][2] #both strand+ strand- and flag 0&16 start<end
            match_end = sm_main_chr_dic[key][3]
            start = sm_main_chr_best_alignment_dict[key]["sstart"] # -strand start>send
            send= sm_main_chr_best_alignment_dict[key]["send"]
            if strand == "-":
                deletion_size = abs(send - match_end)
            else:
                deletion_size = abs(match_start - send)
            # Generate all possible coordinate pairs
            # sm_main_chr_dic[read.qname] = [read.query_sequence, cigar]
            #sm_main_chr_dic[read.qname] = [read.query_sequence, s_sequence, match_start_index, match_end_index, cigar]
            # Construct value list
            value_list = [sm_main_chr_dic[key][0], sm_main_chr_dic[key][4]]

            # Classify and store in the corresponding dictionary
            if 0 <= deletion_size < 200:
                small_deletion[key] = value_list
                deletion_type_length[key] = {
                    "deletion_type": "sm_small_indels",  # small_indels_dic
                    "deletion_length": deletion_size
                }
            elif 200 <= deletion_size < 400:
                medium_deletion[key] = value_list
                deletion_type_length[key] = {
                    "deletion_type": "sm_medium_indels",  # medium_indels_dic
                    "deletion_length": deletion_size
                }
            elif deletion_size >= 400:
                large_deletion[key] = value_list
                deletion_type_length[key] = {
                    "deletion_type": "sm_large_indels",  # large_indels_dic
                    "deletion_length": deletion_size
                }
            else:
                print(deletion_size)
                deletion_type_length[key] = {
                    "deletion_type": 'sm_other',  # other_dic
                    "deletion_length": deletion_size
                }
                #logging.debug(f"Key: {key}, Deletion size: {deletion_size}, Values: {value_list}")

    return small_deletion, medium_deletion, large_deletion


def categorize_ms(ms_main_chr_best_alignment_dict, ms_main_chr_dic,main_chr,deletion_type_length,strand):
    # ms_main_chr_best_alignment_dict blast seq
    # ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
    small_deletion = {}
    medium_deletion = {}
    large_deletion = {}
    translocation = {}
    ki_events = {}
    ms_translocation_type={}
    for key, alignment_info in ms_main_chr_best_alignment_dict.items():
        if key in ms_main_chr_dic:
            sseqid = alignment_info["sseqid"]  # get align chrom number
            sstart = alignment_info["sstart"]
            send = alignment_info["send"]
            value_list = [ms_main_chr_dic[key][0], ms_main_chr_dic[key][4]]

            if sseqid == main_chr:
                # match_end- sstart
                # Calculate the missing size
                match_start = ms_main_chr_dic[key][2]
                match_end = ms_main_chr_dic[key][3]
                sstart = ms_main_chr_best_alignment_dict[key]["sstart"]
                send = ms_main_chr_best_alignment_dict[key]["send"]
                if strand == "-":
                    deletion_size = abs(match_start - sstart)
                else:
                    deletion_size = abs(sstart - match_end)
                if 0 <= deletion_size < 200:
                    small_deletion[key] = value_list
                    deletion_type_length[key] = {
                        "deletion_type": "ms_small_indels",  # small_indels_dic
                        "deletion_length": deletion_size
                    }
                elif 200 <= deletion_size < 400:
                    medium_deletion[key] = value_list
                    deletion_type_length[key] = {
                        "deletion_type": "ms_medium_indels",  # small_medium_dic
                        "deletion_length": deletion_size
                    }
                elif deletion_size >= 400:
                    large_deletion[key] = value_list
                    deletion_type_length[key] = {
                        "deletion_type": "ms_large_indels",  # small_large_dic
                        "deletion_length": deletion_size
                    }
                else:
                    deletion_type_length[key] = {
                        "deletion_type": "ms_other",  # other_dic
                        "deletion_length": deletion_size
                    }
                    print(f"Unexpected deletion size: {deletion_size}")

            elif sseqid == "chr_transgene":
                ki_events[key] = value_list  # regard as KI event

            else:
                ms_translocation_type[key] = {"ref_name": sseqid,
                                                     "chr_seq": main_chr,
                                              "reference_start": sstart,
                                              "reference_end": send
                                              }  # both strand + - read.reference_start<reference_end
                translocation[key] = value_list  # 认为是易位（translocation）

    return small_deletion, medium_deletion, large_deletion, translocation, ki_events,ms_translocation_type

def extract_prefix(sample_name):
    """
    Extract the part before the first digit from the sample name
    as 'AH2-OT1' -> 'AH2'
    """
    match = re.match(r"^[^\d]*\d", sample_name)  # Match the characters from the beginning to the first digit
    if match:
        return match.group(0)
    else:
        return sample_name  # If no match is found, return the original original name


def add_columns_from_dicts(excel_file, deletion_type_length_type, ms_translocation_type, sm_translocation_type,
                           output_file):
    # read Excel file
    df = pd.read_excel(excel_file, index_col=1)  # Use the second column as the index

    # Handle the deletion_type_length_type dictionary and add the deletion_type and deletion_length columns
    df['deletion_type'] = df.index.map(lambda key: deletion_type_length_type.get(key, {}).get('deletion_type', None))
    df['deletion_length'] = df.index.map(
        lambda key: deletion_type_length_type.get(key, {}).get('deletion_length', None))

    # Handle the ms_translocation_type dictionary and add the translocation_chr column
    df['main_chr'] = df.index.map(lambda key: ms_translocation_type.get(key, {}).get('chr_seq', None))
    df['translocation_chr'] = df.index.map(lambda key: ms_translocation_type.get(key, {}).get('ref_name', None))
    df['reference_start'] = df.index.map(lambda key: ms_translocation_type.get(key, {}).get('reference_start', None))
    df['reference_end'] = df.index.map(lambda key: ms_translocation_type.get(key, {}).get('reference_end', None))
    # Handle the sm_translocation_type dictionary and make sure to convert it to Series
    df['main_chr'] = df['main_chr'].combine_first(
        pd.Series(df.index.map(lambda key: sm_translocation_type.get(key, {}).get('chr_seq', None)), index=df.index)
    )
    df['translocation_chr'] = df['translocation_chr'].combine_first(
        pd.Series(df.index.map(lambda key: sm_translocation_type.get(key, {}).get('ref_name', None)), index=df.index)
    )
    df['reference_start'] = df['reference_start'].combine_first(
        pd.Series(df.index.map(lambda key: sm_translocation_type.get(key, {}).get('reference_start', None)),
                  index=df.index)
    )
    df['reference_end'] = df['reference_end'].combine_first(
        pd.Series(df.index.map(lambda key: sm_translocation_type.get(key, {}).get('reference_end', None)),
                  index=df.index)
    )

    df.to_excel(output_file)
    print(f"✅ data have save to {output_file}")

def save_dicts_to_excel(output_file,output_evaluate_sample_dir, **dicts):
    """
    Save multiple dictionaries to one Excel file, with the name of each dictionary as "Category"

    parameters：
    - output_file: str, output Excel filepath
    - **dicts: key parameters，each dict key is sequence ID，value is [sequence, cigar]
    """
    all_data = []  # Store all the data of the dictionaries
    category_counts = {}  # Used for counting the quantity of each category

    for category, data_dict in dicts.items():
        count = 0  # Count the number of current categories
        for key, values in data_dict.items():
            if len(values) >= 2:  # Make sure there is sufficient data
                sequence, cigar = values[0], values[1]
                all_data.append([category, key, sequence, cigar])
                count += 1  # Update the count of the current category
        category_counts[category] = count  # The number of saved categories

    df = pd.DataFrame(all_data, columns=["Category", "Key", "Sequence", "CIGAR"])

    df.to_excel(output_file, index=False)
    print(f"✅ data have save to {output_file}")

    # Print the quantity for each category
    print("📊 the quantity for each category：")
    for category, count in category_counts.items():
        print(f"  {category}: {count} pieces of data")
    # Calculate the total number of categories
    total_count = sum(category_counts.values())

    # Calculate the proportion of each category
    category_proportions = {category: count / total_count for category, count in category_counts.items()}

    # Print the proportion of each category
    print("📊 the proportion of each category：")
    for category, proportion in category_proportions.items():
        print(f"  {category}: {proportion * 100:.2f}%")
    # Draw a stacked bar chart
    categories = list(category_proportions.keys())
    values = list(category_proportions.values())
    # extract sample
    sample_name = os.path.basename(output_evaluate_sample_dir)

    plt.bar(categories, values, color='skyblue')


    plt.title(f"{sample_name} Proportion of Each Category in the Total Data", fontsize=16)
    plt.xlabel("Categories", fontsize=14)
    plt.ylabel("Proportion", fontsize=14)


    plt.xticks(rotation=45, ha='right', fontsize=12)  # ha='right' 使标签右对齐，避免偏移


    for i, (category, proportion) in enumerate(category_proportions.items()):
        plt.text(i, proportion + 0.01, f"{proportion * 100:.2f}%", ha='center', va='bottom', fontsize=6)

    plt.text(0.95, 0.95, f"Total Count: {total_count}", ha='right', va='top', fontsize=10,
             transform=plt.gca().transAxes)
    plt.ylim(0, 1.05)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.savefig(os.path.join(output_evaluate_sample_dir,"category_proportions_freq.png"))
    plt.show()


def plot_deletion_length_distribution(excel_file, output_evaluate_sample_dir):
    df = pd.read_excel(excel_file)

    # Convert to a numeric value and remove NaN
    deletion_lengths = pd.to_numeric(df['deletion_length'], errors='coerce').dropna()

    # Filter out particularly small or extreme values, such as those less than 1
    deletion_lengths = deletion_lengths[deletion_lengths > 0]

    # Perform a log10 transformation on the data to avoid stretching large values
    log_deletion_lengths = np.log10(deletion_lengths)

    # extract sample
    sample_name = os.path.basename(output_evaluate_sample_dir)

    # 1. Draw logarithmic X-axis and logarithmic Y-axis histograms
    plt.figure(figsize=(10, 6))
    plt.hist(log_deletion_lengths, bins=50, color='skyblue', edgecolor='black')
    plt.yscale("log")  
    plt.title(f"{sample_name} Deletion Length Frequency Distribution (Log X)", fontsize=16)
    plt.xlabel("Log10(Deletion Length)", fontsize=14)
    plt.ylabel("Frequency (log scale)", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_evaluate_sample_dir, "deletion_length_distribution_freq_logx.png"))
    plt.show()
    # Filtering data: Only retain data points where deletion_length is ≤ 1000
    filtered_deletion_lengths = deletion_lengths[deletion_lengths <= 1000]
    plt.figure(figsize=(10, 6))
    plt.hist(filtered_deletion_lengths, bins=50, color='salmon', edgecolor='black')
    plt.yscale("log")  
    plt.xlim([0, 1000])  # Limit the X-axis range to within 1k
    plt.title(f"{sample_name} Deletion Length Frequency Distribution (≤ 1k)", fontsize=16)
    plt.xlabel("Deletion Length", fontsize=14)
    plt.ylabel("Frequency (log scale)", fontsize=14)

    # Set the X-axis scale, marking one every 100
    plt.xticks(np.arange(0, 1001, 100))

    plt.tight_layout()
    plt.savefig(os.path.join(output_evaluate_sample_dir, "deletion_length_distribution_freq_1k_logy.png"))
    plt.show()
def read_sample_primer_distance_excel(file_path):
    # real Excel file
    df = pd.read_excel(file_path, dtype=str)  # Read in string format to prevent data loss

    # Build a dictionary, with the first column being "key" and the last two columns in dictionary form
    data_dict = {
        row["sample"]: {"primer_plus_5bp": row["primer_plus_5bp"],
                        "primer_to_cut_plus20bp": row["primer_to_cut_plus20bp"]}
        for _, row in df.iterrows()
    }

    return data_dict


def process_sample(sample, cm_dir, output_evaluate_dir,excel_file):
    # sample = 'Asn2-T2-T3'
    # cm_dir = '/data5/wangxin/20241001_wcx'  # TODO: /data5/wangxin/20241001_wcx 后侧没有/
    # output_evaluate_dir ="/data5/wangxin/20241001_wcx/shuyu/20250307"


    sample_data=get_sample_info_from_excel(excel_file,sample)
    strand = sample_data["strand"]
    primer_to_cut_plus20bp=int(sample_data["primer_to_cut_plus20bp"])
    primer_plus_5bp_sequence=sample_data["primer_plus_5bp"].upper()
    chr_seq = sample_data["chr_seq"]
    cut_site = sample_data["cut_site"] # not use in T group
    adapter_sequence ="ACCACGCGTGCTCTACA"
##############  This part code only needs to be run once.
    ######### step1: merge R1 and R2 files to one merged fq file
    flash_stitch(cm_dir, sample)

    ######### step2: mapping using bwa (from merged fq file to bam)
    output_dir = cm_dir  + sample + '/bwa/'
    cmd = "mkdir {}".format(output_dir)
    print(cmd)
    os.system(cmd)

    ######## step3: UMI normalization
    output_dir = cm_dir  + sample + '/consolidate/'
    cmd = "mkdir {}".format(output_dir)
    print(cmd)
    os.system(cmd)
    fq_file = cm_dir  + sample + '/flash/' + sample + '.merge.fastq'
    consolidated_fq_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated.fq'
    cmd = '/home/wangxin/miniconda3_1/envs/PEM-Q/bin/python consolidate_batch.py {} {}  {} 14 0.7 1 2 {} {} {}'.format(
        sample, \
        fq_file, \
        consolidated_fq_file, \
        adapter_sequence, \
        primer_to_cut_plus20bp, \
        primer_plus_5bp_sequence,\
    )
    print(cmd)
    os.system(cmd)
    ######## step4: remove adapter using trimmomatic
    output_dir = cm_dir  + sample + '/consolidate/'
    consolidated_fq_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated.fq'
    consolidated_fq_gz_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated.fq.gz'
    consolidated_fq_tri_gz_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri.fq.gz'
    adapter_sequence_fa_file = '/data5/wangxin/20241001_wcx/PEM-seq-sequences-AAVS1/adapter_sequence.fa'

    cmd = 'gzip -f {}'.format(consolidated_fq_file)
    print(cmd)
    os.system(cmd)

    cmd = 'java -jar /data2/wangxin/biosoft/Trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar SE \
                                                                                                -phred33 \
                                                                                                {} \
                                                                                                {} \
                                                                                                ILLUMINACLIP:{}:2:12:10'.format(
        consolidated_fq_gz_file, \
        consolidated_fq_tri_gz_file, \
        adapter_sequence_fa_file)
    print(cmd)
    os.system(cmd)
    cmd = 'gunzip -f {}'.format(consolidated_fq_tri_gz_file)
    print(cmd)
    os.system(cmd)
#########filtered by sequence length
    consolidated_fq_tri_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri.fq'
    consolidated_fq_tri_filtered_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_filtered.fq'
    filter_fastq_by_length(consolidated_fq_tri_file, consolidated_fq_tri_filtered_file, primer_to_cut_plus20bp)
    ####### step5: mapping from consolidated_tri_fq to ref genome (add transgene) using bwa
    output_dir = cm_dir  + sample + '/consolidate/'
    output_sam = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_change.sam'
    genome_name= group_genome_mapping[sample]
    genome_add_transgene = os.path.join("/data5/shuyu/result/modify_genome_db/", genome_name, f"{genome_name}_transgene.fa")

    cmd = 'bwa mem  {} {} -t 16 > {}'.format(
        genome_add_transgene, consolidated_fq_tri_filtered_file, output_sam)
    print(cmd)
    os.system(cmd)

    filter_and_reverse_output_sam = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_filter_and_reverse.sam'

    filter_and_reverse_sam(output_sam, filter_and_reverse_output_sam)


    output_filter_and_reverse_bam = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_filter_and_reverse.bam'
    cmd = 'samtools sort {} -@ 15 -o {} -O BAM'.format(filter_and_reverse_output_sam, output_filter_and_reverse_bam)
    print(cmd)
    os.system(cmd)

    ######## step6 : from consolidated_tri_fq to fa
    output_dir = cm_dir  + sample + '/consolidate/'  # TODO: 这个位置的输入 /data5/wangxin/20241001_wcx/  /cleandata/20241227/这个就得变成cleandata/20241227/

############# Classify for different situations
    ######## step7: counts of different editing events (KI, small indels, large deletions, translocations, nojunctions, other editing)
    output_dir = cm_dir  + sample + '/consolidate/'
#######
    output_filter_and_reverse_bam = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_filter_and_reverse.bam'
    consolidated_fq_tri_filtered_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_filtered.fq'
##############
    ####### generate blastn db

    db= os.path.join("/data5/shuyu/result/modify_genome_db/", genome_name,f"{genome_name}_transgene_db")
    #db = '/data5/wangxin/20241001_wcx/PEM-seq-sequences-AAVS1/bwa_mask_KIplasmid_EF1A_addtransgene383bp_hg38.fa_blastndb'
    sm_ki_dic, sm_translocation_dic,small_sm_deletion, medium_sm_deletion, large_sm_deletion,small_indels_dic, wt_dic,small_ms_deletion, medium_ms_deletion, large_ms_deletion, ms_ki_dic, ms_translocation_dic, other_dic,ms_other_dict,ms_translocation_type,sm_translocation_type,deletion_type_length= get_edit_count(output_filter_and_reverse_bam, chr_seq, consolidated_fq_tri_filtered_file, output_dir, sample, db,strand)
    # call save_dicts_to_excel  passes in all the dictionaries

    output_evaluate_sample_dir = os.path.join(output_evaluate_dir,sample)
    cmd = "mkdir {}".format(output_evaluate_sample_dir)
    print(cmd)
    os.system(cmd)

    save_dicts_to_excel(
        os.path.join(output_evaluate_dir,sample,"output.xlsx"),
        output_evaluate_sample_dir,
        sm_ki_dic=sm_ki_dic,
        sm_translocation_dic=sm_translocation_dic,
        small_sm_deletion=small_sm_deletion,
        medium_sm_deletion=medium_sm_deletion,
        large_sm_deletion=large_sm_deletion,
        small_indels_dic=small_indels_dic,
        wt_dic=wt_dic,
        small_ms_deletion=small_ms_deletion,
        medium_ms_deletion=medium_ms_deletion,
        large_ms_deletion=large_ms_deletion,
        ms_ki_dic=ms_ki_dic,
        ms_translocation_dic=ms_translocation_dic,
        other_dic=other_dic,
        ms_other_dict=ms_other_dict
    )
    add_columns_from_dicts(os.path.join(output_evaluate_dir,sample,"output.xlsx"), deletion_type_length, ms_translocation_type, sm_translocation_type,
                           os.path.join(output_evaluate_dir,sample,"add_output.xlsx"))
    # call function and pass in Excel file path
    excel_file = os.path.join(output_evaluate_dir,sample,"add_output.xlsx")
    plot_deletion_length_distribution(excel_file,output_evaluate_sample_dir)
