# -*- coding: utf-8 -*-
#!usr/bin/python3
import regex as re2
import os
import pysam
import re
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd

from common_analysis_utils import (
    filter_and_reverse_sam,
    flash_stitch,
    blastn,
    filter_fastq_by_length,
    get_sample_info_from_excel,
    group_genome_mapping
)

def get_edit_count(bam, chr_seq, consolidated_fq_tri_filtered_file, output_dir, sample, db,target_group,cut_site):
    # Initialize the dictionary
    KI={}
    reverseKI={}
    offtargetKI={}
    ms_dic = {}  # Store a dictionary of type MS

    # A dictionary for storing sequences
    raw_fa_dic = {}
    print("consolidated_fq_tri_filtered_file",consolidated_fq_tri_filtered_file)
    # Read FASTQ 
    with open(consolidated_fq_tri_filtered_file, "r") as raw_fq:
        for record in SeqIO.parse(raw_fq, "fastq"):
            if record.seq is None:
                print(f"empty sequence: {record.id}")
            else:
                raw_fa_dic[record.id] = str(record.seq)  # record ID and seq

    # read BAM file
    bam_file = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    for read in bam_file:
        if read.qname in raw_fa_dic:
            if read.reference_name != "chr_transgene":
                continue
            if "AP" in sample:
                unwanted_sequence = "cggtggGAGCTGGACGGCGACGTAAACGGGGCACTTTTCGGGGAAATG"  # Replace it with your target sequence
                max_mismatches = 1  # The maximum allowable misfit
                pattern = f"({unwanted_sequence}){{e<={max_mismatches}}}"  # A maximum of two mismatches are allowed

                if re.search(pattern, read.query_sequence, re2.IGNORECASE):
                    print("plasmid")
                    continue  # If a match that meets the requirements is found, skip this read
            cigar = str(read.cigarstring)
            letter = re.findall('\D', cigar)
            letters = "".join(re.findall(r'\D', cigar))
            if set(letters) != {'M', 'S'}:  # only allow M and S
                continue

            if not letters.startswith('M'):  # M It must be before S
                continue
            number = re.findall('\d+', cigar)
            number = list(map(int, number))  # Convert the numeric part to an integer


            condition2 = read.query_sequence.count('N') <= len(read.query_sequence) * 0.05  # The quantity of N does not exceed 5%

            if condition2:
                # deal MS type
                s_length = number[letter.index('S')]
                s_sequence = read.query_sequence[-s_length:]
                match_start_index = read.blocks[0][0]
                match_end_index = read.blocks[0][1]
                ms_dic[read.qname] = [read.query_sequence, s_sequence, match_start_index, match_end_index, cigar]

    print('MS type align numbers:', len(ms_dic))
    # Write the S part sequence of the MS type to the FASTA file
    os.makedirs(os.path.join(output_dir, 'blastn'), exist_ok=True)
    if len(ms_dic)!=0:
        ms_s_fa = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_primeradd5bp_noadapter_consolidated_tri_MS_S.fa')
        with open(ms_s_fa, 'w') as f:
            for key in ms_dic.keys():
                s_seq = ms_dic[key][1]
                f.write(f'>{key}\n{s_seq}\n')

        # define blastn output file
        ms_s_fa_output = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_noadapter_consolidated_tri_MS_S_fa_to_hg38fa_blastn_result.txt')

        # using blastn to align
        blastn(ms_s_fa,
               db, ms_s_fa_output)

        # Extract the best alignment information from the blastn comparison results
        ms_main_chr_best_alignment_dict, ms_other_dict = extract_best_alignment_qend(ms_s_fa_output,ms_dic)

        # Category deletion type (classify based on the comparison results)
        #ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
        KI,reverseKI,offtargetKI= categorize_ms(ms_main_chr_best_alignment_dict, ms_dic,chr_seq,target_group,cut_site)

    return  KI,reverseKI,offtargetKI

def extract_best_alignment_qend(blastn_result,ms_dic):
    """
    Extract the best alignment for each query sequence in the BLASTN alignment results (based on qend == s_length).
    If there is no matching of qend == s_length, put it in the 'other' dictionary.

    parameters：
    - blastn_result: BLASTN compares the file path of the result
    - s_length: The sequence length used for comparison, filter the alignment of qend == s_length

    output：
    - best_alignment_dict: Store the best alignment information for each query sequence
    - other_dict: Store the query sequence that does not match qend == s_length
    """
    # read BLASTN align result
    try:
        df = pd.read_csv(blastn_result, delimiter='\t', header=None)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        # Group by qseqid
        grouped = df.groupby('qseqid')

        # A dictionary for storing the results
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

                # save in dict
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
        print(f"read file {blastn_result} have mistake: {e}")
        return {}, {}


def categorize_ms(ms_main_chr_best_alignment_dict, ms_main_chr_dic,main_chr,target_group,cut_site):
    # ms_main_chr_best_alignment_dict blast seq
    # ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
    KI={}
    reverseKI={}
    offtargetKI={}

    for key, alignment_info in ms_main_chr_best_alignment_dict.items():
        if key in ms_main_chr_dic:
            sseqid = alignment_info["sseqid"]  # Obtain the chromosome number that has been matched
            start = ms_main_chr_best_alignment_dict[key]["sstart"] # -strand start>send
            send= ms_main_chr_best_alignment_dict[key]["send"]
            # Construct value list
            value_list = [ms_main_chr_dic[key][0], ms_main_chr_dic[key][4],sseqid,start,send] # read.query_sequence,cigar

            if sseqid == main_chr:
                if target_group == "A":
                    # sstart<send
                    if start<send: # same to + strand because reverse in
                        if abs(start-cut_site)<1000:
                            KI[key] = value_list
                        else:
                            offtargetKI[key] = value_list
                    else:
                        reverseKI[key] = value_list
                else:
                    if start>send:
                        if abs(start-cut_site)<1000:
                            KI[key] = value_list
                        else:
                            offtargetKI[key] = value_list
                    else:
                        reverseKI[key] = value_list
            else:
                offtargetKI[key] = value_list  # regard as  KI event


    return KI, reverseKI, offtargetKI

def save_dicts_to_excel(output_file,output_evaluate_sample_dir, **dicts):
    """
    Save multiple dictionaries to one Excel file, with each dictionary named "Category".

    parameter：
    - output_file: str, output Excel file path
    - **dicts: keyword parameter，each dict key is sequence ID，value is [sequence, cigar]
    """
    all_data = []  # Store all the data of the dictionaries
    category_counts = {}  # Used for counting the quantity of each category

    # Traverse the dictionary
    for category, data_dict in dicts.items():
        count = 0  # Count the number of current categories
        for key, values in data_dict.items():
            if len(values) >= 2:  # Make sure there is sufficient data
                sequence, cigar,chr_seq,start,send = values[0], values[1],values[2],values[3],values[4]
                all_data.append([category, key, sequence, cigar,chr_seq,start,send])
                count += 1  # Update the count of the current category
        category_counts[category] = count  # The number of saved categories

    # create DataFrame
    df = pd.DataFrame(all_data, columns=["Category", "Key", "Sequence", "CIGAR", "chr_seq","Start", "Send"])

    # save to Excel
    df.to_excel(output_file, index=False)
    print(f"✅ data have save to {output_file}")

    # Print the quantity of each category
    print("📊 Statistics of quantities in various categories：")
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
    categories = list(category_proportions.keys())
    values = list(category_proportions.values())
    sample_name = os.path.basename(output_evaluate_sample_dir)
    plt.bar(categories, values, color='skyblue')
    plt.title(f"{sample_name} Proportion of Each Category in the Total Data", fontsize=16)
    plt.xlabel("Categories", fontsize=14)
    plt.ylabel("Proportion", fontsize=14)

    plt.xticks(rotation=45, ha='right', fontsize=12)  # ha='right' 使标签右对齐，避免偏移
    for i, (category, proportion) in enumerate(category_proportions.items()):
        plt.text(i, proportion + 0.01, f"{proportion * 100:.2f}%", ha='center', va='bottom', fontsize=6)
    # Add the total number of categories in the upper right corner
    plt.text(0.95, 0.95, f"Total Count: {total_count}", ha='right', va='top', fontsize=10,
             transform=plt.gca().transAxes)

    plt.ylim(0, 1.05)
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.savefig(os.path.join(output_evaluate_sample_dir,"category_proportions_freq.png"))
    plt.show()

def process_sample(sample, cm_dir, output_evaluate_dir,excel_file):
    # sample = 'AH1-KI'
    # cm_dir = '/data5/wangxin/20241001_wcx'  # TODO: /data5/wangxin/20241001_wcx 后侧没有/
    # clean_data_dir= '/cleandata/20241227/'
    sample_data=get_sample_info_from_excel(excel_file,sample)
    primer_to_cut_plus20bp=int(sample_data["primer_to_cut_plus20bp"])
    primer_plus_5bp_sequence=sample_data["primer_plus_5bp"].upper()
    chr_seq = sample_data["chr_seq"]
    cut_site = int(sample_data["cut_site"])
    adapter_sequence ="ACCACGCGTGCTCTACA"
#############  This part code only needs to be run once.
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


    #fq_file = cm_dir + clean_data_dir + sample + '/flash/' + sample +".extendedFrags.fastq"
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
    output_dir = cm_dir  + sample + '/consolidate/'
    consolidated_fq_tri_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri.fq'
    consolidated_fq_tri_filtered_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_filtered.fq'
    filter_fastq_by_length(consolidated_fq_tri_file, consolidated_fq_tri_filtered_file, primer_to_cut_plus20bp)
    ####### step5: mapping from consolidated_tri_fq to ref genome (add transgene) using bwa
    output_dir = cm_dir  + sample + '/consolidate/'


    output_sam = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_change.sam'

    sample_without_number_name = group_genome_mapping[sample]
    genome_add_transgene = os.path.join("/data5/shuyu/result/modify_genome_db/", sample_without_number_name, f"only_{sample_without_number_name}_transgene.fa")

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


# ############# Classify for different situations
#     ######## step7: counts of different editing events (KI, small indels, large deletions, translocations, nojunctions, other editing)
    output_dir = cm_dir  + sample + '/consolidate/'
    ####### generate blastn db
    output_filter_and_reverse_bam = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_filter_and_reverse.bam'
    consolidated_fq_tri_filtered_file = output_dir + sample + '.aln_filtered_primeradd5bp_noadapter_consolidated_tri_filtered.fq'
    db= os.path.join("/data5/shuyu/result/modify_genome_db/", sample_without_number_name,f"{sample_without_number_name}_transgene_db")
    target_group=sample[0] # AP1-KI -> A
    KI,reverseKI,offtargetKI= get_edit_count(output_filter_and_reverse_bam, chr_seq, consolidated_fq_tri_filtered_file, output_dir, sample, db,target_group,cut_site)
    # call function
    # call save_dicts_to_excel function And pass in all the dictionaries
    output_evaluate_sample_dir = os.path.join(output_evaluate_dir,sample)
    cmd = "mkdir {}".format(output_evaluate_sample_dir)
    print(cmd)
    os.system(cmd)

    save_dicts_to_excel(
        os.path.join(output_evaluate_dir,sample,"output.xlsx"),
        output_evaluate_sample_dir,
        KI=KI,
        reverseKI=reverseKI,
        offtargetKI=offtargetKI
    )
