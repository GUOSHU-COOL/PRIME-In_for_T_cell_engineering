# encoding:utf-8

from __future__ import division
import HTSeq
import sys
import os
from collections import defaultdict
from collections import Counter
import pysam
import re
import pandas as pd
from pandas.core.sample import sample
import pickle
from collections import defaultdict
import re
from Bio import SeqIO

# umi_start = -12
# umi_end = 0
# umi_count_threshold = 2
def save_pickle(dict, file_path):
    # 保存 filtered_umi_counter 到 pickle 文件
    with open(file_path, 'wb') as f:
        pickle.dump(dict, f)
    print(f"Filtered UMI counter saved to {file_path}")
def load_pickle(file_path="filtered_umi_counter.pkl"):
    """从 pickle 文件读取 filtered_umi_counter"""
    with open(file_path, 'rb') as f:
        filtered_umi_counter = pickle.load(f)
    print(f"Filtered UMI counter loaded from {file_path}")
    return filtered_umi_counter

def filter_raw_bam(bam, primer_sequence, adapter_sequence,
                   output_bam):  # contains 20bp primer sequences and 20bp adapter sequences.

    i_bam = pysam.AlignmentFile(bam, 'rb')
    o_bam = pysam.AlignmentFile(output_bam, 'wb', template=i_bam)

    for read in i_bam:
        # print(read.query_sequence)
        condition1 = re.search(primer_sequence, str(read.query_sequence))
        # condition2 = re.search(adapter_sequence, str(read.query_sequence))
        # if condition1 and condition2:
        if condition1:
            # print(read.query_sequence)
            o_bam.write(read)
    print('get filtered bam, done!')

def check_umi(umi):
    """检查 UMI 是否长度为 12，且第 4 位和第 9 位都是 'A' 或 'T'"""
    if len(umi) == 12 and umi[3] in {'A', 'T'} and umi[8] in {'A', 'T'}:
        return True
    return False
def check_umi(umi):
    """检查 UMI 是否符合 NNNTNNNNTNNN 模式"""
    if len(umi) == 12 and umi[3] == 'T' and umi[8] == 'T':
        return True
    return False

def read_umis_fq(fastq_file, umi_count_threshold, primer_plus_5bp_sequence, adapter_sequence):
    umi_counter = defaultdict(int)
    umi_reads = defaultdict(list)  # 存储 UMI 对应的 read.query_sequence

    linker_sequence = "TGATAG"

    # 读取 FASTQ 文件中的每一条记录
    for record in SeqIO.parse(fastq_file, "fastq"):
        seq = str(record.seq)  # 获取序列
        if re.search(primer_plus_5bp_sequence, seq, re.IGNORECASE):  # ignore
            match = list(re.finditer(linker_sequence, seq))
            match1 = list(re.finditer(adapter_sequence, seq))
            if match1:
                match_end = match1[0].end()  # 第一个 match 位点

                umi_end = match_end + 12
                umi = seq[match_end: umi_end]  # 提取 UMI 序列

                if check_umi(umi):
                    umi_counter[umi] += 1
                    umi_reads[umi].append(record)  # 存储对应的 SeqRecord 对象
                else:
                    continue

    # 过滤 UMI 计数，生成 filtered_umi_counter
    filtered_umi_counter = {umi: count for umi, count in umi_counter.items() if count >= umi_count_threshold}

    # 仅保留 filtered_umi_counter 里的 UMI
    filtered_umi_reads = {umi: umi_reads[umi] for umi in filtered_umi_counter}

    print('get umi counts ' + str(len(umi_counter)) + ' done!')
    print('get filtered umi counts ' + str(len(filtered_umi_counter)) + ' done!')

    return filtered_umi_counter, filtered_umi_reads




def consolidate_position(bases, quals, min_qual, min_freq):
    length = len(list(bases))

    final_correct_seq = ''
    final_correct_qual = []
    num = {}
    qual = {}

    num['A'] = num['C'] = num['G'] = num['T'] = num['N'] = 0
    qual['A'] = qual['C'] = qual['G'] = qual['T'] = qual['N'] = 0

    zipped = zip(bases, quals)
    for i, (tuple1, tuple2) in enumerate(zipped):
        counter = Counter(tuple1)
        tuple1_length = len(tuple1)
        most_common_base = counter.most_common(1)[0][0]
        index = tuple1.index(most_common_base)
        s_qual = tuple2[index]
        s_num = counter[most_common_base]

        freq = round(float(s_num) / tuple1_length, 2)
        if freq > min_freq:
            x = True
            y = most_common_base
            z = s_qual
            # return True, most_common_base, qual[most_common_base]
        else:
            x = False
            y = 'N'
            z = 0
        final_correct_seq = final_correct_seq + y
        final_correct_qual.append(z)

    # print(final_correct_seq + '\t' + final_correct_qual)
    return final_correct_seq, final_correct_qual


def consolidate(min_qual, min_freq, filtered_umi_counter, filtered_umi_reads,primer_to_cut_plus20bp):
    num_input_reads = 0
    num_consolidated_reads = 0
    dic = {}

    for umi, umi_data in filtered_umi_reads.items():  # umi_data 现在是 SeqRecord 列表
        read_ids = [record.id for record in umi_data]  # 提取 read ID
        read_sequences = [str(record.seq) for record in umi_data]  # 提取序列
        read_qualities = [record.letter_annotations["phred_quality"] for record in umi_data]  # 提取质量分数

        num_input_reads += len(read_sequences)
        num_consolidated_reads += 1

        # 进行碱基合并
        read_bases_list = list(zip(*read_sequences))  # 读取所有序列的碱基
        read_quals_list = list(zip(*read_qualities))  # 读取所有序列的质量分数

        # 调用 consolidate_position 进行整合
        final_correct_seq, final_correct_qual = consolidate_position(read_bases_list, read_quals_list, min_qual,
                                                                     min_freq)

        # 计算 UMI 过滤后的计数
        molecular_id_length = filtered_umi_counter[umi]

        # 质量分数转换回 ASCII 格式
        quality = ''.join(chr(value + 33) for value in final_correct_qual)

        # 存储 UMI 处理后的结果
        dic[umi] = [final_correct_seq, quality, molecular_id_length, read_ids]  # 保留 read_id 信息

    print(f'get consolidated sequence, done! final umi count: {len(dic)}')
    return dic



def consolidate_bam(dic_bam, consolidated_fq_file):
    with open(consolidated_fq_file, 'w') as fq_out:
        for umi, values in dic_bam.items():
            final_correct_seq, quality, molecular_id_length, read_ids = values

            # 生成 read_id，保留 UMI 和 read 数量信息
            fa_read_id = f"{umi}_:{molecular_id_length}_:{len(read_ids)}"

            # 写入 FASTQ 格式
            fq_out.write(f"@{fa_read_id}\n")
            fq_out.write(f"{final_correct_seq}\n")
            fq_out.write("+\n")
            fq_out.write(f"{quality}\n")

        print(f'write consolidated sequences to {consolidated_fq_file}, done!')

    return


def get_sample_primer_and_adapter_sequence(samples_info_file):
    df = pd.read_csv(samples_info_file, sep='\t')
    # print(df)
    # sample  cut_site        chr_seq start   end     strand  primer  adapter loci    trangene20bp_seq        linker_add30bp  primer_add5bp_seq
    key = 'sample'
    primer = 'primer'
    adapter = 'adapter'
    strand = 'strand'
    cutsite = 'cut_site'
    chr_seq = 'chr_seq'
    start = 'start'
    end = 'end'
    transgene_20bp = 'trangene20bp_seq'
    linker_add30bp = 'linker_add30bp_seq'
    primer_add5bp = 'primer_add5bp_seq'
    downstream_genome20bp = 'downstream_genome20bp_seq'
    p_to_d20_distance = 'p_to_d20_distance'

    samples_primer_dic = df.set_index(key)[primer].to_dict()
    samples_adapter_dic = df.set_index(key)[adapter].to_dict()
    samples_strand_dic = df.set_index(key)[strand].to_dict()
    samples_cutsite_dic = df.set_index(key)[cutsite].to_dict()
    samples_chrseq_dic = df.set_index(key)[chr_seq].to_dict()
    samples_start_dic = df.set_index(key)[start].to_dict()
    samples_end_dic = df.set_index(key)[end].to_dict()
    samples_transgene_20bp_dic = df.set_index(key)[transgene_20bp].to_dict()
    samples_add30bp_dic = df.set_index(key)[linker_add30bp].to_dict()
    samples_primer_add5bp_dic = df.set_index(key)[primer_add5bp].to_dict()
    samples_downstream_genome20bp_dic = df.set_index(key)[downstream_genome20bp].to_dict()
    samples_p_to_d20_distance_dic = df.set_index(key)[p_to_d20_distance].to_dict()

    print('get samples" information, done!')
    return samples_primer_dic, samples_adapter_dic, \
        samples_strand_dic, samples_cutsite_dic, \
        samples_chrseq_dic, \
        samples_start_dic, samples_end_dic, \
        samples_transgene_20bp_dic, \
        samples_add30bp_dic, \
        samples_primer_add5bp_dic, \
        samples_downstream_genome20bp_dic, \
        samples_p_to_d20_distance_dic


def part_reverse_complement_based_on_strand(samples_seq_dic, samples_strand_dic):
    samples_seq_tiaozheng_dic = {}
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    for key in samples_seq_dic.keys():
        if key in samples_strand_dic.keys():

            dna = samples_seq_dic[key]
            if samples_strand_dic[key] == '+':
                primer_seq = dna
                samples_seq_tiaozheng_dic[key] = primer_seq
            else:
                reversed_dna = dna[::-1]
                rc_dna = ''
                for base in reversed_dna:
                    rc_dna += complement[base]
                    primer_seq = rc_dna
                samples_seq_tiaozheng_dic[key] = primer_seq

    return samples_seq_tiaozheng_dic


def read_excel(file_path):
    # 读取 Excel 文件
    df = pd.read_excel(file_path, dtype=str)  # 以字符串格式读取，防止数据丢失

    # 构建字典，第一列为 key，后两列为字典形式
    data_dict = {
        row["sample"]: {"primer_plus_5bp": row["primer_plus_5bp"],
                        "primer_to_cut_plus20bp": row["primer_to_cut_plus20bp"]}
        for _, row in df.iterrows()
    }

    return data_dict




def main():
    if len(sys.argv) < 4:
        print('Usage: python consolidate.py sample \
                                            fq_file \
                                            consolidated_fq_file \
                                            min_qual \
                                            min_freq \
                                            N_freq \
                                            umi_count_threshold \
                                            adapter_sequence\
                                            primer_to_cut_plus20bp\
                                            primer_plus_5bp_sequence\    ')
        sys.exit()


    sample, fq_file, consolidated_fq_file = sys.argv[1:4]
    min_qual = int(sys.argv[4])
    min_freq = float(sys.argv[5])
    N_freq = float(sys.argv[6])
    umi_count_threshold = int(sys.argv[7])
    adapter_sequence = sys.argv[8]
    primer_to_cut_plus20bp = int(sys.argv[9])
    primer_plus_5bp_sequence =sys.argv[10].upper()

    filtered_umi_counter, filtered_umi_reads = read_umis_fq(fq_file, \
                                                            umi_count_threshold, \
                                                            primer_plus_5bp_sequence, \
                                                            adapter_sequence)
    # filtered_umi_counter = load_pickle("/data5/wangxin/20241001_wcx/shuyu/temp/filtered_umi_counter.pkl")
    # filtered_umi_reads = load_pickle("/data5/wangxin/20241001_wcx/shuyu/temp/filtered_umi_reads.pkl")
    dic_bam = consolidate(min_qual, \
                          min_freq, \
                          filtered_umi_counter, \
                          filtered_umi_reads, \
                          primer_to_cut_plus20bp,\
                          )
    consolidate_bam(dic_bam, \
                    consolidated_fq_file, \
                    )


if __name__ == '__main__':
    main()
