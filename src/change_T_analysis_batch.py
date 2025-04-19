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
    # 初始化字典
    sm_main_chr_dic = {}  # 存储主染色体上的SM类型比对
    ms_dic = {}  # 存储MS类型的比对
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

    # 存储序列的字典
    raw_fa_dic = {}
    print("consolidated_fq_tri_filtered_file",consolidated_fq_tri_filtered_file)
    # 读取 FASTQ 文件
    with open(consolidated_fq_tri_filtered_file, "r") as raw_fq:
        for record in SeqIO.parse(raw_fq, "fastq"):
            if record.seq is None:
                print(f"empty序列: {record.id}")
            else:
                raw_fa_dic[record.id] = str(record.seq)  # 记录 ID 和序列

    # 读取BAM文件
    bam_file = pysam.AlignmentFile(bam, 'rb', check_sq=False)
    for read in bam_file:
        if read.qname in raw_fa_dic:
            cigar = str(read.cigarstring)
            letter = re.findall('\D', cigar)
            number = re.findall('\d+', cigar)
            number = list(map(int, number))  # 将数字部分转换为整数

            # 过滤包含H的CIGAR字符串
            if 'H' in letter:
                print("H in CIGAR")
                continue


            condition2 = read.query_sequence.count('N') <= len(read.query_sequence) * 0.05  # N的数量不超过5%

            if  condition2:
                # 判断CIGAR字符串类型
                if 'S' in letter and 'M' in letter:
                    if letter.count('S') == 1 and letter.count('M') == 1:
                        if letter.index('S') < letter.index('M'):
                            # 提取S部分的序列
                            s_length = number[letter.index('S')]
                            s_sequence = read.query_sequence[:s_length]
                            match_start_index = read.blocks[0][0] # 获取sam序列中match序列的起始 # 这个位点加一才是匹配上的 这个起始是0
                            match_end_index = read.blocks[0][1]

                            # 获取M部分的比对染色体
                            ref_name = read.reference_name  # 比对染色体名称

                            # 判断M部分的比对染色体
                            if ref_name == "chr_transgene":
                                # 如果是transgene，归类到KI
                                #sm_ki_dic[read.qname] = [read.query_sequence, match_start_index, match_end_index, cigar,ref_name]
                                sm_ki_dic[read.qname] = [read.query_sequence, cigar] # for unify
                            elif ref_name != chr_seq:
                                # 如果是其他染色体，归类到Translocation
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
                                # 如果是主染色体，归类到SM
                                #sm_main_chr_dic[read.qname] = [read.query_sequence, cigar]
                                sm_main_chr_dic[read.qname] = [read.query_sequence, s_sequence, match_start_index,match_end_index, cigar]
                        else:
                            # 处理MS类型
                            s_length = number[letter.index('S')]
                            s_sequence = read.query_sequence[-s_length:]
                            match_start_index = read.blocks[0][0]
                            match_end_index = read.blocks[0][1]
                            ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
                            #ms_dic[read.qname] = [read.query_sequence,cigar]
                elif 'M' in letter and ('D' in letter or 'I' in letter):
                    # 严格判断是否为 M...M 形式
                    if len(letter) >= 3 and letter[0] == 'M' and letter[-1] == 'M':
                        # 只要中间包含 D 或 I，就归入 small_indels_dic
                        if any(mid in letter[1:-1] for mid in ['D', 'I', 'M']):
                            small_indels_dic[read.qname] = [read.query_sequence, cigar]
                            match = re.findall(r'(\d+)D', cigar)
                            deletion_length = sum(map(int, match))
                            # 存储到字典
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
    print('SM类型比对数量:', len(sm_main_chr_dic))
    print('KI类型比对数量:', len(sm_ki_dic))
    print('Translocation类型比对数量:', len(sm_translocation_dic))
    print('Small Indels类型比对数量:', len(small_indels_dic))
    print('WT类型比对数量:', len(wt_dic))
    print('MS类型比对数量:', len(ms_dic))
    print('Other类型比对数量:', len(other_dic))

    # 将SM类型的S部分序列写入FASTA文件
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
        # 调用blastn进行比对（仅对主染色体上的SM类型）
    # 将MS类型的S部分序列写入FASTA文件
    os.makedirs(os.path.join(output_dir, 'blastn'), exist_ok=True)

    # 将MS类型的S部分序列写入FASTA文件
    if len(ms_dic)!=0:
        ms_s_fa = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_primeradd5bp_noadapter_consolidated_tri_MS_S.fa')
        with open(ms_s_fa, 'w') as f:
            for key in ms_dic.keys():
                s_seq = ms_dic[key][1]
                f.write(f'>{key}\n{s_seq}\n')

        # 定义blastn输出文件
        ms_s_fa_output = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_noadapter_consolidated_tri_MS_S_fa_to_hg38fa_blastn_result.txt')

        # 调用blastn进行比对（仅对主染色体上的MS类型）
        blastn(ms_s_fa, db, ms_s_fa_output)

        # 提取blastn比对结果中最佳对齐信息
        ms_main_chr_best_alignment_dict, ms_other_dict = extract_best_alignment_qend(ms_s_fa_output,ms_dic)

        # 分类删除类型（根据比对的结果进行分类）
        #ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
        small_ms_deletion, medium_ms_deletion, large_ms_deletion, ms_translocation_dic,ms_ki_dic,ms_translocation_type= (
            categorize_ms(ms_main_chr_best_alignment_dict, ms_dic,chr_seq,deletion_type_length,strand))

    return sm_ki_dic, sm_translocation_dic,small_sm_deletion, medium_sm_deletion, large_sm_deletion,small_indels_dic, wt_dic,small_ms_deletion, medium_ms_deletion, large_ms_deletion, ms_ki_dic, ms_translocation_dic, other_dic,ms_other_dict,ms_translocation_type,sm_translocation_type,deletion_type_length

def extract_best_alignment_qstart(blastn_result):
    """
    提取BLASTN比对结果中每个查询序列的最佳比对（根据 qstart == 0）。

    参数：
    - blastn_result: BLASTN 比对结果文件路径

    返回：
    - best_alignment_dict: 存储每个查询序列的最佳比对信息
    """
    # 读取BLASTN比对结果
    try:
        df = pd.read_csv(blastn_result, delimiter='\t', header=None)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        # 按 qseqid 分组
        grouped = df.groupby('qseqid')

        # 存储结果的字典
        best_alignment_dict = {}

        # 遍历每组，提取 qstart == 0 的行
        for name, group in grouped:
            group = group[group['qstart'] == 1]  # 仅保留 qstart == 1 的行

            if not group.empty:
                # 按照 序列比对百分比最大值的行
                max_row = group.loc[group['pident'].idxmax()]

                # 解析所需字段
                qseqid = max_row['qseqid']  # 作为 key
                s_sseqid = max_row['sseqid']
                s_sstart = int(max_row['sstart'])
                s_send = int(max_row['send'])

                # 存入字典
                best_alignment_dict[qseqid] = {
                    'sseqid': s_sseqid,
                    'sstart': s_sstart,
                    'send': s_send
                }

        return best_alignment_dict
    except Exception as e:
        print(f"读取文件 {blastn_result} 时发生错误: {e}")
        return {}

def extract_best_alignment_qend(blastn_result,ms_dic):
    """
    提取BLASTN比对结果中每个查询序列的最佳比对（根据 qend == s_length）。
    如果没有匹配上 qend == s_length 的比对，则将其放入 'other' 字典中。

    参数：
    - blastn_result: BLASTN 比对结果文件路径
    - s_length: 用于比较的序列长度，筛选 qend == s_length 的比对

    返回：
    - best_alignment_dict: 存储每个查询序列的最佳比对信息
    - other_dict: 存储未匹配 qend == s_length 的查询序列
    """
    # 读取BLASTN比对结果
    try:
        df = pd.read_csv(blastn_result, delimiter='\t', header=None)
        df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                      'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

        # 按 qseqid 分组
        grouped = df.groupby('qseqid')

        # 存储结果的字典
        best_alignment_dict = {}
        other_dict = {}

        # 遍历每组，提取 qend == s_length 的行
        for name, group in grouped:
            group_id=group['qseqid'].iloc[0]
            s_length=len(ms_dic[group_id][1])
            origin_seq=ms_dic[group_id][0]
            cigar=ms_dic[group_id][4]
            qseqid = group['qseqid'].iloc[0]
            group = group[group['qend'] == s_length]  # 所有qseqid都是一个样

            if not group.empty:
                # 按照 序列比对百分比 最大值的行
                max_row = group.loc[group['pident'].idxmax()]

                # 解析所需字段
                qseqid = max_row['qseqid']  # 作为 key
                s_sseqid = max_row['sseqid']
                s_sstart = int(max_row['sstart'])
                s_send = int(max_row['send'])

                # 存入字典
                best_alignment_dict[qseqid] = {
                    'sseqid': s_sseqid,
                    'sstart': s_sstart,
                    'send': s_send
                }
            else:
                # 如果没有匹配 qend == s_length，则放入 'other' 字典
                other_dict[qseqid] = [origin_seq,cigar]

        return best_alignment_dict, other_dict
    except Exception as e:
        print(f"读取文件 {blastn_result} 时发生错误: {e}")
        return {}, {}

def categorize_sm(sm_main_chr_best_alignment_dict, sm_main_chr_dic,deletion_type_length,strand):
    # sm_main_chr_best_alignment_dict blast seq
    small_deletion = {}
    medium_deletion = {}
    large_deletion = {}

    for key in sm_main_chr_best_alignment_dict:
        if key in sm_main_chr_dic:
            # 计算缺失大小
            match_start = sm_main_chr_dic[key][2] #both strand+ strand- and flag 0&16 start<end
            match_end = sm_main_chr_dic[key][3]
            start = sm_main_chr_best_alignment_dict[key]["sstart"] # -strand start>send
            send= sm_main_chr_best_alignment_dict[key]["send"]
            if strand == "-":
                deletion_size = abs(send - match_end)
            else:
                deletion_size = abs(match_start - send)
            # 生成所有可能的坐标对
            # sm_main_chr_dic[read.qname] = [read.query_sequence, cigar]
            #sm_main_chr_dic[read.qname] = [read.query_sequence, s_sequence, match_start_index, match_end_index, cigar]
            # 构造值列表
            value_list = [sm_main_chr_dic[key][0], sm_main_chr_dic[key][4]]

            # 分类存入相应字典
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
                # 使用logging记录到日志文件
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
            sseqid = alignment_info["sseqid"]  # 获取比对到的染色体编号
            sstart = alignment_info["sstart"]
            send = alignment_info["send"]
            # 构造值列表
            value_list = [ms_main_chr_dic[key][0], ms_main_chr_dic[key][4]]

            if sseqid == main_chr:
                # match_end- sstart
                # 计算缺失大小
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
                ki_events[key] = value_list  # 认为是 KI 事件

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
    从样本名中提取第一个数字之前的部分
    例如 'AH2-OT1' -> 'AH2'
    """
    match = re.match(r"^[^\d]*\d", sample_name)  # 匹配从头开始到第一个数字之前的字符
    if match:
        return match.group(0)
    else:
        return sample_name  # 如果没有找到匹配，返回原样本名


def add_columns_from_dicts(excel_file, deletion_type_length_type, ms_translocation_type, sm_translocation_type,
                           output_file):
    # 读取 Excel 文件
    df = pd.read_excel(excel_file, index_col=1)  # 使用第二列作为索引

    # 处理 deletion_type_length_type 字典，添加 deletion_type 和 deletion_length 列
    df['deletion_type'] = df.index.map(lambda key: deletion_type_length_type.get(key, {}).get('deletion_type', None))
    df['deletion_length'] = df.index.map(
        lambda key: deletion_type_length_type.get(key, {}).get('deletion_length', None))

    # 处理 ms_translocation_type 字典，添加 translocation_chr 列
    df['main_chr'] = df.index.map(lambda key: ms_translocation_type.get(key, {}).get('chr_seq', None))
    df['translocation_chr'] = df.index.map(lambda key: ms_translocation_type.get(key, {}).get('ref_name', None))
    df['reference_start'] = df.index.map(lambda key: ms_translocation_type.get(key, {}).get('reference_start', None))
    df['reference_end'] = df.index.map(lambda key: ms_translocation_type.get(key, {}).get('reference_end', None))
    # 处理 sm_translocation_type 字典，并确保转换为 Series
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

    # 保存为新的 Excel 文件
    df.to_excel(output_file)
    print(f"✅ 数据已保存至 {output_file}")

def save_dicts_to_excel(output_file,output_evaluate_sample_dir, **dicts):
    """
    将多个字典保存到一个 Excel 文件中，每个字典的名字作为 "Category"（分类）。

    参数：
    - output_file: str, 输出 Excel 文件路径
    - **dicts: 关键字参数，每个字典的 key 是序列 ID，value 是 [sequence, cigar]
    """
    all_data = []  # 存储所有字典的数据
    category_counts = {}  # 用于统计每个类别的数量

    # 遍历字典
    for category, data_dict in dicts.items():
        count = 0  # 统计当前类别的数量
        for key, values in data_dict.items():
            if len(values) >= 2:  # 确保有足够的数据
                sequence, cigar = values[0], values[1]
                all_data.append([category, key, sequence, cigar])
                count += 1  # 更新当前类别的计数
        category_counts[category] = count  # 保存类别数量

    # 创建 DataFrame
    df = pd.DataFrame(all_data, columns=["Category", "Key", "Sequence", "CIGAR"])

    # 保存为 Excel
    df.to_excel(output_file, index=False)
    print(f"✅ 数据已保存至 {output_file}")

    # 打印每个类别的数量
    print("📊 各类别数量统计：")
    for category, count in category_counts.items():
        print(f"  {category}: {count} 条数据")
    # 计算类别总数
    total_count = sum(category_counts.values())

    # 计算每个类别的比例
    category_proportions = {category: count / total_count for category, count in category_counts.items()}

    # 打印每个类别的比例
    print("📊 各类别比例统计：")
    for category, proportion in category_proportions.items():
        print(f"  {category}: {proportion * 100:.2f}%")
    # 绘制堆叠柱状图
    categories = list(category_proportions.keys())
    values = list(category_proportions.values())
    # 提取 sample
    sample_name = os.path.basename(output_evaluate_sample_dir)
    # 绘制堆叠柱状图
    plt.bar(categories, values, color='skyblue')

    # 设置图表标题和标签
    plt.title(f"{sample_name} Proportion of Each Category in the Total Data", fontsize=16)
    plt.xlabel("Categories", fontsize=14)
    plt.ylabel("Proportion", fontsize=14)

    # 将 x 轴标签旋转 45 度，避免堆叠
    plt.xticks(rotation=45, ha='right', fontsize=12)  # ha='right' 使标签右对齐，避免偏移

    # 添加比例文本标签
    for i, (category, proportion) in enumerate(category_proportions.items()):
        plt.text(i, proportion + 0.01, f"{proportion * 100:.2f}%", ha='center', va='bottom', fontsize=6)
    # 在右上角添加类别总数
    plt.text(0.95, 0.95, f"Total Count: {total_count}", ha='right', va='top', fontsize=10,
             transform=plt.gca().transAxes)

    # 设置 y 轴上限以提供更多空间
    plt.ylim(0, 1.05)
    # 调整布局，避免标签被裁剪
    plt.tight_layout()
    # 如果标签仍然超出，增加调整边距
    plt.subplots_adjust(top=0.85)
    # 保存图像
    plt.savefig(os.path.join(output_evaluate_sample_dir,"category_proportions_freq.png"))
    # 显示图表
    plt.show()


def plot_deletion_length_distribution(excel_file, output_evaluate_sample_dir):
    df = pd.read_excel(excel_file)

    # 转换为数值，去除 NaN
    deletion_lengths = pd.to_numeric(df['deletion_length'], errors='coerce').dropna()

    # 过滤掉特别小或极端的值，比如小于1的值
    deletion_lengths = deletion_lengths[deletion_lengths > 0]

    # 对数据取 log10 变换，避免大数值拉伸
    log_deletion_lengths = np.log10(deletion_lengths)

    # 提取 sample
    sample_name = os.path.basename(output_evaluate_sample_dir)

    # 1. 绘制对数 x 轴 + 对数 y 轴 直方图
    plt.figure(figsize=(10, 6))
    plt.hist(log_deletion_lengths, bins=50, color='skyblue', edgecolor='black')
    plt.yscale("log")  # 对数 y 轴
    plt.title(f"{sample_name} Deletion Length Frequency Distribution (Log X)", fontsize=16)
    plt.xlabel("Log10(Deletion Length)", fontsize=14)
    plt.ylabel("Frequency (log scale)", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_evaluate_sample_dir, "deletion_length_distribution_freq_logx.png"))
    plt.show()
    # 过滤数据：只保留 deletion_length ≤ 1000 的数据点
    filtered_deletion_lengths = deletion_lengths[deletion_lengths <= 1000]
    plt.figure(figsize=(10, 6))
    plt.hist(filtered_deletion_lengths, bins=50, color='salmon', edgecolor='black')
    plt.yscale("log")  # 对数 y 轴
    plt.xlim([0, 1000])  # 限制 X 轴范围在 1k 内
    plt.title(f"{sample_name} Deletion Length Frequency Distribution (≤ 1k)", fontsize=16)
    plt.xlabel("Deletion Length", fontsize=14)
    plt.ylabel("Frequency (log scale)", fontsize=14)

    # 设置 X 轴刻度，每隔 100 标一个
    plt.xticks(np.arange(0, 1001, 100))

    plt.tight_layout()
    plt.savefig(os.path.join(output_evaluate_sample_dir, "deletion_length_distribution_freq_1k_logy.png"))
    plt.show()
def read_sample_primer_distance_excel(file_path):
    # 读取 Excel 文件
    df = pd.read_excel(file_path, dtype=str)  # 以字符串格式读取，防止数据丢失

    # 构建字典，第一列为 key，后两列为字典形式
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
    # TODO: 这里的比对基因组需要根据不同的样本进行调整


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
    # 调用函数
    # 调用 save_dicts_to_excel 函数并将所有字典传入
    # 调用 save_dicts_to_excel 函数并将每个字典传入对应的位置

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
    # 调用函数并传入 Excel 文件路径
    excel_file = os.path.join(output_evaluate_dir,sample,"add_output.xlsx")
    plot_deletion_length_distribution(excel_file,output_evaluate_sample_dir)