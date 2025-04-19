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
    # 初始化字典
    KI={}
    reverseKI={}
    offtargetKI={}
    ms_dic = {}  # 存储MS类型的字典

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
            if read.reference_name != "chr_transgene":
                continue
            if "AP" in sample:
                unwanted_sequence = "cggtggGAGCTGGACGGCGACGTAAACGGGGCACTTTTCGGGGAAATG"  # 替换为你的目标序列
                max_mismatches = 1  # 允许的最大错配数
                pattern = f"({unwanted_sequence}){{e<={max_mismatches}}}"  # 允许最大2个错配

                if re.search(pattern, read.query_sequence, re2.IGNORECASE):
                    print("plasmid")
                    continue  # 如果找到符合要求的匹配，跳过该 read
            cigar = str(read.cigarstring)
            letter = re.findall('\D', cigar)
            letters = "".join(re.findall(r'\D', cigar))
            if set(letters) != {'M', 'S'}:  # 只允许 M 和 S
                continue

            if not letters.startswith('M'):  # M 必须在 S 之前
                continue
            number = re.findall('\d+', cigar)
            number = list(map(int, number))  # 将数字部分转换为整数

            condition2 = read.query_sequence.count('N') <= len(read.query_sequence) * 0.05  # N的数量不超过5%

            if condition2:
                # 处理MS类型
                s_length = number[letter.index('S')]
                s_sequence = read.query_sequence[-s_length:]
                match_start_index = read.blocks[0][0]
                match_end_index = read.blocks[0][1]
                ms_dic[read.qname] = [read.query_sequence, s_sequence, match_start_index, match_end_index, cigar]

    print('MS类型比对数量:', len(ms_dic))
    # 将MS类型的S部分序列写入FASTA文件
    os.makedirs(os.path.join(output_dir, 'blastn'), exist_ok=True)
    if len(ms_dic)!=0:
        ms_s_fa = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_primeradd5bp_noadapter_consolidated_tri_MS_S.fa')
        with open(ms_s_fa, 'w') as f:
            for key in ms_dic.keys():
                s_seq = ms_dic[key][1]
                f.write(f'>{key}\n{s_seq}\n')

        # 定义blastn输出文件
        ms_s_fa_output = os.path.join(output_dir, 'blastn', f'{sample}.aln_filtered_noadapter_consolidated_tri_MS_S_fa_to_hg38fa_blastn_result.txt')

        # 调用blastn进行比对
        blastn(ms_s_fa,
               db, ms_s_fa_output)

        # 提取blastn比对结果中最佳对齐信息
        ms_main_chr_best_alignment_dict, ms_other_dict = extract_best_alignment_qend(ms_s_fa_output,ms_dic)

        # 分类删除类型（根据比对的结果进行分类）
        #ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
        KI,reverseKI,offtargetKI= categorize_ms(ms_main_chr_best_alignment_dict, ms_dic,chr_seq,target_group,cut_site)

    return  KI,reverseKI,offtargetKI

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


def categorize_ms(ms_main_chr_best_alignment_dict, ms_main_chr_dic,main_chr,target_group,cut_site):
    # ms_main_chr_best_alignment_dict blast seq
    # ms_dic[read.qname] = [read.query_sequence,s_sequence, match_start_index, match_end_index, cigar]
    KI={}
    reverseKI={}
    offtargetKI={}

    for key, alignment_info in ms_main_chr_best_alignment_dict.items():
        if key in ms_main_chr_dic:
            sseqid = alignment_info["sseqid"]  # 获取比对到的染色体编号
            start = ms_main_chr_best_alignment_dict[key]["sstart"] # -strand start>send
            send= ms_main_chr_best_alignment_dict[key]["send"]
            # 构造值列表
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
                offtargetKI[key] = value_list  # 认为是 KI 事件


    return KI, reverseKI, offtargetKI

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
                sequence, cigar,chr_seq,start,send = values[0], values[1],values[2],values[3],values[4]
                all_data.append([category, key, sequence, cigar,chr_seq,start,send])
                count += 1  # 更新当前类别的计数
        category_counts[category] = count  # 保存类别数量

    # 创建 DataFrame
    df = pd.DataFrame(all_data, columns=["Category", "Key", "Sequence", "CIGAR", "chr_seq","Start", "Send"])

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
    # TODO: 这里的比对基因组需要根据不同的样本进行调整


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
        KI=KI,
        reverseKI=reverseKI,
        offtargetKI=offtargetKI
    )
