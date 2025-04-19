import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
from collections import Counter
import re
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
def process_sample_blastn(sample_raw_name,sample, output_dir, target_seq, rna_number):
    if sample_raw_name == "Csn-KI":
        sample_raw_name = "CH3-KI"
    os.makedirs(output_dir, exist_ok=True)

    file_path = os.path.join("/data5/wangxin/20241001_wcx/shuyu/20250307/", sample, "output.xlsx")
    df = pd.read_excel(file_path)

    filtered_df = df[df.iloc[:, 0] == 'offtargetKI']

    result_dict = {}
    for _, row in filtered_df.iterrows():
        chr_seq = row["chr_seq"]
        start, send = row["Start"], row["Send"]
        result_dict.setdefault(chr_seq, []).append({"start": start, "send": send})# for same start and send

    #genome_file = "/data5/wangxin/20241001_wcx/shuyu/genome/hg38_remove.fa"
    sample_name= extract_prefix(sample_raw_name) # for example 'AH2-OT1' -> 'AH2'
    genome_file = os.path.join("/data5/wangxin/20241001_wcx/shuyu/genome/", sample_name, f"{sample_name}_transgene.fa")

    # 找出 `filtered_df` 里 `chr_seq` 但 `genome_file` 里没有的染色体
    df_chr_set = set(filtered_df["chr_seq"])
    fasta_chr_set = set(record.id for record in SeqIO.parse(genome_file, "fasta"))
    missing_chr = df_chr_set - fasta_chr_set

    print(f"这些 chr_seq 在 filtered_df 里有，但 genome_file 里没有: {missing_chr}")

    sequences = []
    empty_count = 0
    for record in SeqIO.parse(genome_file, "fasta"):
        if record.id in result_dict:
            for entry in result_dict[record.id]:
                start, send=entry["start"], entry["send"]
                min_pos = max(0, min(start, send) - 100)
                max_pos = max(start, send) + 100
                seq_length = len(record.seq)
                if min_pos >= seq_length or max_pos >= seq_length:
                    print(f"染色体 {record.id}: start={start}, send={send} 超出范围，染色体长度={seq_length}")
                extracted_seq = record.seq[max(0, send - 1 - 100): start + 100] if start > send else record.seq[max(0,
                                                                                                                    start - 1 - 100): send + 100]
                if len(extracted_seq) == 0:
                    empty_count += 1
                    print(f"空序列: {record.id}, start={start}, send={send}")
                sequences.append(SeqRecord(Seq(str(extracted_seq)), id=f"{record.id}_{start}_{send}", description=""))
    print(f"空序列数量: {empty_count}")
    extracted_fasta = os.path.join(output_dir, f"extracted_sequences_{rna_number}.fasta")
    SeqIO.write(sequences, extracted_fasta, "fasta")
    sequences = list(SeqIO.parse(extracted_fasta, "fasta"))
    print(f"Filtered DataFrame 行数: {len(filtered_df)}")
    print(f"提取的 sequences 数量: {len(sequences)}")
    window_size = 23
    target_seq = target_seq.upper()
    match_results = []
    matched_ids = set()
    window_seq_list = []
    for seq_record in sequences:
        seq = str(seq_record.seq).upper()

        # 遍历每个窗口进行比对
        for i in range(len(seq) - window_size + 1):
            window_seq = seq[i:i + window_size]
            # 计算窗口序列与目标序列的匹配情况
            if window_seq.endswith("GG"):
                matches = sum(1 for a, b in zip(window_seq[:-3], target_seq[:-3]) if a == b)
                mismatch_count = window_size - matches-3  # 错配数

                if mismatch_count <6:  # 允许错配数小于6
                    print("match")
                    match_results.append({
                        "sseqid": seq_record.id,
                        "target_seq": target_seq,
                        "window_seq": window_seq,
                        "Identity (%)": (matches / window_size) * 100,
                        "Mismatch": mismatch_count,
                        "Start": i,
                        "End": i + window_size
                    })
                    matched_ids.add(seq_record.id)
                    window_seq_list.append(window_seq)

            elif window_seq.startswith("CC"):
                reverse_complement_seq = str(Seq(window_seq).reverse_complement())

                reverse_matches = sum(1 for a, b in zip(reverse_complement_seq[3:], target_seq[3:]) if a == b)
                reverse_mismatch_count = window_size - reverse_matches - 3  # 错配数（排除结尾的 "GG"）

                if reverse_mismatch_count <6:  # 允许错配数小于6
                    print("reverse match")
                    match_results.append({
                        "sseqid": seq_record.id,
                        "target_seq": target_seq,
                        "window_seq": reverse_complement_seq,
                        "Identity (%)": (reverse_matches / window_size) * 100,
                        "Mismatch": reverse_mismatch_count,
                        "Start": i,
                        "End": i + window_size
                    })
                    matched_ids.add(seq_record.id)
                    window_seq_list.append(reverse_complement_seq)
            else:
                pass




    filtered_output = os.path.join(output_dir, f"filtered_results_{rna_number}.json")
    with open(filtered_output, "w") as f:
        for result in match_results:
            f.write(json.dumps(result) + "\n")

    total_sequences = len(sequences)
    sgRNA_caused = len(matched_ids)
    random = total_sequences - sgRNA_caused
    all_UMIs = df.shape[0]
    summary_file = os.path.join(output_dir, "summary.xlsx")
    pd.DataFrame({"Sample": [sample], "sgRNA_caused": [sgRNA_caused], "random": [random],"all_UMIs":[all_UMIs]}).to_excel(summary_file,
                                                                                                    index=False)
    # 统计 window_seq_rc 出现次数
    rc_counts = Counter(window_seq_list)

    # 创建 DataFrame 并输出到 Excel
    rc_df = pd.DataFrame(rc_counts.items(), columns=["window_seq", "Count"])
    rc_df["Normalized Count"] = (rc_df["Count"] / all_UMIs) * 1000  # 归一化
    excel_output_path = os.path.join(output_dir, f"window_seq_counts_{rna_number}.xlsx")
    rc_df.to_excel(excel_output_path, index=False)
    print(f"window_seq_rc 统计已保存到: {excel_output_path}")

    # 仅筛选第一列为 'offtargetKI' 的行
    filtered_df = df[df.iloc[:, 0] == 'offtargetKI']

    # 仅对 filtered_df 进行匹配，并将结果填充回原 DataFrame
    df.loc[df.iloc[:, 0] == 'offtargetKI', "Category"] = filtered_df.apply(
        lambda row: "sgRNA_caused" if f"{row['chr_seq']}_{row['Start']}_{row['Send']}" in matched_ids else "random",
        axis=1
    )
    new_file_path = os.path.join(output_dir,"new_output.xlsx")
    df.to_excel(new_file_path, index=False)

    return summary_file


file_path = "/data5/wangxin/20241001_wcx/shuyu/20250307/sgRNA.xlsx"
df = pd.read_excel(file_path)
filtered_df = df[df.iloc[:, 0].astype(str).str.contains("KI")]

for _, row in filtered_df.iterrows():
    sample_raw_name = row["sample"]
    sample = row["Rename"]
    sgRNAs = [row["sgRNA1"], row["sgRNA2"], row["sgRNA3"]]
    for i, sgRNA in enumerate(sgRNAs, start=1):
        print(sample,i)
        if pd.notna(sgRNA):
            process_sample_blastn(sample_raw_name,sample, os.path.join("/data5/wangxin/20241001_wcx/shuyu/sgRNA/", sample, str(i)),
                                  sgRNA, i)
