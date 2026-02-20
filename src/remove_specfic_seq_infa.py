import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio import SeqIO

# 使用示例
input_file = '/data5/wangxin/20241001_wcx/shuyu/genome/hg38.fa'  # 输入FASTA文件
#sequence_to_remove = 'GTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGGGAGGGGTCGGCAATTGAACCGGTGCCTAGAGAAGGTGGCGCGGGGTAAACTGGGAAAGTGATGTCGTGTACTGGCTCCGCCTTTTTCCCGAGGGTGGGGGAGAACCGTATATAAGTGCAGTAGTCGCCGTGAACGTTCTTTTTCGCAACGGGTTTGCCGCCAGAACACAGGTAAGTGCCGTGTGTGGTTCCCGCGGGCCTGGCCTCTTTACGGGTTATGGCCCTTGCGTGCCTTGAATTACTTCCACCTGGCTGCAGTACGTGATTCTTGATCCCGAGCTTCGGGTTGGAAGTGGGTGGGAGAGTTCGAGGCCTTGCGCTTAAGGAGCCCCTTCGCCTCGTGCTTGAGTTGAGGCCTGGCCTGGGCGCTGGGGCCGCCGCGTGCGAATCTGGTGGCACCTTCGCGCCTGTCTCGCTGCTTTCGATAAGTCTCTAGCCATTTAAAATTTTTGATGACCTGCTGCGACGCTTTTTTTCTGGCAAGATAGTCTTGTAAATGCGGGCCAAGATCTGCACACTGGTATTTCGGTTTTTGGGGCCGCGGGCGGCGACGGGGCCCGTGCGTCCCAGCGCACATGTTCGGCGAGGCGGGGCCTGCGAGCGCGGCCACCGAGAATCGGACGGGGGTAGTCTCAAGCTGGCCGGCCTGCTCTGGTGCCTGGCCTCGCGCCGCCGTGTATCGCCCCGCCCTGGGCGGCAAGGCTGGCCCGGTCGGCACCAGTTGCGTGAGCGGAAAGATGGCCGCTTCCCGGCCCTGCTGCAGGGAGCTCAAAATGGAGGACGCGGCGCTCGGGAGAGCGGGCGGGTGAGTCACCCACACAAAGGAAAAGGGCCTTTCCGTCCTCAGCCGTCGCTTCATGTGACTCCACGGAGTACCGGGCGCCGTCCAGGCACCTCGATTAGTTCTCGAGCTTTTGGAGTACGTCGTCTTTAGGTTGGGGGGAGGGGTTTTATGCGATGGAGTTTCCCCACACTGAGTGGGTGGAGACTGAAGTTAGGCCAGCTTGGCACTTGATGTAATTCTCCTTGGAATTTGCCCTTTTTGAGTTTGGATCTTGGTTCATTCTCAAGCCTCAGACAGTGGTTCAAAGTTTTTTTCTTCCATTTCAGGTGTCGTGA'  # 需要删除的序列片段
sequence_to_remove = "GTCAGTATCAATTCTGGAAGAATTTCCAGACA"

#fasta_file = "/data5/wangxin/20241001_wcx/shuyu/genome/EEF1A1_genome.fa"
#
# # 读取FASTA文件中的所有序列
# sequences = SeqIO.parse(fasta_file, "fasta")
#
# # 假设你只需要获取第一个序列作为查询序列
# sequence_to_remove = str(next(sequences).seq)
output_file = "/data5/wangxin/20241001_wcx/shuyu/genome/hg38_remove.fa"  # 输出的FASTA文件

def modify_sequence(original_sequence, start, end):
    """
    逐个读取碱基，并在指定范围（start 到 end）删除。
    """
    modified_sequence = []
    for i, base in enumerate(original_sequence, start=1):  # 1-based index
        if start <= i <= end:
            modified_sequence.append('')  # 替换为空，相当于删除
        else:
            modified_sequence.append(base)
    return ''.join(modified_sequence)
records = []  # 用于存储更新后的记录
for record in SeqIO.parse(input_file, "fasta"):
    if record.id == "chr6":  # 只处理chr6的基因组序列
        original_sequence = str(record.seq)  # 获取原始序列
        updated_sequence = original_sequence[:73520048 - 1] + original_sequence[73521223:] # blastn get indice
        if len(updated_sequence) != len(original_sequence):
            print(f"remove")
        print(len(updated_sequence), len(original_sequence))
        record.seq = Seq(updated_sequence)
        records.append(record)

    else:
        pass

        original_sequence = str(record.seq)
        updated_2_sequence =original_sequence.replace(sequence_to_remove, "")  # 删除特定序列
        print(len(updated_2_sequence), len(original_sequence))
        # 更新记录中的序列
        record.seq = Seq(updated_2_sequence)
        records.append(record)

        # 检查是否进行了更改
        if sequence_to_remove in original_sequence and sequence_to_remove not in updated_2_sequence:
            print(f"Sequence removed successfully from {record.id}")
        else:
            print(f"No change for sequence {record.id}")

# 写回更新后的FASTA文件
SeqIO.write(records, output_file, "fasta")

print(f"Updated file saved as {output_file}")



