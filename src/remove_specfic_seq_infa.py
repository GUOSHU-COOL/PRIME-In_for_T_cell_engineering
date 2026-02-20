import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio import SeqIO

# example
input_file = '/data5/wangxin/20241001_wcx/shuyu/genome/hg38.fa'  # 输入FASTA文件
#sequence_to_remove = 'GTGCCCGTCAGTGGGCAGAGCGCACATCGCCCACAGTCCCCGAGAAGTTGGGGGGAGGGGTCGGCAATTGAACCGGTGCCTAGAGAAGGTGGCGCGGGGTAAACTGGGAAAGTGATGTCGTGTACTGGCTCCGCCTTTTTCCCGAGGGTGGGGGAGAACCGTATATAAGTGCAGTAGTCGCCGTGAACGTTCTTTTTCGCAACGGGTTTGCCGCCAGAACACAGGTAAGTGCCGTGTGTGGTTCCCGCGGGCCTGGCCTCTTTACGGGTTATGGCCCTTGCGTGCCTTGAATTACTTCCACCTGGCTGCAGTACGTGATTCTTGATCCCGAGCTTCGGGTTGGAAGTGGGTGGGAGAGTTCGAGGCCTTGCGCTTAAGGAGCCCCTTCGCCTCGTGCTTGAGTTGAGGCCTGGCCTGGGCGCTGGGGCCGCCGCGTGCGAATCTGGTGGCACCTTCGCGCCTGTCTCGCTGCTTTCGATAAGTCTCTAGCCATTTAAAATTTTTGATGACCTGCTGCGACGCTTTTTTTCTGGCAAGATAGTCTTGTAAATGCGGGCCAAGATCTGCACACTGGTATTTCGGTTTTTGGGGCCGCGGGCGGCGACGGGGCCCGTGCGTCCCAGCGCACATGTTCGGCGAGGCGGGGCCTGCGAGCGCGGCCACCGAGAATCGGACGGGGGTAGTCTCAAGCTGGCCGGCCTGCTCTGGTGCCTGGCCTCGCGCCGCCGTGTATCGCCCCGCCCTGGGCGGCAAGGCTGGCCCGGTCGGCACCAGTTGCGTGAGCGGAAAGATGGCCGCTTCCCGGCCCTGCTGCAGGGAGCTCAAAATGGAGGACGCGGCGCTCGGGAGAGCGGGCGGGTGAGTCACCCACACAAAGGAAAAGGGCCTTTCCGTCCTCAGCCGTCGCTTCATGTGACTCCACGGAGTACCGGGCGCCGTCCAGGCACCTCGATTAGTTCTCGAGCTTTTGGAGTACGTCGTCTTTAGGTTGGGGGGAGGGGTTTTATGCGATGGAGTTTCCCCACACTGAGTGGGTGGAGACTGAAGTTAGGCCAGCTTGGCACTTGATGTAATTCTCCTTGGAATTTGCCCTTTTTGAGTTTGGATCTTGGTTCATTCTCAAGCCTCAGACAGTGGTTCAAAGTTTTTTTCTTCCATTTCAGGTGTCGTGA'  # 需要删除的序列片段
sequence_to_remove = "GTCAGTATCAATTCTGGAAGAATTTCCAGACA"

output_file = "/data5/wangxin/20241001_wcx/shuyu/genome/hg38_remove.fa"  # output fa file

def modify_sequence(original_sequence, start, end):
    """
    Read the bases one by one and delete them within the specified range (start to end).
    """
    modified_sequence = []
    for i, base in enumerate(original_sequence, start=1):  # 1-based index
        if start <= i <= end:
            modified_sequence.append('')  # Replace it with blank, which is equivalent to deletion
        else:
            modified_sequence.append(base)
    return ''.join(modified_sequence)
records = []  # Used for storing updated records
for record in SeqIO.parse(input_file, "fasta"):
    if record.id == "chr6":  # have known target chrom
        original_sequence = str(record.seq)  # Obtain the original sequence
        updated_sequence = original_sequence[:73520048 - 1] + original_sequence[73521223:] # blastn get indice
        if len(updated_sequence) != len(original_sequence):
            print(f"remove")
        print(len(updated_sequence), len(original_sequence))
        record.seq = Seq(updated_sequence)
        records.append(record)

    else:
        pass

        original_sequence = str(record.seq)
        updated_2_sequence =original_sequence.replace(sequence_to_remove, "")  # Delete a specific sequence
        print(len(updated_2_sequence), len(original_sequence))
        # Update the sequence in the record
        record.seq = Seq(updated_2_sequence)
        records.append(record)

        # Check if any changes have been made
        if sequence_to_remove in original_sequence and sequence_to_remove not in updated_2_sequence:
            print(f"Sequence removed successfully from {record.id}")
        else:
            print(f"No change for sequence {record.id}")

# Write back to the updated FASTA file
SeqIO.write(records, output_file, "fasta")

print(f"Updated file saved as {output_file}")




