import os
import pandas as pd


def create_db(fa, db):
    """使用 makeblastdb 生成 blastn 参考数据库"""
    cmd = f"makeblastdb -dbtype nucl -in {fa} -input_type fasta -parse_seqids -out {db}"
    os.system(cmd)
    print(f'Created BLAST database for {fa}, done!')


def process_transgenes(reference_content,excel_file, output_dir):
    """
    遍历所有样本：
    1. 读取 Excel 文件，获取 transgene 序列
    2. 为每个样本在其文件夹中创建新的 .fa 文件，包含参考 .fa 文件和转基因序列
    3. 为每个样本的 .fa 文件生成 BWA 索引
    """
    # 读取 Excel 文件
    df = pd.read_excel(excel_file, usecols=[0, 1])  # 只读取第一列和第二列
    data_dict = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))  # 创建样本名到转基因序列的字典

    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)

    # 遍历所有样本
    for sample, transgene in data_dict.items():
        #have_run=["AP1-KI","AH1-KI"]
        # if sample in have_run:
        #     continue

        sample_dir = os.path.join(output_dir, sample)
        os.makedirs(sample_dir, exist_ok=True)  # 为每个样本创建目录

        # 新的 fa 文件路径
        new_fa = os.path.join(sample_dir, f"only_{sample}_transgene.fa")

        # 创建新的 .fa 文件，并写入参考 .fa 文件内容和该样本的转基因序列
        with open(new_fa, "w") as fa_file:
            #fa_file.write(reference_content)  # 写入参考 .fa 文件内容
            fa_file.write(f">chr_transgene\n{transgene}\n")  # 追加转基因序列
            print(f"Generated {new_fa}")

        # 为新的 .fa 文件生成 BWA 索引
        cmd_bwa = f"bwa index {new_fa}"
        print(f"Running: {cmd_bwa}")
        os.system(cmd_bwa)
        # 为新的 .fa 文件生成 BLAST 数据库
        db_name = os.path.join(sample_dir, f"only_{sample}_transgene_db")
        create_db(new_fa, db_name)
    print("All processing completed.")


# 调用示例
excel_file = "./group_name_transgene_KI.xlsx"
output_dir = "./genome"
reference_fa = "./hg38_remove.fa"  # 替换为你的参考 .fa 文件路径

# 读取参考 .fa 文件内容
with open(reference_fa, "r") as ref_file:
    reference_content = ref_file.read()
process_transgenes(reference_content,excel_file, output_dir)