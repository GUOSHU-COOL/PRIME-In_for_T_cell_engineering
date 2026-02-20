import os
import pandas as pd


def create_db(fa, db):
    """Generate the blastn reference database using makeblastdb"""
    cmd = f"makeblastdb -dbtype nucl -in {fa} -input_type fasta -parse_seqids -out {db}"
    os.system(cmd)
    print(f'Created BLAST database for {fa}, done!')


def process_transgenes(reference_content,excel_file, output_dir):
    """
    Traverse all samples
    1. Read the Excel file to obtain the transgene sequence
    2. Create a new.fa file for each sample in its folder, including the reference.fa file and the transgenic sequence
    3. Generate a BWA index for the.fa file of each sample
    """
    # read Excel file
    df = pd.read_excel(excel_file, usecols=[0, 1])  # Only read the first column and the second column
    data_dict = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))  # Create a dictionary from sample names to transgenic sequences

    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Traverse all samples
    for sample, transgene in data_dict.items():
        #have_run=["AP1-KI","AH1-KI"]
        # if sample in have_run:
        #     continue

        sample_dir = os.path.join(output_dir, sample)
        os.makedirs(sample_dir, exist_ok=True)  # Create a directory for each sample

        # The new path of the fa file
        new_fa = os.path.join(sample_dir, f"only_{sample}_transgene.fa")

        # Create a new.fa file and write the content of the reference.fa file and the transgenic sequence of this sample
        with open(new_fa, "w") as fa_file:
            #fa_file.write(reference_content)  # Write the content of the reference.fa file
            fa_file.write(f">chr_transgene\n{transgene}\n")  # Additional transgenic sequence
            print(f"Generated {new_fa}")

        # Generate a BWA index for the new.fa file
        cmd_bwa = f"bwa index {new_fa}"
        print(f"Running: {cmd_bwa}")
        os.system(cmd_bwa)
        # Generate a BLAST database for the new.fa file
        db_name = os.path.join(sample_dir, f"only_{sample}_transgene_db")
        create_db(new_fa, db_name)
    print("All processing completed.")


# example
excel_file = "./group_name_transgene_KI.xlsx"
output_dir = "./genome"
reference_fa = "./hg38_remove.fa"  # Replace it with the path of your reference.fa file

# Read the content of the reference.fa file
with open(reference_fa, "r") as ref_file:
    reference_content = ref_file.read()

process_transgenes(reference_content,excel_file, output_dir)
