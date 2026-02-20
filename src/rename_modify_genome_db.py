import os
import shutil

def copy_and_rename_folder(src_dir: str, dst_root: str, old_str: str, new_str: str):
    if not os.path.exists(src_dir):
        raise FileNotFoundError(f"Source directory not found: {src_dir}")

    for root, dirs, files in os.walk(src_dir):
        # Build the target path (and replace CP2 in the path)

        dst_path = os.path.join(dst_root,new_str)
        os.makedirs(dst_path, exist_ok=True)
        for file in files:
            src_file = os.path.join(root, file)
            new_file_name = file.replace(old_str, new_str)
            dst_file = os.path.join(dst_path, new_file_name)

            shutil.copy2(src_file, dst_file)
            print(f"Copied: {src_file} → {dst_file}")

# example
genome_directory = "/data5/wangxin/20241001_wcx/shuyu/genome/"
dst_directory = "/data5/shuyu/result/modify_genome_db/"
change_name_db=[["CP2","CCR5-P-T"],["CH2","CCR5-H-T"],["AP2","AAVS1-P-T"],["AH1","AAVS1-H-T"],["AP-KI","AAVS1-P-KI"],["AH-KI","AAVS1-H-KI"],["CP-KI","CCR5-P-KI"],["CH-KI","CCR5-H-KI"]]
for i in change_name_db:
    src_directory = os.path.join(genome_directory, i[0])
    # call function
    copy_and_rename_folder(src_directory, dst_directory, old_str=i[0], new_str=i[1])

