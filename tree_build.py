# Author : Peter(周辉鑫)
# Time   : 2023/08/31  16:54

import os
import re
from Bio import SeqIO
from tqdm import tqdm  # 可视化模块
import subprocess  # 命令行调用的模块
from Bio import pairwise2


# 一、序列合并
# 合并后的序列地址， 可自己更改
original_sequences = '/tree_build/original_sequences/original_sequences.fasta'

# 读取并合并所有序列文件
def read_files_in_directory(path):
    len = 0
    # 遍历文件夹中的根目录，子文件夹以及所有文件
    for root, dirs, file in os.walk(path):
        for file_name in file:
            len += 1
            # 拼接文件路径
            file_path = os.path.join(root, file_name)
            # 检查文件类型， 这里以fasta文件为例
            if file_path.endswith("fasta"):  # 这里根据所下载的序列格式进行更改
                with open(file_path, "r") as f:
                    content = f.readlines()
                with open(original_sequences, "a+") as o:
                    o.writelines(content)

    print("成功合并%s个文件！" % len)


path = input("请输入需要合并的序列的所在文件夹地址：")
read_files_in_directory(path)

# 二、序列同名筛选
# 创建一个排除列表与一个接受列表
lst_remove = []
lst = []

# 读取数据
sequences = list(SeqIO.parse(original_sequences, 'fasta'))

# 可视化
for i in tqdm(sequences):
    lst_id = [s.id for s in sequences]
    if i.id not in lst_remove:  # 排除已筛选过的重名序列
        if lst_id.count(i.id) >= 2:  # 未筛选过的重名序列
            lst_remove.append(i.id)
            longest_sequences = max([s for s in sequences if s.id == i.id], key=lambda x: len(x.seq))
            lst.append(longest_sequences)
        else:  # 非重名序列
            lst.append(i)

nameselect_sequences = '/tree_build/name_selection/nameselect_sequences.fasta'
# 写入到新的fasta文件中， 可根据自己的路径进行修改
with open(nameselect_sequences, "w") as f:
    SeqIO.write(lst, nameselect_sequences, 'fasta')
    print("同名筛选完成！！")

# 三、CD-HIT去重

print("三、CD-HIT去重")
c = input("设置CD-HIT的阈值（0-1）：")
n = input("设置n值：")
cd_hit_sequences = '/tree_build/cd-hit/cd_hit_sequences.fasta'
# 定义要运行的命令
command = "cd-hit -i {0} -o cd_py_H1N1_original_data.fasta -c {1} -n {2}" .format(cd_hit_sequences, c, n)
# 使用subprocess运行命令
subprocess.run(command, shell=True)
print("CD-HIT去重运行成功！")


# 四、mafft序列比对

print("四、mafft序列比对")
mafft_sequences = '/tree_build/mafft/mafft_sequences.fasta'
# 定义要运行的命令
command = "mafft {0} > {1}" .format(cd_hit_sequences, mafft_sequences)
# 使用subprocess运行命令
subprocess.run(command, shell=True)
print("mafft序列比对运行成功！")

# 五、读取需要截取蛋白序列的模板文件并截取模板文件， 以H3N2的HA1为例子
template = "/tree_build/template/template.txt"  # 模板文件保存于此
with open(template, "r") as f:
    content = f.read()
    first_ten_characters = content[:10]  # 提取前十个作为对比
    last_ten_characters = content[-10:] # 提取后十个作为对比

# 读取FASTA文件
sequences = list(SeqIO.parse(mafft_sequences, "fasta"))

# 找到大多数序列中的gap位置
gap_position = None
threshold = len(sequences) * 0.9  # 大多数的阈值，可以根据需要进行调整

for i in range(len(sequences[0])):
    gap_count = sum(1 for seq in sequences if seq[i] == "-")
    if gap_count > threshold:
        gap_position = i
        break

# 删除所有序列的gap位置的碱基
sequences = [seq[:gap_position] + seq[gap_position + 1:] for seq in sequences]

# 存储匹配结果
matched_sequences = []

# 比对序列与模板，相似度大于95%的保留
for seq in sequences:
    alignment = pairwise2.align.localxx(template, seq, one_alignment_only=True)[0]
    similarity = (alignment[2] + alignment[3]) / len(template)
    if similarity >= 0.95:
        # 截取与模板长度相同的部分
        start = alignment[3]
        end = alignment[3] + len(template)
        matched_sequences.append(seq[start:end])

# matched_sequences 中包含了与模板长度相同的相似度大于95%的序列片段

cut_sequences = "/tree_build/cut_sequences/cut_sequences.fasta"
# 可以将结果保存到文件
with open(cut_sequences, "w") as output_file:
    for i, seq in enumerate(matched_sequences):
        SeqIO.write(seq, output_file, "fasta")

# 六、删除大于10%以上的gaps列
# 创建一个接收列表
lst = []
# 读取数据
sequences = list(SeqIO.parse(cut_sequences, 'fasta'))
# 可视化
for i in tqdm(sequences):
    lst_split = []
    if "-" in i.seq or "X" in i.seq:
        total = i.seq.count("X") + i.seq.count("-") # 计算未知碱基长度
        if total/len(i.seq) < 0.1:  # 计算未知碱基占总序列的长度比例
            lst.append(i)
    else:  # 不含未知碱基的序列
        lst.append(i)

delete_gaps_sequences = "/tree_build/delete_gaps_sequences/delete_gaps_sequences.fasta"
# 写入文件
with open(delete_gaps_sequences,"w"):
    SeqIO.write(lst, delete_gaps_sequences, 'fasta')

# 七、未知碱基填充
# 获取比对数据
with open(delete_gaps_sequences) as f:
    seqs_align = {}
    line = 1
    tem_name = tem_seq = ''
    while line:
        line = f.readline()
        if '>' in line:
            seqs_align[tem_name] = tem_seq
            tem_name = line.split(' | ')[0]+"_"+line.split(' | ')[1].split("/")[1]+"_"+line.split(' | ')[2]

            tem_seq = ''
        else:
            tem_seq += line.strip('\n')

    seqs_align[tem_name] = tem_seq
    seqs_align.pop('')

t = 0
seq_filled = {}  # 填充之后的序列

for i in tqdm(seqs_align.keys()):
    seq_to_fill = seqs_align[i]
    similar_seqs = {}  # 最相似的序列
    similar_site = list(range(20))  # 相似位点的数量
    for j in similar_site:
        similar_seqs[j] = ''

    # 遍历所有序列，寻找最相似的序列
    for seq in seqs_align.values():
        if seq != seq_to_fill:
            count = 0
            for aa in range(len(seq_to_fill)):
                if seq[aa] == seq_to_fill[aa]:
                    count += 1
            if count >= min(similar_site):
                index = similar_site.index(min(similar_site))
                similar_site[index] = count
                similar_seqs[index] = seq

    # 替换X和-
    tem_seq = ''
    for aa in range(len(seq_to_fill)):
        if seq_to_fill[aa] not in ['-', 'X']:
            tem_seq += seq_to_fill[aa]
        else:
            frequent_aa = []
            for seq in similar_seqs.values():
                if len(seq) > aa and seq[aa] not in ['-', 'X']:
                    frequent_aa.append(seq[aa])
            if frequent_aa:
                tem_seq += max(frequent_aa, key=frequent_aa.count)
            else:
                # 使用序列中出现频率最高的氨基酸替换
                all_aa = [seq[aa] for seq in seqs_align.values() if len(seq) > aa and seq[aa] not in ['-', 'X']]
                if all_aa:
                    tem_seq += max(set(all_aa), key=all_aa.count)
                else:
                    tem_seq += seq_to_fill[aa]

    seq_filled[i] = tem_seq or seq_to_fill  # 如果 tem_seq 为空，则直接使用原始序列

filled_sequences = "/tree_build/filled_sequences/filled_sequences.fasta"  # 新文件路径

if seq_filled:
    with open(filled_sequences, "w") as f:
        for seq_name, seq in seq_filled.items():
            f.write(f"{seq_name}\n{seq}\n")

    print("填充后的序列已写入文件:", filled_sequences)
else:
    print("无需填充，输出文件为空。")

# 八、fasttree建树
print("八、fasttree建树")
tree = '/tree_build/fasttree/rough_tree.nwk'
# 定义要运行的命令
command = "fasttree -lg {0} > {1}" .format(filled_sequences, tree)
# 使用subprocess运行命令
subprocess.run(command, shell=True)
print("fasttree建树成功！")

# 九、时间校准
fasta_file = filled_sequences
date_csv_file = "/tree_build/treetime/date.csv"


# 打开 FASTA 和 date.csv 文件
with open(fasta_file, "r") as f_fasta, open(date_csv_file, "w") as f_date:
    # 写入 date.csv 列标题
    f_date.write("name,date\n")

    # 解析 FASTA 文件
    sequences = SeqIO.parse(f_fasta, "fasta")

    # 遍历每个序列
    for sequence in sequences:
        sequence_name = sequence.id

        # 从序列描述中查找日期信息
        description = sequence.description
        match = re.search(r"\|\s*(\d{4}-\d{2}-\d{2})\s*\|", description)
        if match:
            sequence_date = match.group(1)
        else:
            # 如果无法从序列描述中提取日期，则跳过该序列
            print(f"无法识别序列 {sequence_name} 的日期，已跳过。")
            continue

        # 将序列名称和日期信息写入 date.csv 文件
        f_date.write(f"{sequence_name},{sequence_date}\n")

print("date.csv 文件已生成！")

final_tree = '/tree_build/final_tree/tree.nwk'

# 定义要运行的命令
command = "treetime --tree {0} –aln {1} --dates {2}" .format(tree, filled_sequences, date_csv_file)
# 使用subprocess运行命令
subprocess.run(command, shell=True)
print("时间校准成功！")

# 十、打开figtree
command = "figtree"
# 使用subprocess运行命令
subprocess.run(command, shell=True)




