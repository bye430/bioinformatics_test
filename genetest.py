from Bio import Entrez, SeqIO
import pandas as pd

# 设置 NCBI Entrez 的电子邮件和 API 密钥
Entrez.email = ""  # 使用您提供的实际邮箱
Entrez.api_key = ""  # 使用您提供的实际 NCBI API 密钥

def fetch_gene_sequence(organism, gene_name):
    """从 NCBI 数据库中查询特定物种和基因的序列"""
    query = f"{gene_name}[Gene] AND {organism}[Organism]"
    try:
        with Entrez.esearch(db="nucleotide", term=query, retmax=1) as handle:
            record = Entrez.read(handle)
        if record["IdList"]:
            gene_id = record["IdList"][0]
            with Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text") as handle:
                gene_seq_record = SeqIO.read(handle, "genbank")
            return gene_seq_record
        else:
            print(f"未找到 {organism} 的 {gene_name} 基因序列")
            return None
    except Exception as e:
        print(f"查询发生错误: {e}")
        return None

def find_fuzzy_matches(seq1, seq2):
    """查找模糊匹配并返回匹配的下标"""
    matches = []
    longest_match = ("", -1, -1)  # (匹配子序列, 西蓝花下标, 卷心菜下标)
    len1, len2 = len(seq1), len(seq2)
    
    # 逐个子序列进行匹配
    for start in range(len1):
        for end in range(start + 1, len1 + 1):
            subseq = seq1[start:end]
            if subseq in seq2:
                match_index = seq2.index(subseq)
                matches.append((start, match_index, subseq))
                
                # 更新最长匹配
                if len(subseq) > len(longest_match[0]):
                    longest_match = (subseq, start, match_index)
    
    return matches, longest_match

# 查询西蓝花和卷心菜的 rbcL 基因序列
broccoli_gene = fetch_gene_sequence("Brassica oleracea var. italica", "rbcL")
cabbage_gene = fetch_gene_sequence("Brassica oleracea var. capitata", "rbcL")

# 确保两者都被成功获取
if broccoli_gene and cabbage_gene:
    # 输出源序列和序列长度
    print(f"西蓝花基因序列：{broccoli_gene.seq}")
    print(f"卷心菜基因序列：{cabbage_gene.seq}")
    print(f"西蓝花序列长度：{len(broccoli_gene.seq)}")
    print(f"卷心菜序列长度：{len(cabbage_gene.seq)}")

    # 查找模糊匹配
    matches, longest_match = find_fuzzy_matches(str(broccoli_gene.seq), str(cabbage_gene.seq))

    # 输出匹配信息
    match_data = []
    for (broccoli_index, cabbage_index, matched_subseq) in matches:
        match_data.append({
            "西蓝花下标": broccoli_index,
            "卷心菜下标": cabbage_index,
            "匹配子序列": matched_subseq
        })

    match_df = pd.DataFrame(match_data)

    # 输出比对表格
    print(match_df)

    # 输出最长匹配信息
    longest_subseq, longest_broccoli_index, longest_cabbage_index = longest_match
    longest_length = len(longest_subseq)
    print(f"\n最长匹配子序列: {longest_subseq}")
    print(f"西蓝花下标: {longest_broccoli_index}, 卷心菜下标: {longest_cabbage_index}, 匹配长度: {longest_length}")

    # 计算匹配比例
    shorter_length = min(len(broccoli_gene.seq), len(cabbage_gene.seq))
    matching_ratio = longest_length / shorter_length if shorter_length > 0 else 0
    print(f"匹配比例: {matching_ratio:.2f}")

    # 根据匹配比例分析是否为同一个物种
    ratio_threshold = 0.8  # 设置比例阈值
    if matching_ratio >= ratio_threshold:
        print("这两个物种可能属于同一个物种。")
    else:
        print("这两个物种可能不属于同一个物种。")

    # 保存到 CSV 文件
    match_df.to_csv("fuzzy_matches.csv", index=False, encoding='utf-8-sig')
    print("\n模糊匹配结果已保存到 fuzzy_matches.csv 文件。")

else:
    print("无法完成基因比对。")

