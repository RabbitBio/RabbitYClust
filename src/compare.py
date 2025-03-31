import argparse
from collections import defaultdict
from itertools import combinations

# 读取聚类文件，构建：序列ID -> 聚类编号 的映射
def parse_cluster_file(filepath):
    seq_to_cluster = {}
    with open(filepath) as f:
        cluster_id = None
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                cluster_id = int(line.split()[1])
            elif line:
                seq_id = line.split()[1].lstrip('>')  # 去掉前导的 >
                seq_to_cluster[seq_id] = cluster_id
    return seq_to_cluster

# 构造每个聚类内部的所有 pair
def cluster_pairs(seq_to_cluster):
    cluster_map = defaultdict(list)
    for seq_id, cid in seq_to_cluster.items():
        cluster_map[cid].append(seq_id)

    pairs = set()
    for seqs in cluster_map.values():
        if len(seqs) > 1:
            for a, b in combinations(seqs, 2):
                pairs.add(tuple(sorted((a, b))))
    return pairs

def main():
    parser = argparse.ArgumentParser(description="Precise comparison of two clustering results (pair-level, only flag truly mis-clustered sequences)")
    parser.add_argument("my_file", help="Your clustering result")
    parser.add_argument("cdhit_file", help="CD-HIT clustering result")
    args = parser.parse_args()

    # 读取两个聚类映射
    my_map = parse_cluster_file(args.my_file)
    cdhit_map = parse_cluster_file(args.cdhit_file)

    # 构建 pair 集合
    my_pairs = cluster_pairs(my_map)
    cdhit_pairs = cluster_pairs(cdhit_map)

    # 构建反向聚类（聚类编号 -> 序列集合）
    reverse_my_map = defaultdict(set)
    reverse_cdhit_map = defaultdict(set)
    for seq, cid in my_map.items():
        reverse_my_map[cid].add(seq)
    for seq, cid in cdhit_map.items():
        reverse_cdhit_map[cid].add(seq)

    # 初始化错聚序列集合
    only_in_my_seqs = set()
    only_in_cdhit_seqs = set()

    # CD-HIT 聚在一起，但你没聚的 -> 看是否有被孤立的
    for a, b in cdhit_pairs - my_pairs:
        if my_map.get(a) != my_map.get(b):
            if len(reverse_my_map.get(my_map.get(a), [])) == 1:
                only_in_cdhit_seqs.add(a)
            if len(reverse_my_map.get(my_map.get(b), [])) == 1:
                only_in_cdhit_seqs.add(b)

    # 你聚在一起，但 CD-HIT 没聚的 -> 看是否有被孤立的
    for a, b in my_pairs - cdhit_pairs:
        if cdhit_map.get(a) != cdhit_map.get(b):
            if len(reverse_cdhit_map.get(cdhit_map.get(a), [])) == 1:
                only_in_my_seqs.add(a)
            if len(reverse_cdhit_map.get(cdhit_map.get(b), [])) == 1:
                only_in_my_seqs.add(b)

    # 总序列数
    all_seqs = set(my_map.keys()) | set(cdhit_map.keys())
    total_seqs = len(all_seqs)

    # 输出结果
    print(f"Total sequences: {total_seqs}")
    print(f"Unique sequences only in your clusters: {len(only_in_my_seqs)} ({len(only_in_my_seqs)/total_seqs:.4f})")
    print(f"Unique sequences only in CD-HIT clusters: {len(only_in_cdhit_seqs)} ({len(only_in_cdhit_seqs)/total_seqs:.4f})")

    # 可选：输出错聚的序列ID到文件
    with open("wrong_in_my.txt", "w") as f:
        for seq in sorted(only_in_my_seqs):
            f.write(seq + "\n")

    with open("wrong_in_cdhit.txt", "w") as f:
        for seq in sorted(only_in_cdhit_seqs):
            f.write(seq + "\n")

    print("\nError sequence lists written to: wrong_in_my.txt and wrong_in_cdhit.txt")

if __name__ == "__main__":
    main()
