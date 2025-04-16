import re
# 读cdhit的结果
def read_clstr(filepath):
#    clusters_dict = {}
    reversed_clusters_dict = {}
    key = None
    current_cluster = []
#    pattern = re.compile(r'>(\w+[\S.]*)\b')
    pattern = re.compile(r'>(.*?)\.\.\.')
    with open(filepath,'r') as file:
        for line in file.readlines():
            line = line.strip()
            if line.startswith('>Cluster'):
                # 保存当前 cluster 的数据
                if current_cluster is not None:
#                    clusters_dict[key] = current_cluster
                    for seq in current_cluster:
                        reversed_clusters_dict[seq] = key
                current_cluster = []
                key = None
            else:
                match = pattern.search(line)
                if match:
                    if line.endswith('*'):
                        key = match.group(1)
                        reversed_clusters_dict[key] = key
                    else:
                        current_cluster.append(match.group(1))
    if current_cluster is not None:
#        clusters_dict[key] = current_cluster
        for seq in current_cluster:
            reversed_clusters_dict[seq] = key
    print(f'read {len(reversed_clusters_dict)} sequences from cd-hit', file=sys.stderr)
#    return clusters_dict, reversed_clusters_dict
    return reversed_clusters_dict

# 读取yclust分组结果
from collections import defaultdict
def read_yclust_res(filepath):
    reversed_result = {}
    result = defaultdict(list)  # 自动初始化空列表作为值
    with open(filepath, 'r') as file:
        for line in file.readlines():
            line = line.strip().split()
            if line:
                reversed_result[line[0][1:]] = line[1][1:]
#                result[line[1][1:].append(line[0][1:])  # 自动处理键不存在的情况
    print(f'read {len(reversed_result)} sequences from yclust results', file=sys.stderr)
#    return dict(result),reversed_result 
    return reversed_result

from collections import defaultdict
def parse_yclust_res(filepath):
    result = defaultdict(list)  # 自动初始化空列表作为值
    with open(filepath, 'r') as file:
        for line in file.readlines():
            line = line.strip().split()
            if line:
                result[line[1][1:]].append(line[0][1:])  # 自动处理键不存在的情况
    parse_clust_res(result)

def parse_clust_res(result):
    print(f'read {len(dict(result))} clusters from yclust results')
    lt1w=0
    lt2k=0
    lt1k=0
    lt100=0
    lt10=0
    lt1=0
    only1=0
    for rep, seqs in dict(result).items():
        length = len(seqs)
        if length > 10000:
            lt1w+=1
        elif length > 2000:
            lt2k+=1
        elif length > 1000:
            lt1k+=1
        elif length > 100:
            lt100+=1
        elif length > 10:
            lt10+=1
        elif length > 1:
            lt1+=1
        else:
            only1+=1
    print(f'large than 10000 {lt1w}')
    print(f'large than 2000 {lt2k}')
    print(f'large than 1000 {lt1k}')
    print(f'large than 100 {lt100}')
    print(f'large than 10 {lt10}')
    print(f'large than 1 {lt1}')
    print(f'only 1 {only1}')

from collections import defaultdict
# 求被拆开分到不同组的类的个数
def cross_group_consistency(dictA, dictB):
    commons = set(dictA.keys()).intersection(dictB.keys())
    print(f'common seqs: {len(commons)}')
    
    clusterA = [dictA[seq] for seq in commons]
    clusterB = [dictB[seq] for seq in commons]
 
    cluster_to_groups = defaultdict(set)
    for group, cluster in zip(clusterA, clusterB):
        cluster_to_groups[cluster].add(group)
    
    inconsistent_clusters = sum(1 for groups in cluster_to_groups.values() if len(groups) > 1)
    total_clusters = len(cluster_to_groups)
    print(f'sepearted cluster numbers: {inconsistent_clusters}')
    return inconsistent_clusters / total_clusters

def identify_disjoint_elements(dictA, dictB):
    count = 0
    exclusive_elements = 0
    seperater_ancestor = set()
    for seq, rep in dictA.items():
        if not dictB.get(seq):
            exclusive_elements += 1
            continue
        if dictB[seq] != dictB[rep]:
            count += 1
            seperater_ancestor.add(rep)
    print(f'{count} sequences seperated from its ancestor', file=sys.stderr)
    print(f'{len(seperater_ancestor)} ancestor seperated from its childs', file=sys.stderr)

def find_different_groups(dictA, dictB):
    different_res = 0
    commons = set(dictA.keys()).intersection(dictB.keys())
    print(f'common seqs: {len(commons)}')
    for seq, rep in dictA.items():
        if seq not in commons or rep not in commons:
            continue
        if dictB[seq] != dictB[rep]:
            different_res += 1
    print(f'different: {different_res}')

# 从dict转成list
def get_list_from_dict(dict):
    result = [value for key, value in sorted(dict.items())]
    return result

import sys
import struct
num_arguments=len(sys.argv)-1
if num_arguments==1 :
    yclust_path = sys.argv[1]
    parse_yclust_res(yclust_path)
else:
    cdhit_path = sys.argv[1]
    yclust_path = sys.argv[2]
    cdhit_dict = read_clstr(cdhit_path)
    yclust_dict = read_yclust_res(yclust_path)
    find_different_groups(cdhit_dict, yclust_dict)



from collections import defaultdict
def check_cross_group_clustering(data, pre_groups, clusters):
    """
    检查是否有不同组的元素被分到了同一个类中。
    
    参数:
    - data: 数据集的元素列表，例如 [1, 2, 3, 4, 5]。
    - pre_groups: 数据的预分组结果，列表长度与 data 一致，例如 [0, 0, 1, 1, 2]。
    - clusters: 数据的聚类结果，列表长度与 data 一致，例如 [0, 1, 1, 2, 0]。
    
    输出:
    - has_error: 是否存在不同组分到同一个类的情况。
    - errors: 字典，键为类，值为属于该类的多个组。
    """
    # 建立类 -> 组的映射
    cluster_to_groups = defaultdict(set)
    for item, group, cluster in zip(data, pre_groups, clusters):
        cluster_to_groups[cluster].add(group)
    
    # 检查是否有多个组被分到同一个类
    has_error = False
    errors = {}
    for cluster, groups in cluster_to_groups.items():
        if len(groups) > 1:
            has_error = True
            errors[cluster] = groups
    
    return has_error, errors

