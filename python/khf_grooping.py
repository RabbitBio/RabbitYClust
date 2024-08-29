class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, u):
        if self.parent[u] != u:
            self.parent[u] = self.find(self.parent[u])
        return self.parent[u]

    def union(self, u, v):
        root_u = self.find(u)
        root_v = self.find(v)
        if root_u != root_v:
            if self.rank[root_u] > self.rank[root_v]:
                self.parent[root_v] = root_u
            elif self.rank[root_u] < self.rank[root_v]:
                self.parent[root_u] = root_v
            else:
                self.parent[root_v] = root_u
                self.rank[root_u] += 1

def group_sequences(sequences):
    num_sequences = len(sequences)
    sketch_len = len(sequences[0])
    
    # Initialize Union-Find for num_sequences elements
    uf = UnionFind(num_sequences)

    # Group by each column
    for col in range(sketch_len):
        value_to_ids = {}
        for seq_id in range(num_sequences):
            value = sequences[seq_id][col]
            if value not in value_to_ids:
                value_to_ids[value] = []
            value_to_ids[value].append(seq_id)
        
        # Union groups with the same value
        for ids in value_to_ids.values():
            for i in range(1, len(ids)):
                uf.union(ids[0], ids[i])
    
    # Collect the final groups
    groups = {}
    for seq_id in range(num_sequences):
        root = uf.find(seq_id)
        if root not in groups:
            groups[root] = []
        groups[root].append(seq_id)
    
    return groups

def read_sequences_from_binary(file_path, m):
    sequences = []

    with open(file_path, 'rb') as file:
        data = file.read(m * 8)
        while data:
            # read m uint64 value（each uint64 is 8bytes）
            if len(data) < 8:
                print(f"incomplete data:{data}")
                continue
            # struct.unpack: tranfrom binary data to uint64 tuple
            sequence = struct.unpack(f'{m}Q', data)  # 'Q' 表示 unsigned long long (uint64)
            sequences.append(sequence)
            data = file.read(m * 8)

    return sequences

def read_sequence_info(file_path):
    seq_dict = {}
    with open(file_path, 'r') as file:
        seq_counter = 0
        while True:
            identifier_line = file.readline().strip()
            sequence_line = file.readline().strip()
            if not identifier_line or not sequence_line:
                break
            name, size = identifier_line.split()
            seq_list = [name[1:], size, sequence_line]
            seq_dict.update({seq_counter : seq_list})
            seq_counter += 1
    return seq_dict

def select_center_sequences(seq_dict, groups):
    center_sequences_list = []
    group_dict = {}
    for root, group in groups.items():
        max_size = -1
        max_seq_id = None
        entries = []
        for seq_id in group:
            entries.append(seq_dict[seq_id][0])
            if max_size == int(seq_dict[seq_id][1]):
                max_seq_id = int(seq_id) if int(seq_id) < max_seq_id else max_seq_id
            elif max_size < (int)(seq_dict[seq_id][1]):
                max_size = (int)(seq_dict[seq_id][1])
                max_seq_id = seq_id
# wirte element with the max size of the group result to the FASTQ file
# method1: use if statement to compare to choose the max one in each iteration
# method2: put all elements'size and id in a list and use max(list, key=lambda x : x[1]) to get the max one
#        group_list = [seq_dict[key] for key in group if key in seq_dict]
#        max_size = max(group_list, key = lambda x : x[1])
#        max_size_counter = sum(1 for element in group_list if element[1] == max_size)
#        if max_size_counter > 1:
#            
        seq_name = seq_dict[max_seq_id][0]
        center_sequences_list.append((seq_name))
        group_dict[seq_name] = entries

        #for i in range(0, len(seq), 80):
        #    file.write(seq[i:i+80] + '\n')
    return center_sequences_list, group_dict

#import itertools
def read_fasta_from_clust(file_path):
#    seq_dict = {}
    with open(file_path, 'r') as file:
        content = file.read()
    entries = content.strip().split('>')
    # reserve sequences names
    entries = [entry[:entry.find('\n')] for entry in entries if entry]
#    t_items = itertools.islice(seq_dict.items(), 3)
#    for key, value in seq_dict.items():
#        print(f"{key}:{value}")
    entries = [entry.split()[0] for entry in entries if entry]
    seq_list = [f'>{entry}' for entry in entries]
    return seq_list

import re
def read_clust_clstr(file_path):
    clusters_dict={}
    current_cluster = []
    key = None
    # 正则表达式用于匹配以 > 开头并以 .1 结尾的条目
    pattern = re.compile(r'>(\S+)\.\.\.')
    with open(file_path,'r') as file:
        for line in file.readlines()[1:]:
            line = line.strip()
            if line.startswith('>Cluster'):
                # 保存当前 cluster 的数据
                if current_cluster is not None:
                    clusters_dict[key] = current_cluster
                #print(f"cluster{key} : {current_cluster}")
                # 开始新的 cluster
                current_cluster = []
                key = None
            else:
                match = pattern.search(line)
                if match:
                    current_cluster.append(match.group(1))
                    if line.endswith('*'):
                        key = match.group(1)
    if current_cluster is not None:
        clusters_dict[key] = current_cluster
    return clusters_dict


def find_unique_sequence(group_list, clusters_dict):
    clust_list = list(clusters_dict.keys())
    print(f"read {len(clust_list)} center sequence name from cd-hit clust results")
    group_set = set(group_list)
    clust_set = set(clust_list)
    unique_elements = group_set - clust_set
    #unique_elements = clust_set - group_set
    print(f"{len(unique_elements)} sequences should be in a cluster but in different group")
    for seq in unique_elements:
        if seq in clust_list:
            print(f"{seq} : {clust_list[seq]}")
        else:
            print(f"{seq} : NULL")

def check_group_result(group_dict, clusters_dict):
    group_rev_dict = {}
    clust_rev_dict = {}
    for key, items in group_dict.items():
        for item in items:
            group_rev_dict[item] = key
    for key, items in clusters_dict.items():
        for item in items:
            clust_rev_dict[item] = key
    error_group_list = []
    for seq, center_seq in group_rev_dict.items():
        if center_seq != clust_rev_dict[seq]:
            error_group_list.append(key)
            print(f'right center seq of {seq} is {clust_rev_dict[seq]} not {center_seq}')
    print(len(error_group_list))
    for item in error_group_list:
        print(f'{item}')
    
import random
import struct
import sys

if len(sys.argv) < 4:
    print("usage: python3 test.py binary_input_file seq_FASTA_input cd-hit_clust_input")
    sys.exit(1)
# sequences = [[random.randint(1, 2**31-1) for _ in range(15)] for _ in range(100000000)]

# for row in sequences:
#       print(" ".join(map(str, row)))

input_file_path = sys.argv[1]
m = 15
seq_input = read_sequences_from_binary(input_file_path, m)
#for seq in seq_input:
#    print(f"{seq}")
print(f"length of sequences read from {input_file_path}:{len(seq_input)}")
groups = group_sequences(seq_input)
print(f"grouping results: {len(groups)}")
#
FASTA_file_path = sys.argv[2]
sequences = read_sequence_info(FASTA_file_path)
print(f"read sequences' name, size and content from {FASTA_file_path}")
center_sequences_list, group_dict = select_center_sequences(sequences, groups)
print(f"select center sequences of grooping results")

clust_file_path = sys.argv[3]
clusters_dict = read_clust_clstr(clust_file_path)
#check_group_result(group_dict, clusters_dict)
find_unique_sequence(center_sequences_list, clusters_dict)
# Print the resulting groups
#for root, group in groups.items():
#    print(f"Group {root}: {group}")
#print(f"groupnumber: {len(groups)}")
