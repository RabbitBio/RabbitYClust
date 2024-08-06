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
    seq_list = []
    with open(file_path, 'r') as file:
        while True:
            identifier_line = file.readline().strip()
            sequence_line = file.readline().strip()
            if not identifier_line or not sequence_line:
                break
            name, size = identifier_line.split()
            seq_list.append((name, size, sequence_line))
    return seq_list

def select_center_sequences(seq_list, groups):
    center_sequences_list = []
    for root, group in groups.items():
        for seq_id in group:
            max_size = 0
            max_seq_id = 0
#max_size = max_size if seq_list[seq_id][1] > max_size else seq_list[seq_id][1]
            if max_size < (int)(seq_list[seq_id][1]):
                max_size = (int)(seq_list[seq_id][1])
                max_seq_id = seq_id
# wirte element with the max size of the group result to the FASTQ file
# method1: use if statement to compare to choose the max one in each iteration
# method2: put all elements'size and id in a list and use max(list, key=lambda x : x[1]) to get the max one
# we use methon1 at here
        seq_name = seq_list[max_seq_id][0]
        #seq = seq_list[max_seq_id][2]
        center_sequences_list.append((seq_name))
        #for i in range(0, len(seq), 80):
        #    file.write(seq[i:i+80] + '\n')
    return center_sequences_list

import itertools
def read_fasta_from_clust(file_path):
    seq_dict = {}
    with open(file_path, 'r') as file:
        content = file.read()
    entries = content.strip().split('>')
    # reserve sequences names
    entries = [entry[:entry.find('\n')] for entry in entries if entry]
    for entry in entries:
        keys = entry.split('\x01')
        if len(keys) > 1: # put clusts with multiple sequences into seq_dict
            keys = [f'>{key.split()[0]}' for key in keys]
            value = keys[0]
            seq_dict.update({key : value for key in keys})
#    t_items = itertools.islice(seq_dict.items(), 3)
#    for key, value in seq_dict.items():
#        print(f"{key}:{value}")
    entries = [entry.split()[0] for entry in entries if entry]
    seq_list = [f'>{entry}' for entry in entries]
    return seq_list, seq_dict

def find_unique_sequence(group_list, clust_list, clust_dict):
    group_set = set(group_list)
    clust_set = set(clust_list)
    #unique_elements = group_set - clust_set
    unique_elements = clust_set - group_set
    print(len(unique_elements))
    #print(f"{len(unique_elements)} sequences should be in a cluster but in different group")
   # for seq in unique_elements:
   #     if seq in clust_dict:
   #         print(f"{seq} : {clust_dict[seq]}")
   #     else:
   #         print(f"{seq} : NULL")
import random
import struct
import sys

if len(sys.argv) < 4:
    print("usage: python3 test.py binary_input_file FASTA_input clust_input")
    sys.exit(1)
# sequences = [[random.randint(1, 2**31-1) for _ in range(15)] for _ in range(100000000)]

# for row in sequences:
#       print(" ".join(map(str, row)))

input_file_path = sys.argv[1]
m = 15
seq_input = read_sequences_from_binary(input_file_path, m)
print(f"length of sequences read from {input_file_path}:{len(seq_input)}")
groups = group_sequences(seq_input)
print(f"grouping results: {len(groups)}")
#
FASTA_file_path = sys.argv[2]
sequences = read_sequence_info(FASTA_file_path)
print(f"read sequences' name, size and content from {FASTA_file_path}")
center_sequences_list = select_center_sequences(sequences, groups)
print(f"select center sequences of grooping results")

clust_file_path = sys.argv[3]
clust_list, clust_dict = read_fasta_from_clust(clust_file_path)
print(f"read {len(clust_list)} center sequence name from cd-hit clust results")

print(len(clust_dict))
find_unique_sequence(center_sequences_list, clust_list, clust_dict)
# Print the resulting groups
#for root, group in groups.items():
#    print(f"Group {root}: {group}")
#print(f"groupnumber: {len(groups)}")
