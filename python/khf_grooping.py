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

# Ensure there are exactly 100 sequences each with 10 values
import random
import struct
import sys

if len(sys.argv) < 2:
    print("usage: python3 test.py binary_input_file")
    sys.exit(1)
# sequences = [[random.randint(1, 2**31-1) for _ in range(15)] for _ in range(100000000)]

# for row in sequences:
#       print(" ".join(map(str, row)))

file_path = sys.argv[1]
m = 15
sequences = read_sequences_from_binary(file_path, m)
print(f"length of sequences read from file:{len(sequences)}")
groups = group_sequences(sequences)

# Print the resulting groups
for root, group in groups.items():
    print(f"Group {root}: {group}")
print(f"groupnumber: {len(groups)}")
