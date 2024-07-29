#This testing scripts for KMV-based method.
#This is only used for proof-of-concept.
#Disjoint-set data structure is the key to reduce the compute complexity.

import random
import sys

class UnionFind:
    def __init__(self, size):
        self.parent = list(range(size))
        self.rank = [1] * size

    def find(self, p):
        if self.parent[p] != p:
            self.parent[p] = self.find(self.parent[p])  # Path compression
        return self.parent[p]

    def union(self, p, q):
        rootP = self.find(p)
        rootQ = self.find(q)
        if rootP != rootQ:
            # Union by rank
            if self.rank[rootP] > self.rank[rootQ]:
                self.parent[rootQ] = rootP
            elif self.rank[rootP] < self.rank[rootQ]:
                self.parent[rootP] = rootQ
            else:
                self.parent[rootQ] = rootP
                self.rank[rootP] += 1

def find_disjoint_groups(sets):
    element_to_set = {}
    uf = UnionFind(len(sets))

    for i, s in enumerate(sets):
        for elem in s:
            if elem in element_to_set:
                uf.union(i, element_to_set[elem])
            element_to_set[elem] = i

    groups = {}
    for i, s in enumerate(sets):
        root = uf.find(i)
        if root not in groups:
            groups[root] = []
        groups[root].append(s)

    return list(groups.values())

# Generate 10000 sets with 5 random 32-bit integers each
def generate_random_sets(num_sets, set_size):
    return [{random.randint(0, 1000000) for _ in range(set_size)} for _ in range(num_sets)]

# Example usage
sets = generate_random_sets(10**4, 10)

# Find disjoint groups
disjoint_groups = find_disjoint_groups(sets)

print(f"number of groups: {len(disjoint_groups)}", file=sys.stderr)
# Display the results
for idx, group in enumerate(disjoint_groups[:5]):  # Displaying only the first 5 groups for brevity
    print(f"Group {idx + 1}: {group}")

