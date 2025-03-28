#ifndef __UNIFIND_H__
#define __UNIFIND_H__

#include <vector>
#include <unordered_map>
#include <chrono>
#include <iostream>
#include <numeric> //iota
#include <algorithm>
#include <variant>

using namespace std;

struct Data {
	vector<uint64_t> value;
	int id;
	Data(int id, vector<uint64_t> hash) : id(id), value(hash) {}
	Data() {}
};

struct GroupNode {
	int id;
	int root;
	GroupNode(int id, int root_id) : id(id), root(root_id) {}
	GroupNode() {}
};

class UnionFind {
private:
//	vector<int> parent;
	vector<int> rank;

public:	
	vector<int> parent;

    UnionFind(int n) {
        parent.resize(n);
		rank.resize(n,0);
        iota(parent.begin(), parent.end(), 0); // 初始化为自身
        // iota(rank.begin(), rank.end(), 0); 
		//cout << "unionfind generation ends." << endl;
    }
	//拷贝
	UnionFind(const UnionFind&) = default;
	UnionFind& operator=(const UnionFind&) = default;
    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]); // 路径压缩
        }
        return parent[x];
    }


    void unite(int x, int y) {
	    int rootX = find(x);
	    int rootY = find(y);
	    if (rootX != rootY) {
	        if (rank[rootX] > rank[rootY]) {
	            parent[rootY] = rootX;
	        } else if (rank[rootX] < rank[rootY]) {
	            parent[rootX] = rootY;
	        } else {
	            parent[rootY] = rootX;
	            rank[rootX]++;
	        }
	    }
	}

	int countSetsSize() {
		int count = 0;
		for (int i = 0; i < parent.size(); i++) {
			if (parent[i] == i) count++;
		}
		return count;
	}

	void updateParent(const vector<int> & new_parents) {
		for(int i=0; i < parent.size(); i++) {
			parent[i] = new_parents[i];
		}
	}

	void findRoot (vector<int>& root) {
		for(int i=0; i < parent.size(); i++) {
			int root_id = find(parent[i]);
			root[i] = root_id;
		}
	}

	void findRoot (vector<GroupNode>& root) {
		for(int i=0; i < parent.size(); i++) {
			int root_id = find(parent[i]);
			root[i].id = i;
			root[i].root = root_id;
		}
	}

	void findRoot (vector<Data>& root) {
		for(int i=0; i < parent.size(); i++) {
			int root_id = find(parent[i]);
			root[i].id = i;
			root[i].value[0] = root_id;
		}
	}

};
#endif
