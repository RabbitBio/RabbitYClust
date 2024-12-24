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

struct Data1 {
	uint64_t value;
	int id;

	// 构造函数用于 r=1 的情况
	Data1(int id, uint64_t hash1) : id(id), value(hash1){}

	// 打印结构体数据，方便调试
	void print() const {
		cout << "Hash: " << value << ", ID: " << id << endl;
	}
	Data1() {}
};

struct Data {
	uint64_t value1, value2;
	int id;
	Data(int id, uint64_t hash1, uint64_t hash2) : id(id), value1(hash1), value2(hash2) {}
	Data() {}
};

class UnionFind {
private:
	vector<int> parent;
	vector<int> rank;

public:
	
    UnionFind(int n) {
        parent.resize(n);
		rank.resize(n);
        iota(parent.begin(), parent.end(), 0); // 初始化为自身
        iota(rank.begin(), rank.end(), 0); 
		cout << "unionfind generation ends." << endl;
    }

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

	void findRoot (vector<Data>& root) {
		for(int i=0; i < parent.size(); i++) {
			int root_id = find(parent[i]);
			root[i].id = i;
			root[i].value1 = root_id;
		}
	}

};
#endif
