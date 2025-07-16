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
    vector<int> lastround;
    unordered_map<int, int> groups_cnt;
    int unite_condition = 10000000;

    UnionFind(int n) {
        parent.resize(n);
        lastround.resize(n);
		rank.resize(n, 0);
        iota(parent.begin(), parent.end(), 0); // 初始化为自身
        iota(lastround.begin(), lastround.end(), 0); // 初始化为自身
        //iota(rank.begin(), rank.end(), 0); 
		//cout << "unionfind generation ends." << endl;
        for(int i = 0; i < n; i++){
            groups_cnt[i] = 1;
        }
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
            // check UniteCondition
            if (groups_cnt[rootX] + groups_cnt[rootY] > unite_condition){
            return;
        }
	        if (rank[rootX] > rank[rootY]) {
                groups_cnt[rootX] += groups_cnt[rootY];
                groups_cnt.erase(rootY);
	            parent[rootY] = rootX;
	        } else if (rank[rootX] < rank[rootY]) {
                groups_cnt[rootY] += groups_cnt[rootX];
                groups_cnt.erase(rootX);
	            parent[rootX] = rootY;
	        } else {
                groups_cnt[rootX] += groups_cnt[rootY];
                groups_cnt.erase(rootY);
	            parent[rootY] = rootX;
	            rank[rootX]++;
	        }
	    }
	}

    // 感觉应该在聚类完成后
    // 因为聚类还要更新一次groups_cnt
    void countGroupsSizeofSeqs(vector<int>& seqs) {
        for(int i = 0; i < parent.size(); i++){
            int root_i = find(i);
            seqs[i] = groups_cnt[root_i];
        }
    }

    void updateGroupSizeCnt(unordered_map<int, vector<int>>& map_after_cluster){
        // 遍历删掉不再存在的root
        for(auto& [root, size] : groups_cnt){
            if(!map_after_cluster.count(root)){
                groups_cnt.erase(root);
            }
        }
        // 增加新的root
        for(auto& [root, groups] : map_after_cluster){
            groups_cnt[root] = groups.size();
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
