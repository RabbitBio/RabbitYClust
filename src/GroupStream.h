#ifndef __GROUPSTREAM_H__
#define __GROUPSTREAM_H__

#include "unionfind.h"
#include "cluster.h"

extern unordered_map<uint64_t, string> fa_map;

class GroupStream {
public:
	UnionFind uf;
	int items;
	int M;
	int R = 1;
	int L = 1;
	int cluster_condition = 200;
	int second_condition = 100000;
	int num_threads = 8;
	bool slide = true;
	bool cluster_on = false;
	bool temp_output_on = false;
	bool second_group = false;
	vector<Data> hash_vec;
	string folder_name = "nr-15/";

//	vector<GroupNode> id_root_map;
	vector<int> id_root_map;
	// 存储seq-id到root-id的映射
	vector<uint64_t> seq_ids; 
	// 存储序列在hash-vec里的顺序和读入顺序的映射

	cluster cluster_cdhit;

	GroupStream(int n, int m, int r, int l) : uf(n), items(n), R(r), L(l), M(m) {
		resize(items);
	}

	GroupStream(int n) : uf(n), items(n) { resize(items); }

	void setIDs(const vector<uint64_t>& seq_ids) {
		this->seq_ids = seq_ids;
	}
	

	void setClusterOn() { cluster_on = true; }
	void setClusterCondition(int conditon) { cluster_condition = conditon; }
	void setSlideOff() { slide = false; }
	void setSecondGroup() { second_group = true; }
	void setM(int m) { M = m; }
	void setR(int r) { R = r; }
	void setL(int l) { L = l; }
	void setNumThreads(int threads) { num_threads = threads; }
	void resize(int n) {
		hash_vec.resize(items);
		id_root_map.resize(items, -1);
		for(auto& data : hash_vec)
			data.value.resize(L * R);
	}
	void Group(vector<vector<uint64_t>>& hashes, unordered_map<int, vector<int>>& group_map);
	// grouping sequences with m hash-functions

	void Unite(vector<Data> dataList, UnionFind& uf);
 	// use unionfind to unite group results by per column
	void Sort(vector<Data>& dataList);
	// sort hash-vec by value 1.used for fast unite; 2.used for constructing GroupResMap
	void GroupByCol(vector<Data>& hash_vec, UnionFind& uf);
	// grouping by column

	void fillHashVec(const vector<vector<uint64_t>>& vec, vector<Data>& hash_vec, int m);
	// construct a sorted struct Data(hash-vec) for a column of hash-funtions(vec)

	void countGroupSize(UnionFind& uf,int m,vector<vector<uint64_t>>& hashes);
	void countGroupSizeBySort(UnionFind& uf);
	
	void getGroupMap(UnionFind& uf,unordered_map<int, vector<int>>& group_map);
	//unordered_map<int, vector<int>> getGroupMap();

	void clusterEachGroup(vector<int>& seq_ids);

	void tempOutput(vector<vector<int>>& cluster_sequences);

	void SecondGroup(vector<int>& group_seqs, int m);

	void SecondUpdate(vector<int>& group_seqs);
};
#endif
