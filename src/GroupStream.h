#ifndef __GROUPSTREAM_H__
#define __GROUPSTREAM_H__

#include "unionfind.h"
#include "libcdhit/cdhit.h"
#include <math.h>

extern void initOptions();
extern void setOptionsClusterThd(float cluster_thd);

struct Task {
	vector<vector<int>> task_cluster;
	int required_threads;
	Task() = default;
	Task(const vector<vector<int>>& _cluster, int _threads)
 		: task_cluster(_cluster), required_threads(_threads) {}
};
class GroupStream {
	UnionFind uf;
	ClusterWS ws;
	int items;
	int M = 15;
	int R = 1;
	int L = 1;
	int K = 8; 
	int cluster_condition = 500000;
	int num_threads = 8;
	bool slide = true;
	bool cluster_on = false;
	bool threadPool_on = false;
	bool output_on = false;
	string res_file = "";
	vector<Data> hash_vec;
	vector<string> names;
	unordered_map<uint64_t, string> fa_map;


	// count time
	unordered_map<int, int> cdhit_cnt;
	unordered_map<int, int> build_cnt;
	unordered_map<int, int> update_cnt;
	unordered_map<int, int> tasks_cnt;
	void init_cnt(){
		tasks_cnt[10000]=0;
		tasks_cnt[50000]=0;
		tasks_cnt[100000]=0;
		tasks_cnt[500000]=0;
		tasks_cnt[1000000]=0;
		tasks_cnt[5000000]=0;
		tasks_cnt[10000000]=0;

		cdhit_cnt[10000]=0;
		cdhit_cnt[50000]=0;
		cdhit_cnt[100000]=0;
		cdhit_cnt[500000]=0;
		cdhit_cnt[1000000]=0;
		cdhit_cnt[5000000]=0;
		cdhit_cnt[10000000]=0;

		build_cnt[10000]=0;
		build_cnt[50000]=0;
		build_cnt[100000]=0;
		build_cnt[500000]=0;
		build_cnt[1000000]=0;
		build_cnt[5000000]=0;
		build_cnt[10000000]=0;

		update_cnt[10000]=0;
		update_cnt[50000]=0;
		update_cnt[100000]=0;
		update_cnt[500000]=0;
		update_cnt[1000000]=0;
		update_cnt[5000000]=0;
		update_cnt[10000000]=0;
	}
	// 存储seq-id到root-id的映射
	vector<int> id_root_map;
	// 存储序列在hash-vec里的顺序和读入顺序的映射
	vector<uint64_t> seq_ids; 


 	// use unionfind to unite group results by per column
	void Unite(const vector<Data>& dataList, UnionFind& uf);

	// sort hash-vec by value 1.used for fast unite; 2.used for constructing GroupResMap
	void Sort(vector<Data>& dataList);

	// grouping by column
	void GroupByCol(vector<Data>& hash_vec, UnionFind& uf);

	// construct a sorted struct Data(hash-vec) for a column of hash-funtions(vec)
	void fillHashVec(const vector<vector<uint64_t>>& vec, vector<Data>& hash_vec, int m);
	void fillHashVec(string sketch_filename, vector<Data>& hash_vec, int m);
	
	void checkEdges(vector<Data>& hash_vec, UnionFind& cur_uf); 
	void unite_by_edges(UnionFind& col_uf);

	void countGroupSize(int m, UnionFind& uf);
	void countGroupSizeBySort(UnionFind& uf);
	
	void get_group_res(UnionFind& uf,unordered_map<int, vector<int>>& group_map);

	void clusterEachGroup(vector<int>& seq_ids);
	void clusterEachGroup(vector<int>& seq_ids, int needed_threads);

	void Cluster(vector<vector<int>>& cluster_sequences);

    void build_connected_components(vector<int>& group_seqs, int needed_threads);

	void cut_edges(vector<vector<int>>& sequences_collisions, int huge_groups_cnt);

	void tempOutput(vector<vector<int>>& cluster_sequences);

	void outputClstr();

public:
	GroupStream(int n, int m, int r, int l) : uf(n), items(n), R(r), L(l), M(m) {
		valid_items = n;
		resize(items);
		initOptions();
	}

	GroupStream(int n) : uf(n), items(n) {
		resize(items);
		initOptions();
    }

	void setClusterThd(float cluster_thd) {
		setOptionsClusterThd(cluster_thd);
	}
	void setIDs(const vector<uint64_t>& seq_ids) {
		this->seq_ids = seq_ids;
	}
	
	int valid_items;
	bool rep_only_group = false;
	// 分组时只考虑代表序列
	bool final_cluster_on = false;
	// 最后一轮全聚类,cluster-condition=1
	bool second_group = false;
	
	void setOutput(string res_file_name) {
		output_on = true; 
		res_file = res_file_name;
	}

	void setNames(vector<string>&& sequences_names) {
		names = std::move(sequences_names);
	}
	void setNames(vector<string>& sequences_names) {
		names = sequences_names;
	}
	
	void setSequences(unordered_map<uint64_t, string>&& sequences_aa) {
		fa_map = std::move(sequences_aa);
	}
	void setSequences(unordered_map<uint64_t, string>& sequences_aa) {
		fa_map = sequences_aa;
	}

	void setFinalClusterOn() { final_cluster_on = true; }
	void setThreadPool() { threadPool_on = true; }
	void setClusterOn() { cluster_on = true; }
	void setClusterCondition(int conditon) { cluster_condition = conditon; }
	void setSlideOff(){ slide = false; }
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
	void Group(
		vector<vector<uint64_t>>& hashes, 
		unordered_map<int, vector<int>>& group_map
		);
	void Group(string sketch_filename, unordered_map<int, vector<int>>& group_map);
	// grouping sequences with m hash-functions

   
};
#endif
