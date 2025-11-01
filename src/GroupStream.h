#ifndef __GROUPSTREAM_H__
#define __GROUPSTREAM_H__

#include "unionfind.h"
#include "libcdhit/cdhit.h"
#include <math.h>
#include "ConfigData.h"
#include "SharedData.h"
#include <unordered_map>

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
public:
	struct Config{
		int items;
		int M = 15;
		int R = 1;
		int L = 1;
		int num_threads = 10;
		bool cluster_on = false;
		bool final_cluster_on = false;
		int cluster_condition = 500000;
		bool output_on = true;
		string res_file = "";
	};

	explicit GroupStream(const Config& cfg);

/* unused
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
*/

	void resize(int n) {
		hash_vec.resize(gs_config.items);
		id_root_map.resize(gs_config.items, -1);
		for(auto& data : hash_vec)
			data.value.resize(gs_config.L * gs_config.R);
	}

	//void Group(
	//	vector<vector<uint64_t>>& hashes, 
	//	unordered_map<int, vector<int>>& group_map
	//	);
	void Group(string sketch_filename, const ProteinData& proteindata);

private:
	Config gs_config;
	UnionFind uf;
	
	vector<Data> hash_vec;

	// 存储seq-id到root-id的映射
	vector<int> id_root_map;
	// 存储序列在hash-vec里的顺序和读入顺序的映射


 	// use unionfind to unite group results by per column
	void Unite(const vector<Data>& dataList, UnionFind& uf);

	// sort hash-vec by value 1.used for fast unite; 2.used for constructing GroupResMap
	void Sort(vector<Data>& dataList);

	// grouping by column
	void GroupByCol(vector<Data>& hash_vec, const unordered_map<uint64_t, string>& fa_map);

	// construct a sorted struct Data(hash-vec) for a column of hash-funtions(vec)
	void fillHashVec(const vector<vector<uint64_t>>& vec, vector<Data>& hash_vec, int m);
	void fillHashVec(string sketch_filename, vector<Data>& hash_vec, int m);
	
	void checkEdges(vector<Data>& hash_vec, UnionFind& cur_uf, const unordered_map<uint64_t, string>& fa_map); 
	void unite_by_edges(UnionFind& col_uf);

	void countGroupSize(int m, UnionFind& uf, const unordered_map<uint64_t, string>& fa_map);
	void countGroupSizeBySort(UnionFind& uf);
	
	void get_group_res(UnionFind& uf,unordered_map<int, vector<int>>& group_map);

	void clusterEachGroup(vector<int>& group_seqs, int needed_threads, const unordered_map<uint64_t, string>& fa_map);

	void Cluster(vector<vector<int>>& cluster_sequences, const unordered_map<uint64_t, string>& fa_map);

    void build_connected_components(vector<int>& group_seqs, int needed_threads, const unordered_map<uint64_t, string>& fa_map);
    void build_connected_components_st(vector<int>& group_seqs, int needed_threads, const unordered_map<uint64_t, string>& fa_map);
    void build_connected_components_st_reuse(vector<int>& group_seqs, int needed_threads, const unordered_map<uint64_t, string>& fa_map, ClusterWS& ws);

	void cut_edges(vector<vector<int>>& sequences_collisions, int huge_groups_cnt, const unordered_map<uint64_t, string>& fa_map);

	void outputClstr(const vector<string>& names, const unordered_map<uint64_t, string>& fa_map);



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
   
};
#endif
