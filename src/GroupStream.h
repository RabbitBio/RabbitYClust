#ifndef __GROUPSTREAM_H__
#define __GROUPSTREAM_H__

#include "unionfind.h"
#include "libcdhit/cdhit.h"
#include <math.h>
#include "ConfigData.h"
#include "SharedData.h"
#include <unordered_map>
#include "ProteinAAStore.hpp"

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
		int cluster_condition = -1;
		float similarity = 0.9;
		bool output_on = true;
		string res_file = "";
        string names_path = "";
	};

	explicit GroupStream(const Config& cfg);

	void resize(int n) {
		hash_vec.resize(gs_config.items);
		for(auto& data : hash_vec)
			data.value.resize(gs_config.L * gs_config.R);
	}

	void Group( const ProteinSketchData& sketchdata, const ProteinData& proteindata);
	void Group(string sketch_filename, const ProteinData& proteindata);
	void Group(string sketch_filename, ProteinAAStore& store);

private:
	int greedy_condition = 1;
	Config gs_config;
	UnionFind uf;
	double tau = 0.5; // TODO 根据用户输入的similarity—threshold计算tau
	double ed_thres = 0.5;
	
	vector<Data> hash_vec;
	vector<sharedData> seq_vec;

	void initVec(int n);
	void reset(vector<pair<uint32_t, uint32_t>>& minhash_collisions, vector<sharedData>& seq_hash_vec);


 	// use unionfind to unite group results by per column
	void Unite(const vector<Data>& dataList, UnionFind& uf);

	void Sort(vector<Data>& dataList);

	// grouping by column
	void GroupByCol(vector<Data>& hash_vec, ProteinAAStore& store);
	void GroupByCol(vector<Data>& hash_vec, const vector<string>& fa_map);
	void GroupByCol(vector<pair<uint32_t, uint32_t>>& minhash_collisions, ProteinAAStore& store);

	void fillHashVec(const ProteinSketchData& sketchdata, vector<Data>& hash_vec, int m);
	void fillHashVec(string sketch_filename, vector<Data>& hash_vec, int m);

	void fillHashVecAndDetectMinHash(string sketch_filename, vector<sharedData>& seq_vec, int m, vector<pair<uint32_t, uint32_t>>& minhash_collisions);
	
	void checkEdges(vector<Data>& hash_vec, ProteinAAStore& store, vector<vector<pair<int, int>>>& minhash_collisions); 
	void checkEdges(vector<Data>& hash_vec, const vector<string>& fa_map, vector<vector<pair<int, int>>>& minhash_collisions); 
	
	void uniteByEdges(
		vector<pair<uint32_t, uint32_t>>& minhash_collisions,
		int start_update_idx
		);
	void uniteByEdges(vector<vector<pair<int, int>>>& minhash_collisions);
	void uniteByEdges(UnionFind& col_uf);
	void uniteByEdges(vector<int>& id_root_map);

	void countGroupSize(int m, UnionFind& uf, const vector<string>& fa_map, vector<int>& id_root_map);
	void countGroupSizeBySort(vector<pair<uint32_t, uint32_t>>& need_to_clutser, int cluster_condition);
	
	void getGroupRes(UnionFind& uf,unordered_map<int, vector<int>>& group_map, vector<int>& id_root_map);

	void clusterEachGroup(vector<int>& group_seqs, int needed_threads, const vector<string>& fa_map);
	//void clusterEachGroup(vector<int>& group_seqs, int needed_threads, const unordered_map<uint64_t, string>& fa_map);

	void ClusterFinally(
		vector<pair<uint32_t, uint32_t>>& need_to_cluster, 
		ProteinAAStore& store, 
		bool is_cluster, 
		int start_idx,
		int end_idx
		);
	void Cluster(vector<pair<uint32_t, uint32_t>>& need_to_clutser, ProteinAAStore& store);
	void Cluster(vector<vector<int>>& cluster_sequences, const vector<string>& fa_map);
	//void Cluster(vector<vector<int>>& cluster_sequences, const unordered_map<uint64_t, string>& fa_map);

	void callLib_cdhit(
			vector<pair<uint32_t, uint32_t>>& minhash_collisions,
			ProteinAAStore& store,
			int cdhit_thres// 从哪个下标开始做断边
			);
	void buildConnectedComponentsByLib_cdhit(
			int needed_threads, 
			vector<pair<uint32_t, uint32_t>>& task,
			ProteinAAStore& store
			);
	
    void buildConnectedComponentsByLib_cdhit(
		int needed_threads, 
		uint32_t start_idx,
		uint32_t end_idx,
		ProteinAAStore& store
		);
    vector<uint64_t> buildConnectedComponents(
		int needed_threads, 
		uint32_t start_idx,
		uint32_t end_idx,
		ProteinAAStore& store
		);
    vector<uint64_t> buildConnectedComponents(vector<pair<int,int>>& group_seqs, int needed_threads, const vector<string>& fa_map, ClusterWS& ws);
    //vector<uint64_t> buildConnectedComponents(vector<int>& group_seqs, int needed_threads, const unordered_map<uint64_t, string>& fa_map, ClusterWS& ws);

    vector<uint64_t> buildConnectedComponents_st(vector<pair<int, int>>& group_seqs, ProteinAAStore& store, int use_wt);
    vector<uint64_t> buildConnectedComponents_st(vector<pair<int, int>>& group_seqs, const vector<string>& fa_map, int use_wt);
    //vector<uint64_t> buildConnectedComponents_st( vector<pair<int, int>>& group_seqs, const unordered_map<uint64_t, string>& fa_map, int use_wt);
    //vector<uint64_t> buildConnectedComponents_st( vector<int>& group_seqs, const unordered_map<uint64_t, string>& fa_map, int use_wt);
    //void buildConnectedComponents_st( vector<int>& group_seqs, const unordered_map<uint64_t, string>& fa_map, int use_wt);

	void cutEdges(
		vector<pair<uint32_t, uint32_t>>& minhash_collisions, 
		ProteinAAStore& store,
		int start_idx,
		int end_idx
		);
	void cutEdges(vector<vector<pair<int, int>>>& sequences_collisions, ProteinAAStore& store);
	void cutEdges(vector<vector<pair<int, int>>>& sequences_collisions, int huge_groups_cnt, const vector<string>& fa_map);
	//void cutEdges(vector<vector<pair<int, int>>>& sequences_collisions, int huge_groups_cnt, const unordered_map<uint64_t, string>& fa_map);
	//void cutEdges(vector<vector<int>>& sequences_collisions, int huge_groups_cnt, const unordered_map<uint64_t, string>& fa_map);

	void outputClstr(const vector<string>& names);
    void outputClstr(ProteinAAStore& store);



	// count time
	//unordered_map<int, int> cdhit_cnt;
	//unordered_map<int, int> build_cnt;
	//unordered_map<int, int> update_cnt;
	//unordered_map<int, int> tasks_cnt;
	//void init_cnt(){
	//	tasks_cnt[10000]=0;
	//	tasks_cnt[50000]=0;
	//	tasks_cnt[100000]=0;
	//	tasks_cnt[500000]=0;
	//	tasks_cnt[1000000]=0;
	//	tasks_cnt[5000000]=0;
	//	tasks_cnt[10000000]=0;

	//	cdhit_cnt[10000]=0;
	//	cdhit_cnt[50000]=0;
	//	cdhit_cnt[100000]=0;
	//	cdhit_cnt[500000]=0;
	//	cdhit_cnt[1000000]=0;
	//	cdhit_cnt[5000000]=0;
	//	cdhit_cnt[10000000]=0;

	//	build_cnt[10000]=0;
	//	build_cnt[50000]=0;
	//	build_cnt[100000]=0;
	//	build_cnt[500000]=0;
	//	build_cnt[1000000]=0;
	//	build_cnt[5000000]=0;
	//	build_cnt[10000000]=0;

	//	update_cnt[10000]=0;
	//	update_cnt[50000]=0;
	//	update_cnt[100000]=0;
	//	update_cnt[500000]=0;
	//	update_cnt[1000000]=0;
	//	update_cnt[5000000]=0;
	//	update_cnt[10000000]=0;
	//}
   
};
#endif
