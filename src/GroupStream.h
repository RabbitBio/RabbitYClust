#ifndef __GROUPSTREAM_H__
#define __GROUPSTREAM_H__

#include "unionfind.h"
#include "cluster.h"

extern unordered_map<uint64_t, string> fa_map;
extern vector<string> names;
extern void initOptions();

struct Task {
	vector<vector<int>> task_cluster;
	int required_threads;
	Task() = default;
	Task(const vector<vector<int>>& _cluster, int _threads)
 		: task_cluster(_cluster), required_threads(_threads) {}
};
class GroupStream {
public:
	UnionFind uf;
	int items;
	int M;
	int R = 1;
	int L = 1;
	int cluster_condition = 2000;
	int num_threads = 8;
	bool slide = true;
	bool cluster_on = false;
	bool threadPool_on = false;
	bool temp_output_on = false;
	bool output_on = false;
	bool small_data_merge_cluster_on = true;
	string res_file = "";
	vector<Data> hash_vec;
	string folder_name = "nr-15/";

	//	vector<GroupNode> id_root_map;
	vector<int> id_root_map;
	// 存储seq-id到root-id的映射
	vector<uint64_t> seq_ids; 
	// 存储序列在hash-vec里的顺序和读入顺序的映射

	GroupStream(int n, int m, int r, int l) : uf(n), items(n), R(r), L(l), M(m) {
		valid_items = n;
		resize(items);
		initOptions();
	}

	GroupStream(int n) : uf(n), items(n) { resize(items);
		initOptions();
       	}

	void setIDs(const vector<uint64_t>& seq_ids) {
		this->seq_ids = seq_ids;
	}
	
	int valid_items;
	int total_clusters;
	int redundant_seqs;
	bool rep_only_group = false;
	// 分组时只考虑代表序列
	bool rep_only_cluster = false;
	// 聚类时只考虑代表序列
	bool final_cluster_on = false;
	// 最后一轮全聚类,cluster-condition=1
	bool second_group = false;
	
	vector<bool> valid_seqs;
	void setRepGroupOn(){
		rep_only_group = true;
		valid_seqs = vector<bool>(items, true);
	}
	void setRepGroupAndClusterOn(){
		rep_only_group = true;
		rep_only_cluster = true;
		valid_seqs = vector<bool>(items, true);
	}
	void setValidStatus(vector<int>& group_seqs);

	void setOutput(string res_file_name) {
		output_on = true; 
		res_file = res_file_name;
	}

	void setSecondGroup() {
		if(R == 1){
			second_group = true;
		}
	}
	void setSmallDataMergeClusterOff() { small_data_merge_cluster_on = false; }
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
	void Group(vector<vector<uint64_t>>& hashes, unordered_map<int, vector<int>>& group_map);
	// grouping sequences with m hash-functions

	void Unite(const vector<Data>& dataList, UnionFind& uf);
 	// use unionfind to unite group results by per column
	void Sort(vector<Data>& dataList);
	// sort hash-vec by value 1.used for fast unite; 2.used for constructing GroupResMap
	void GroupByCol(vector<Data>& hash_vec, UnionFind& uf);
	// grouping by column

	void fillHashVec(const vector<vector<uint64_t>>& vec, vector<Data>& hash_vec, int m);
	// construct a sorted struct Data(hash-vec) for a column of hash-funtions(vec)

	void countGroupSize(int m, UnionFind& uf);
	void countGroupSizeBySort(UnionFind& uf);
	
	void getGroupMap(UnionFind& uf,unordered_map<int, vector<int>>& group_map);
	//unordered_map<int, vector<int>> getGroupMap();

	void clusterEachGroup(vector<int>& seq_ids);
	void clusterEachGroup(vector<int>& seq_ids, int neededThread);

	void Cluster(int m, vector<vector<int>>& cluster_sequences);

	void tempOutput(vector<vector<int>>& cluster_sequences);

	void outputClstr();
};
#endif
