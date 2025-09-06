#ifndef __GROUPSTREAM_H__
#define __GROUPSTREAM_H__

#include "unionfind.h"
#include "libcdhit/cdhit.h"
#include <math.h>

extern unordered_map<uint64_t, string> fa_map;
extern vector<string> names;
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
	UnionFind uf;
	int items;
	int M = 15;
	int R = 1;
	int L = 1;
	int K = 8; 
	int cluster_condition = 1000000;
	int num_threads = 8;
	bool slide = true;
	bool cluster_on = false;
	bool threadPool_on = false;
	bool temp_output_on = false;
	bool output_on = false;
	bool small_data_merge_cluster_on = true;
    bool break_unite = false;
    bool break_directly = false;
	string res_file = "";
	vector<Data> hash_vec;
	string folder_name = "nr-15/";

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

	// 利用jaccard相似度来进一步筛选
	float S = 0.9;
	float distance = 1 - S;
	int sharedKmersThd = 2;
/*
	void computeSharedKmersThd() {
	    float exp_neg_kd = std::exp(-k * d);
	    sharedKmersThd = static_cast<int>(M * exp_neg_kd);
    	// // 根据公式：j = e^(-k * d) / (2 - e^(-k * d))
	    // float jaccard = exp_neg_kd / (2 - exp_neg_kd);
	    // // 然后计算共享 k-mers: w = j * (2n) / (1 + j)
	    // sharedKmersThd = jaccard * (2.0 * n) / (1.0 + jaccard);
	}
*/
	bool checkJaccard(int node1, int node2);
	//	vector<GroupNode> id_root_map;
	vector<int> id_root_map;
	// 存储seq-id到root-id的映射
	vector<uint64_t> seq_ids; 
	// 存储序列在hash-vec里的顺序和读入顺序的映射
	static vector<int> GroupSizeCnt;

	GroupStream(int n, int m, int r, int l) : uf(n), items(n), R(r), L(l), M(m) {
		valid_items = n;
		resize(items);
		initOptions();
	}

	GroupStream(int n) : uf(n), items(n) {
		resize(items);
		initOptions();
    }

	void setSimilarityThd(float similarity_thd) {
		S = similarity_thd;
		distance = 1 - S;
	}

	void setClusterThd(float cluster_thd) {
		setOptionsClusterThd(cluster_thd);
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
        cerr << "break_unite: "<< break_unite << endl;
        cerr << "break_directly: "<< break_directly << endl;
		hash_vec.resize(items);
		id_root_map.resize(items, -1);
		GroupSizeCnt.resize(items, 1);
		for(auto& data : hash_vec)
			data.value.resize(L * R);
	}
	// 这个hashes还是得改成全局的
	void Group(vector<vector<uint64_t>>& hashes, unordered_map<int, vector<int>>& group_map);
	// grouping sequences with m hash-functions

	void checkEdges(vector<Data>& hash_vec, UnionFind& cur_uf);
	void unite_by_edges(UnionFind& col_uf);
	void Unite(const vector<Data>& dataList, UnionFind& uf);
 	// use unionfind to unite group results by per column
	void Sort(vector<Data>& dataList);
	// sort hash-vec by value 1.used for fast unite; 2.used for constructing GroupResMap
	void GroupByCol(vector<Data>& hash_vec, UnionFind& uf);
	// grouping by column

	void fillHashVec(const vector<vector<uint64_t>>& vec, vector<Data>& hash_vec, int m);
	// construct a sorted struct Data(hash-vec) for a column of hash-funtions(vec)

	void countGroupSize(UnionFind& uf);
	void countGroupSizeBySort(UnionFind& uf);
	
	void getGroupMap(UnionFind& uf,unordered_map<int, vector<int>>& group_map);
	//unordered_map<int, vector<int>> getGroupMap();

	void clusterEachGroup(vector<int>& seq_ids);
	void clusterEachGroup(vector<int>& seq_ids, int neededThread);
	void clusterEachGroup_st(vector<int>& group_seqs);

	void Cluster(vector<vector<int>>& cluster_sequences);

	void tempOutput(vector<vector<int>>& cluster_sequences);

	void outputClstr();
	void checkUnite(unordered_map<int, vector<int>>& map, UnionFind& uf);
    
    static bool compareByHashAndGroupSize(const Data &a, const Data &b) {
        if(a.value != b.value)
    	    return a.value < b.value;
        return GroupSizeCnt[a.id] < GroupSizeCnt[b.id];
    }
};
#endif
