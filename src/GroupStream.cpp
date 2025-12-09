#include "GroupStream.h"
#include <queue>
#include <omp.h>
#include <thread>
#include <algorithm>
#include <atomic>
#include <unordered_set>
#include <fstream>

#include "ips2ra.hpp"
#include "ips4o.hpp"
// 静态成员变量的定义

int round_cnt=0;

struct minheapcompare {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

bool compareByHash(const Data &a, const Data &b) {
	return a.value < b.value;
}

GroupStream::GroupStream(const Config& cfg) : gs_config(cfg), uf(cfg.items){
	resize(gs_config.items);
	initOptions();
	if(cfg.similarity == 0.9) {
		tau = 0.36;
	}else {
		tau = 0.05;
	}
}

/*
 *
void GroupStream::tempOutput(vector<vector<int>>& cluster_sequences) {	
	unordered_map<int, vector<int>> map_after_cluster;
	priority_queue<pair<int,int>, vector<pair<int, int>>, minheapcompare> minHeap;
	for(int i = 0; i < cluster_sequences.size(); i++) {
		for(int x : cluster_sequences[i]) {
			map_after_cluster[id_root_map[x]].push_back(x);
		}
	}
	for(auto &[root_id, seqs] : map_after_cluster){
		if (seqs.size() > 10000){
			minHeap.push({seqs.size(), root_id});
		}
	}
	ofstream log("output.log");
	cerr.rdbuf(log.rdbuf());
	while(!minHeap.empty()){
		int rootid = minHeap.top().second;
		string filename = to_string(rootid) + ".fa";
		ofstream ofile(filename);
		for(int id : map_after_cluster[rootid]){
			ofile << ">" << id << endl;
			ofile << fa_map[id] << endl;
		}
		ofile.close();
		cerr << "cluster: " << rootid << " contains " << minHeap.top().first << " seqs stored in: " << filename << endl; 
		minHeap.pop();
	}
	cerr << endl;
}
*/

void ips4o_sort_single_thread(vector<Data>& hash_vec, int r = 1) {
    auto cmp = [r](const Data& a, const Data& b) {
        for (int i = 0; i < r; i++) {
            if (a.value[i] < b.value[i])
                return true;
            else if (a.value[i] > b.value[i])
                return false;
        }
        return false;
    };
    ips4o::sort(hash_vec.begin(), hash_vec.end(), cmp);
}

void ips4o_sort_multi_thread(vector<Data>& hash_vec, int r = 1) {
    auto cmp = [r](const Data& a, const Data& b) {
        for (int i = 0; i < r; i++) {
            if (a.value[i] < b.value[i])
                return true;
            else if (a.value[i] > b.value[i])
                return false;
        }
        return false;
    };
    ips4o::parallel::sort(hash_vec.begin(), hash_vec.end(), cmp);
}

void ips2ra_sort_single_thread(vector<Data>& hash_vec, const int r = 1) {
    ips2ra::sort(hash_vec.begin(), hash_vec.end(), [](const Data& r) { return r.value[0]; });
}

void ips2ra_sort_multi_thread(vector<Data>& hash_vec, const int r = 1) {
    ips2ra::parallel::sort(hash_vec.begin(), hash_vec.end(), [](const Data& r) { return r.value[0]; }, 1);
}


void GroupStream::Sort(vector<Data>& dataList){
	auto start_time = chrono::high_resolution_clock::now();
//#ifdef USE_PARALLEL
	if (gs_config.num_threads == 1) ips2ra_sort_single_thread(dataList, dataList.size());
    else ips2ra_sort_multi_thread(dataList, dataList.size());

//	if (gs_config.num_threads == 1) ips4o_sort_single_thread(dataList);
//    else ips4o_sort_multi_thread(dataList);
	// r = 1 用ips2ra
//#else
//	sort(dataList.begin(), dataList.end(), [](const Data& a, const Data& b){
//		return a.value < b.value;
//		});
//		sort(dataList.begin(), dataList.end(), compareByHash);
//#endif
	auto end_time = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
	cerr << "sort time (seconds): " << duration  << endl;
}

void GroupStream::checkEdges(
	vector<Data>& hash_vec, 
	UnionFind& cur_uf,
	const unordered_map<uint64_t, string>& fa_map
	) {
	cur_uf.findRoot(id_root_map);
	unordered_map<int, vector<int>> map;

	for(int i = 0; i < gs_config.items; i++) {
		map[id_root_map[i]].push_back(i);
	}

	vector<vector<int>> first_hit_sequences;
    int huge_groups_cnt = 0;
	for(auto &[key, seqs] : map){
		if(seqs.size() > 1) {
			first_hit_sequences.emplace_back(seqs);
		}
        if(seqs.size() > 100000) {
            huge_groups_cnt++;
        }
	}
	cerr << "MinHash collisions : " << first_hit_sequences.size() << endl;
	sort(first_hit_sequences.begin(), first_hit_sequences.end(),
		[](const vector<int>& a, const vector<int>& b){
			return a.size() > b.size();
		});

	cerr << "    Top 10 largest collisions: " ;
	for(int i = 0; i < std::min(10, (int)first_hit_sequences.size()); i++){
		cerr << first_hit_sequences[i].size() << " ";
	}
	cerr << endl;


	auto start_time = chrono::high_resolution_clock::now();
	cutEdges(first_hit_sequences, huge_groups_cnt, fa_map);
	cur_uf.updateParent(id_root_map);

	auto end_time = chrono::high_resolution_clock::now();
	auto duration_cc = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
	cerr << "	Break the bad edges time(seconds): " << duration_cc  << endl;

	unordered_map<int, vector<int>> groups_after_filter; // rootid:[seq0, seq1...]
	priority_queue<int, vector<int>, greater<int>> minHeap;

	for(int i = 0; i < gs_config.items; i++) {
		groups_after_filter[id_root_map[i]].push_back(i);
	}

	for(auto &[root_id, seqs] : groups_after_filter){
		minHeap.push(seqs.size());
		if (minHeap.size() > 10){
			minHeap.pop();
		}

	}
	cerr << "After break the bad edges, connected components number are: " << groups_after_filter.size() << endl;
	cerr << "    Top 10 largest connected components: ";
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;
}


void GroupStream::Unite(const vector<Data>& dataList, UnionFind& this_uf) {
	vector<uint64_t> cur_value = dataList[0].value;
	for (int i = 1; i < dataList.size(); i++) {
		auto thisone = dataList[i];
    	auto lastone = dataList[i-1];
		if(thisone.value == lastone.value) {
		    this_uf.unite(thisone.id, lastone.id);
		}
    }
}


void GroupStream::uniteByEdges(UnionFind& cur_uf) {
	for(int i = 0; i < gs_config.items; i++) {
		uf.unite(uf.find(i), cur_uf.find(i));
	}
}
								
void GroupStream::GroupByCol(
	vector<Data>& hash_vec, 
	const unordered_map<uint64_t, string>& fa_map
	) {
	Sort(hash_vec);
	if(gs_config.cluster_on) {
		UnionFind col_uf(gs_config.items);
		Unite(hash_vec, col_uf);
		checkEdges(hash_vec, col_uf, fa_map);
		uniteByEdges(col_uf);
	}else {
		Unite(hash_vec, uf);
	}
	int groups_size = uf.countSetsSize();
	//cerr << "---------------------------------------------------" << endl;
	cout << "Group Size after merging:" << groups_size << endl;
}

void GroupStream::fillHashVec(const ProteinSketchData& sketchdata, vector<Data>& hash_vec, int m) {
// 1. use std::transform
//	int index = 0;
//	transform(vec.begin(), vec.end(), hash_vec.begin(), 
//		[&index](int value) {
//		return Data{value, (index++)};
//		});
// 2. iteration to construct pair
    auto& tmp_hashes = sketchdata.hashes[m];
	for (int i = 0; i < gs_config.items; i++) {
		hash_vec[i].id = i;
        hash_vec[i].value[0] = tmp_hashes[i];
		//copy(vec[i].begin() + m * gs_config.L, vec[i].begin() + m * gs_config.L + gs_config.R * gs_config.L, hash_vec[i].value.begin());
	}
	cerr << gs_config.items << " valid items in round " << m << endl;
}


void GroupStream::fillHashVec(string sketch_filename, vector<Data>& hash_vec, int m) {
	int valid_items = 0;
	std::ifstream ifs(sketch_filename, ios::binary);
	if(!ifs){
		cerr << "Error opening file!" << endl;
		return;
	}
	ifs.seekg(sizeof(ProteinSketchData::Config)  + m * gs_config.items * sizeof(uint64_t), std::ios::beg);
	vector<uint64_t> read_hashes(gs_config.items);
	ifs.read(reinterpret_cast<char*>(read_hashes.data()), gs_config.items * sizeof(uint64_t));
	cout << read_hashes.size() << endl;
	for (int i = 0; i < gs_config.items; i++) {
		hash_vec[i].id = i;
		hash_vec[i].value[0] = read_hashes[i];
		valid_items++;
	}
	cerr << valid_items << " valid items in round " << m << endl;
}

void GroupStream::getGroupRes(UnionFind& uf, unordered_map<int, vector<int>>& group_map) {
	uf.findRoot(id_root_map);
	for(int i = 0; i < gs_config.items; i++) {
		int id = i;
		int root_id = id_root_map[id];
		group_map[root_id].push_back(id);
	}
}

void GroupStream::cutEdges(
	vector<vector<int>>& sequences_collisions, 
	int huge_groups_cnt, // 需要多线程libcdhit的组的个数
	const unordered_map<uint64_t, string>& fa_map
	) {
	int avail_threads = gs_config.num_threads;
 	omp_set_num_threads(avail_threads);
	int mt_seqs = 0;
	cout << "avail threads: " << omp_get_max_threads() << endl;
    cerr << "Collisions number processed in multi-thread: " << huge_groups_cnt << endl;
	cerr << "Collisions number processed in single-thread: " << (sequences_collisions.size() - huge_groups_cnt) << endl;

	// huge collisions in multiple thread libcdhit
	auto start_huge_time = chrono::high_resolution_clock::now();
	ClusterWS ws_for_hugegroup;
	for(int i = 0; i < huge_groups_cnt; i++) {
		buildConnectedComponents(sequences_collisions[i], avail_threads, fa_map, ws_for_hugegroup);
		mt_seqs += sequences_collisions[i].size();
	}
	auto end_huge_time = chrono::high_resolution_clock::now();
	auto duration_huge = chrono::duration_cast<chrono::seconds>(end_huge_time - start_huge_time).count();
	cerr << "multi-thread libcdhit (use all threads once): " << duration_huge << endl;
	
	// small collisions in single thread libcdhit
	auto start_small_time = chrono::high_resolution_clock::now();
    #pragma omp parallel num_threads(avail_threads) 
	{
		ClusterWS ws;
		#pragma omp for schedule(dynamic)
    	for(int i = huge_groups_cnt; i < sequences_collisions.size(); i++) {
			buildConnectedComponents(sequences_collisions[i], 1, fa_map, ws);
    	}
	}
	auto end_small_time = chrono::high_resolution_clock::now();
	auto duration_small = chrono::duration_cast<chrono::seconds>(end_small_time - start_small_time).count();
	cerr << "single-thread libcdhit (use only 1 threads each group): " << duration_small << endl;

	cerr << "Seqs number processed in mt_libcdhit: " << mt_seqs << endl;
	cerr << "Seqs number processed in st_libcdhit: " << gs_config.items - mt_seqs << endl;
}

void GroupStream::Cluster(
		vector<vector<int>>& cluster_sequences,
		const unordered_map<uint64_t, string>& fa_map
	) {
    init_cnt();
	std::atomic<int> thread_pool;
 	int TOTAL_THREADS;
 	TOTAL_THREADS = gs_config.num_threads;
    thread_pool = TOTAL_THREADS;
 	omp_set_num_threads(TOTAL_THREADS);
 	omp_set_nested(1);
 	vector<Task> tasks;


	vector<vector<int>> temp_cluster_sequences;
	vector<int>temp_temp_cluster_sequences;
	int count=0;
    int max_threads = omp_get_max_threads()-1;
	for(int i=0;i<cluster_sequences.size();i++){
		if (cluster_sequences[i].size()>=100000)
		{
			if(cluster_sequences[i].size() >= 10000000){
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 60);
			} else if(cluster_sequences[i].size() >= 1000000) {
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 32);
			} else if(cluster_sequences[i].size() >= 500000){
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 16);
			} else if(cluster_sequences[i].size() >= 100000){
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 8);
			}else{
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 4);
			}
			
		} else {
			if(cluster_sequences[i].size()<10000){
				temp_temp_cluster_sequences.insert(
					temp_temp_cluster_sequences.end(),
					cluster_sequences[i].begin(),
					cluster_sequences[i].end()
				);
				if(temp_temp_cluster_sequences.size()>=10000){
					temp_cluster_sequences.emplace_back(temp_temp_cluster_sequences);
					temp_temp_cluster_sequences.clear();
					count++;
					if(count >=1){
						tasks.emplace_back(temp_cluster_sequences,1);
						count=0;
						temp_cluster_sequences.clear();
					}
				}
			}else{
				temp_cluster_sequences.emplace_back(cluster_sequences[i]);
				count++;
				if(count >=1){
					tasks.emplace_back(temp_cluster_sequences,1);
					count=0;
					temp_cluster_sequences.clear();
				}
			}
		}
	}

//		} else {
//			tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 1);
//		}
/**
 * gyj old version
		if (cluster_sequences[i].size()>100000)
		{
			tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]},16);
		}
		else
		{
			if(cluster_sequences[i].size()<100000){
				temp_temp_cluster_sequences.insert(
						temp_temp_cluster_sequences.end(),
						cluster_sequences[i].begin(),
						cluster_sequences[i].end()
						);
				if(temp_temp_cluster_sequences.size()>=100000){
					temp_cluster_sequences.emplace_back(temp_temp_cluster_sequences);
					temp_temp_cluster_sequences.clear();
					count++;
					if(count >=1){
						tasks.emplace_back(temp_cluster_sequences,1);
						count=0;
						temp_cluster_sequences.clear();
					}
				}
			} else{
				temp_cluster_sequences.emplace_back(cluster_sequences[i]);
				count++;
				if(count >=1){
						tasks.emplace_back(temp_cluster_sequences,1);
						count=0;
						temp_cluster_sequences.clear();
				}
			}
		}
**/
	if(!temp_temp_cluster_sequences.empty()){
		temp_cluster_sequences.emplace_back(temp_temp_cluster_sequences);
	}
	if (!temp_cluster_sequences.empty())
	{
		tasks.emplace_back(temp_cluster_sequences, 1);
	}
    std::sort(tasks.begin(), tasks.end(), [](const Task& a, const Task& b) {
                return a.required_threads < b.required_threads;
                    });
	cerr<<"--------------------------"<<endl;
	cerr<<"task size      "<<tasks.size()<<endl;
	auto timestart = chrono::high_resolution_clock::now();
#pragma omp parallel
{
#pragma omp single
{
	for (auto& task : tasks) {
#pragma omp task firstprivate(task)
{
		// // 等待足够的线程资源
		while (true) {
			int available = thread_pool.load(std::memory_order_relaxed);
			if (available >= task.required_threads) {
				int prev = thread_pool.fetch_sub(task.required_threads, std::memory_order_acquire);
				if (prev >= task.required_threads) break;
				thread_pool.fetch_add(task.required_threads, std::memory_order_release);
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}

		// 执行任务

		for(int i=0;i<task.task_cluster.size();i++){
			clusterEachGroup(task.task_cluster[i],task.required_threads, fa_map);

		}
		// 释放线程资源
		thread_pool.fetch_add(task.required_threads, std::memory_order_release);
}
	}
	
#pragma omp taskwait
	//printf("All tasks complete.\n");
}
}
	auto timeend = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(timeend - timestart).count();
	cerr << "cdhit cluster time: " << duration << endl;
    // 打印时间
	/*
    if(tasks_cnt[10000000] > 0){
        int x = 10000000;
        cerr << "大于10,000,000: " ;
        cerr << "任务个数: " << tasks_cnt[x] << endl;
        cerr << "聚类总时间: " << cdhit_cnt[x] << " ";
        cerr << "聚类平均时间: " << (cdhit_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类前准备总时间: " << build_cnt[x] << " ";
        cerr << "聚类前准备平均时间: " << (build_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类后更新总时间: " << update_cnt[x] << " ";
        cerr << "聚类后平均时间: " << (update_cnt[x] / tasks_cnt[x]) << endl;
    }
    if(tasks_cnt[5000000] > 0){
        int x = 5000000;
        cerr << "大于5,000,000: " ;
        cerr << "任务个数: " << tasks_cnt[x] << endl;
        cerr << "聚类总时间: " << cdhit_cnt[x] << " ";
        cerr << "聚类平均时间: " << (cdhit_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类前准备总时间: " << build_cnt[x] << " ";
        cerr << "聚类前准备平均时间: " << (build_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类后更新总时间: " << update_cnt[x] << " ";
        cerr << "聚类后平均时间: " << (update_cnt[x] / tasks_cnt[x]) << endl;
    }
    if(tasks_cnt[1000000] > 0){
        int x = 1000000;
        cerr << "大于1,000,000: " ;
        cerr << "任务个数: " << tasks_cnt[x] << endl;
        cerr << "聚类总时间: " << cdhit_cnt[x] << " ";
        cerr << "聚类平均时间: " << (cdhit_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类前准备总时间: " << build_cnt[x] << " ";
        cerr << "聚类前准备平均时间: " << (build_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类后更新总时间: " << update_cnt[x] << " ";
        cerr << "聚类后平均时间: " << (update_cnt[x] / tasks_cnt[x]) << endl;
    }
    if(tasks_cnt[500000] > 0){
        int x = 500000;
        cerr << "大于500,000: " ;
        cerr << "任务个数: " << tasks_cnt[x] << endl;
        cerr << "聚类总时间: " << cdhit_cnt[x] << " ";
        cerr << "聚类平均时间: " << (cdhit_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类前准备总时间: " << build_cnt[x] << " ";
        cerr << "聚类前准备平均时间: " << (build_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类后更新总时间: " << update_cnt[x] << " ";
        cerr << "聚类后平均时间: " << (update_cnt[x] / tasks_cnt[x]) << endl;
    }
    if(tasks_cnt[100000] > 0){
        int x = 100000;
        cerr << "大于100,000: " ;
        cerr << "任务个数: " << tasks_cnt[x] << endl;
        cerr << "聚类总时间: " << cdhit_cnt[x] << " ";
        cerr << "聚类平均时间: " << (cdhit_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类前准备总时间: " << build_cnt[x] << " ";
        cerr << "聚类前准备平均时间: " << (build_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类后更新总时间: " << update_cnt[x] << " ";
        cerr << "聚类后平均时间: " << (update_cnt[x] / tasks_cnt[x]) << endl;
    }
    if(tasks_cnt[50000] > 0){
        int x = 50000;
        cerr << "大于50,000: " ;
        cerr << "任务个数: " << tasks_cnt[x] << endl;
        cerr << "聚类总时间: " << cdhit_cnt[x] << " ";
        cerr << "聚类平均时间: " << (cdhit_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类前准备总时间: " << build_cnt[x] << " ";
        cerr << "聚类前准备平均时间: " << (build_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类后更新总时间: " << update_cnt[x] << " ";
        cerr << "聚类后平均时间: " << (update_cnt[x] / tasks_cnt[x]) << endl;
    }
    if(tasks_cnt[10000] > 0){
        int x = 10000;
        cerr << "大于10,000: " ;
        cerr << "任务个数: " << tasks_cnt[x] << endl;
        cerr << "聚类总时间: " << cdhit_cnt[x] << " ";
        cerr << "聚类平均时间: " << (cdhit_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类前准备总时间: " << build_cnt[x] << " ";
        cerr << "聚类前准备平均时间: " << (build_cnt[x] / tasks_cnt[x]) << endl;
        cerr << "聚类后更新总时间: " << update_cnt[x] << " ";
        cerr << "聚类后平均时间: " << (update_cnt[x] / tasks_cnt[x]) << endl;
    }
	*/
}

void GroupStream::countGroupSize(int m, UnionFind& uf, const unordered_map<uint64_t, string>& fa_map) {
// FIXME:用结构体GroupNode存储id-root的映射还是用hash_vec继续存
// 用GroupNode增加内存但是如果排序的话要搬移的数据少
	uf.findRoot(id_root_map);
// FIXME:用map来统计还是排序后统计
// 1.用map来统计分组结果 增加内存 只遍历一次
	unordered_map<int, vector<int>> map;

	for(int i = 0; i < gs_config.items; i++) {
		map[id_root_map[i]].push_back(i);
	}
	//统计超过cluster—condition的组的个数
	vector<vector<int>> cluster_sequences;
	for(auto &[key, seqs] : map){
		if(seqs.size() > gs_config.cluster_condition) {
			cluster_sequences.emplace_back(seqs);
		}
	}
	//cout << round_cnt << " " << m << endl;

	if(gs_config.cluster_on && (gs_config.final_cluster_on && round_cnt == gs_config.M - gs_config.R || cluster_sequences.size() > 0)) {
	    cerr << "---------------------------------------------------" << endl;
		if(gs_config.final_cluster_on && round_cnt == gs_config.M - gs_config.R){
			cerr << "Use cdhit in Clustering: groups size larger than " << gs_config.cluster_condition << " number of sequences: " << cluster_sequences.size() << endl;
		}else if(cluster_sequences.size() > 0){
			cerr << "start rescure: groups larger than" << gs_config.cluster_condition << " : " << cluster_sequences.size() << endl;
		}

		sort(cluster_sequences.begin(), cluster_sequences.end(), [](const vector<int>& a, const vector<int>& b){
				return a.size() > b.size();
			});
		cerr << "Before clustering, Top 10 largest group size is: ";
		for(int i = 0; i < std::min(10, (int)cluster_sequences.size()); i++){
			cerr << cluster_sequences[i].size() << " ";
		}
		cerr << endl;

		Cluster(cluster_sequences, fa_map);

		uf.updateParent(id_root_map);

		map.clear();

		for(int i = 0; i < gs_config.items; i++) {
			map[id_root_map[i]].push_back(i);
		}
	}
	cerr << "Top 10 largest group size in this round is: ";

	priority_queue<int, vector<int>, greater<int>> minHeap;
	for(auto &[root_id, seqs] : map){
		minHeap.push(seqs.size());
		if (minHeap.size() > 10){
			 minHeap.pop();
		}
	}
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;
	cerr << "---------------------------------------------------" << endl;

}

void GroupStream::countGroupSizeBySort(UnionFind& uf) {
	priority_queue<int, vector<int>, greater<int>> minHeap;
//	2.直接排序统计分组结果 不增加内存 但多了排序的时间
//	sort(id_root_map.begin(), id_root_map.end(), [](const GroupNode& a, const GroupNode& b){
//		return a.root < b.root;
//		});
//	int group_size = 0;
//	int cur_root = id_root_map[0].root;//  第一个root-id
//    for (const auto& p : id_root_map) {
//    	if(cur_root == p.root) {
//			group_size++;
//		}else{
//			minHeap.push(group_size);
//			if (minHeap.size() > 10){
//				 minHeap.pop();
//			}
//			cur_root = p.root;
//			group_size = 1;
//		}
//    }
//	if(group_size >= 1) {
//		minHeap.push(group_size);
//		if (minHeap.size() > 10){
//			 minHeap.pop();
//		}
//	}
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;
}

void GroupStream::Group(
    const ProteinSketchData& sketchdata,
	const ProteinData& proteindata
	) {
	for(int m=0; m < gs_config.M - gs_config.R+1; m++){
		cerr << "round "<<  m << endl;
		fillHashVec(sketchdata, hash_vec, m);
		GroupByCol(hash_vec, proteindata.sequence_map);
		if(m == gs_config.M - gs_config.R && gs_config.final_cluster_on) {
			gs_config.cluster_condition = 1;
		}
		countGroupSize(m, uf, proteindata.sequence_map);
		round_cnt++;
	}

	//getGroupRes(uf, group_map);
	if(gs_config.output_on) {
		outputClstr(proteindata.names, proteindata.sequence_map);
	}
}

void GroupStream::Group(
	string sketch_filename,
	const ProteinData& proteindata
	) {
	for(int m=0; m < gs_config.M-gs_config.R+1; m++){
		cerr << "round "<<  m << endl;
		fillHashVec(sketch_filename, hash_vec, m);
		GroupByCol(hash_vec, proteindata.sequence_map);
		return;
		if(m == gs_config.M-gs_config.R && gs_config.final_cluster_on) {
			gs_config.cluster_condition = 1;
		}
		countGroupSize(m, uf, proteindata.sequence_map);
		round_cnt++;
	}

	if(gs_config.output_on) {
		outputClstr(proteindata.names, proteindata.sequence_map);
	}
}

void GroupStream::buildConnectedComponents(
	vector<int>& group_seqs, 
	int needed_threads,
	const unordered_map<uint64_t, string>& fa_map,
	ClusterWS& ws
	) {
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		sequences.emplace_back(group_seqs[i], fa_map.at(group_seqs[i]).c_str());
	}
	if(needed_threads > 1) {
		cluster_sequences(sequences, 5, tau, needed_threads); 
	}else {
		if(group_seqs.size() < 1000){
			cluster_sequences_st_less10(sequences, 5, tau); 
		}else{
			cluster_sequences_st_reuse(sequences, 5, tau, ws); 
		}
	}
	for(int i = 0; i < group_seqs.size(); i++)
	{
		id_root_map[sequences[i].seq_id] = sequences[i].new_root_id;
	}
}

void GroupStream::clusterEachGroup(
	vector<int>& group_seqs,
	int needed_threads,
	const unordered_map<uint64_t, string>& fa_map
	) {
	auto start_time_build = chrono::high_resolution_clock::now();
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		sequences.emplace_back(group_seqs[i], fa_map.at(group_seqs[i]).c_str());
	}
	auto end_time_build = chrono::high_resolution_clock::now();
    auto duration_build = chrono::duration_cast<chrono::seconds>(end_time_build - start_time_build).count();

	//读取FAI获取data

	//auto start_time = chrono::high_resolution_clock::now();
		cluster cluster_cdhit;
		cluster_cdhit.cdhit_cluster(sequences, id_root_map, needed_threads);
	//auto end_time = chrono::high_resolution_clock::now();
	//auto duration_cdhit = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();

    //序列输出ID...
	/*
    if(group_seqs.size() >= 10000000){ // > 10,000,000
        tasks_cnt[10000000]++;
        build_cnt[10000000]+=duration_build;
        cdhit_cnt[10000000]+=duration_cdhit;
    }else if(group_seqs.size() >= 5000000){ // >5,000,000
        tasks_cnt[5000000]++;
        build_cnt[5000000]+=duration_build;
        cdhit_cnt[5000000]+=duration_cdhit;
    }else if(group_seqs.size() >= 1000000){ // >1,000,000
        tasks_cnt[1000000]++;
        build_cnt[1000000]+=duration_build;
        cdhit_cnt[1000000]+=duration_cdhit;
    }else if(group_seqs.size() >= 5000000){ // > 500,000
        tasks_cnt[500000]++;
        build_cnt[500000]+=duration_build;
        cdhit_cnt[500000]+=duration_cdhit;
    }else if(group_seqs.size() >= 100000){ // > 100,000
        tasks_cnt[100000]++;
        build_cnt[100000]+=duration_build;
        cdhit_cnt[100000]+=duration_cdhit;
    }else if(group_seqs.size() >= 50000){ // > 50,000
        tasks_cnt[50000]++;
        build_cnt[50000]+=duration_build;
        cdhit_cnt[50000]+=duration_cdhit;
    }else if(group_seqs.size() >= 10000){ // > 10,000
        tasks_cnt[10000]++;
        build_cnt[10000]+=duration_build;
        cdhit_cnt[10000]+=duration_cdhit;
    }else{
        cerr << "<10,000" << endl;
    }
	*/
}

void GroupStream::outputClstr(
	const vector<string>& names,
	const unordered_map<uint64_t, string>& fa_map
) {
	cerr << "Total Clusters: " << uf.countSetsSize() << endl;
	cerr << "cluster result stored: " << gs_config.res_file << endl;
	ofstream seq_id(gs_config.res_file);
	streambuf* origin_cout = cout.rdbuf();
	cout.rdbuf(seq_id.rdbuf());

	uf.findRoot(id_root_map);
	for(int i = 0; i < gs_config.items; i++) {
		cout << ">" << names[i] << " " << ">" << names[id_root_map[i]] << endl;
	}
	cout.rdbuf(origin_cout);
}
