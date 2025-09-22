#include "GroupStream.h"
#include <queue>
#include <omp.h>
#include <thread>
#include <algorithm>
#include <atomic>
#include <unordered_set>
#include <fstream>

// 静态成员变量的定义

int round_cnt=0;
int max_size=0;

struct minheapcompare {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

bool compareByHash(const Data &a, const Data &b) {
	return a.value < b.value;
}

//vector<int> GroupStream::GroupSizeCnt;
//bool compareByHashAndGroupSize(const Data &a, const Data &b) {
//    if(a.value != b.value)
//	    return a.value < b.value;
//    return GroupSizeCnt[a.id] < GroupSizeCnt[b.id];
//}

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
		string filename = folder_name + to_string(rootid) + ".fa";
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

void GroupStream::Sort(vector<Data>& dataList){

#ifdef parallel
	size_t n = dataList.size();
	int num_threads = omp_get_max_threads();
	size_t chunk_size = n / num_threads;
	cout << num_threads << " threads are used in sorting" << endl;

#pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        size_t start = tid * chunk_size;
        size_t end = (tid == num_threads - 1) ? n : start + chunk_size;
        sort(dataList.begin()+start, dataList.begin()+end, compareByHash);
    }

//#pragma omp parallel for num_threads(num_threads)
	for (int i = 1; i < num_threads; ++i) {
		inplace_merge(dataList.begin(), dataList.begin() + i * chunk_size,
			dataList.begin() + ((i == num_threads - 1) ? n : (i + 1) * chunk_size),
			compareByHash);
	}
#else
//	sort(dataList.begin(), dataList.end(), [](const Data& a, const Data& b){
//		return a.value < b.value;
//		});
		sort(dataList.begin(), dataList.end(), compareByHash);
#endif
}

void GroupStream::checkEdges(vector<Data>& hash_vec, UnionFind& cur_uf) {
	cur_uf.findRoot(id_root_map);
	unordered_map<int, vector<int>> map;

	for(int i = 0; i < items; i++) {
		map[id_root_map[i]].push_back(i);
	}

	int cnt_1000_100 = 0;
	int cnt_100_10 = 0;
	int cnt_10_1 = 0;

	vector<vector<int>> first_hit_sequences;
    int huge_groups_cnt = 0;
	for(auto &[key, seqs] : map){
		if(seqs.size() > 1) {
			first_hit_sequences.emplace_back(seqs);
		}
        if(seqs.size() > 10000) {
            huge_groups_cnt++;
        }
        
		// 统计完全图的大小
		int size = seqs.size();
		//if(size > 1000){
		//	cnt_1000_100++;
		//}else if (size > 10){
		//	cnt_100_10++;
		//}else if (size > 1){
		//	cnt_10_1++;
		//}
	}
	cerr << "MinHash collisions : " << first_hit_sequences.size() << endl;
	//cerr << "    groups size count: " << endl;
	//cerr << "       1000-100: " << cnt_1000_100 << endl;
	//cerr << "       100-10: " << cnt_100_10 << endl;
	//cerr << "       10-1: " << cnt_10_1 << endl;

	sort(first_hit_sequences.begin(), first_hit_sequences.end(),
		[](const vector<int>& a, const vector<int>& b){
			return a.size() > b.size();
		});

	//cerr << "    Top 10 largest groups: " ;
	cerr << "    Top 10 largest collisions: " ;
	for(int i = 0; i < std::min(10, (int)first_hit_sequences.size()); i++){
		cerr << first_hit_sequences[i].size() << " ";
	}
	cerr << endl;

    cut_edges(first_hit_sequences, huge_groups_cnt);
	cur_uf.updateParent(id_root_map);

	unordered_map<int, vector<int>> groups_after_filter; // rootid:[seq0, seq1...]
	priority_queue<int, vector<int>, greater<int>> minHeap;

	for(int i = 0; i < items; i++) {
		groups_after_filter[id_root_map[i]].push_back(i);
	}

	for(auto &[root_id, seqs] : groups_after_filter){
		minHeap.push(seqs.size());
		if (minHeap.size() > 10){
			minHeap.pop();
		}

	}
	cerr << "After break the bad edges, connected components number are: " << groups_after_filter.size() << endl;
	cerr << "    Top 10 largest connected components: " << endl;
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;
}


void GroupStream::Unite(const vector<Data>& dataList, UnionFind& this_uf) {
	vector<uint64_t> cur_value = dataList[0].value;
	for (int i = 1; i < valid_items; i++) {
		auto thisone = dataList[i];
    	auto lastone = dataList[i-1];
		if(thisone.value == lastone.value) {
		    this_uf.unite(thisone.id, lastone.id);
		}
    }
}


void GroupStream::unite_by_edges(UnionFind& cur_uf) {
	for(int i = 0; i < items; i++) {
		uf.unite(uf.find(i), cur_uf.find(i));
	}
}
								
void GroupStream::GroupByCol(vector<Data>& hash_vec, UnionFind& uf) {
	Sort(hash_vec);
	UnionFind col_uf(items);
	Unite(hash_vec, col_uf);
	checkEdges(hash_vec, col_uf);
	unite_by_edges(col_uf);
	int groups_size = uf.countSetsSize();
	cerr << "---------------------------------------------------" << endl;
	cout << "Group Size after merging:" << groups_size << endl;
}

void GroupStream::fillHashVec(const vector<vector<uint64_t>>& vec, vector<Data>& hash_vec, int m) {
// 1. use std::transform
//	int index = 0;
//	transform(vec.begin(), vec.end(), hash_vec.begin(), 
//		[&index](int value) {
//		return Data{value, (index++)};
//		});
// 2. iteration to construct pair
//    for (int i = 0; i < vec.size(); ++i) {
//        hash_vec[i].id = i;
//        hash_vec[i].value = vec[i];
//    }
//
//	if(rep_only_group){
//		valid_items = 0;
//		for (int i = 0; i < items; i++) {
//			if(valid_seqs[seq_ids[i]]){
//				hash_vec[valid_items].id = seq_ids[i];
//				copy(vec[i].begin() + m * L, vec[i].begin() + m * L + R * L, hash_vec[valid_items].value.begin());
//				valid_items++;
//			}
//		}
//		cerr << valid_items << " valid items in round " << m << endl;
//	}else{
//		for (int i = 0; i < items; i++) {
//			//hash_vec[i].id = i;
//			hash_vec[i].id = seq_ids[i];
//			copy(vec[i].begin() + m * L, vec[i].begin() + m * L + R * L, hash_vec[i].value.begin());
//		}
//		cerr << valid_items << " valid items in round " << m << endl;
//	}


	valid_items = 0;
	std::ifstream ifs(sketch_file_name, ios::binary);
	if(!ifs){
		cerr << "Error opening file!" << endl;
		return;
	}
	ifs.seekg(m * items * sizeof(uint64_t), std::ios::beg);
	vector<uint64_t> read_hashes(items);
	ifs.read(reinterpret_cast<char*>(read_hashes.data()), items * sizeof(uint64_t));
	cout << read_hashes.size() << endl;
	for (int i = 0; i < items; i++) {
		hash_vec[i].id = i;
		hash_vec[i].value[0] = read_hashes[i];
		valid_items++;
	}
	cerr << valid_items << " valid items in round " << m << endl;
}

void GroupStream::get_group_res(UnionFind& uf, unordered_map<int, vector<int>>& group_map) {
	uf.findRoot(id_root_map);
	for(int i = 0; i < items; i++) {
		int id = i;
		int root_id = id_root_map[id];
		group_map[root_id].push_back(id);
	}
}

void GroupStream::cut_edges(vector<vector<int>>& sequences_collisions, int huge_groups_cnt) {
	int avail_threads = num_threads > 1 ? num_threads - 1 : 1;
 	omp_set_num_threads(avail_threads);

	// huge collisions
    cerr << "Collisions size large than 10000: " << huge_groups_cnt << endl;
    cerr << "Collisions size small than 10000: " << (sequences_collisions.size() - huge_groups_cnt) << endl;
	auto start_huge_time = chrono::high_resolution_clock::now();
    for(int i = 0; i < huge_groups_cnt; i++) {
    	build_connected_components(sequences_collisions[i], avail_threads);
    }
	auto end_huge_time = chrono::high_resolution_clock::now();
	auto duration_huge = chrono::duration_cast<chrono::seconds>(end_huge_time - start_huge_time).count();
	cerr << "Cut edges time(seconds): " << endl;
	cerr << "    huge collision groups (use all threads once)" << duration_huge << endl;

	// small collisions
	auto start_small_time = chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(int i = huge_groups_cnt; i < sequences_collisions.size(); i++) {
    	build_connected_components(sequences_collisions[i], 1);
    }
	auto end_small_time = chrono::high_resolution_clock::now();
	auto duration_small = chrono::duration_cast<chrono::seconds>(end_small_time - start_small_time).count();
	cerr << "    small collision groups (use only 1 threads each group)" << duration_small << endl;

}
void GroupStream::Cluster(vector<vector<int>>& cluster_sequences) {
    init_cnt();
	std::atomic<int> thread_pool;
 	int TOTAL_THREADS;
 	TOTAL_THREADS=num_threads;
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
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 40);
			} else if(cluster_sequences[i].size() >= 1000000) {
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 32);
			} else if(cluster_sequences[i].size() >= 500000){
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 16);
			} else{
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 8);
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
			clusterEachGroup(task.task_cluster[i],task.required_threads);

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

void GroupStream::countGroupSize(int m, UnionFind& uf) {
// FIXME:用结构体GroupNode存储id-root的映射还是用hash_vec继续存
// 用GroupNode增加内存但是如果排序的话要搬移的数据少
	uf.findRoot(id_root_map);
// FIXME:用map来统计还是排序后统计
// 1.用map来统计分组结果 增加内存 只遍历一次
	unordered_map<int, vector<int>> map;

	for(int i = 0; i < items; i++) {
		map[id_root_map[i]].push_back(i);
	}
	//统计超过cluster—condition的组的个数
	vector<vector<int>> cluster_sequences;
	for(auto &[key, seqs] : map){
		if(seqs.size() > cluster_condition) {
			cluster_sequences.emplace_back(seqs);
		}
	}
	cout << round_cnt << " " << m << endl;

	if(cluster_on && (round_cnt == M-R || cluster_sequences.size() > 0)) {
	    cerr << "---------------------------------------------------" << endl;
		if(round_cnt == M-R){
			cerr << "Final Cluster: groups larger than" << cluster_condition << " : " << cluster_sequences.size() << endl;
		}else if(cluster_sequences.size() > 0){
			cerr << "start rescure: groups larger than" << cluster_condition << " : " << cluster_sequences.size() << endl;
		}

		sort(cluster_sequences.begin(), cluster_sequences.end(), [](const vector<int>& a, const vector<int>& b){
				return a.size() > b.size();
			});
		cerr << "Before clustering, Top 10 largest group size is: ";
		for(int i = 0; i < std::min(10, (int)cluster_sequences.size()); i++){
			cerr << cluster_sequences[i].size() << " ";
		}
		cerr << endl;

		if(threadPool_on){
            max_size=0;
			Cluster(cluster_sequences);
		}else{
			for(int i = 0; i < cluster_sequences.size(); i++) {
					clusterEachGroup(cluster_sequences[i]);
			}
		}

		uf.updateParent(id_root_map);

		map.clear();

		for(int i = 0; i < items; i++) {
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

void GroupStream::Group(vector<vector<uint64_t>>& hashes, unordered_map<int, vector<int>>& group_map) {
	cerr << "==========Group Parameters==========" << endl;
	cerr << "cluster:" << cluster_on << endl;
	cerr << "cluster-condition: " << cluster_condition << endl;
	cerr << "clustering all sequences in last round: " << final_cluster_on << endl;
	cerr << "==========Group Parameters==========" << endl;
	if(slide) {
		for(int m=0; m < M-R+1; m++){
			cerr << "round "<<  m << endl;
			fillHashVec(hashes, hash_vec, m);
			GroupByCol(hash_vec, uf);
			if(m == M-R && final_cluster_on) {
				cluster_condition = 1;
			}
			countGroupSize(m, uf);
            round_cnt++;
		}
	}

	get_group_res(uf, group_map);
	if(output_on) {
		outputClstr();
	}
}

void GroupStream::build_connected_components(vector<int>& group_seqs, int needed_threads)
{
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());
	}
	if(needed_threads > 1) {
		cluster_sequences(sequences, id_root_map, 5, 0.05, needed_threads); 
	}else {
		cluster_sequences_st(sequences, id_root_map, 5, 0.05); 
	}
}

void GroupStream::clusterEachGroup(vector<int>& group_seqs){
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());
	}
	cluster cluster_cdhit;
	cluster_cdhit.cdhit_cluster(sequences, id_root_map, 1);
}

void GroupStream::clusterEachGroup(vector<int>& group_seqs,int needed_threads) {
	auto start_time_build = chrono::high_resolution_clock::now();
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());
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
    if(group_seqs.size() == max_size && group_seqs.size() > 10000000) {
            ostringstream oss;
            ofstream seqs("max_group_seqs");
            streambuf* origin_cout = cout.rdbuf();
            cout.rdbuf(seqs.rdbuf());

            uf.findRoot(id_root_map);
            for (const auto& seq : group_seqs) {
                oss << ">" << names[seq] << "\n";
            }
            cout << oss.str();
            cout.rdbuf(origin_cout);
    }

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

//void GroupStream::setValidStatus(vector<int>& group_seqs){
//	for(int seq : group_seqs){
//		if(seq != id_root_map[seq]){
//			valid_seqs[seq] = false;
//		}else{
//			valid_seqs[seq] = true;
//		}
//	}
//}

void GroupStream::outputClstr() {
	cerr << "cluster result stored: " << res_file << endl;
	ofstream seq_id(res_file);
	streambuf* origin_cout = cout.rdbuf();
	cout.rdbuf(seq_id.rdbuf());

	uf.findRoot(id_root_map);
	for(int i = 0; i < items; i++) {
		cout << ">" << names[i] << " " << ">" << names[id_root_map[i]] << endl;
	}
	cout.rdbuf(origin_cout);
}
