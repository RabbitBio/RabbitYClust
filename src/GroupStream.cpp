#include "GroupStream.h"
#include <queue>
#include <omp.h>
#include <thread>
#include <algorithm>
#include <atomic>

struct minheapcompare {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

bool compareByHash(const Data &a, const Data &b) {
	return a.value < b.value;
}

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
	cout << num_threads << " threads are working" << endl;
	
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
		//sort(dataList.begin(), dataList.end(), compareByHash);
		sort(dataList.begin(), dataList.begin() + valid_items, compareByHash);
#endif
}

void GroupStream::Unite(const vector<Data>& dataList, UnionFind& uf) {
#ifdef VERBOSE
	int count = 1;
#endif
	vector<uint64_t> cur_value = dataList[0].value;
	int cur_head = dataList[0].id;
	for (int i = 0; i < valid_items; i++) {
		auto data = dataList[i];
		if(data.value == cur_value) {
			uf.unite(data.id, cur_head);
		}else{
#ifdef VERBOSE
			count++;
#endif
			cur_value = data.value;
			cur_head = data.id;
		}
	}
#ifdef VERBOSE
	cerr << "groups number in this col is: " << count << endl;
#endif
}

void GroupStream::GroupByCol(vector<Data>& hash_vec, UnionFind& uf) {
#ifdef TIMING
	auto sort_start_time = chrono::high_resolution_clock::now();
	Sort(hash_vec);
	auto sort_end_time = chrono::high_resolution_clock::now();
	auto sort_duration = chrono::duration_cast<chrono::seconds>(sort_end_time - sort_start_time).count();
	cout << "sort time: "<< sort_duration << endl;
	auto union_start_time = chrono::high_resolution_clock::now();
	Unite(hash_vec, uf);
	auto union_end_time = chrono::high_resolution_clock::now();
	auto union_duration = chrono::duration_cast<chrono::seconds>(union_end_time - union_start_time).count();
	cout << "union time: " << " is " << union_duration << endl;

	int groups_size = uf.countSetsSize();
	cout << "Group Size is " << groups_size << endl;
#else
	auto timestart = chrono::high_resolution_clock::now();
	Sort(hash_vec);
	Unite(hash_vec, uf);
	auto timeend = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(timeend - timestart).count();
	cerr << "sort and unite time: " << duration << endl;
#endif
#ifdef VERBOSE
	int groups_size = uf.countSetsSize();
	cerr << "Group Size is " << groups_size << endl;
#endif
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
	if(rep_only_group){
		valid_items = 0;
		for (int i = 0; i < items; i++) {
			if(valid_seqs[seq_ids[i]]){
				hash_vec[valid_items].id = seq_ids[i];
				copy(vec[i].begin() + m * L, vec[i].begin() + m * L + R * L, hash_vec[valid_items].value.begin());
				valid_items++;
			}
		}
		cerr << valid_items << " valid items in round " << m << endl;
	}else{
		for (int i = 0; i < items; i++) {
			//hash_vec[i].id = i;
			hash_vec[i].id = seq_ids[i];
			copy(vec[i].begin() + m * L, vec[i].begin() + m * L + R * L, hash_vec[i].value.begin());
		}
		cerr << valid_items << " valid items in round " << m << endl;
	}

// 3. memory alignment and memcpy

}

void GroupStream::getGroupMap(UnionFind& uf, unordered_map<int, vector<int>>& group_map) {
	uf.findRoot(id_root_map);
//  1.用vector<GroupNode>存储id-root-map时候要对其进行kj
//	sort(id_root_map.begin(), id_root_map.end(), [](const GroupNode& a, const GroupNode& b){
//		return a.root < b.root; 
//		});
//    for (const auto& p : id_root_map) {
//    	group_map[p.root].push_back(p.id); 
//    }

	for(int i = 0; i < items; i++) {
		int id = seq_ids[i];
		int root_id = id_root_map[id];
		group_map[root_id].push_back(id);
	}
}

void GroupStream::Cluster(int m, vector<vector<int>>& cluster_sequences) {
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
//	if(m == M-R){
//	for(int i=0;i<cluster_sequences.size();i++){
//		if (cluster_sequences[i].size()>=100000)
//		{
//			if(cluster_sequences[i].size() >= 10000000){
//				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 60);
//			} else if(cluster_sequences[i].size() >= 1000000) {
//				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 32);
//			} else if(cluster_sequences[i].size() >= 500000){
//				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 16;
//			} else{
//				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 8);
//			}
//			
//		} else {
//			if(cluster_sequences[i].size()<10000){
//				temp_temp_cluster_sequences.insert(
//					temp_temp_cluster_sequences.end(),
//					cluster_sequences[i].begin(),
//					cluster_sequences[i].end()
//				);
//				if(temp_temp_cluster_sequences.size()>=10000){
//					temp_cluster_sequences.emplace_back(temp_temp_cluster_sequences);
//					temp_temp_cluster_sequences.clear();
//					count++;
//					if(count >=1){
//						tasks.emplace_back(temp_cluster_sequences,1);
//						count=0;
//						temp_cluster_sequences.clear();
//					}
//				}
//			}else{
//				temp_cluster_sequences.emplace_back(cluster_sequences[i]);
//				count++;
//				if(count >=1){
//					tasks.emplace_back(temp_cluster_sequences,1);
//					count=0;
//					temp_cluster_sequences.clear();
//				}
//			}
//		}
//	}
//	}else{
	for(int i=0;i<cluster_sequences.size();i++){
		if (cluster_sequences[i].size()>=100000)
		{
			if(cluster_sequences[i].size() >= 10000000){
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 60);
			} else if(cluster_sequences[i].size() >= 1000000) {
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 32);
			} else if(cluster_sequences[i].size() >= 500000){
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 16);
			} else{
				tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 8);
			}

		} else {
			tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]}, 1);
		}
	}
//	}
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
		// cerr<<"thread_pool     "<<thread_pool<<endl;
}
	}
	
#pragma omp taskwait
	printf("All tasks complete.\n");
}
}
	auto timeend = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(timeend - timestart).count();
	cerr << "cdhit cluster time: " << duration << endl;
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

	if(cluster_on) {
		int seqs_count = 0;
		vector<vector<int>> cluster_sequences;
		for(auto &[key, seqs] : map){
			if(seqs.size() > cluster_condition) {
				cluster_sequences.emplace_back(seqs);
			}
		}
		cerr << "Groups larger than " << cluster_condition << " : " << cluster_sequences.size() << endl;

		sort(cluster_sequences.begin(), cluster_sequences.end(), [](const vector<int>& a, const vector<int>& b){
		return a.size() > b.size();
		});
		cerr << "Top 10 largest group size is: ";
		for(int i = 0; i < std::min(10, (int)cluster_sequences.size()); i++){
			cerr << cluster_sequences[i].size() << " ";
		}
		cerr << endl;

		total_clusters = 0;
		redundant_seqs = 0;

		if(threadPool_on){
			Cluster(m, cluster_sequences);
		}else{
			//#pragma omp parallel for num_threads(num_threads)
			for(int i = 0; i < cluster_sequences.size(); i++) {
					//cerr << i << " is doing cluster " << cluster_sequences[i].size() << " sequences" << endl;
					clusterEachGroup(cluster_sequences[i]);
			}
		}
//		if(rep_on){
//			cerr << "聚类了 " << total_clusters << " 个类" << endl;
//			cerr << "共去除了冗余序列 " << redundant_seqs - total_clusters << "条" << endl;
//		}
		uf.updateParent(id_root_map);

		unordered_map<int, vector<int>> map_after_cluster;
		priority_queue<int, vector<int>, greater<int>> minHeap;

		for(int i = 0; i < items; i++) {
			map_after_cluster[id_root_map[i]].push_back(i);
		}
		int largethan1w = 0;
		for(auto &[root_id, seqs] : map_after_cluster){
			if(seqs.size() > cluster_condition)
				largethan1w++;
			minHeap.push(seqs.size());
			if (minHeap.size() > 10){
				 minHeap.pop();
			}

		}
		cerr << "After clustering, clusters size large than " << cluster_condition << " : "<< largethan1w << endl;
		while(!minHeap.empty()){
//			if(temp_output_on){
//				int rootid = minHeap.top().second;
//				string filename = folder_name + to_string(rootid) + ".fa";
//				ofstream ofile(filename);
//				for(int id : map_after_cluster[rootid]){
//					ofile << ">" << id << endl;
//					ofile << fa_map[id] << endl;
//				}
//				ofile.close();
//				cerr << "cluster: " << rootid << " contains " << minHeap.top().first << " seqs stored in: " << filename << endl; 
//			}
			cerr << minHeap.top() << " ";
			minHeap.pop();
		}
		cerr << endl;

	}else {
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
	}
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
	cerr << "only use reps in grouping: " << rep_only_group << endl;
	cerr << "only use reps in clustering: " << rep_only_cluster << endl;
	cerr << "clustering all sequences in last round: " << final_cluster_on << endl;
	cerr << "==========Group Parameters==========" << endl;
	if(slide) {
		for(int m=0; m < M-R+1; m++){
			cerr << "round "<<  m << endl;
			fillHashVec(hashes, hash_vec, m);
			GroupByCol(hash_vec, uf);
			if(m == M-R && final_cluster_on) {
				cluster_condition = 1;
				rep_only_cluster = false;
			}
			countGroupSize(m, uf);
		}
	}else{
		for(int m=0; m < M / R; m++){
			cerr << "round "<<  m << endl;
			fillHashVec(hashes, hash_vec, m * R);
			GroupByCol(hash_vec, uf);
			countGroupSize(m, uf);
		}
		if(M % R != 0){
			cerr << "round "<<  M / R;
			setR(M % R);
			fillHashVec(hashes, hash_vec, (M/R) * R);
			GroupByCol(hash_vec, uf);
			countGroupSize(M, uf);
		}
	}
#ifdef TIMING
	auto transpose_start_time = chrono::high_resolution_clock::now();
	getGroupMap(uf, group_map);
	auto transpose_end_time = chrono::high_resolution_clock::now();
	auto transpose_duration = chrono::duration_cast<chrono::seconds>(transpose_end_time - transpose_start_time).count();
	cout << "transpose time is " << transpose_duration << endl;
#else
	getGroupMap(uf, group_map);
#endif
	if(output_on) {
		outputClstr();
	}
}

void GroupStream::clusterEachGroup(vector<int>& group_seqs){
	vector<Sequence_new> sequences;
	if(rep_only_cluster && group_seqs.size() > cluster_condition ){
		for(int i = 0; i < group_seqs.size(); i++) {
			if(valid_seqs[group_seqs[i]]){
				sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());
			}
		}
	}else{
		for(int i = 0; i < group_seqs.size(); i++) {
			sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());
		}
	}
	//读取FAI获取data

	cluster cluster_cdhit;
	cluster_cdhit.cdhit_cluster(sequences, id_root_map);

	//从中挑选出代表序列作为以后分组和聚类的唯一代表
	if(rep_only_group){
		redundant_seqs += sequences.size();
		setValidStatus(group_seqs);
	}
}
void GroupStream::clusterEachGroup(vector<int>& group_seqs,int neededThread) {
	vector<Sequence_new> sequences;
	if(rep_only_cluster && group_seqs.size() > cluster_condition){
		for(int i = 0; i < group_seqs.size(); i++) {
			if(valid_seqs[group_seqs[i]]){
				sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());
			}
		}
	}else{
		for(int i = 0; i < group_seqs.size(); i++) {
			sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());
		}
	}

	//读取FAI获取data

	cluster cluster_cdhit;
	cluster_cdhit.cdhit_cluster(sequences, id_root_map, neededThread);

	if(rep_only_group){
//		redundant_seqs += sequences.size();
		setValidStatus(group_seqs);
	}

}

void GroupStream::setValidStatus(vector<int>& group_seqs){
	for(int seq : group_seqs){
		if(seq != id_root_map[seq]){
			valid_seqs[seq] = false;
		}else{
			valid_seqs[seq] = true;
//			total_clusters++;
		}

	}
}

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
