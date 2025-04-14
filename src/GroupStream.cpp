#include "GroupStream.h"
#include <queue>
#include <omp.h>
#include <thread>
#include <algorithm>

struct minheapcompare {
	bool operator()(const pair<int, int>& a, const pair<int, int>& b) {
		return a.first > b.first;
	}
};

bool compareByHash(const Data& a, const Data& b) {
	return a.value < b.value;
}

bool compareById(const Data& a, const Data& b) {
	return a.id < b.id;
}
double Gettime_ms(){
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return (double)tv.tv_sec *1000+ (double) tv.tv_usec/1000;
	
}
void GroupStream::tempOutput(vector<vector<int>>& cluster_sequences) {
	unordered_map<int, vector<int>> map_after_cluster;
	priority_queue<pair<int, int>, vector<pair<int, int>>, minheapcompare> minHeap;
	for (int i = 0; i < cluster_sequences.size(); i++) {
		for (int x : cluster_sequences[i]) {
			map_after_cluster[id_root_map[x]].push_back(x);
		}
	}
	for (auto& [root_id, seqs] : map_after_cluster) {
		if (seqs.size() > 10000) {
			minHeap.push({ seqs.size(), root_id });
		}
	}
	ofstream log("output.log");
	cerr.rdbuf(log.rdbuf());
	while (!minHeap.empty()) {
		int rootid = minHeap.top().second;
		string filename = folder_name + to_string(rootid) + ".fa";
		ofstream ofile(filename);
		for (int id : map_after_cluster[rootid]) {
			ofile << ">" << id << endl;
			ofile << fa_map[id] << endl;
		}
		ofile.close();
		cerr << "cluster: " << rootid << " contains " << minHeap.top().first << " seqs stored in: " << filename << endl;
		minHeap.pop();
	}
	cerr << endl;
}

void GroupStream::Sort(vector<Data>& dataList) {

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
		sort(dataList.begin() + start, dataList.begin() + end, compareByHash);
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

void GroupStream::Unite(vector<Data> dataList, UnionFind& uf) {
#ifdef VERBOSE
	int count = 1;
#endif
	vector<uint64_t> cur_value = dataList[0].value;
	int cur_head = dataList[0].id;
	for (const auto& data : dataList) {
		if (data.value == cur_value) {
			uf.unite(data.id, cur_head);
		}
		else {
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
	cout << "sort time: " << sort_duration << endl;
	auto union_start_time = chrono::high_resolution_clock::now();
	Unite(hash_vec, uf);
	auto union_end_time = chrono::high_resolution_clock::now();
	auto union_duration = chrono::duration_cast<chrono::seconds>(union_end_time - union_start_time).count();
	cout << "union time: " << " is " << union_duration << endl;

	int groups_size = uf.countSetsSize();
	cout << "Group Size is " << groups_size << endl;
#else
	Sort(hash_vec);
	Unite(hash_vec, uf);
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
	//	cerr << m << " " << m*L << " " << m * L + R << endl; 
	for (int i = 0; i < items; i++) {
		//hash_vec[i].id = i;
		hash_vec[i].id = seq_ids[i];
		copy(vec[i].begin() + m * L, vec[i].begin() + m * L + R * L, hash_vec[i].value.begin());
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

	for (int i = 0; i < items; i++) {
		int id = seq_ids[i];
		int root_id = id_root_map[id];
		group_map[root_id].push_back(id);
	}
}
void GroupStream::use_thread_pool(vector<vector<int>>& cluster_sequences){
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
		for(int i=0;i<cluster_sequences.size();i++){
			if (cluster_sequences[i].size()>100000)
			{
				if(cluster_sequences[i].size()>1000000){
					tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]},16);
				}
				else{
					tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]},4);
				}
				
			}
			// else if(cluster_sequences[i].size()<1000)
			else
			{
					if(cluster_sequences[i].size()<10000){
						temp_temp_cluster_sequences.insert(
							temp_temp_cluster_sequences.end(),
							cluster_sequences[i].begin(),
							cluster_sequences[i].end()
						);
						// cerr<<"temp_temp_cluster_sequences size    "<<temp_temp_cluster_sequences.size()<<endl;
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
					}
					else{
						temp_cluster_sequences.emplace_back(cluster_sequences[i]);
						count++;
						if(count >=1){
							tasks.emplace_back(temp_cluster_sequences,1);
							count=0;
							temp_cluster_sequences.clear();
						}
					}
			}
			// else{
			// 	tasks.emplace_back(std::vector<vector<int>>{cluster_sequences[i]},1);
				
			// }
			
		}
		if(!temp_temp_cluster_sequences.empty()){
			temp_cluster_sequences.emplace_back(temp_temp_cluster_sequences);
		}
		if (!temp_cluster_sequences.empty())
		{
			tasks.emplace_back(temp_cluster_sequences, 1);
		}
		cerr<<"--------------------------"<<endl;
		cerr<<"task size      "<<tasks.size()<<endl;
		// cerr<<"block_condition      "<<block_condition<<endl;
		double t0 =Gettime_ms();
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
		double t1=Gettime_ms();
		cerr<<"Total time   "<<t1-t0<<endl;


}
void GroupStream::countGroupSize(UnionFind& uf, int m, vector<vector<uint64_t>>& hashes) {
	// FIXME:用结构体GroupNode存储id-root的映射还是用hash_vec继续存
	// 用GroupNode增加内存但是如果排序的话要搬移的数据少
	uf.findRoot(id_root_map);
	// FIXME:用map来统计还是排序后统计
	// 1.用map来统计分组结果 增加内存 只遍历一次
	unordered_map<int, vector<int>> map;
	
	for (int i = 0; i < items; i++) {
		map[id_root_map[i]].push_back(i);
	}

	if (second_group) {
		vector<vector<int>> cluster_sequences;
		// 用于记录大小超过指定阈值的序列集合
		for (auto& [key, seqs] : map) {
			if (seqs.size() > second_condition) {
				cluster_sequences.emplace_back(seqs);
			}
		}

		int total_num = 0;
		for (int i = 0;i < items;i++) {
			if (id_root_map[i] == i) total_num++;
		}
		cerr << " >>> Total groups in this round BEFORE second group: " << total_num << endl;

		cerr << "Groups larger than " << second_condition << " : " << cluster_sequences.size() << endl;

		sort(cluster_sequences.begin(), cluster_sequences.end(), [](const vector<int>& a, const vector<int>& b) {
			return a.size() > b.size();
			});

		cerr << "Top 10 largest group size is: ";
		for (int i = 0; i < std::min(10, (int)cluster_sequences.size()); i++) {
			cerr << cluster_sequences[i].size() << " ";
		}
		cerr << endl;
		if (m + 1 < M) {
		#pragma omp parallel for num_threads(num_threads)
			for (int i = 0;i < cluster_sequences.size();i++) {
				SecondUpdate(cluster_sequences[i]);
			}

			uf.updateParent(id_root_map);

		#pragma omp parallel for num_threads(num_threads)
			for (int i = 0; i < cluster_sequences.size(); i++) {
				SecondGroup(cluster_sequences[i], m, hashes);
			}
			uf.findRoot(id_root_map);

			// 下面是输出
			total_num = 0;
			for (int i = 0;i < items;i++) {
				if (id_root_map[i] == i) total_num++;
			}
			cerr << " >>> Total groups in this round AFTER second group: " << total_num << endl;
			unordered_map<int, vector<int>> map_after_cluster;
			priority_queue<int, vector<int>, greater<int>> minHeap;

			for (int i = 0; i < items; i++) {
				map_after_cluster[id_root_map[i]].push_back(i);
			}
			int largethan1w = 0;
			for (auto& [root_id, seqs] : map_after_cluster) {
				if (seqs.size() > second_condition)
					largethan1w++;
				minHeap.push(seqs.size());
				if (minHeap.size() > 10) {
					minHeap.pop();
				}
			}
			cerr << "After second grouping, groups size large than " << second_condition << " : " << largethan1w << endl;
			while (!minHeap.empty()) {
				cerr << minHeap.top() << " ";
				minHeap.pop();
			}
			cerr << endl;
		}
		// else {
		// 	#pragma omp parallel for num_threads(num_threads)
		// 	for(int i = 0; i < cluster_sequences.size(); i++) {
		// 		// cerr << i << " is doing cluster " << cluster_sequences[i].size() << " sequences" << endl;
		// 		clusterEachGroup(cluster_sequences[i]);
		// 	}
		// 	uf.updateParent(id_root_map);
		// }
		uf.findRoot(id_root_map);
		map.clear();
		for (int i = 0; i < items; i++) {
			map[id_root_map[i]].push_back(i);
		}
	}

	if (cluster_on) {
		vector<vector<int>> cluster_sequences;
		for (auto& [key, seqs] : map) {
			if (seqs.size()>cluster_condition) {
				cluster_sequences.emplace_back(seqs);
			}
		}
		// vector<vector<int>> temp_cluster_sequences;
		// int count=0;
		
		// 用于记录大小超过指定阈值的序列集合
		// for (auto &[key, seqs] : map)
		// {

		// 	if (seqs.size() > cluster_condition)
		// 	{  
				
		// 		if (seqs.size() > 100000)
		// 		{	vector<std::vector<int>> block_seqs;
		// 			block_seqs.push_back(seqs);
		// 			tasks.emplace_back(block_seqs, 4);
		// 			cluster_sequences.emplace_back(seqs);
		// 		}
		// 		else if(seqs.size()<200){
		// 			cluster_sequences.emplace_back(seqs);
		// 			temp_cluster_sequences.emplace_back(seqs);
		// 			count ++;
		// 			if(count>=50){
		// 			tasks.emplace_back(temp_cluster_sequences, 1);
		// 			temp_cluster_sequences.clear();
		// 			count=0;
		// 			}
		// 		}

		
		// 		else
		// 		{	
		// 			vector<std::vector<int>> block_seqs;
		// 			block_seqs.push_back(seqs);
		// 			tasks.emplace_back(block_seqs, 1);
		// 			cluster_sequences.emplace_back(seqs);
		// 		}

		// 	}


		// 	// cerr<<"seqs size"<<cluster_sequences.back().size()<<endl;
		
		// }
		// int block_condition=count/TOTAL_THREADS;
		// int count_1=0;
		// vector<vector<int>> temp_temp_cluster_sequences;
		// for(int i=0;i<temp_cluster_sequences.size();i++){
		// 	temp_temp_cluster_sequences.emplace_back(temp_cluster_sequences[i]);
		// 	count_1++;
		// 	if(count_1>block_condition){
		// 		tasks.emplace_back(temp_temp_cluster_sequences, 1);
		// 		temp_temp_cluster_sequences.clear();
		// 		count_1=0;
		// 	}
		// }
		// if (!temp_cluster_sequences.empty())
		// {
		// 	tasks.emplace_back(temp_cluster_sequences, 1);
		// }
		// cerr<<"--------------------------"<<endl;
		// cerr<<"task size      "<<tasks.size()<<endl;
		// // cerr<<"block_condition      "<<block_condition<<endl;
		int total_num = 0;
		for (int i = 0;i < items;i++) {
			if (id_root_map[i] == i) total_num++;
		}
		cerr << " >>> Total groups in this round before cluster: " << total_num << endl;

		cerr << "Groups larger than " << cluster_condition << " : " << cluster_sequences.size() << endl;

		sort(cluster_sequences.begin(), cluster_sequences.end(), [](const vector<int>& a, const vector<int>& b) {
			return a.size() > b.size();
			});

		cerr << "Top 10 largest group size is: ";
		for (int i = 0; i < std::min(10, (int)cluster_sequences.size()); i++) {
			cerr << cluster_sequences[i].size() << " ";
		}
		cerr << endl;
		if(thread_pool){
			use_thread_pool(cluster_sequences);
		}

		else{
			#pragma omp parallel for num_threads(num_threads)
		for (int i = 0; i < cluster_sequences.size(); i++) {
			// cerr << i << " is doing cluster " << cluster_sequences[i].size() << " sequences" << endl;
			no_thread_clusterEachGroup(cluster_sequences[i]);
		}
		}
	
		uf.updateParent(id_root_map);

		unordered_map<int, vector<int>> map_after_cluster;
		priority_queue<int, vector<int>, greater<int>> minHeap;

		total_num = 0;
		for (int i = 0;i < items;i++) {
			if (id_root_map[i] == i) total_num++;
		}
		cerr << " >>> Total groups in this round after cluster: " << total_num << endl;

		for (int i = 0; i < items; i++) {
			map_after_cluster[id_root_map[i]].push_back(i);
		}
		int largethan1w = 0;
		for (auto& [root_id, seqs] : map_after_cluster) {
			if (seqs.size() > cluster_condition)
				largethan1w++;
			minHeap.push(seqs.size());
			if (minHeap.size() > 10) {
				minHeap.pop();
			}
		}
		cerr << "After clustering, clusters size large than " << cluster_condition << " : " << largethan1w << endl;
		while (!minHeap.empty()) {
			cerr << minHeap.top() << " ";
			minHeap.pop();
		}
		cerr << endl << "==================" << endl;
		cerr << endl;

	}
	else {
		vector<vector<int>> cluster_sequences;
		for (auto& [key, seqs] : map) {
			cluster_sequences.emplace_back(seqs);
		}

		sort(cluster_sequences.begin(), cluster_sequences.end(), [](const vector<int>& a, const vector<int>& b) {
			return a.size() > b.size();
			});

		cerr << ">>> The number of groups in this step: " << cluster_sequences.size() << endl;
		cerr << ">>> Top 10 largest group size is: " << endl;
		for (int i = 0; i < std::min(10, (int)cluster_sequences.size()); i++) {
			cerr << cluster_sequences[i].size() << " ";
		}
		cerr << endl << "==================" << endl << endl;
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
	while (!minHeap.empty()) {
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl << endl;
}

void GroupStream::Group(vector<vector<uint64_t>>& hashes, unordered_map<int, vector<int>>& group_map) {
	if (slide) {
		for (int m = 0; m < M - R + 1; m++) {
			uf.findRoot(temp_id_root_map);
			cerr << "ROUND " << m << endl;
			fillHashVec(hashes, hash_vec, m * L);
			GroupByCol(hash_vec, uf);
			if (m == M - R) {
				temp_output_on = true;
				cluster_condition = 1;
			}
			if (m <= 10) second_condition = 10000;
			else second_condition = 100000;
			countGroupSize(uf, m, hashes);
		}
	}
	else {
		for (int m = 0; m < M / R; m++) {
			cerr << "round " << m << endl;
			fillHashVec(hashes, hash_vec, m * R * L);
			GroupByCol(hash_vec, uf);
			countGroupSize(uf, m, hashes);
		}
		if (M % R != 0) {
			cerr << "round " << M / R;
			setR(M % R);
			fillHashVec(hashes, hash_vec, (M / R) * R * L);
			GroupByCol(hash_vec, uf);
			// 此处会存在问题
			countGroupSize(uf, 0, hashes);
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
}

void GroupStream::SecondUpdate(vector<int>& group_seqs) {
	for (int i = 0;i < group_seqs.size();i++) {
		int id = group_seqs[i];
		id_root_map[id] = temp_id_root_map[id];
	}
}

void GroupStream::SecondGroup(vector<int>& group_seqs, int m, const vector<vector<uint64_t>>& hashes) {
	vector<Data> temp_hash_vec(group_seqs.size());
	for (int i = 0;i < group_seqs.size();i++) {
		int id = group_seqs[i];
		temp_hash_vec[i].id = id;
		temp_hash_vec[i].value.resize(2);    // 因为你 copy 了 2 个元素 // 似乎这一句不需要
		copy(hashes[id].begin() + m, hashes[id].begin() + m + 2, temp_hash_vec[i].value.begin());
	}
	Sort(temp_hash_vec);
	Unite(temp_hash_vec, uf);
	return;
}

void GroupStream::clusterEachGroup(vector<int>& group_seqs,int need_thread) {
	vector<Sequence_new> sequences;
	// cerr<<"group_seqssize"<<group_seqs.size()<<endl;
	// cerr<<"need_thread"<<need_thread<<endl;
	for (int i = 0; i < group_seqs.size(); i++) {

		sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());

	}
	
	//读取FAI获取data
	cluster_cdhit.cdhit_cluster(sequences, id_root_map,need_thread);
	
}


void GroupStream::no_thread_clusterEachGroup(vector<int>& group_seqs) {
	vector<Sequence_new> sequences;
	// cerr<<"group_seqssize"<<group_seqs.size()<<endl;
	// cerr<<"need_thread"<<need_thread<<endl;
	for (int i = 0; i < group_seqs.size(); i++) {

		sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());

	}
	
	//读取FAI获取data
	cluster_cdhit.no_thread_cdhit_cluster(sequences, id_root_map);
	
}
