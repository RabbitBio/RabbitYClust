#include "GroupStream.h"
#include <queue>
#include <omp.h>
#include <thread>
#include <algorithm>

struct minheapcompare {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

bool compareByHash(const Data &a, const Data &b) {
	return a.value < b.value;
}

bool compareById(const Data &a,const Data & b){
	return a.id < b.id;
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

	for(int i = 0; i < items; i++) {
		int id = seq_ids[i];
		int root_id = id_root_map[id];
		group_map[root_id].push_back(id);
	}
}

void GroupStream::countGroupSize(UnionFind& uf,int m,vector<vector<uint64_t>>& hashes) {
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
		vector<vector<int>> cluster_sequences;
		// 用于记录大小超过指定阈值的序列集合
		for (auto& [key, seqs] : map) {
			if(seqs.size() > cluster_condition) {
				cluster_sequences.emplace_back(seqs);
			}
		}

		int total_num = 0;
		for (int i = 0;i < items;i++) {
			if (id_root_map[i] == i) total_num++;
		}
		cerr << " >>> Total groups in this round before cluster: " << total_num << endl;

		cerr << "Groups larger than " << cluster_condition << " : " << cluster_sequences.size() << endl;

		sort(cluster_sequences.begin(), cluster_sequences.end(), [](const vector<int>& a, const vector<int>& b){
		return a.size() > b.size();
			});

		cerr << "Top 10 largest group size is: ";
		for(int i = 0; i < std::min(10, (int)cluster_sequences.size()); i++){
			cerr << cluster_sequences[i].size() << " ";
		}
		cerr << endl;

		if (second_group) {
			if(m+1<=M){
				// 这里仍旧考虑两种方案，第一种是 (m+1)%M,另一种是m+1>M时直接return
				fillHashVec(hashes, hash_vec, (m + 1) * L);
				sort(hash_vec.begin(), hash_vec.end(), compareById);
			#pragma omp parallel for num_threads(num_threads)
				for (int i = 0; i < cluster_sequences.size(); i++) {
					SecondGroup(cluster_sequences[i], m, hashes);
				}
				uf.findRoot(id_root_map);
				total_num = 0;
				for (int i = 0;i < items;i++) {
					if (id_root_map[i] == i) total_num++;
				}
				cerr << " >>> Total groups in this round after second group: " << total_num << endl;
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
				cerr << "After second grouping, groups size large than " << cluster_condition << " : "<< largethan1w << endl;
				while(!minHeap.empty()){
					cerr << minHeap.top() << " ";
					minHeap.pop();
				}
				cerr << endl;
			}
			else {
				#pragma omp parallel for num_threads(num_threads)
				for(int i = 0; i < cluster_sequences.size(); i++) {
					// cerr << i << " is doing cluster " << cluster_sequences[i].size() << " sequences" << endl;
					clusterEachGroup(cluster_sequences[i]);
				}
				uf.updateParent(id_root_map);
			}
		}


		#pragma omp parallel for num_threads(num_threads)
		for(int i = 0; i < cluster_sequences.size(); i++) {
			// cerr << i << " is doing cluster " << cluster_sequences[i].size() << " sequences" << endl;
			clusterEachGroup(cluster_sequences[i]);
		}
		uf.updateParent(id_root_map);

		unordered_map<int, vector<int>> map_after_cluster;
		priority_queue<int, vector<int>, greater<int>> minHeap;

		total_num = 0;
		for (int i = 0;i < items;i++) {
			if (id_root_map[i] == i) total_num++;
		}
		cerr << " >>> Total groups in this round after cluster: " << total_num << endl;

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
			cerr << minHeap.top() << " ";
			minHeap.pop();
		}
		cerr << endl << "==================" << endl;
		cerr << endl;

	}
	else {
		vector<vector<int>> cluster_sequences;
		for(auto &[key, seqs] : map){
			cluster_sequences.emplace_back(seqs);
		}

		sort(cluster_sequences.begin(), cluster_sequences.end(), [](const vector<int>& a, const vector<int>& b){
		return a.size() > b.size();
			});

		cerr << "Information for this group step:" << endl;
		cerr << ">>> The number of groups in this step: " << cluster_sequences.size() << endl;
		cerr << ">>> Top 10 largest group size is: ";
		for(int i = 0; i < std::min(10, (int)cluster_sequences.size()); i++){
			cerr << cluster_sequences[i].size() << " ";
		}
		cerr << endl<< "==================" << endl;
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
	if(slide) {
		for(int m=0; m < M-R+1; m++){
			cerr << "round "<<  m << endl;
			fillHashVec(hashes, hash_vec, m * L);
			GroupByCol(hash_vec, uf);
			if(m == M-R) {
				temp_output_on = true;
			}
			countGroupSize(uf,m,hashes);
		}
	}else{
		for(int m=0; m < M / R; m++){
			cerr << "round "<<  m << endl;
			fillHashVec(hashes, hash_vec, m * R * L);
			GroupByCol(hash_vec, uf);
			countGroupSize(uf,m,hashes);
		}
		if(M % R != 0){
			cerr << "round "<<  M / R;
			setR(M % R);
			fillHashVec(hashes, hash_vec, (M/R) * R * L );
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


void GroupStream::SecondGroup(vector<int>& group_seqs, int m, vector<vector<uint64_t>>& hashes) {
	vector<Data> temp_hash_vec;
	for (int i = 0;i < group_seqs.size();i++) {
		id_root_map[group_seqs[i]] = group_seqs[i];
		temp_hash_vec.emplace_back(hash_vec[i]);
	}
	uf.updateParent(id_root_map);
	Sort(temp_hash_vec);
	Unite(temp_hash_vec, uf);
	return;
}

void GroupStream::clusterEachGroup(vector<int>& group_seqs) {
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]].c_str());
	}
	//读取FAI获取data
	cluster_cdhit.cdhit_cluster(sequences, id_root_map);
}
