#include "GroupStream.h"
#include <queue>

bool compareByHash(const Data &a, const Data &b) {
	return a.value < b.value;
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

void GroupStream::countGroupSize(UnionFind& uf) {
// FIXME:用结构体GroupNode存储id-root的映射还是用hash_vec继续存
// 用GroupNode增加内存但是如果排序的话要搬移的数据少
	uf.findRoot(id_root_map);
	priority_queue<int, vector<int>, greater<int>> minHeap;
// FIXME:用map来统计还是排序后统计
// 1.用map来统计分组结果 增加内存 只遍历一次
	unordered_map<int, vector<int>> map;
//	for(auto &p : id_root_map){
//    	map[p.root].push_back(p.id); 
//	}
	for(int i = 0; i < items; i++) {
		map[id_root_map[i]].push_back(i);
	}

	if(cluster_on) {
		for(auto &[key, seqs] : map){
			if(seqs.size() > 1) {
				clusterEachGroup(seqs);
			}
		}
		uf.updateParent(id_root_map);

		unordered_map<int, vector<int>> map_after_cluster;
		for(int i = 0; i < items; i++) {
			map_after_cluster[id_root_map[i]].push_back(i);
		}
		for(auto &[root_id, seqs] : map_after_cluster){
			minHeap.push(seqs.size());
			if (minHeap.size() > 10){
				 minHeap.pop();
			}
		}
	}else {

		for(auto &[root_id, seqs] : map){
			minHeap.push(seqs.size());
			if (minHeap.size() > 10){
				 minHeap.pop();
			}
		}
	}
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
			countGroupSize(uf);
		}
	}else{
		for(int m=0; m < M / R; m++){
			cerr << "round "<<  m << endl;
			fillHashVec(hashes, hash_vec, m * R * L);
			GroupByCol(hash_vec, uf);
			countGroupSize(uf);
		}
		if(M % R != 0){
			cerr << "round "<<  M / R;
			setR(M % R);
			fillHashVec(hashes, hash_vec, (M/R) * R * L );
			GroupByCol(hash_vec, uf);
			countGroupSize(uf);
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

void GroupStream::clusterEachGroup(vector<int>& group_seqs){
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		sequences.emplace_back(group_seqs[i], fa_map[group_seqs[i]]);
	}
	//读取FAI获取data
	cluster_cdhit.cdhit_cluster(sequences, id_root_map);
}
