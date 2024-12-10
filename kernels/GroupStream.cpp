#include "GroupStream.h"

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
	uint64_t cur_value = dataList[0].value;
	int cur_head = dataList[0].id;
	for (const auto& data : dataList) {
		if(data.value == cur_value){
			uf.unite(data.id, cur_head);
		}else{
			cur_value = data.value;
			cur_head = data.id;
		}
	}
}

void GroupStream::GroupByCol(vector<Data>& hash_vec, UnionFind& uf) {
#ifdef timing
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
	
    for (int i = 0; i < items; ++i) {
		hash_vec[i].id = i;
		hash_vec[i].value = vec[i][m];
	}
// 3. memory alignment and memcpy

}

void GroupStream::getGroupMap(UnionFind& uf, unordered_map<int, vector<int>>& group_map) {
#ifdef timing
	auto construct_start = chrono::high_resolution_clock::now();
	uf.findRoot(hash_vec);
	auto construct_end = chrono::high_resolution_clock::now();
	auto construct_duration = chrono::duration_cast<chrono::seconds>(construct_end - construct_start).count();
	cout << "construct vector<id, root> time needed is " << construct_duration << endl;

	auto sort_start = chrono::high_resolution_clock::now();
	Sort(hash_vec);
	auto sort_end = chrono::high_resolution_clock::now();
	auto sort_duration = chrono::duration_cast<chrono::seconds>(construct_end - construct_start).count();
	cout << "sort vector<id, root> time needed is " << sort_duration << endl;

	auto map_start = chrono::high_resolution_clock::now();
    for (const auto& p : hash_vec) {
    	group_map[p.value].push_back(p.id); 
    }

	cout << "size of hashmap is " << group_map.size() << endl;
	auto map_end = chrono::high_resolution_clock::now();
	auto map_duration = chrono::duration_cast<chrono::seconds>(map_end - map_start).count();
	cout << "time need for constructing GroupResMap is " << map_duration << endl;
#else
	uf.findRoot(hash_vec);
	Sort(hash_vec);
    for (const auto& p : hash_vec) {
    	group_map[p.value].push_back(p.id); 
    }
#endif
}

void GroupStream::Group(vector<vector<uint64_t>>& hashes, unordered_map<int, vector<int>>& group_map) {
	for(int m=0; m < M; m++){
		cerr << m << " sort and unite" << endl;
		fillHashVec(hashes, hash_vec, m);
		GroupByCol(hash_vec, uf);
	}
	
#ifdef timing
	auto transpose_start_time = chrono::high_resolution_clock::now();
	getGroupMap(uf, group_map);
	auto transpose_end_time = chrono::high_resolution_clock::now();
	auto transpose_duration = chrono::duration_cast<chrono::seconds>(transpose_end_time - transpose_start_time).count();
	cout << "transpose time is " << transpose_duration << endl;
#else
	getGroupMap(uf, group_map);
#endif
}
