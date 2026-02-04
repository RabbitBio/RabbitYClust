#include "GroupStream.h"
#include <queue>
#include <omp.h>
#include <thread>
#include <algorithm>
#include <atomic>
#include <unordered_set>
#include <fstream>
#include <functional>
#include <numeric>

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
	//resize(gs_config.items);
	initVec(gs_config.items);
	initOptions();
    tau = 0.5;
}
void GroupStream::reset(
	vector<pair<uint32_t, uint32_t>>& minhash_collisions,
	vector<sharedData>& seq_hash_vec
	) {
	minhash_collisions.clear();
	for(int i = 0; i < gs_config.items; i++) {
		seq_vec[i].seq_id = i;
		seq_vec[i].group_id = 0;
		seq_vec[i].value = 0;
	}
}


void GroupStream::initVec(int n) {
	seq_vec.reserve(n);
	for(int i = 0; i < n; i++) {
		seq_vec.emplace_back(i, 0, 0);
	}
}

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
    if(dataList[0].value.size() > 1){
        if (gs_config.num_threads == 1) ips4o_sort_single_thread(dataList, gs_config.R);
        else ips4o_sort_multi_thread(dataList, gs_config.R);
    }else{
        if (gs_config.num_threads == 1) ips2ra_sort_single_thread(dataList, dataList.size());
        else ips2ra_sort_multi_thread(dataList, dataList.size());
    }
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
	ProteinAAStore& store,
	vector<vector<pair<int, int>>>& minhash_collisions
	) {

	std::cout << "Traversing hash_vec to collect minhash info" << std::endl;
	auto start_collect = chrono::high_resolution_clock::now();
    int huge_groups_cnt = 0;
	for(int i = 0; i < gs_config.items; i++) {
		auto thisone = hash_vec[i];
		vector<pair<int,int>> collision;
		if(i == gs_config.items - 1 || hash_vec[i+1].value != thisone.value) continue;
		while(i < gs_config.items) {
			auto nextone = hash_vec[i];
			if(thisone.value != nextone.value) break;
			collision.emplace_back(nextone.id, thisone.id);
			i++;
		}
		if(collision.size() > 5000) huge_groups_cnt++;
		sort(collision.begin(), collision.end());
		minhash_collisions.emplace_back(collision);
		i -= 1;
	}
	auto end_collect = chrono::high_resolution_clock::now();
	auto duration_collect = chrono::duration_cast<chrono::seconds>(end_collect - start_collect).count();

	cerr << "MinHash collisions : " << minhash_collisions.size() << endl;
	sort(minhash_collisions.begin(), minhash_collisions.end(),
		[](const vector<pair<int,int>>& a, const vector<pair<int,int>>& b){
			return a.size() > b.size();
		});
	cerr << "    Top 10 largest collisions: " ;
	for(int i = 0; i < std::min(10, (int)minhash_collisions.size()); i++){
		cerr << minhash_collisions[i].size() << " ";
	}
	cerr << endl;

	cerr << endl;
	auto start_time = chrono::high_resolution_clock::now();
	cutEdges(minhash_collisions, store);
	auto end_time = chrono::high_resolution_clock::now();
	auto duration_cc = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
	cerr << "	Break the bad edges time(seconds): " << duration_cc  << endl;
}

void GroupStream::checkEdges(
	vector<Data>& hash_vec, 
	const vector<string>& fa_map,
	//const unordered_map<uint64_t, string>& fa_map,
	vector<vector<pair<int, int>>>& minhash_collisions
	) {
	std::cout << "Traversing hash_vec to collect minhash info" << std::endl;
	auto start_collect = chrono::high_resolution_clock::now();
    int huge_groups_cnt = 0;
	for(int i = 0; i < gs_config.items; i++) {
		auto thisone = hash_vec[i];
		vector<pair<int,int>> collision;
		if(i == gs_config.items - 1 || hash_vec[i+1].value != thisone.value) continue;
		while(i < gs_config.items) {
			auto nextone = hash_vec[i];
			if(thisone.value != nextone.value) break;
			collision.emplace_back(nextone.id, thisone.id);
			i++;
		}
		if(collision.size() > 5000) huge_groups_cnt++;
		sort(collision.begin(), collision.end());
		minhash_collisions.emplace_back(collision);
		i -= 1;
	}
	auto end_collect = chrono::high_resolution_clock::now();
	auto duration_collect = chrono::duration_cast<chrono::seconds>(end_collect - start_collect).count();

	cerr << "MinHash collisions : " << minhash_collisions.size() << endl;
	sort(minhash_collisions.begin(), minhash_collisions.end(),
		[](const vector<pair<int,int>>& a, const vector<pair<int,int>>& b){
			return a.size() > b.size();
		});
	cerr << "    Top 10 largest collisions: " ;
	for(int i = 0; i < std::min(10, (int)minhash_collisions.size()); i++){
		cerr << minhash_collisions[i].size() << " ";
	}
	cerr << endl;

	cerr << endl;
	auto start_time = chrono::high_resolution_clock::now();
	cutEdges(minhash_collisions, huge_groups_cnt, fa_map);

	auto end_time = chrono::high_resolution_clock::now();
	auto duration_cc = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
	cerr << "	Break the bad edges time(seconds): " << duration_cc  << endl;

	//unordered_map<int, vector<int>> groups_after_filter; // rootid:[seq0, seq1...]
	//priority_queue<int, vector<int>, greater<int>> minHeap;
	//for(int i = 0; i < gs_config.items; i++) {
	//	groups_after_filter[id_root_map[i]].push_back(i);
	//}

	//for(auto &[root_id, seqs] : groups_after_filter){
	//	minHeap.push(seqs.size());
	//	if (minHeap.size() > 10){
	//		minHeap.pop();
	//	}
	//}
	//cerr << "After break the bad edges, connected components number are: " << groups_after_filter.size() << endl;
	//cerr << "    Top 10 largest connected components: ";
	//while(!minHeap.empty()){
	//	cerr << minHeap.top() << " ";
	//	minHeap.pop();
	//}
	//cerr << endl;
}


void GroupStream::Unite(const vector<Data>& dataList, UnionFind& this_uf) {
	std::cout << "Uniting Unionfind." << std::endl;
	auto start = chrono::high_resolution_clock::now();
	vector<uint64_t> cur_value = dataList[0].value;
	for (int i = 1; i < dataList.size(); i++) {
		auto thisone = dataList[i];
    	auto lastone = dataList[i-1];
		if(thisone.value == lastone.value) {
		    this_uf.unite(thisone.id, lastone.id);
		}
    }
	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
}

// 不用额外的id_root_map版 只用minhashcollisions 记录所有序列的id
void GroupStream::uniteByEdges(vector<vector<pair<int, int>>>& minhash_collisions) {
	priority_queue<int, vector<int>, greater<int>> minHeap;
	int total_cc_cnt = 0;
	for(auto& vec : minhash_collisions) {
		int n = vec.size();
		int cc_size = 1;
		uf.unite(vec[0].first, vec[0].second);
		for(int i = 1; i < n; i++) {
			auto& thisone = vec[i];
			auto& lastone = vec[i-1];
			uf.unite(thisone.first, thisone.second);
			if(thisone.second != lastone.second){
				total_cc_cnt++;
				minHeap.push(cc_size);
				if(minHeap.size() > 10) minHeap.pop();
				cc_size = 1;
				continue;
			}
			cc_size++;
		}
		minHeap.push(cc_size);
		if(minHeap.size() > 10) minHeap.pop();
		total_cc_cnt++;
	}
	cerr << "After break the bad edges, connected components number are: " << total_cc_cnt << endl;
	cerr << "    Top 10 largest connected components: ";
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;

}

void GroupStream::uniteByEdges(vector<int>& id_root_map) {
	for(int i = 0; i < gs_config.items; i++) {
		uf.unite(uf.find(i), uf.find(id_root_map[i]));
	}
}

void GroupStream::uniteByEdges(UnionFind& cur_uf) {
	for(int i = 0; i < gs_config.items; i++) {
		uf.unite(uf.find(i), cur_uf.find(i));
	}
}

void GroupStream::countGroupSizeBySort(
	vector<pair<uint32_t, uint32_t>>& need_to_clutser,
	int cluster_condition
	){
	int groups_size = uf.countSetsSize();
	cout << "Group Size after merging:" << groups_size << endl;
	cerr << "Group Size after merging:" << groups_size << endl;

	// output top10 group
	// collect group-id
	for(int i = 0; i < gs_config.items; i++) {
		seq_vec[i].seq_id = i;
		seq_vec[i].group_id = uf.find(i);
	}

	ips2ra::parallel::sort(seq_vec.begin(), seq_vec.end(), [](const sharedData& r){return r.group_id;}, gs_config.num_threads);

	priority_queue<int, vector<int>, greater<int>> minHeap;
	int last_group_id = seq_vec[0].group_id;
	int last_idx = 0;
    int total_groups = 0;
    int only_one = 0;
	for(int i = 0; i < gs_config.items; i++) {
		if(seq_vec[i].group_id != last_group_id){
            total_groups++;
            if(last_idx == i - 1) only_one++;
			if(i - last_idx > cluster_condition) need_to_clutser.emplace_back(last_idx, i - last_idx);
			minHeap.push(i - last_idx);
			if (minHeap.size() > 10){
				minHeap.pop();
			}
			last_group_id = seq_vec[i].group_id;
			last_idx = i;
			
		}
	}
	if(gs_config.items - last_idx > cluster_condition) need_to_clutser.emplace_back(last_idx, gs_config.items - last_idx);
    total_groups++;
    if(last_idx == gs_config.items - 1) only_one++;
    cerr << "Number of total_groups: " << total_groups << endl;
    cerr << "Number of groups only having 1 sequence: " << only_one << endl;

	//对need_to_cluster按照group大小排序
	ips2ra::parallel::sort(need_to_clutser.begin(), need_to_clutser.end(), [](const pair<uint32_t, uint32_t>& r){return r.second;}, gs_config.num_threads);

	minHeap.push(gs_config.items - last_idx);
	if (minHeap.size() > 10){
		minHeap.pop();
	}
	cerr << "Top 10 largest groups size in this round is: ";
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;
}

//可以直接改成unite seq_id and group_id
void GroupStream::uniteByEdges(
	vector<pair<uint32_t, uint32_t>>& minhash_collisions,
	int start_update_idx
	){
	//priority_queue<int, vector<int>, greater<int>> minHeap;


	auto time_start = chrono::high_resolution_clock::now();
	if(gs_config.cluster_on){
		int cc_cnt = 0;
		for(int ptr = start_update_idx; ptr < minhash_collisions.size(); ptr++) {
			int size = minhash_collisions[ptr].second;
			int start_idx = minhash_collisions[ptr].first;
			int end_idx = start_idx + size;

			int last_group_id = seq_vec[start_idx].group_id;
			int last_idx = start_idx;
			cc_cnt++;

			for(int i = start_idx; i < end_idx; i++) {
				int seq_id = seq_vec[i].seq_id;
				int this_group_id = seq_vec[i].group_id;
				// TODO 统计cc个数还有top10都可以挪到libcc里面
				//if(this_group_id != last_group_id){ // 为了统计Top10 connected_components
				//	minHeap.push(i - last_idx);
				//	if (minHeap.size() > 10){
				//		minHeap.pop();
				//	}
				//	last_idx = i;
				//	last_group_id = this_group_id;
				//	cc_cnt++;
				//}
				uf.unite(seq_id, this_group_id);
			}
			//minHeap.push(end_idx - last_idx);
			//if (minHeap.size() > 10){
			//	minHeap.pop();
			//}
		}
		//cerr << "Number of connected components created in libcc: " << cc_cnt << endl;
		//cerr << "Sum of connected components (including one node): " << cc_cnt + start_update_idx << endl;
		//cerr << "Top 10 largest connected components size in this round is: ";
		//while(!minHeap.empty()){
		//cerr << minHeap.top() << " ";
		//minHeap.pop();
		//}
		//cerr << endl;
	}else{
		for(int ptr = start_update_idx; ptr < minhash_collisions.size(); ptr++) {
			int size = minhash_collisions[ptr].second;
			int start_idx = minhash_collisions[ptr].first;
			int end_idx = start_idx + size;
			int head_seq_id = seq_vec[start_idx].seq_id;
			for(int i = start_idx; i < end_idx; i++) {
				int seq_id = seq_vec[i].seq_id;
				uf.unite(seq_id, head_seq_id);
			}
		}
	}

}

void GroupStream::GroupByCol(
	vector<pair<uint32_t, uint32_t>>& minhash_collisions,
	ProteinAAStore& store
	) {
	cout << "into GroupByCol " << endl;
	int collision_cnt = minhash_collisions.size();
	cerr << "Number of different MinHash collisions : " << minhash_collisions.size() << endl;
	cerr << "    Top 10 largest collisions: " ;
	for(int i = 0; i < 10; i++) {
		cerr << minhash_collisions[collision_cnt - 1 - i].second << " ";
	}
	cerr << endl;

	// minhash_collisions升序排序 找到第一个size(p.first)>1的位置
	int start_libcc_idx = minhash_collisions.size();
	int start_greedy_idx = minhash_collisions.size();
	auto it_large1 = std::upper_bound(
			minhash_collisions.begin(), minhash_collisions.end(), 1u,
			[](uint32_t key, const auto& p) { return key < p.second; }
			);
	auto it_large_greedy_condition = std::upper_bound(
			minhash_collisions.begin(), minhash_collisions.end(), greedy_condition,
			[](uint32_t key, const auto& p) { return key < p.second; }
			);
	start_libcc_idx = it_large1 - minhash_collisions.begin();
	start_greedy_idx = it_large_greedy_condition - minhash_collisions.begin();
	cerr << "Number of MinHash collisions containing more than 1 sequence: " << collision_cnt - start_libcc_idx << endl;
	cerr << "Number of MinHash collisions processed in libcc: " << start_greedy_idx - start_libcc_idx << endl;
	cerr << "Number of MinHash collisions processed in greedy (groups size large than ";
	cerr << greedy_condition << " ): " << collision_cnt - start_greedy_idx << endl;

	if(gs_config.M > 1 && gs_config.cluster_on && (start_greedy_idx < collision_cnt || start_libcc_idx < collision_cnt)) {
		cerr << "Already in breaking bad edges" << endl;
		auto start_time = chrono::high_resolution_clock::now();
		setOptionsSkipAlign(true);
		if(start_greedy_idx < collision_cnt)
			ClusterLargeThanRescueCondition(minhash_collisions, store, false, start_greedy_idx, collision_cnt);
		if(start_libcc_idx < start_greedy_idx)
			cutEdges(minhash_collisions, store, start_libcc_idx, start_greedy_idx);
		auto end_time = chrono::high_resolution_clock::now();
		auto duration_cc = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
		cerr << "Break the bad edges time(seconds): " << duration_cc  << endl;
	}

	// temp check! make sure the size of minhash_collisions is 1; must be 1!!! in test greedy no-sw libcdhit
	cerr << "check! size of minhash_collisions(must be 1): " << minhash_collisions.size() << endl; 
	auto start_unite_time = chrono::high_resolution_clock::now();
	uniteByEdges(minhash_collisions, start_libcc_idx);
	auto end_unite_time = chrono::high_resolution_clock::now();
	auto duration_unite = chrono::duration_cast<chrono::seconds>(end_unite_time - start_unite_time).count();
	cerr << "Time of update unionfind: " << duration_unite << endl;
}

void GroupStream::GroupByCol(
	vector<Data>& hash_vec, 
	ProteinAAStore& store
	) {
	Sort(hash_vec);
	if(gs_config.cluster_on) {
		vector<vector<pair<int, int>>> minhash_collisions;
		checkEdges(hash_vec, store, minhash_collisions);
		uniteByEdges(minhash_collisions);
	}else {
		Unite(hash_vec, uf);
	}
	int groups_size = uf.countSetsSize();
	cout << "Group Size after merging:" << groups_size << endl;
}
								
void GroupStream::GroupByCol(
	vector<Data>& hash_vec, 
	const vector<string>& fa_map
	//const unordered_map<uint64_t, string>& fa_map
	) {
	Sort(hash_vec);
	if(gs_config.cluster_on) {
		//UnionFind col_uf(gs_config.items);
		//Unite(hash_vec, col_uf);
		//checkEdges(hash_vec, col_uf, fa_map);
		//uniteByEdges(col_uf);
		vector<vector<pair<int, int>>> minhash_collisions;
		checkEdges(hash_vec, fa_map, minhash_collisions);
		//uniteByEdges(id_root_map);
		uniteByEdges(minhash_collisions);
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
    for(int r = 0; r < gs_config.R; r++) {
        auto& tmp_hashes = sketchdata.hashes[m+r];
        for (int i = 0; i < gs_config.items; i++) {
            hash_vec[i].id = i;
            hash_vec[i].value[r] = tmp_hashes[i];
            //copy(vec[i].begin() + m * gs_config.L, vec[i].begin() + m * gs_config.L + gs_config.R * gs_config.L, hash_vec[i].value.begin());
        }
    }
	cerr << gs_config.items << " valid items in round " << m << endl;
}


void GroupStream::fillHashVec(string sketch_filename, vector<Data>& hash_vec, int m) {
	int valid_items = 0;
	std::ifstream ifs(sketch_filename, ios::binary);
	std::cout << "start reading hash function from f" << m << " to  f" << (m+gs_config.R) << std::endl;
	if(!ifs){
		cerr << "Error opening file!" << endl;
		return;
	}
    for(int r = 0; r < gs_config.R; r++){
        valid_items = 0;
        vector<uint64_t> read_hashes(gs_config.items);
		ifs.seekg(sizeof(ProteinSketchData::Config)  + static_cast<long long>(m+r) *static_cast<long long>(gs_config.items) * sizeof(uint64_t), std::ios::beg);
		ifs.read(reinterpret_cast<char*>(read_hashes.data()), static_cast<long long>(gs_config.items) * sizeof(uint64_t));
        for (int i = 0; i < gs_config.items; i++) {
            hash_vec[i].id = i;
            hash_vec[i].value[r] = read_hashes[i];
            valid_items++;
        }
	}
	cerr << valid_items << " valid items in round " << m << endl;
}

void GroupStream::fillHashVecAndDetectMinHash(
	string sketch_filename, 
	vector<sharedData>& seq_hash_vec, 
	int m,
	vector<pair<uint32_t, uint32_t>>& minhash_collisions // first collisions_size; second offset in seq_hash_vec
	) {
	std::ifstream ifs(sketch_filename, ios::binary);
	std::cout << "start reading hash function from f" << m << " to  f" << (m+gs_config.R) << std::endl;
	if(!ifs){
		cerr << "Error opening file!" << endl;
		return;
	}
	// 后面几轮 把seq_vec minhash_collisions清空
	if(m > 0) {
		reset(minhash_collisions, seq_vec);
	}

    for(int r = 0; r < gs_config.R; r++){
        vector<uint64_t> read_hashes(gs_config.items);
		ifs.seekg(sizeof(ProteinSketchData::Config)  + static_cast<long long>(m+r) *static_cast<long long>(gs_config.items) * sizeof(uint64_t), std::ios::beg);
		ifs.read(reinterpret_cast<char*>(read_hashes.data()), static_cast<long long>(gs_config.items) * sizeof(uint64_t));
		cout << "successfully read!" <<endl;

		// R > 1 时做二次排序的前提
		if(r > 0) ips2ra::parallel::sort(seq_hash_vec.begin(), seq_hash_vec.end(), [](const sharedData& r) { return r.seq_id; }, gs_config.num_threads);

        for (int i = 0; i < gs_config.items; i++) {
            seq_hash_vec[i].value = read_hashes[i];
        }
		cout << "successfully fill hash value!" <<endl;

		if(r > 0) {
			//std::sort(seq_hash_vec.begin(), seq_hash_vec.end(), [](const sharedData& a, const sharedData& b) {
			//		return a.group_id != b.group_id ? (a.group_id < b.group_id) : (a.value < b.value);
			//});
			ips4o::parallel::sort(seq_hash_vec.begin(), seq_hash_vec.end(), 
					[](const sharedData& a, const sharedData& b) {
					return a.group_id != b.group_id ? (a.group_id < b.group_id) : (a.value == b.value ? a.seq_id < b.seq_id : a.value < b.value); },
					gs_config.num_threads);
			cout << "successfully sort seq_hash_vec to detect the same MinHash!" <<endl;
		}else {
			// 稳定排序版单线程
			//stable_sort(seq_hash_vec.begin(), seq_hash_vec.end(), [](const sharedData& a, const sharedData& b) { return a.value < b.value; });
			// 快速 但是结果不稳定版
			//ips2ra::parallel::sort(seq_hash_vec.begin(), seq_hash_vec.end(), [](const sharedData& r) { return r.value; }, gs_config.num_threads);
			ips4o::parallel::sort(seq_hash_vec.begin(), seq_hash_vec.end(), 
					[](const sharedData& a, const sharedData& b) { return a.value == b.value ? a.seq_id < b.seq_id : a.value < b.value; }, 
					gs_config.num_threads);
			cout << "successfully sort seq_hash_vec to detect the same MinHash!" <<endl;
		}

		// update group_id
		if(r == gs_config.R - 1) {
			int update_group_id = seq_hash_vec[0].seq_id;
			int this_gid = seq_hash_vec[0].group_id;
			uint64_t this_value = seq_hash_vec[0].value;

			int last_i = 0;
			for(int i = 0; i < gs_config.items; i++) {
				if(seq_hash_vec[i].group_id != this_gid 
						|| seq_hash_vec[i].value != this_value) {
					this_gid = seq_hash_vec[i].group_id;
					this_value = seq_hash_vec[i].value;
					minhash_collisions.emplace_back(last_i, i - last_i); // first是size second是偏移
					last_i = i;
					update_group_id = seq_hash_vec[i].seq_id;
				}
				seq_hash_vec[i].group_id = update_group_id;
			}
			minhash_collisions.emplace_back(last_i, gs_config.items - last_i); // 最后一组的信息
			cout << "successfully collect MinHash collisions info!" <<endl;
			ips2ra::parallel::sort(minhash_collisions.begin(), minhash_collisions.end(), [](const pair<uint32_t, uint32_t>& c) {return c.second; }, gs_config.num_threads);
			cout << "successfully sort MinHash collisions by size!" <<endl;
		}else {
			int update_group_id = seq_hash_vec[0].seq_id;
			int this_gid = seq_hash_vec[0].group_id;
			uint64_t this_value = seq_hash_vec[0].value;
			for(int i = 0; i < gs_config.items; i++) {
				if(seq_hash_vec[i].group_id != this_gid 
						|| seq_hash_vec[i].value != this_value) {
					this_gid = seq_hash_vec[i].group_id;
					this_value = seq_hash_vec[i].value;
					update_group_id = seq_hash_vec[i].seq_id;
				}
				seq_hash_vec[i].group_id = update_group_id;
			}
		}
	}
}

void GroupStream::getGroupRes(UnionFind& uf, unordered_map<int, vector<int>>& group_map, vector<int>& id_root_map) {
	uf.findRoot(id_root_map);
	for(int i = 0; i < gs_config.items; i++) {
		int id = i;
		int root_id = id_root_map[id];
		group_map[root_id].push_back(id);
	}
}

constexpr int aa_value(char c) {
    switch (c) {
        case 'A': return 0;  case 'C': return 1;  case 'D': return 2;  case 'E': return 3;
        case 'F': return 4;  case 'G': return 5;  case 'H': return 6;  case 'I': return 7;
        case 'K': return 8;  case 'L': return 9;  case 'M': return 10; case 'N': return 11;
        case 'P': return 12; case 'Q': return 13; case 'R': return 14; case 'S': return 15;
        case 'T': return 16; case 'V': return 17; case 'W': return 18; case 'Y': return 19;
        case 'X': return 20;
        default:  return 0;
    }
}

constexpr std::array<int,256> aa_map = []{
    std::array<int,256> a{};                 // 全部置 0
    for (int i = 0; i < 256; ++i)
        a[static_cast<unsigned char>(i)] = aa_value(static_cast<char>(i));
    return a;
}();

bool computeGlobalUniqueKmers(
    const vector<int>& seqs, 
    int k,
	const vector<string>& fa_map
	//const unordered_map<uint64_t, string>& fa_map
    ) {
    //unordered_map<int, int> unique_map;
    unordered_set<int> unique;
    int total_kmers = 0;
    int NAA = k;
    for (const auto& seq_id : seqs) {
        const char* seq = fa_map.at(seq_id).c_str();
        int len = fa_map.at(seq_id).size();
        if (len < k) continue;
        total_kmers += len - k + 1;
        
        // EncodeWords 
        int aan_no = len - NAA + 1;
        unsigned char k, k1;
        for (int j = 0; j < aan_no; j++) {
            const char* word = seq + j;
            int encode = 0;
            for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) {
                encode += aa_map[(unsigned char)word[k]] * NAAN_array[k1];  // 修改为使用 aa_map
            }
            unique.insert(encode);
            //unique_map[encode]++;
        }
    }
    return unique.size() * 32 <= total_kmers;
    // TODO 使用组合数

}
// 测试一下这个kernel
//bool computeGlobalUniqueKmers(
//    const vector<int>& seqs, 
//    int k,
//	const unordered_map<uint64_t, string>& fa_map
//    ) {
//    unordered_set<uint64_t> unique;
//    uint64_t base = 1;
//    for (int i = 1; i < k; ++i) base *= 131;  // 用131更稳
//
//    uint64_t total_kmers = 0;
//    for (const auto& seq_id : seqs) {
//        const char* seq = fa_map.at(seq_id).c_str();
//        int len = fa_map.at(seq_id).size();
//        total_kmers += len - k + 1;
//        
//        if (len < k) continue;
//        uint64_t h = 0;
//        for (int i = 0; i < k; ++i)
//            h = h * 131 + (seq[i] - 'A' + 1);
//
//        unique.insert(h);
//        for (size_t i = k; i < len; ++i) {
//            h = h * 131 
//                - (seq[i-k] - 'A' + 1) * base 
//                + (seq[i]   - 'A' + 1);
//            unique.insert(h);
//        }
//    }
//    // 以4%作为是否使用wordtable的标准
//    return unique.size() * 25 <= total_kmers;
//}
// 进度条函数
void print_progress(int current, int total, int bar_width = 50) {
	if (total == 0) return;
	
	double progress = (double)current / total;
	int pos = (int)(bar_width * progress);
	
	std::cout << "\r[";
	for (int i = 0; i < bar_width; ++i) {
		if (i < pos) std::cout << "=";
		else if (i == pos) std::cout << ">";
		else std::cout << " ";
	}
	std::cout << "] " << int(progress * 100.0) << "% (" << current << "/" << total << ")";
	std::cout.flush();
	
	if (current == total) std::cout << std::endl;
}

void GroupStream::cutEdges(
	vector<pair<uint32_t, uint32_t>>& minhash_collisions,
	ProteinAAStore& store,
	int start_idx, // 从哪个下标开始做断边
	int end_idx // start_idx < end_idx, minhash_collisions升序排序
	){
	// minhash_collisions升序排序 找到第一个size(p.first)>1的位置
	auto it_multi = std::upper_bound(
			minhash_collisions.begin() + start_idx, minhash_collisions.begin() + end_idx, 5000u,
			[](uint32_t key, const auto& p) { return key < p.second; }
			);
	size_t start_multi_idx = it_multi == minhash_collisions.begin() + end_idx ? end_idx : it_multi - minhash_collisions.begin();
	//start_multi_idx = it_multi == minhash_collisions.end() ? collision_cnt : it_multi - minhash_collisions.begin();
	cerr << "Number of MinHash collision processed in multi-threading: " << end_idx - start_multi_idx << endl;
	cerr << "Number of MinHash collision processed in single-threading: " << start_multi_idx - start_idx << endl;
    uint64_t validated_edges = 0;
    uint64_t cross_edge_sum = 0;
    uint64_t high_cj_edge_sum = 0;
	uint64_t filter_ed_edge_sum = 0;
	uint64_t pass_ed_edge_sum = 0;
	uint64_t last_round_jump_cnt = 0;
	uint64_t this_round_jump_cnt = 0;
	uint64_t need_edlib_edge = 0;

	int avail_threads = gs_config.num_threads;
 	omp_set_num_threads(avail_threads);

	if(start_multi_idx < end_idx) {
		auto start_huge_time = chrono::high_resolution_clock::now();
		for(int i = end_idx - 1 ; i >= start_multi_idx; i--) {
		//for(int i = start_multi_idx; i < minhash_collisions.size(); i++) {
			// 显示进度条（每10组更新一次，或最后一个）
			cerr << "doing task " << i-start_multi_idx << "of " << end_idx - start_multi_idx << endl;
			//if (i % 10 == 0 || i == end_idx - 1) {
			//	print_progress(i + 1, end_idx - start_idx);
			//}
			vector<uint64_t>  edge_stat = buildConnectedComponents(
					avail_threads, //threads
					minhash_collisions[i].first, // start_pos
					minhash_collisions[i].second + minhash_collisions[i].first, // end_pos
					store);
			validated_edges += edge_stat[0];
			cross_edge_sum += edge_stat[1];
			high_cj_edge_sum += edge_stat[2];
			filter_ed_edge_sum += edge_stat[3];
			pass_ed_edge_sum += edge_stat[4];
			last_round_jump_cnt += edge_stat[5];
			this_round_jump_cnt += edge_stat[6];
			need_edlib_edge += edge_stat[7];
		}
		auto end_huge_time = chrono::high_resolution_clock::now();
		auto duration_huge = chrono::duration_cast<chrono::seconds>(end_huge_time - start_huge_time).count();
		cerr << "Time of multi-thread libcdhit (use all threads once): " << duration_huge << endl;
	}

	int total_small_tasks = start_multi_idx - start_idx;
	int one_step = total_small_tasks / 100;
	auto start_small_time = chrono::high_resolution_clock::now();
    #pragma omp parallel num_threads(avail_threads) 
	{
        int tid = omp_get_thread_num();
		#pragma omp for schedule(dynamic,1) reduction(+:cross_edge_sum,validated_edges,high_cj_edge_sum,pass_ed_edge_sum,filter_ed_edge_sum,last_round_jump_cnt,this_round_jump_cnt,need_edlib_edge)
    	for(int i = start_multi_idx-1; i >= start_idx; i--) {
    	//for(int i = start_check_idx; i < start_multi_idx; i++) {
			vector<uint64_t>  edge_stat = buildConnectedComponents(
					1, //threads
					minhash_collisions[i].first, // start_pos
					minhash_collisions[i].second + minhash_collisions[i].first, // end_pos
					store);
			validated_edges += edge_stat[0];
			cross_edge_sum += edge_stat[1];
			high_cj_edge_sum += edge_stat[2];
			filter_ed_edge_sum += edge_stat[3];
			pass_ed_edge_sum += edge_stat[4];
			last_round_jump_cnt += edge_stat[5];
			this_round_jump_cnt += edge_stat[6];
			need_edlib_edge += edge_stat[7];

			if (i % one_step == 0 || i == start_multi_idx - 1) {
				#pragma omp critical
				print_progress(i - start_idx + 1, total_small_tasks);
			}
    	}
	}
	auto end_small_time = chrono::high_resolution_clock::now();
	auto duration_small = chrono::duration_cast<chrono::seconds>(end_small_time - start_small_time).count();
	cerr << "Time of single-thread libcdhit (use only 1 threads each group): " << duration_small << endl;

    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cerr << "Number of validated edges build among CC: " << validated_edges << endl;
    cerr << "Number of edges cross origin group among CC: " << cross_edge_sum << endl;
    cerr << "Number of edges get CJ >=0.6 among CC: " << high_cj_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6 and pass EDLIB: " << pass_ed_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6 and no pass EDLIB: " << filter_ed_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6: " << pass_ed_edge_sum + filter_ed_edge_sum << endl;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

}
void GroupStream::cutEdges(
	vector<vector<pair<int,int>>>& sequences_collisions, 
	ProteinAAStore& store
	){
    uint64_t validated_edges = 0;
    uint64_t cross_edge_sum = 0;
    uint64_t high_cj_edge_sum = 0;
	uint64_t filter_ed_edge_sum = 0;
	uint64_t pass_ed_edge_sum = 0;
    uint64_t minhash_edge_sum = 0;

	int avail_threads = gs_config.num_threads;
 	omp_set_num_threads(avail_threads);
	auto start_small_time = chrono::high_resolution_clock::now();
    #pragma omp parallel num_threads(avail_threads) 
	{
        int tid = omp_get_thread_num();
		#pragma omp for schedule(runtime) reduction(+:cross_edge_sum,minhash_edge_sum,validated_edges,high_cj_edge_sum,pass_ed_edge_sum,filter_ed_edge_sum)
    	for(int i = 0; i < sequences_collisions.size(); i++) {
			vector<uint64_t>  edge_stat = buildConnectedComponents_st(sequences_collisions[i], store, 1);
			validated_edges += edge_stat[0];
			cross_edge_sum += edge_stat[1];
			high_cj_edge_sum += edge_stat[2];
			filter_ed_edge_sum += edge_stat[3];
			pass_ed_edge_sum += edge_stat[4];
            minhash_edge_sum += (sequences_collisions[i].size() * (sequences_collisions[i].size()-1)) >> 1;
    	}
	}
	auto end_small_time = chrono::high_resolution_clock::now();
	auto duration_small = chrono::duration_cast<chrono::seconds>(end_small_time - start_small_time).count();
	cerr << "Time of single-thread libcdhit (use only 1 threads each group): " << duration_small << endl;

    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cerr << "Number of validated edges build among CC: " << validated_edges << endl;
    cerr << "Number of edges cross origin group among CC: " << cross_edge_sum << endl;
    cerr << "Number of edges get CJ >=0.6 among CC: " << high_cj_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6 and pass EDLIB: " << pass_ed_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6 and no pass EDLIB: " << filter_ed_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6: " << pass_ed_edge_sum + filter_ed_edge_sum << endl;
    cerr << "Number of total edges in MinHash collisions: " << minhash_edge_sum << endl;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

}

void GroupStream::cutEdges(
	vector<vector<pair<int,int>>>& sequences_collisions, 
	//vector<vector<int>>& sequences_collisions, 
	int huge_groups_cnt, // 需要多线程libcdhit的组的个数
	const vector<string>& fa_map
	//const unordered_map<uint64_t, string>& fa_map
	) {

	int avail_threads = gs_config.num_threads;
 	omp_set_num_threads(avail_threads);
	int mt_seqs = 0;
	cout << "avail threads: " << omp_get_max_threads() << endl;
    cerr << "Collisions number processed in multi-thread: " << huge_groups_cnt << endl;
	cerr << "Collisions number processed in single-thread: " << (sequences_collisions.size() - huge_groups_cnt) << endl;
    // edge stat
    uint64_t validated_edges = 0;
    uint64_t cross_edge_sum = 0;
    uint64_t high_cj_edge_sum = 0;
	uint64_t filter_ed_edge_sum = 0;
	uint64_t pass_ed_edge_sum = 0;
    uint64_t minhash_edge_sum = 0;

	cerr << "Huge task in multi-threading..." <<endl;
	// huge collisions in multiple thread libcdhit
	auto start_huge_time = chrono::high_resolution_clock::now();
	ClusterWS ws_for_hugegroup;
	for(int i = 0; i < huge_groups_cnt; i++) {
		// 显示进度条（每10组更新一次，或最后一个）
		if (i % 10 == 0 || i == huge_groups_cnt - 1) {
			print_progress(i + 1, huge_groups_cnt);
		}
        vector<uint64_t> edge_stat = buildConnectedComponents(sequences_collisions[i], avail_threads, fa_map, ws_for_hugegroup);
        validated_edges += edge_stat[0];
        cross_edge_sum += edge_stat[1];
        high_cj_edge_sum += edge_stat[2];
		filter_ed_edge_sum += edge_stat[3];
		pass_ed_edge_sum += edge_stat[4];
        minhash_edge_sum += (sequences_collisions[i].size() * (sequences_collisions[i].size()-1)) >> 1;
		mt_seqs += sequences_collisions[i].size();
	}
	auto end_huge_time = chrono::high_resolution_clock::now();
	auto duration_huge = chrono::duration_cast<chrono::seconds>(end_huge_time - start_huge_time).count();
	cerr << "Time of multi-thread libcdhit (use all threads once): " << duration_huge << endl;

    vector<int> tasks(sequences_collisions.size()); 
	auto start_build = chrono::high_resolution_clock::now();
    // task_id in sequences_collisions, type of libcdhit
    // eg. 0:less10; >1:word_table 
    // ---------------
    // 统计小于100的任务个数
    //for(int i = huge_groups_cnt; i < sequences_collisions.size(); i++)
    //{
    //    if(sequences_collisions[i].size() <= 100){
    //        cerr << "Number of groups(size <= 100): " <<  sequences_collisions.size()-i << endl;
    //        break;
    //    }
    //}
    // ---------------
    // 统计wordtable 根据word table的分布来决定使用wt or direct
    //int use_direct = 0;
	//#pragma omp parallel for schedule(dynamic) reduction(+:use_direct)
    //for(int i = huge_groups_cnt; i < sequences_collisions.size(); i++)
    //{
    //    // 统计不同的kmers
    //    if(sequences_collisions[i].size() <= 100 || computeGlobalUniqueKmers(sequences_collisions[i], 5, fa_map)) {
    //        tasks[i] = 0;
    //        use_direct++;
    //    }else {
    //        tasks[i] = 1;
    //    }
    //}
	//auto end_build = chrono::high_resolution_clock::now();
	//auto duration_build = chrono::duration_cast<chrono::seconds>(end_build - start_build).count();
    //cerr << "Number of groups use direct: " << use_direct << endl;
	//cerr << "Time of computing kmers density: " << duration_build << endl;
    // ---------------
    // 统计最耗时的任务
    //vector<int> per_thread_max_times(avail_threads, 0);
    //vector<int> use_type(avail_threads, 0);
    //vector<int> group_ids(avail_threads, 0);
    //vector<uint64_t> containing_aas(avail_threads, 0);
	cerr << "Small task in single-threading..." <<endl;
	auto start_small_time = chrono::high_resolution_clock::now();
	int small_groups_cnt = sequences_collisions.size() - huge_groups_cnt;
	int sequences_collisions_cnt = sequences_collisions.size();
    #pragma omp parallel num_threads(avail_threads) 
	{
        int tid = omp_get_thread_num();
		#pragma omp for schedule(runtime) reduction(+:cross_edge_sum,minhash_edge_sum,validated_edges,high_cj_edge_sum,pass_ed_edge_sum,filter_ed_edge_sum)
    	for(int i = huge_groups_cnt; i < sequences_collisions.size(); i++) {
	        //auto start1 = chrono::high_resolution_clock::now();
			vector<uint64_t>  edge_stat = buildConnectedComponents_st(sequences_collisions[i], fa_map, 1);
			validated_edges += edge_stat[0];
			cross_edge_sum += edge_stat[1];
			high_cj_edge_sum += edge_stat[2];
			filter_ed_edge_sum += edge_stat[3];
			pass_ed_edge_sum += edge_stat[4];
            minhash_edge_sum += (sequences_collisions[i].size() * (sequences_collisions[i].size()-1)) >> 1;
            //auto end1 = chrono::high_resolution_clock::now();
            //auto duration1 = chrono::duration_cast<chrono::seconds>(end1 - start1).count();
            //if(duration1 >= per_thread_max_times[tid])
            //{
            //    per_thread_max_times[tid] = duration1;
            //    use_type[tid] = tasks[i];
            //    group_ids[tid] = i;
            //    containing_aas[tid] = total_aas;
            //}
    	}
	}
	auto end_small_time = chrono::high_resolution_clock::now();
	auto duration_small = chrono::duration_cast<chrono::seconds>(end_small_time - start_small_time).count();
	cerr << "Time of single-thread libcdhit (use only 1 threads each group): " << duration_small << endl;
//    for(int i = 0; i < avail_threads; i++){
//        int group_id = group_ids[i];
//        cerr << "times: " << per_thread_max_times[i] << endl;
//        if(use_type[i] == 1) cerr << "type: wt" << endl;
//        else cerr << "type: direct" << endl;
//        cerr << "number of seqs processed: " << sequences_collisions[group_id].size()  << endl;
//        cerr << "number of AAs processed: " << containing_aas[i] << endl;
//        //cerr << "fa name: " << i << endl;
//        //string filename = "fa" + to_string(i);
//	    //ofstream ofs(filename);
//        //for(int j = 0; j < sequences_collisions[group_id].size(); j++) {
//        //    ofs << ">seq" << j << "\n";
//		//    ofs << fa_map.at(sequences_collisions[group_id][j]).c_str() << "\n";
//        //}
//		//ofs.close();
//	}


//	// small collisions in single thread libcdhit
//	auto start_small_time = chrono::high_resolution_clock::now();
//    #pragma omp parallel num_threads(avail_threads) 
//	{
//		#pragma omp for schedule(dynamic)
//    	for(int i = huge_groups_cnt; i < sequences_collisions.size(); i++) {
//			buildConnectedComponents_st(sequences_collisions[i], fa_map, 1);
//    	}
//	}
//	auto end_small_time = chrono::high_resolution_clock::now();
//	auto duration_small = chrono::duration_cast<chrono::seconds>(end_small_time - start_small_time).count();
//	cerr << "Time of single-thread libcdhit (use only 1 threads each group): " << duration_small << endl;

	cerr << "Seqs number processed in mt_libcdhit: " << mt_seqs << endl;
	cerr << "Seqs number processed in st_libcdhit: " << gs_config.items - mt_seqs << endl;

    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cerr << "Number of validated edges build among CC: " << validated_edges << endl;
    cerr << "Number of edges cross origin group among CC: " << cross_edge_sum << endl;
    cerr << "Number of edges get CJ >=0.6 among CC: " << high_cj_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6 and pass EDLIB: " << pass_ed_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6 and no pass EDLIB: " << filter_ed_edge_sum << endl;
    cerr << "Number of edges  CJ < 0.6: " << pass_ed_edge_sum + filter_ed_edge_sum << endl;
    cerr << "Number of total edges in MinHash collisions: " << minhash_edge_sum << endl;
    cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
}

void GroupStream::callLib_cdhit(
	vector<pair<uint32_t, uint32_t>>& minhash_collisions,
	ProteinAAStore& store,
	int cdhit_thres
	){
	bool output_max_group = false;
	if(output_max_group) {
	for(int i = 0; i < 3; i++) {
		int start_pos = minhash_collisions[minhash_collisions.size()-1-i].first;
		int end_pos = minhash_collisions[minhash_collisions.size()-1-i].second + start_pos;
		string filename = "file_" + to_string(i);
		ofstream ofs(filename);
		int end = minhash_collisions[i].second + minhash_collisions[i].first;
		for(int j = start_pos; j < end_pos; j++){
			int seq_id = seq_vec[j].seq_id;
			ofs << ">sequence_" << seq_id << "\n";
			string s = store.get(seq_id);
			ofs<< s << "\n";
		}
		ofs.close();
	}
	return;
	}
	int collision_cnt = minhash_collisions.size();
	auto it_rescure = std::upper_bound(
			minhash_collisions.begin(), minhash_collisions.end(), 500000u,
			[](uint32_t key, const auto& p) { return key < p.second; }
			);
	size_t start_check_idx = it_rescure == minhash_collisions.end() ? collision_cnt : it_rescure - minhash_collisions.begin();
	// minhash_collisions升序排序 找到第一个size(p.first)>1的位置
	auto it_multi = std::upper_bound(
			minhash_collisions.begin(), minhash_collisions.end(), 100000u,
			[](uint32_t key, const auto& p) { return key < p.second; }
			);
	size_t start_multi_idx = it_multi == minhash_collisions.end() ? collision_cnt : it_multi - minhash_collisions.begin();
	if(start_multi_idx < start_check_idx) start_multi_idx = start_check_idx;
	cerr << "Number of groups processed in multi-threading cd-hit: " << collision_cnt - start_multi_idx << endl;

	int avail_threads = gs_config.num_threads;
 	omp_set_num_threads(avail_threads);
	if(it_multi != minhash_collisions.end()) {
		auto start_huge_time = chrono::high_resolution_clock::now();
		for(int i = start_multi_idx; i < minhash_collisions.size(); i++) {
			// 显示进度条（每10组更新一次，或最后一个）
			cerr << "doing task " << i-start_multi_idx << "of " << minhash_collisions.size() - start_multi_idx << endl;
			//调用cdhit
			buildConnectedComponentsByLib_cdhit(
					avail_threads, //threads
					minhash_collisions[i].first, // start_pos
					minhash_collisions[i].second + minhash_collisions[i].first, // end_pos
					store);
		}
		auto end_huge_time = chrono::high_resolution_clock::now();
		auto duration_huge = chrono::duration_cast<chrono::seconds>(end_huge_time - start_huge_time).count();
		cerr << "Time of multi-thread libcdhit (use all threads once): " << duration_huge << endl;
	}

	if(start_multi_idx <= start_check_idx) return;
	cerr << "Number of groups processed in single-threading cd-hit: " << start_multi_idx - start_check_idx << endl;
	int total_small_tasks = start_multi_idx - start_check_idx;
	int one_step = total_small_tasks / 100;
	auto start_small_time = chrono::high_resolution_clock::now();
    #pragma omp parallel num_threads(avail_threads) 
	{
        int tid = omp_get_thread_num();
		#pragma omp for schedule(runtime) 
    	for(int i = start_check_idx; i < start_multi_idx; i++) {
			//调用cdhit
			buildConnectedComponentsByLib_cdhit(
					1, //threads
					minhash_collisions[i].first, // start_pos
					minhash_collisions[i].second + minhash_collisions[i].first, // end_pos
					store);
			if (i % one_step == 0 || i == start_multi_idx - 1) {
				#pragma omp critical
				print_progress(i - start_check_idx + 1, total_small_tasks);
			}
    	}
	}
	auto end_small_time = chrono::high_resolution_clock::now();
	auto duration_small = chrono::duration_cast<chrono::seconds>(end_small_time - start_small_time).count();
	cerr << "Time of single-thread libcdhit (use only 1 threads each group): " << duration_small << endl;

}
	

void GroupStream::Cluster(
	vector<pair<uint32_t, uint32_t>>& need_to_clutser,
	ProteinAAStore& store
	){
	ips2ra::parallel::sort(need_to_clutser.begin(), need_to_clutser.end(), [](const pair<uint32_t, uint32_t>& r){return r.second;}, gs_config.num_threads);
	cluster cluster_cdhit;
	setOptionsSkipAlign(true);
	for(auto& pair : need_to_clutser){
		int start_idx = pair.first;
		int end_idx = pair.second + start_idx;
		auto start_1 = chrono::high_resolution_clock::now();
		vector<Sequence_new> sequences;
		vector<string> seqs;
		seqs.reserve(end_idx-start_idx);
		for(int i = start_idx; i < end_idx; i++) {
			int seq_id = seq_vec[i].seq_id;
			seqs.emplace_back(store.get(seq_id));
			sequences.emplace_back(seq_id, uf.find(seq_id), seqs[i-start_idx].c_str());
		}
		auto end_1 = chrono::high_resolution_clock::now();
		auto duration_1 = chrono::duration_cast<chrono::seconds>(end_1 - start_1).count();
		cerr << "Time of collect " << pair.second << " sequences' groups info: " << duration_1 << endl;

		auto start_2 = chrono::high_resolution_clock::now();
		cluster_cdhit.cdhit_cluster(sequences, gs_config.num_threads);
		auto end_2 = chrono::high_resolution_clock::now();
		auto duration_2 = chrono::duration_cast<chrono::seconds>(end_2 - start_2).count();
		cerr << "Time of cluster " << pair.second << " sequences in cdhit: " << duration_2 << endl;

		for(int i = start_idx, j = 0; i < end_idx; i++, j++) {
			seq_vec[i].seq_id = sequences[j].seq_id;
			seq_vec[i].group_id = sequences[j].new_root_id;
			uf.updateOneParent(seq_vec[i].seq_id, seq_vec[i].group_id);
		}
		ips2ra::parallel::sort(seq_vec.begin() + start_idx, seq_vec.begin() + end_idx, [](const sharedData& r) { return r.group_id; }, gs_config.num_threads);
	}

	//count group size after cd-hit
	int groups_size = uf.countSetsSize();
	cout << "Group Size after using cd-hit to cluster:" << groups_size << endl;
	
	auto start_1 = chrono::high_resolution_clock::now();
	priority_queue<int, vector<int>, greater<int>> minHeap;
	int last_group_id = seq_vec[0].group_id;
	int last_idx = 0;
	for(int i = 0; i < gs_config.items; i++) {
		if(seq_vec[i].group_id != last_group_id){
			minHeap.push(i - last_idx);
			if (minHeap.size() > 10){
				minHeap.pop();
			}
			last_group_id = seq_vec[i].group_id;
			last_idx = i;
			
		}
	}
	minHeap.push(gs_config.items - last_idx);
	if (minHeap.size() > 10){
		minHeap.pop();
	}
	cerr << "Top 10 largest groups size in this round is: ";
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;
	auto end_1 = chrono::high_resolution_clock::now();
	auto duration_1 = chrono::duration_cast<chrono::seconds>(end_1 - start_1).count();
	cerr << "Time of calculate top groups: " << duration_1 << endl;

}
void GroupStream::ClusterLargeThanRescueCondition(
		vector<pair<uint32_t, uint32_t>>& need_to_cluster,
		ProteinAAStore& store,
		bool is_cluster,
		int start_idx,
		int end_idx
	) {

	vector<vector<pair<uint32_t, uint32_t>>> tasks;
	int group_size_ge_condition = 0;
	int group_size_gt1 = 0;
	//for(int i = need_to_cluster.size() - 1; i >= 0 && need_to_cluster[i].second > 1; i--){
	// 遍历minhash_collisions 去掉group_size < 1的组
	for(int i = end_idx - 1; i >= start_idx && need_to_cluster[i].second > 1; i--){
		uint32_t start_idx = need_to_cluster[i].first;
		uint32_t continue_size = need_to_cluster[i].second;
		vector<pair<uint32_t, uint32_t>> one_task;
		// 第一个pair是序列条数 和线程数 先预设再修改
		one_task.emplace_back(continue_size, 1);
		one_task.emplace_back(start_idx, start_idx + continue_size);
		while(continue_size < 10000u && i > 0 && need_to_cluster[i-1].second > 1)
		{
			i--;
			continue_size += need_to_cluster[i].second;
			one_task.emplace_back(need_to_cluster[i].first, need_to_cluster[i].first + need_to_cluster[i].second);
		}
		//根据最终序列条数修改
		one_task[0].first = continue_size;
		if(is_cluster){
			if(continue_size >= 500000u) one_task[0].second = gs_config.num_threads;
        	else if(continue_size >= 20000u) one_task[0].second = 16;
        	else if(continue_size >= 10000u) one_task[0].second = 8;
		}else{
			if(continue_size >= 500000u) one_task[0].second = 24;
        	else if(continue_size >= 100000u) one_task[0].second = 8;
        	else if(continue_size >= 50000u) one_task[0].second = 4;
        	else if(continue_size >= 20000u) one_task[0].second = 2;
		}
		tasks.emplace_back(one_task);
		group_size_gt1 = i;
	}
	cerr << "Number of groups larger than 1: " << need_to_cluster.size() - group_size_gt1 << endl;
	ips2ra::parallel::sort(tasks.begin(), tasks.end(), [](const vector<pair<uint32_t, uint32_t>>& r){return r[0].first;}, gs_config.num_threads);
	cerr << "Number of groups need to cluster in cdhit(after collection small groups into groups>10,000): " << tasks.size() << endl;

	std::atomic<int> thread_pool;
 	int TOTAL_THREADS;
 	TOTAL_THREADS = gs_config.num_threads;
    thread_pool = TOTAL_THREADS;
 	omp_set_num_threads(TOTAL_THREADS);
 	omp_set_nested(1);

	struct ResourceManager {
		int available_threads;
		std::mutex mtx;
		std::condition_variable cv;

		ResourceManager(int total) : available_threads(total) {std::cerr << "Total threads in thread pool: " << available_threads << std::endl;}
	} rm(TOTAL_THREADS);

	auto timestart = chrono::high_resolution_clock::now();
#pragma omp parallel
{
#pragma omp single
{
	for (int i = tasks.size()-1; i >= 0; i--) 
	//for (int i = 0; i < cluster_tasks.size(); i++) 
	{
		//std::cerr << "Current i: " << i << std::endl;
        //if (i % 100 == 0 || i == 0) {
		//	std::cerr << "Before print progress. i: " << i << std::endl;
		//	print_progress(tasks.size()-i, tasks.size());
        //}
        auto& task = tasks[i];
		int required_threads = task[0].second;
#pragma omp task firstprivate(task, required_threads, i)
{
		// 1
		//while (true) {
		//	int available = thread_pool.load(std::memory_order_relaxed);
		//	if (available >= required_threads) {
		//		int prev = thread_pool.fetch_sub(required_threads, std::memory_order_acquire);
		//		if (prev >= required_threads) break;
		//		thread_pool.fetch_add(required_threads, std::memory_order_release);
		//	}
		//	std::this_thread::sleep_for(std::chrono::milliseconds(1));
		//	std::cerr << "task: " << i << " is waitting for resources" << std::endl;
		//}
		
		{
			std::unique_lock<std::mutex> lock(rm.mtx);
			rm.cv.wait(lock, [&] { return rm.available_threads >= required_threads; });
			rm.available_threads -= required_threads;
		}
		
		std::cerr << "launching task: " << i;
		std::cerr << " with " << required_threads;
		std::cerr << " threads, containing " << task[0].first << "sequences." << endl;

		if((tasks.size() - i - 1) % 100 == 0)
			print_progress(tasks.size()-i, tasks.size());

		buildConnectedComponentsByLib_cdhit(required_threads, task, store);
		// 1 
		//thread_pool.fetch_add(required_threads, std::memory_order_release);
		//std::cerr << "finishing task " << i << " finished." << std::endl;

		// 2
		{
			std::lock_guard<std::mutex> lock(rm.mtx);
			rm.available_threads += required_threads;
		}
		rm.cv.notify_all();
		std::cerr << "finishing task " << i << " finished." << std::endl;
}
	}
#pragma omp taskwait
}
}

	auto timeend = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(timeend - timestart).count();
	cerr << "Time of cdhit clustering finally: " << duration << endl;


	if(!is_cluster) {
		//按照vector<vector<pair<>>的形式来组合后不需要全局更新了
		//ips2ra::parallel::sort(seq_vec.begin(), seq_vec.end(), [](const sharedData& r) { return r.group_id; }, gs_config.num_threads);
		//重新整合minhash_collisions
		//need_to_cluster.clear();
		//need_to_cluster.emplace_back(0, seq_vec.size());
		return;
	}

	//最后一轮聚类可以直接写覆盖掉uf
	// 只在finally clustering更新并查集
	// 这么写是为了不去更新那些group_size=1的
	for(auto& task : tasks) {
		for(int i = 1; i < task.size(); i++) {
			for(int j = task[i].first; j < task[i].second; j++) {
				uf.updateOneParent(seq_vec[j].seq_id, seq_vec[j].group_id);
			}
		}
	}
	//count group size after cd-hit
	ips2ra::parallel::sort(seq_vec.begin(), seq_vec.end(), [](const sharedData& r) { return r.group_id; }, gs_config.num_threads);
    int total_clusters = 0;
	priority_queue<int, vector<int>, greater<int>> minHeap;
	int last_group_id = seq_vec[0].group_id;
	int last_idx = 0;
    int only_one = 0;
	for(int i = 0; i < gs_config.items; i++) {
		if(seq_vec[i].group_id != last_group_id){
            total_clusters++;
            if(last_idx == i - 1) only_one++;
			minHeap.push(i - last_idx);
			if (minHeap.size() > 10){
				minHeap.pop();
			}
			last_group_id = seq_vec[i].group_id;
			last_idx = i;
			
		}
	}
    total_clusters++;
    if(last_idx == gs_config.items - 1) only_one++;
	minHeap.push(gs_config.items - last_idx);
	if (minHeap.size() > 10){
		minHeap.pop();
	}
    cerr << "Number of Total Clusters: " << total_clusters << endl;
    cerr << "Number of clusters only having 1 sequence: " << only_one << endl;
	cerr << "Top 10 largest clusters size in this round is: ";
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;

    cerr << "Using uf to check the results " << endl;
	int groups_size = uf.countSetsSize();
	cerr << "Number of clusters after using cd-hit to cluster:" << groups_size << endl;
	for(int i = 0; i < gs_config.items; i++) {
		seq_vec[i].seq_id = i;
		seq_vec[i].group_id = uf.find(i);
    }
    cout << "finish collect!!!" << endl;

    ips2ra::sort(seq_vec.begin(), seq_vec.end(), [](const sharedData& r) { return r.group_id; });
    cout << " finiash sort !!!" << endl;
	last_group_id = seq_vec[0].group_id;
	last_idx = 0;
	for(int i = 0; i < gs_config.items; i++) {
		if(seq_vec[i].group_id != last_group_id){
			minHeap.push(i - last_idx);
			if (minHeap.size() > 10){
				minHeap.pop();
			}
			last_group_id = seq_vec[i].group_id;
			last_idx = i;		
		}
	}
    cout << "finiash update!!!!" << endl;
	minHeap.push(gs_config.items - last_idx);
	if (minHeap.size() > 10){
		minHeap.pop();
	}
	cerr << "Top 10 largest groups size in this round is: ";
	while(!minHeap.empty()){
		cerr << minHeap.top() << " ";
		minHeap.pop();
	}
	cerr << endl;
}
	
void GroupStream::Cluster(
		vector<vector<int>>& cluster_sequences,
		const vector<string>& fa_map
		//const unordered_map<uint64_t, string>& fa_map
	) {
    //init_cnt();
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

void GroupStream::countGroupSize(
	int m, 
	UnionFind& uf, 
	const vector<string>& fa_map, 
	vector<int>& id_root_map) {
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

void GroupStream::Group(
    const ProteinSketchData& sketchdata,
	const ProteinData& proteindata
	) {
	// 用于验证的变量：保存上一轮的roots
	vector<int> prev_round_roots;
	vector<int> cur_round_roots;
	
	//for(int m=0; m < gs_config.M - gs_config.R+1; m++){
	for(int m=0; m < gs_config.M; m++){
		cerr << "round "<<  m << endl;
		fillHashVec(sketchdata, hash_vec, m+gs_config.R);
		//fillHashVec(sketchdata, hash_vec, m);
		GroupByCol(hash_vec, proteindata.sequence_map);
		if(m == gs_config.M - gs_config.R && gs_config.final_cluster_on) {
			gs_config.cluster_condition = 1;
		}
		//countGroupSize(m, uf, proteindata.sequence_map);
		
		// 验证并查集合并结果的正确性
		saveCurrentRoots(cur_round_roots);
		validateUnionFind(round_cnt, prev_round_roots);
		prev_round_roots = cur_round_roots;
		
		round_cnt++;
	}

    outputClstr(proteindata.names);
}

void GroupStream::Group(
	string sketch_filename,
	const ProteinData& proteindata
	) {
    cerr << "tau in libcdhit: " << tau << endl;
	
	// 用于验证的变量：保存上一轮的roots
	vector<int> prev_round_roots;
	vector<int> cur_round_roots;
	
	//for(int m=0; m < gs_config.M - gs_config.R+1; m++){
	for(int m=0; m < gs_config.M; m++){
		cerr << "round "<<  m << endl;
		//fillHashVec(sketchdata, hash_vec, m);
		fillHashVec(sketch_filename, hash_vec, m*gs_config.R);
		GroupByCol(hash_vec, proteindata.sequence_map);
		if(m == gs_config.M-gs_config.R && gs_config.final_cluster_on) {
			gs_config.cluster_condition = 1;
		}
		//countGroupSize(m, uf, proteindata.sequence_map);
		
		// 验证并查集合并结果的正确性
		saveCurrentRoots(cur_round_roots);
		validateUnionFind(round_cnt, prev_round_roots);
		prev_round_roots = cur_round_roots;
		
		round_cnt++;
	}

	//if(gs_config.output_on) {
	//	outputClstr(proteindata.names, proteindata.sequence_map);
	//}
}

void GroupStream::Group(
	string sketch_filename,
	ProteinAAStore& store
	) {
    cerr << "tau in libcdhit: " << tau << endl;
	vector<pair<uint32_t, uint32_t>> minhash_collisions;
	
	// 用于验证的变量：保存上一轮的roots
	vector<int> prev_round_roots;
	vector<int> cur_round_roots;
	
	for(int m=0; m < gs_config.M; m++){
		cerr << "round "<<  m << endl;
		if(m > 0) minhash_collisions.clear();
		fillHashVecAndDetectMinHash(sketch_filename, seq_vec,  m*gs_config.R, minhash_collisions);
		GroupByCol(minhash_collisions, store);

		// 从这里开始挪出GroupByCol
		vector<pair<uint32_t, uint32_t>> need_to_clutser;
		auto start_count = chrono::high_resolution_clock::now();
		if(round_cnt == gs_config.M - 1 && gs_config.final_cluster_on) countGroupSizeBySort(need_to_clutser, 1);
		else if(gs_config.cluster_condition != -1) countGroupSizeBySort(need_to_clutser, gs_config.cluster_condition);
		auto end_count = chrono::high_resolution_clock::now();
		auto duration_count = chrono::duration_cast<chrono::seconds>(end_count - start_count).count();
		cerr << "Time of count group size: " << duration_count << endl;
		if(need_to_clutser.size() > 0) {
			if(round_cnt == gs_config.M-1) {
				cerr << "Already in final clustering" << endl;
				setOptionsSkipAlign(false);
				ClusterLargeThanRescueCondition(need_to_clutser, store, true, 0, need_to_clutser.size());
			} else {
				Cluster(need_to_clutser, store); // 进rescue-mode
			}
		}

		// 验证并查集合并结果的正确性
		saveCurrentRoots(cur_round_roots);
		validateUnionFind(round_cnt, prev_round_roots);
		
		// 保存当前轮的roots作为下一轮的prev_roots
		prev_round_roots = cur_round_roots;

		round_cnt++;
	}

    outputClstr(store);
}

vector<uint64_t> GroupStream::buildConnectedComponents_st(
	vector<pair<int,int>>& group_seqs, 
	ProteinAAStore& store,
	int use_wt
	) {
	vector<Sequence_new> sequences;
	vector<string> seqs;
	seqs.reserve(group_seqs.size());
	for(int i = 0; i < group_seqs.size(); i++) {
		int seq_id = group_seqs[i].first;
		seqs.emplace_back(store.get(seq_id));
		sequences.emplace_back(seq_id, uf.find(seq_id), seqs[i].c_str());
	}

    vector<uint64_t> edge_stat;
    if(use_wt == 1){
		edge_stat = cluster_sequences_st(sequences, 5, tau); 
    }else{
	    cluster_sequences_st_less10(sequences, 5, tau); 
    }

	for(int i = 0; i < group_seqs.size(); i++)
	{
		group_seqs[i].first = sequences[i].seq_id;
		group_seqs[i].second = sequences[i].new_root_id;
	}
	sort(group_seqs.begin(), group_seqs.end(), [](const pair<int, int>& a, const pair<int, int>& b){
		return a.second < b.second;
	});
    return edge_stat;
}

vector<uint64_t> GroupStream::buildConnectedComponents_st(
	vector<pair<int,int>>& group_seqs, 
	const vector<string>& fa_map,
    int use_wt
	) {
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		int seq_id = group_seqs[i].first;
		//int seq_id = group_seqs[i];
		sequences.emplace_back(seq_id, uf.find(seq_id), fa_map.at(seq_id).c_str());
	}

    vector<uint64_t> edge_stat = {0, 0};
    if(use_wt == 1){
		edge_stat = cluster_sequences_st(sequences, 5, tau); 
    }else{
	    cluster_sequences_st_less10(sequences, 5, tau); 
    }

	for(int i = 0; i < group_seqs.size(); i++)
	{
		group_seqs[i].first = sequences[i].seq_id;
		group_seqs[i].second = sequences[i].new_root_id;
		//id_root_map[sequences[i].seq_id] = sequences[i].new_root_id;
	}
	sort(group_seqs.begin(), group_seqs.end(), [](const pair<int, int>& a, const pair<int, int>& b){
		return a.second < b.second;
	});
    return edge_stat;
}

vector<uint64_t> GroupStream::buildConnectedComponents(
	int needed_threads, 
	uint32_t start_idx,
	uint32_t end_idx,
	ProteinAAStore& store
	) {
	// TODO 直接把这个sequence替换成seq_vec传进去
	vector<Sequence_new> sequences;
	vector<string> seqs;
	seqs.reserve(end_idx-start_idx);
	for(int i = start_idx; i < end_idx; i++) {
		int seq_id = seq_vec[i].seq_id;
		seqs.emplace_back(store.get(seq_id));
		sequences.emplace_back(seq_id, uf.find(seq_id), seqs[i-start_idx].c_str());
		//sequences[i-start_idx].length = seqs[i-start_idx].c_str().size();
	}

	vector<uint64_t> edge_stat;
	if(needed_threads == 1) {
		edge_stat = cluster_sequences_new_st(sequences, 5, tau, ed_thres); 
	}else {
		edge_stat = cluster_sequences_new(sequences, 5, tau, ed_thres, needed_threads); 
	}

	for(int i = start_idx, j = 0; i < end_idx; i++, j++) {
		seq_vec[i].seq_id = sequences[j].seq_id;
		seq_vec[i].group_id = sequences[j].new_root_id;
	}

	// TODO 可以挪到libcc中直接对sequences排序
	// 对做了libcc的结果按照group_id排序
    ips2ra::sort(seq_vec.begin() + start_idx, seq_vec.begin() + end_idx, [](const sharedData& r) { return r.group_id; });
	return edge_stat;
}

void GroupStream::buildConnectedComponentsByLib_cdhit(
	int needed_threads, 
	vector<pair<uint32_t, uint32_t>>& task, //每个pair first的起始位置 second是结束位置
	ProteinAAStore& store
	) {
	vector<Sequence_new> sequences;
	vector<string> seqs;
	seqs.reserve(task[0].first);
	int cnt = 0;

	//更新结果到seq_vec
	for(int i = 1; i < task.size(); i++) {
		// j是序列在seq_vec的偏移
		// 记录这个是为了方便把cdhit的结果写回去
		for(int j = task[i].first; j < task[i].second; j++) {
			int seq_id = seq_vec[j].seq_id;
			seqs.emplace_back(store.get(seq_id));
			sequences.emplace_back(seq_id, uf.find(seq_id), seqs[cnt].c_str());
			cnt++;
		}
	}
	
	cluster cluster_cdhit;
	cluster_cdhit.cdhit_cluster(sequences, needed_threads);

	cnt = 0;
	for(int i = 1; i < task.size(); i++) {
		for(int j = task[i].first; j < task[i].second; j++) {
			seq_vec[j].seq_id = sequences[cnt].seq_id;
			seq_vec[j].group_id = sequences[cnt].new_root_id;
			cnt++;
		}
	}
	//for(int i = 0; i < sequences.size(); i++) {
	//	int id_in_seq_vec = sequences[i].seq_id;
	//	seq_vec[id_in_seq_vec].group_id = sequences[i].new_root_id;
	//}
}


void GroupStream::buildConnectedComponentsByLib_cdhit(
	int needed_threads, 
	uint32_t start_idx,
	uint32_t end_idx,
	ProteinAAStore& store
	) {
	vector<Sequence_new> sequences;
	vector<string> seqs;
	seqs.reserve(end_idx-start_idx);
	for(int i = start_idx; i < end_idx; i++) {
		int seq_id = seq_vec[i].seq_id;
		seqs.emplace_back(store.get(seq_id));
		sequences.emplace_back(seq_id, uf.find(seq_id), seqs[i-start_idx].c_str());
	}

	cluster cluster_cdhit;
	cluster_cdhit.cdhit_cluster(sequences, needed_threads);

	for(int i = start_idx, j = 0; i < end_idx; i++, j++) {
		seq_vec[i].seq_id = sequences[j].seq_id;
		seq_vec[i].group_id = sequences[j].new_root_id;
	}
    ips2ra::parallel::sort(seq_vec.begin() + start_idx, seq_vec.begin() + end_idx, [](const sharedData& r) { return r.group_id; }, needed_threads);
}

vector<uint64_t> GroupStream::buildConnectedComponents(
	vector<pair<int,int>>& group_seqs, 
	int needed_threads,
	const vector<string>& fa_map,
	ClusterWS& ws
	) {
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		int seq_id = group_seqs[i].first;
		sequences.emplace_back(seq_id, uf.find(seq_id), fa_map.at(seq_id).c_str());
	}

	vector<uint64_t> edge_stat;
	if(needed_threads > 1) {
		edge_stat = cluster_sequences(sequences, 5, tau, needed_threads); 
	}else {
		if(group_seqs.size() < 100){
			cluster_sequences_st_less10(sequences, 5, tau); 
		}else{
			cluster_sequences_st_reuse(sequences, 5, tau, ws); 
		}
	}
	for(int i = 0; i < group_seqs.size(); i++)
	{
		group_seqs[i].first = sequences[i].seq_id;
		group_seqs[i].second = sequences[i].new_root_id;
		//id_root_map[sequences[i].seq_id] = sequences[i].new_root_id;
	}
	sort(group_seqs.begin(), group_seqs.end(), [](const pair<int, int>& a, const pair<int, int>& b){
		return a.second < b.second;
	});
    return edge_stat;
}

void GroupStream::clusterEachGroup(
	vector<int>& group_seqs,
	int needed_threads,
	const vector<string>& fa_map
	//const unordered_map<uint64_t, string>& fa_map
	) {
	auto start_time_build = chrono::high_resolution_clock::now();
	vector<Sequence_new> sequences;
	for(int i = 0; i < group_seqs.size(); i++) {
		sequences.emplace_back(group_seqs[i], uf.find(group_seqs[i]), fa_map.at(group_seqs[i]).c_str());
	}
	auto end_time_build = chrono::high_resolution_clock::now();
    auto duration_build = chrono::duration_cast<chrono::seconds>(end_time_build - start_time_build).count();

	//读取FAI获取data

	//auto start_time = chrono::high_resolution_clock::now();
		//cluster cluster_cdhit;
		//cluster_cdhit.cdhit_cluster(sequences, id_root_map, needed_threads);
	//auto end_time = chrono::high_resolution_clock::now();
}

void GroupStream::outputClstr(
	const vector<string>& names
) {
	cerr << "Total Clusters: " << uf.countSetsSize() << endl;
	cerr << "cluster result stored: " << gs_config.res_file << endl;
	string out_file_name = "round_" + to_string(round_cnt) + "_" + gs_config.res_file;
	cerr << "resulf of round " << round_cnt << "write to: " << out_file_name << endl;
	ofstream ofs(out_file_name);
	//uf.findRoot(id_root_map);
	for(int i = 0; i < gs_config.items; i++) {
		ofs << ">" << names[i] << " " << ">" << names[uf.find(i)] << endl;
		//ofs << ">" << names[i] << " " << ">" << names[id_root_map[i]] << endl;
	}
	ofs.close();
	//cerr << "Total Clusters: " << uf.countSetsSize() << endl;
	//cerr << "cluster result stored: " << gs_config.res_file << endl;
	//ofstream seq_id(gs_config.res_file);
	//streambuf* origin_cout = cout.rdbuf();
	//cout.rdbuf(seq_id.rdbuf());

	//uf.findRoot(id_root_map);
	//for(int i = 0; i < gs_config.items; i++) {
	//	cout << ">" << names[i] << " " << ">" << names[id_root_map[i]] << endl;
	//}
	//cout.rdbuf(origin_cout);
}
// 保存当前round的所有节点的根节点
void GroupStream::saveCurrentRoots(vector<int>& roots) {
	roots.resize(gs_config.items);
	for(int i = 0; i < gs_config.items; i++) {
		roots[i] = uf.find(i);
	}
}

// 逆向验证：从已合并的uf结果出发，检查每个合并关系是否能追溯到上一轮或当前轮哈希
bool GroupStream::validateUnionFind(int round_num, const vector<int>& prev_roots) {
	cerr << "========== 逆向验证 Round " << round_num << " 合并结果 ==========" << endl;
	
	bool is_valid = true;
	int error_count = 0;
	const int MAX_ERRORS_TO_SHOW = 10;

	// 构建"合法边"的并查集：只包含上一轮和当前轮哈希分组产生的连通关系
	vector<int> valid_parent(gs_config.items);
	iota(valid_parent.begin(), valid_parent.end(), 0);
	
	std::function<int(int)> valid_find = [&](int x) -> int {
		if(valid_parent[x] != x) {
			valid_parent[x] = valid_find(valid_parent[x]);
		}
		return valid_parent[x];
	};
	
	auto valid_unite = [&](int x, int y) {
		int rx = valid_find(x);
		int ry = valid_find(y);
		if(rx != ry) {
			valid_parent[ry] = rx;
		}
	};
	
	// 1) 添加上一轮的合法连通关系
	int prev_edges = 0;
	if(!prev_roots.empty() && prev_roots.size() == gs_config.items) {
		unordered_map<int, vector<int>> prev_groups; // prev_root -> [节点列表]
		for(int i = 0; i < gs_config.items; i++) {
			prev_groups[prev_roots[i]].push_back(i);
		}
		for(auto& [root, nodes] : prev_groups) {
			for(int k = 1; k < nodes.size(); k++) {
				valid_unite(nodes[0], nodes[k]);
				prev_edges++;
			}
		}
		cerr << "  上一轮合法连通边数: " << prev_edges << " (来自 " << prev_groups.size() << " 个组)" << endl;
	}
	
	// 2) 添加当前轮seq_vec哈希分组的合法连通关系
	int cur_edges = 0;
	unordered_map<int, vector<int>> cur_hash_groups; // group_id -> [seq_id列表]
	for(int i = 0; i < gs_config.items; i++) {
		cur_hash_groups[seq_vec[i].group_id].push_back(seq_vec[i].seq_id);
	}
	for(auto& [group_id, seq_ids] : cur_hash_groups) {
		for(int k = 1; k < seq_ids.size(); k++) {
			valid_unite(seq_ids[0], seq_ids[k]);
			cur_edges++;
		}
	}
	cerr << "  当前轮哈希合法连通边数: " << cur_edges << " (来自 " << cur_hash_groups.size() << " 个哈希组)" << endl;
	
	// 3) 逆向验证：对于实际uf中同组的节点，检查它们是否在合法并查集中也同组
	//    如果uf中同组但合法并查集中不同组，说明有"非法合并"
	cerr << "  开始逆向验证..." << endl;
	
	// 按照uf的分组来检查
	unordered_map<int, vector<int>> uf_groups; // uf_root -> [节点列表]
	for(int i = 0; i < gs_config.items; i++) {
		uf_groups[uf.find(i)].push_back(i);
	}
	
	int invalid_merges = 0;
	for(auto& [uf_root, nodes] : uf_groups) {
		if(nodes.size() <= 1) continue;
		
		// 检查这个组内的节点是否都能通过合法边连通
		int first_valid_root = valid_find(nodes[0]);
		for(int k = 1; k < nodes.size(); k++) {
			int cur_valid_root = valid_find(nodes[k]);
			if(cur_valid_root != first_valid_root) {
				// 找到了非法合并：uf中同组，但上一轮和当前轮哈希都无法解释
				if(error_count < MAX_ERRORS_TO_SHOW) {
					cerr << "  非法合并: 节点 " << nodes[0] << " 和 " << nodes[k] 
						 << " 在uf中同组(root=" << uf_root << ")，但无法通过上一轮或当前轮哈希追溯" << endl;
					// 输出更多调试信息
					if(!prev_roots.empty()) {
						cerr << "    上一轮: prev_root[" << nodes[0] << "]=" << prev_roots[nodes[0]] 
							 << ", prev_root[" << nodes[k] << "]=" << prev_roots[nodes[k]] << endl;
					}
					cerr << "    当前轮: seq_vec中 " << nodes[0] << " 的group_id=";
					for(int i = 0; i < gs_config.items; i++) {
						if(seq_vec[i].seq_id == nodes[0]) { cerr << seq_vec[i].group_id; break; }
					}
					cerr << ", " << nodes[k] << " 的group_id=";
					for(int i = 0; i < gs_config.items; i++) {
						if(seq_vec[i].seq_id == nodes[k]) { cerr << seq_vec[i].group_id; break; }
					}
					cerr << endl;
				}
				error_count++;
				invalid_merges++;
				is_valid = false;
			}
		}
	}
	
	// 4) 统计结果
	unordered_set<int> valid_groups, actual_groups;
	for(int i = 0; i < gs_config.items; i++) {
		valid_groups.insert(valid_find(i));
		actual_groups.insert(uf.find(i));
	}
	
	cerr << "  合法连通后组数: " << valid_groups.size() << endl;
	cerr << "  实际uf组数: " << actual_groups.size() << endl;
	
	if(error_count >= MAX_ERRORS_TO_SHOW) {
		cerr << "  ... 共发现 " << invalid_merges << " 个非法合并，只显示前 " << MAX_ERRORS_TO_SHOW << " 个" << endl;
	}
	
	if(is_valid) {
		cerr << "  通过: 所有uf中的合并都可追溯到上一轮或当前轮哈希分组" << endl;
	} else {
		cerr << "  失败: 存在 " << invalid_merges << " 个无法追溯的非法合并" << endl;
	}
	
	cerr << "========== 验证完成，结果: " << (is_valid ? "通过" : "存在错误") << " ==========" << endl;
	return is_valid;
}

void GroupStream::outputClstr(
    ProteinAAStore& store
) {
	cerr << "Total Clusters: " << uf.countSetsSize() << endl;
	auto start_time = chrono::high_resolution_clock::now();
//  test single rounds results
//	string out_file_name = "round_" + to_string(round_cnt) + "_" + gs_config.res_file;
//	cerr << "resulf of round " << round_cnt << "write to: " << out_file_name << endl;
//	ofstream ofs(out_file_name);
//    if(gs_config.names_path != "") {
//        store.load_names(gs_config.names_path);
//        for(int i = 0; i < gs_config.items; i++) {
//            ofs << ">" << store.name(i) << " " << ">" << store.name(uf.find(i)) << endl;
//        }
//    }else{
//        cerr << "The name_path is not provided!" << endl;
//        //输出seq-id
//        for(int i = 0; i < gs_config.items; i++) {
//            ofs << i << " " << uf.find(i) << "\n";
//        }
//    }
//	ofs.close();

	cerr << "cluster result stored: " << gs_config.res_file << endl;
	ofstream ofs(gs_config.res_file);
    if(gs_config.names_path != "") {
		cerr << "load name from " << gs_config.names_path << endl;
        store.load_names(gs_config.names_path);
        for(int i = 0; i < gs_config.items; i++) {
            ofs  << store.name(uf.find(i)) << " " << store.name(i) << endl;
        }
    }else{
        cerr << "The name_path is not provided!" << endl;
        //输出seq-id
        for(int i = 0; i < gs_config.items; i++) {
            ofs << uf.find(i) << " " << i << "\n";
        }
    }
	auto end_time = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(end_time - start_time).count();
	cerr << "output time (seconds): " << duration  << endl;
    ofs.close();
}
