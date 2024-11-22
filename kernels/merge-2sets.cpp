#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cstdint>
#include <chrono>
#include <unordered_set>

using namespace std;
using namespace chrono;

// 并查集类
class UnionFind {
	private:
		vector<uint64_t> parent; // 父节点数组
		vector<uint64_t> rank;   // 秩数组

	public:
		// 初始化
		UnionFind(size_t n) {
			parent.resize(n);
			rank.resize(n, 0);
			for (size_t i = 0; i < n; ++i) {
				parent[i] = i;
			}
		}

		// 查找操作（路径压缩）
		uint64_t find(uint64_t x) {
			if (x != parent[x]) {
				parent[x] = find(parent[x]);
			}
			return parent[x];
		}

		// 合并操作（按秩合并）
		void unite(uint64_t x, uint64_t y) {
			uint64_t rootX = find(x);
			uint64_t rootY = find(y);

			if (rootX != rootY) {
				if (rank[rootX] > rank[rootY]) {
					parent[rootY] = rootX;
				} else if (rank[rootX] < rank[rootY]) {
					parent[rootX] = rootY;
				} else {
					parent[rootY] = rootX;
					rank[rootX]++;
				}
			}
		}

		// 检查两个元素是否在同一个集合
		bool connected(uint64_t x, uint64_t y) {
			return find(x) == find(y);
		}
		// 获取集合的总数
		size_t countSets() {
			unordered_set<uint64_t> unique_roots;
			for (size_t i = 0; i < parent.size(); ++i) {
				unique_roots.insert(find(i));
			}
			return unique_roots.size();
		}

		size_t countSets_fast() {
			size_t count = 0;
			for (size_t i = 0; i < parent.size(); ++i) {
				if(parent[i] == i) count++;
			}
			return count;
		}
};

int main() {
	const size_t NUM_ELEMENTS = 300000000; // 300M for each vector
	const size_t MAX_SET_ID   =  50000000;    // 最大 set_id = hash space
	vector<pair<uint64_t, uint64_t>> data;
	vector<pair<uint64_t, uint64_t>> extra_data;

	// 随机数生成器
	random_device rd;
	mt19937_64 gen(rd());
	uniform_int_distribution<uint64_t> set_dist(0, MAX_SET_ID - 1);
	uniform_int_distribution<uint64_t> elem_dist(0, NUM_ELEMENTS - 1);

	// 1. 数据生成计时
	cout << "Generating 300M random data for the main vector..." << endl;
	auto start_gen = high_resolution_clock::now();

	data.reserve(NUM_ELEMENTS);
	for (size_t i = 0; i < NUM_ELEMENTS; ++i) {
		data.emplace_back(set_dist(gen), i);
	}

	auto end_gen_main = high_resolution_clock::now();
	auto gen_main_duration = duration_cast<seconds>(end_gen_main - start_gen);
	cout << "Main vector data generation complete. Time taken: " << gen_main_duration.count() << " seconds." << endl;

	cout << "Generating 300M random data for the additional vector..." << endl;
	auto start_gen_extra = high_resolution_clock::now();

	extra_data.reserve(NUM_ELEMENTS);
	for (size_t i = 0; i < NUM_ELEMENTS; ++i) {
		extra_data.emplace_back(set_dist(gen), i);
	}

	auto end_gen_extra = high_resolution_clock::now();
	auto gen_extra_duration = duration_cast<seconds>(end_gen_extra - start_gen_extra);
	cout << "Additional vector data generation complete. Time taken: " << gen_extra_duration.count() << " seconds." << endl;

	// 2. 排序计时
	cout << "Sorting both vectors by set_id..." << endl;
	auto start_sort = high_resolution_clock::now();

	sort(data.begin(), data.end(), [](const pair<uint64_t, uint64_t>& a, const pair<uint64_t, uint64_t>& b) {
			return a.first < b.first;
			});

	sort(extra_data.begin(), extra_data.end(), [](const pair<uint64_t, uint64_t>& a, const pair<uint64_t, uint64_t>& b) {
			return a.first < b.first;
			});

	auto end_sort = high_resolution_clock::now();
	auto sort_duration = duration_cast<seconds>(end_sort - start_sort);
	cout << "Sorting complete. Time taken: " << sort_duration.count() << " seconds." << endl;

	// 3. 并查集合并计时
	cout << "Merging elements with the same set_id from both vectors..." << endl;

	UnionFind uf(NUM_ELEMENTS); // 初始化并查集

	auto start_merge = high_resolution_clock::now();
	uint64_t unite_count = 0;

	// 合并主数据的集合
	for (size_t i = 0; i < data.size();) {
		uint64_t currentSetID = data[i].first;

		size_t j = i;
		while (j < data.size() && data[j].first == currentSetID) {
			j++;
		}

		for (size_t k = i + 1; k < j; ++k) {
			uf.unite(data[i].second, data[k].second);
		}
		
		unite_count += j - (i + 1);

		i = j;
	}

	// 合并额外数据的集合
	for (size_t i = 0; i < extra_data.size();) {
		uint64_t currentSetID = extra_data[i].first;

		size_t j = i;
		while (j < extra_data.size() && extra_data[j].first == currentSetID) {
			j++;
		}

		for (size_t k = i + 1; k < j; ++k) {
			uf.unite(extra_data[i].second, extra_data[k].second);
		}

		unite_count += j - (i + 1);

		i = j;
	}

	auto end_merge = high_resolution_clock::now();
	auto merge_duration = duration_cast<seconds>(end_merge - start_merge);
	cout << "Merging complete. Time taken: " << merge_duration.count() << " seconds." << endl;

	//// 4. 统计集合总数
	//cout << "Counting total number of sets..." << endl;
	//auto start_count = high_resolution_clock::now();

	//size_t total_sets = uf.countSets();

	//auto end_count = high_resolution_clock::now();
	//auto count_duration = duration_cast<seconds>(end_count - start_count);
	//cout << "Counting complete. Time taken: " << count_duration.count() << " seconds." << endl;

	//cout << "Total number of sets: " << total_sets << endl;
	// 5. 统计集合总数
	cout << "Counting total number of sets..." << endl;
	auto start_count2 = high_resolution_clock::now();

	size_t total_sets2 = uf.countSets_fast();

	auto end_count2 = high_resolution_clock::now();
	auto count_duration2 = duration_cast<seconds>(end_count2 - start_count2);
	cout << "Fast counting complete. Time taken: " << count_duration2.count() << " seconds." << endl;

	cout << "Total number of sets (fast): " << total_sets2 << endl;

	// 总计时
	auto total_duration = duration_cast<seconds>(end_count2 - start_gen);
	cout << "Total time taken: " << total_duration.count() << " seconds." << endl;
	cout << "Calling unite count: " << unite_count << endl;

	return 0;
}

