#include "KHFMinHash.h"
#include "kseq.h"
#include <zlib.h>

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <chrono>
#include <random>
#include <queue>
#include <thread>
#include <mutex>
#include <atomic>
#include <pthread.h>
#include <filesystem>
#include <unordered_map>
#include "util.h"

KSEQ_INIT(gzFile, gzread)

std::mutex mtx1;
std::mutex mtx2;
std::atomic<int> num_seqs(0);

vector<uint64_t> seq_ids;
vector<vector<uint64_t>> hashes;
vector<string> names;

vector<vector<string>> thr_names;
vector<vector<vector<uint64_t>>> thr_hashes;
vector<vector<uint64_t>> thr_seq_ids;
unordered_map<uint64_t, string> fa_map;

pthread_barrier_t barrier;
struct ThreadData {
	int tid;
	kseq_t* ks;
	int k, m, min_len;
	bool xxhash_flag;
};
//void consumer(int tid, gzFile fp, kseq_t* ks, int k, int m, bool xxhash_flag, int min_len, std::barrier<>& barrier) {
void* consumer(void* arg) {
	ThreadData* data = static_cast<ThreadData*>(arg);
	int tid = data->tid;
	kseq_t* ks = data->ks;
	int k = data->k;
	int m = data->m;
	bool xxhash_flag = data->xxhash_flag;
	int min_len = data->min_len;

	vector<vector<uint64_t>> local_hashes(m);  // 创建 m 列，每列存储一个哈希
	vector<uint64_t> local_seq_ids;
	vector<string> local_names;
		
    while (true) {
        std::string sequence;
        std::string name;
        int seq_id;

		{
				std::lock_guard<std::mutex> lock(mtx1);
				int length = kseq_read(ks);
				if (length < 0) break;
				if (length < min_len) continue;
				sequence = ks->seq.s;//direct copy?
				name = ks->name.s;
				seq_id = num_seqs.fetch_add(1);
		}

		Sketch::KHFMinHash mh;
		mh.setK(k);
		mh.setM(m);

		if (xxhash_flag)
				mh.buildSketch(sequence.c_str());
		else
				mh.buildSketchByNoSeedAAHash(sequence.c_str());

		auto& sketch = mh.getSektch();
		
		local_seq_ids.emplace_back(seq_id);
		local_names.emplace_back(name);
		for(int i=0; i<m; i++){
			local_hashes[i].emplace_back(sketch.hashes[i]);
		}

	}
	//barrier.arrive_and_wait();
	pthread_barrier_wait(&barrier);
	// 线程完成后汇总其结果
	{
		std::lock_guard<std::mutex> lock(mtx2);
		names.insert(names.end(), std::make_move_iterator(local_names.begin()), std::make_move_iterator(local_names.end()));
		seq_ids.insert(seq_ids.end(), std::make_move_iterator(local_seq_ids.begin()), std::make_move_iterator(local_seq_ids.end()));
		for(int i = 0; i < m; i++) {
			hashes[i].insert(hashes[i].end(), std::make_move_iterator(local_hashes[i].begin()), std::make_move_iterator(local_hashes[i].end()));
		}
	}
	return nullptr;
}

void reorder_swap(vector<uint64_t>& matrix, vector<uint64_t>& indices) {
    size_t n = matrix.size();
    if (n != indices.size() || n == 0) return;

    for (size_t i = 0; i < n; ++i) {
		if(i == indices[i]) continue;

        uint64_t curr_pos = i;          // 当前行索引（uint64_t 类型）
        uint64_t next_pos = indices[static_cast<size_t>(curr_pos)];  // 下一个位置

        while (next_pos != curr_pos) {
			// 交换 matrix 的行指针
			swap(matrix[static_cast<size_t>(curr_pos)], matrix[static_cast<size_t>(next_pos)]);
			uint64_t finised_pos = next_pos;
			next_pos = indices[static_cast<size_t>(next_pos)];
			indices[static_cast<size_t>(finised_pos)] = finised_pos;

		}

		indices[i] = i;
	}
}

void compactSeqIds(int threads)
{
	for(int i = 0;  i < threads; i++)
	{
		seq_ids.insert(seq_ids.end(), make_move_iterator(thr_seq_ids[i].begin()), make_move_iterator(thr_seq_ids[i].end()));
	}
	ostringstream oss;
	for(int i = 0; i < seq_ids.size(); i++){
		oss << seq_ids[i] << " ";
	}
	ofstream ofs("seq_ids");
	if(!ofs){
		cerr << "Error opening file!" << endl;
		return;
	}
	ofs << oss.str();
	ofs.close();
}

template <typename T>
void write_to_file(const std::vector<std::vector<T>>& data, const std::string& filename)
{
	std::ofstream ofs(filename, std::ios::binary);
	if (!ofs) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}
	std::cout << "Sketch write to : " << filename << endl;

	for (int i = 0; i < data.size(); i++) {
		auto& row = data[i];
		ofs.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(T));
		std::cout << "Writing " << row.size() * sizeof(T) << " bytes for row" << std::endl;
	}

}

int buildSketches(string filename, string output_filename, int num_threads, int k, int m, bool xxhash_flag, int min_len, bool output_on)
{
	
	gzFile fp1;
	kseq_t *ks1;

	fp1 = gzopen(filename.c_str(),"r");

	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return 0;
	}
	ks1 = kseq_init(fp1);

	// std=c++20
	//std::barrier sync_point(num_threads);
    //std::vector<std::thread> threads;
	pthread_t threads[num_threads];
	pthread_barrier_init(&barrier, nullptr, num_threads);
	hashes.assign(m, {});
    for (int i = 0; i < num_threads; ++i) {
   	    //threads.emplace_back(consumer, i, fp1, ks1, k, m, xxhash_flag, min_len, std::ref(sync_point));
		ThreadData data;
		data.tid = i;
		data.ks = ks1;
		data.k = k;
		data.m = m;
		data.xxhash_flag = xxhash_flag;
		data.min_len = min_len;
		pthread_create(&threads[i], nullptr, consumer, &data);
   	}

    //for (auto& t : threads) {
    //    t.join();
    //}

    for (int i = 0; i < num_threads; ++i) {
		pthread_join(threads[i], nullptr);
	}
	pthread_barrier_destroy(&barrier);
	cout << "Finish building sketch!" << endl;
	// reorder
	if(num_threads > 1)
	{
		for(int i = 0; i < m; i++) {
			vector<uint64_t> tmp_seq_ids(seq_ids);
			reorder_swap(hashes[i], tmp_seq_ids);
		}
	}
	string sketch_file_name  = output_filename;
	if(output_on) {
		write_to_file(hashes, sketch_file_name);
	}
    gzclose(fp1);
    kseq_destroy(ks1);
	return seq_ids.size();
}

int read_fa(string filename, int min_len)
{
	int cnt_seqs = 0;	
	gzFile fp1;
	kseq_t *ks1;

	fp1 = gzopen(filename.c_str(),"r");
	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return 0;
	}
	ks1 = kseq_init(fp1);

    while (true) {
		int len = kseq_read(ks1);
		if(len < 0) break;
		if(len < min_len) continue;
		fa_map[cnt_seqs] = ks1->seq.s;
		names.emplace_back(ks1->name.s);
		cnt_seqs++;
	}
    gzclose(fp1);
    kseq_destroy(ks1);
	return cnt_seqs;
}
