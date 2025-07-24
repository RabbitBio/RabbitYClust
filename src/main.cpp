#include "KHFMinHash.h"
#include "kseq.h"
#include "CLI11.hpp"
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

#include "GroupStream.h"

#define BUFFER_SIZE (1<<20 * sizeof(char))

using namespace std;
using namespace Sketch;

KSEQ_INIT(gzFile, gzread)

std::mutex mtx1;
std::mutex mtx2;
std::atomic<int> num_seqs(0);

vector<vector<uint64_t>> hashes;
vector<uint64_t> seq_ids;
vector<string> names;
// yy add for cluster
//bool cluster_on = false;
unordered_map<uint64_t, uint64_t> fai_map;
unordered_map<uint64_t, string> fa_map;
int64_t pos = 0;

vector<vector<vector<uint64_t>>> all_hashes;
vector<vector<uint64_t>> thr_seq_ids;
struct compare {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

void consumer(int tid, gzFile fp, kseq_t* ks, int k, int m, bool xxhash_flag, int min_len) {
    // 每个线程有自己独立的按列存储的 vector<vector<uint64_t>> 用于存储该线程的结果
	vector<vector<uint64_t>> local_hashes(m);  // 创建 m 列，每列存储一个哈希
	vector<uint64_t> local_seq_ids;  // 创建 m 列，每列存储一个哈希
		
    while (true) {
        std::string sequence;
        int seq_id;

		{
				std::lock_guard<std::mutex> lock(mtx1);
				int length = kseq_read(ks);
				if (length < 0) break;
				if (length < min_len) continue;
				sequence = ks->seq.s;//direct copy?
				seq_id = num_seqs.fetch_add(1);
		}

		KHFMinHash mh;
		mh.setK(k);
		mh.setM(m);

		if (xxhash_flag)
				mh.buildSketch(sequence.c_str());
		else
				mh.buildSketchByNoSeedAAHash(sequence.c_str());

		auto& sketch = mh.getSektch();
		
		local_seq_ids.emplace_back(seq_id);
		for(int i=0; i<m; i++){
			local_hashes[i].emplace_back(sketch.hashes[i]);
		}

	}

	    // 线程完成后汇总其结果
	{
		std::lock_guard<std::mutex> lock(mtx2);
		all_hashes[tid] = std::move(local_hashes);  // 将线程的结果添加到全局容器
		thr_seq_ids[tid] = std::move(local_seq_ids);
	}
}

void reorderRowsSwap(vector<vector<uint64_t>>& matrix, vector<uint64_t>& indices) {
    size_t n = matrix.size();
    if (n != indices.size() || n == 0) return;

    for (size_t i = 0; i < n; ++i) {
		if(i == indices[i]) continue;

        uint64_t curr_pos = i;          // 当前行索引（uint64_t 类型）
        uint64_t next_pos = indices[static_cast<size_t>(curr_pos)];  // 下一个位置

        while (next_pos != curr_pos) {
				// 交换 matrix 的行指针
				swap(matrix[static_cast<size_t>(curr_pos)], matrix[static_cast<size_t>(next_pos)]);
				swap(names[static_cast<size_t>(curr_pos)], names[static_cast<size_t>(next_pos)]);
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
	
void outputAllHashes(int m, int threads)
{
	for(int i = 0; i < m; i++){
		ostringstream oss;
		for(int t = 0; t < threads; t++){
			for(int j = 0; j < all_hashes[t][i].size(); j++){
				oss << all_hashes[t][i][j] << " ";
			}
			oss << "\n";
		}
		string sketch_file_name = "sketch"+ to_string(i);
		ofstream ofs(sketch_file_name);
		if(!ofs){
			cerr << "Error opening file!" << endl;
			return;
		}
		ofs << oss.str();
		ofs.close();
	}
}

void buildSketches(string filename, int num_threads, int k, int m, bool xxhash_flag, int min_len)
{
	
	gzFile fp1;
	kseq_t *ks1;

	fp1 = gzopen(filename.c_str(),"r");

	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return ;
	}

	
	ks1 = kseq_init(fp1);

	all_hashes.resize(num_threads);
	thr_seq_ids.resize(num_threads);
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
   	    threads.emplace_back(consumer, i, fp1, ks1, k, m, xxhash_flag, min_len);
   	}

    for (auto& t : threads) {
        t.join();
    }

	compactSeqIds(num_threads);
	outputAllHashes(m, num_threads);
    gzclose(fp1);
    kseq_destroy(ks1);
}

int main(int argc, char* argv[])
{

	CLI::App app{"yclust v.0.0.1, extremely fast and scalable protein clustering"};
	int num_threads = 1;
	int min_len = 50;
	int k = 8;
	int m = 15;
	int r = 1;
	float similarity = 0.9;
	float cluster_thd = 0.9;
	string filename = "";
	string res_file = "";
	string fa_output_name = "";

	auto option_threads = app.add_option("-t, --threads", num_threads,  "set the thread number, default 1 thread");
	auto option_min_len = app.add_option("--min-length", min_len, "set the filter minimum length (minLen), protein length less than minLen will be ignore, default 50");
	auto option_min_similarity = app.add_option("-s, --min-similarity", similarity, "set the minimum similarity for clustering, default 0.9");
	auto option_kmer_size = app.add_option("-k, --kmer-size", k, "set the kmer size, default 8");
	auto option_m_size = app.add_option("-m, --m-size", m, "set the number of hash functions will be used, default 15");
	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	auto option_r = app.add_option("-r, --r-size", r, "set the number of block");
	auto option_cluster_thd = app.add_option("--c-thd, --cluser-threshold", cluster_thd, "cluster threshold in cdhit");

	option_input->required();

	bool xxhash_flag = false;
	auto option_xxhash = app.add_flag("-x, --xxhash", xxhash_flag, "Default hash is aahash, if this flag is enabled, use xxhash");

	bool block_on = false;
	auto option_block = app.add_flag("-b, --block-on", block_on, "If this flat is enabled, sort in block mode. eg. m cols unites m / r times");
	option_block->needs("-r");

	bool cluster_on = false;
	auto option_cluster = app.add_flag("-c, --cluster", cluster_on, "If this flat is enabled, clustering during the grouping");

	bool rep_group_on = false;
	auto option_rep_group = app.add_flag("--rep-g, --use-representatives-grouping", rep_group_on, "If this flag is enabled, only use representative sequences during the grouping");

	bool rep_cluster_on = false;
	auto option_rep_cluster = app.add_flag("--rep-c, --use-representatives-clustering", rep_cluster_on, "If this flag is enabled, only use representative sequences during the grouping and clustering");
	option_rep_cluster->needs("--rep-g");

	bool reorder_off = false;
	auto option_reorder = app.add_flag("--reorder-off, --reorderSequences", reorder_off, "If this flat is enabled, reorder sketch vector by the sequence index read order");

	bool threadPool_off = false;
	auto option_thead_pool = app.add_flag("--threadpool-off, --open-thread-pool", threadPool_off, "If this flat is enabled, disable thread pool to allocate threads");

	bool final_cluster_off = false;
	auto option_final_cluster = app.add_flag("-f, --disable-final-cluster", final_cluster_off, "If this flat is disabled, cluster all sequences at the last grouping");

	bool top_on = false;
	auto option_output = app.add_option("-o, --output", res_file, "If this flat is enabled, the top 10 largest clusters suquences content will be written into the output files");

	CLI11_PARSE(app, argc, argv);

	if(num_threads < 1)
	{
		cerr << "Invalid thread number: " << num_threads << endl;
		return 1;
	}

	cerr << "==========Paramters==========" << endl;
	cerr << "Threads: " << num_threads << endl;
	cerr << "K: " << k << endl;
	cerr << "M: " << m << endl;
	cerr << "R: " << r << endl;
	cerr << "Min_len: " << min_len << endl;
	cerr << "Similarity Threshold:" << similarity << endl;
	cerr << "Input: " << filename << endl;
	if(rep_group_on) {
		if(rep_cluster_on) {
			//cerr << "only use representatives in clustering and grouping" << endl;
			cerr << "分组和聚类阶段均只考虑代表序列" << endl;
		}else{
			//cerr << "only use representatives in grouping" << endl;
			cerr << "分组阶段只考虑代表序列" << endl;
		}
	}
	if(reorder_off) {
		cerr << "关闭了重新排序" << endl;
	}else{
		cerr << "对所有序列按照输入顺序排序" << endl;
	}
	if(!threadPool_off) {
		cerr << "use thread pool when clustering" << endl;
	}
	if(!final_cluster_off) {
		cerr << "Clustering all groups in the last round" << endl;
	}
	cerr << "==========End Paramters==========" << endl;


	cerr << "Start Building sketches!" << endl;
	auto generation_start = chrono::high_resolution_clock::now();
	
	buildSketches(filename, num_threads, k, m, xxhash_flag, min_len);

	auto generation_end = chrono::high_resolution_clock::now();
	auto generation_duration = chrono::duration_cast<chrono::seconds>(generation_end - generation_start).count();
	return 0;

	cerr << "Sketching time: " << generation_duration << endl;
	cerr << "Total number of seqs: " << num_seqs.load() << endl;

	//cerr << "Generated Map Size: " << fa_map.size() << endl;


	//grouping
	cerr << "Start grouping!" << endl;
	GroupStream gs(num_seqs.load(), m, r, 1);
	gs.setIDs(seq_ids);
	gs.setNumThreads(num_threads);
	if(rep_group_on) {
		if(rep_cluster_on) {
			gs.setRepGroupAndClusterOn();
		}else{
			gs.setRepGroupOn();
		}
	}
	if(!threadPool_off) {
		gs.setThreadPool();
	}
	if(cluster_on) {
		gs.setClusterOn();
		gs.setClusterThd(cluster_thd);
	}
	if(res_file != "") {
		gs.setOutput(res_file);
		top_on = true;
	}
	if(!final_cluster_off) {
		gs.setFinalClusterOn();
	}
	if(block_on) 
		gs.setSlideOff();
	unordered_map<int, vector<int>> group_map;
	gs.Group(hashes, group_map);

	//输出每个seq和他的root
//	cout.rdbuf(origin_cout);
// 打印代表序列
//	int name_pos = filename.find('.');
//	string rep_name = filename.substr(0, name_pos) + ".rep";
//	ofstream rep_file(rep_name);
//	vector<int> rep_ids;

	string folder_name = "test-output";
//	priority_queue<int, std::vector<int>, std::greater<int>> minHeap;
	priority_queue<pair<int,int>, std::vector<pair<int, int>>, compare> minHeap;
	int max_group_Size = 0;
	for(const auto& pair : group_map) {
/**
 * 把含有多个序列的类输出用cdhit聚类
		if(pair.second.size() > 1) {
			string file_name = folder_name + "/" + to_string(pair.first) + ".fa";
			ofstream out_file(file_name);
			for(int id : pair.second) {
				out_file << ">" << id << endl;
				out_file << fa_map[id]  << endl;
			}
			out_file.close();
		}
**/
		// 输出rep
		//rep_ids.emplace_back(pair.first);
        minHeap.push({pair.second.size(), pair.first});
        if (minHeap.size() > 2) {
            minHeap.pop(); // 保持堆的大小为 2
        }
//		打印id-rootid
//		for(const auto& node : pair.second)
//			cout << node << " " << pair.first << endl;
	}

//	输出代表序列
//	std::copy(rep_ids.begin(), rep_ids.end(), std::ostream_iterator<int>(rep_file, "\n"));

	while(!minHeap.empty()){
		if(top_on) {
			string file_name = "output/" + to_string(minHeap.top().second) + ".fa";
			cerr << "results output to: " << file_name << endl;
			ofstream ofile(file_name);
			int rootid = minHeap.top().second;
			for(int i : group_map[rootid]) {
				ofile << ">" << names[i] << endl;
				ofile << fa_map[i] << endl;
			}
			ofile.close();
		}
		cerr << minHeap.top().first << endl;
		minHeap.pop();
	}

	cerr << "Total Group Nums: " << group_map.size() << endl;
	return 0;
	
}
