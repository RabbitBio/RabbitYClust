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
int cnt_seqs=0;

vector<vector<uint64_t>> hashes;
vector<uint64_t> seq_ids;
vector<string> names;
// yy add for cluster
//bool cluster_on = false;
unordered_map<uint64_t, string> fa_map;
int64_t pos = 0;

vector<vector<vector<uint64_t>>> all_hashes;
vector<vector<uint64_t>> thr_seq_ids;
struct compare {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

void read_fa(string filename, int min_len)
{
	
	gzFile fp1;
	kseq_t *ks1;

	fp1 = gzopen(filename.c_str(),"r");
	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return ;
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
}

int main(int argc, char* argv[])
{

	CLI::App app{"yclust v.0.0.1, extremely fast and scalable protein clustering"};
	int num_threads = 1;
	int min_len = 1;
	int k = 8;
	int m = 15;
	int r = 1;
	float similarity = 0.9;
	float cluster_thd = 0.9;
	string filename = "";
	string res_file = "";
	string fa_output_name = "";
	string sketch_file_name = "";

	auto option_threads = app.add_option("-t, --threads", num_threads,  "set the thread number, default 1 thread");
	auto option_min_len = app.add_option("--min-length", min_len, "set the filter minimum length (minLen), protein length less than minLen will be ignore, default 50");
	auto option_min_similarity = app.add_option("-s, --min-similarity", similarity, "set the minimum similarity for clustering, default 0.9");
	auto option_kmer_size = app.add_option("-k, --kmer-size", k, "set the kmer size, default 8");
	auto option_m_size = app.add_option("-m, --m-size", m, "set the number of hash functions will be used, default 15");
	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	auto option_r = app.add_option("-r, --r-size", r, "set the number of block");
	auto option_cluster_thd = app.add_option("--c-thd, --cluser-threshold", cluster_thd, "cluster threshold in cdhit");
	auto option_sketch = app.add_option("-S, --sketch-file", sketch_file_name, "sketch file name, needed");

	option_input->required();
	option_sketch->required();

	bool xxhash_flag = false;
	auto option_xxhash = app.add_flag("-x, --xxhash", xxhash_flag, "Default hash is aahash, if this flag is enabled, use xxhash");

	bool block_on = false;
	auto option_block = app.add_flag("-b, --block-on", block_on, "If this flat is enabled, sort in block mode. eg. m cols unites m / r times");
	option_block->needs("-r");

	bool cluster_on = false;
	auto option_cluster = app.add_flag("-c, --cluster", cluster_on, "If this flat is enabled, clustering during the grouping");

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


	cerr << "Start reading FA files!" << endl;
	auto read_fa_start = chrono::high_resolution_clock::now();
	read_fa(filename, min_len);
	auto read_fa_end = chrono::high_resolution_clock::now();
	auto read_fa_duration = chrono::duration_cast<chrono::seconds>(read_fa_end - read_fa_start).count();

	cerr << "Reading time: " << read_fa_duration << endl;
	cerr << "Total number of seqs: " << cnt_seqs << endl;

	cerr << "Start grouping!" << endl;
	GroupStream gs(cnt_seqs, m, r, 1);
	gs.setIDs(seq_ids);
	gs.setNumThreads(num_threads);
	gs.set_sketch_filename(sketch_file_name);

	if(!threadPool_off) {
		gs.setThreadPool();
	}
	if(cluster_on) {
		gs.setClusterOn();
		gs.setClusterThd(cluster_thd);
	}
	if(res_file != "") {
		gs.setOutput(res_file);
		//top_on = true;
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

//	string folder_name = "test-output";
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
