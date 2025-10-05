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
#include "util.h"

#define BUFFER_SIZE (1<<20 * sizeof(char))

using namespace std;
using namespace Sketch;

extern unordered_map<uint64_t, string> fa_map;

struct compare {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

int main(int argc, char* argv[])
{

	CLI::App app{"yclust v.0.0.1, extremely fast and scalable protein clustering"};

	int num_threads = 1;
	int min_len = 11;
	int k = 8;
	int m = 15;
	int r = 1;
	int cnt_seqs = 0;
	float similarity = 0.9;
	float cluster_thd = 0.9;
	bool xxhash_flag = false;
	bool cluster_on = false;
	bool final_cluster_off = false;
	bool threadPool_off = false;
	string input_filename = "";
	string res_filename = "";
	string sketch_filename = "";

	auto option_threads = app.add_option("-t, --threads", num_threads,  "set the thread number, default 1 thread");
	auto option_min_len = app.add_option("--min-len", min_len, "set the filter minimum length (minLen), protein length less than minLen will be ignore, default 50");
	auto option_min_similarity = app.add_option("-s, --min-similarity", similarity, "set the minimum similarity for clustering, default 0.9");
	auto option_kmer_size = app.add_option("-k, --kmer-size", k, "set the kmer size, default 8");
	auto option_m_size = app.add_option("-m, --m-size", m, "set the number of hash functions will be used, default 15");
	auto option_r = app.add_option("-r, --r-size", r, "set the number of block");
	auto option_cluster_thd = app.add_option("--c-thd, --cluser-threshold", cluster_thd, "cluster threshold in cdhit");
	auto option_input = app.add_option("-i, --input", input_filename, "input file name, fasta or gziped fasta formats");
	auto option_sketch = app.add_option("-S, --skech-file", sketch_filename, "Sketch file name");
	option_input->required();
	option_sketch->required();

	auto subA = app.add_subcommand("sketch", "Building sketches");
	subA->add_flag("-x", xxhash_flag, "enable this flag to use xxhash in building sketches, default aaHash");
	

	auto subB = app.add_subcommand("clust", "clustering by minhash");
	subB->add_flag("-c", cluster_on, "enable clustering to avoid super-huge group");
	subB->add_flag("-f", final_cluster_off, "turn off the final clustering")->needs("-c");
	subB->add_option("-o", res_filename, "output clusters, like: seq_id : rep_seq_id");

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
	cerr << "Min_len: " << min_len << endl;
	cerr << "Input: " << input_filename << endl;
	cerr << "Output: " << sketch_filename << endl;
	cerr << "==========End Paramters==========" << endl;

	if(*subA) {
		cerr << "Start Building sketches!" << endl;
		auto generation_start = chrono::high_resolution_clock::now();
		cnt_seqs = buildSketches(input_filename, sketch_filename, num_threads, k, m, xxhash_flag, min_len, true);

		//test hashes
		//int n = num_seqs.load();
		//vector<uint64_t> read_hashes(n);
		//ifstream ifs(sketch_filename, ios::binary);
		//ifs.seekg(2*n*sizeof(uint64_t), ios::beg);
		//ifs.read(reinterpret_cast<char*>(read_hashes.data()), n * sizeof(uint64_t));
		//for(int i=0; i<10; i++)
		//{
		//	std::cout << read_hashes[i] << endl;
		//	if(read_hashes[i] != hashes[2][i]){
		//		cout << " wrong!" << endl;
		//	}
		//}

		auto generation_end = chrono::high_resolution_clock::now();
		auto generation_duration = chrono::duration_cast<chrono::seconds>(generation_end - generation_start).count();
		cerr << "Sketching time: " << generation_duration << endl;
		cerr << "Total number of seqs: " << cnt_seqs << endl;
	}
	if (*subB && !*subA && option_sketch->count() == 0) {
		throw CLI::RequiredError("--sketchfile-input");
		return 1;
	}
	if (*subB) {
		if(!threadPool_off) {
			cerr << "use thread pool when clustering" << endl;
		}
		if(!final_cluster_off) {
			cerr << "Clustering all groups in the last round" << endl;
		}
		cerr << "Start reading FA files!" << endl;
		auto read_fa_start = chrono::high_resolution_clock::now();
		cnt_seqs = read_fa(input_filename, min_len);
		auto read_fa_end = chrono::high_resolution_clock::now();
		auto read_fa_duration = chrono::duration_cast<chrono::seconds>(read_fa_end - read_fa_start).count();

		cerr << "Reading time: " << read_fa_duration << endl;
		cerr << "Total number of seqs: " << cnt_seqs << endl;

		if(cnt_seqs == 0) {
			cerr << "no sequences in grouping!" << endl;
			return 1;
		}

		cerr << "Start grouping!" << endl;
		GroupStream gs(cnt_seqs, m, r, 1);
		// no use seq_ids
		//gs.setIDs(seq_ids);
		gs.setNumThreads(num_threads);

		if(!threadPool_off) {
			gs.setThreadPool();
		}
		if(cluster_on) {
			gs.setClusterOn();
			gs.setClusterThd(cluster_thd);
		}
		if(res_filename != "") {
			gs.setOutput(res_filename);
		}
		if(!final_cluster_off) {
			gs.setFinalClusterOn();
		}
		unordered_map<int, vector<int>> group_map;
		gs.Group(sketch_filename, group_map);
		//if(!*subA) {
		//	gs.Group(hashes, group_map);
		//}else {
		//	gs.Group(sketch_filename, group_map);
		//}

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

		//while(!minHeap.empty()){
		//	if(top_on) {
		//		string file_name = "output/" + to_string(minHeap.top().second) + ".fa";
		//		cerr << "results output to: " << file_name << endl;
		//		ofstream ofile(file_name);
		//		int rootid = minHeap.top().second;
		//		for(int i : group_map[rootid]) {
		//			ofile << ">" << names[i] << endl;
		//			ofile << fa_map[i] << endl;
		//		}
		//		ofile.close();
		//	}
		//	cerr << minHeap.top().first << endl;
		//	minHeap.pop();
		//}

		cerr << "Total Group Nums: " << group_map.size() << endl;
	}
	
	return 0;
}
