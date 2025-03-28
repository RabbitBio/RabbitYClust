#include "KHFMinHash.h"
#include "kseq.h"
#include "CLI11.hpp"
#include <zlib.h>

#include <iostream>
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
// yy add for cluster
bool cluster_on = false;
bool second_group = false;
unordered_map<uint64_t, uint64_t> fai_map;
unordered_map<uint64_t, string> fa_map;
int64_t pos = 0;

struct compare {
	bool operator()(const pair<int, int> &a, const pair<int, int> &b) {
		return a.first > b.first;
	}
};

void consumer(int tid, gzFile fp, kseq_t* ks, int k, int m, bool xxhash_flag, int min_len) {
    while (true) {
        std::string sequence;
        int seq_id;

		{
				std::lock_guard<std::mutex> lock(mtx1);
				int length = kseq_read(ks);
				if (length < 0) break;
				if (length < min_len) continue;
				sequence = ks->seq.s;//direct copy?
				
//				cout << ks->name.s << " " << seq_id << endl;
		}

		KHFMinHash mh;
		mh.setK(k);
		mh.setM(m);

		if (xxhash_flag)
				mh.buildSketch(sequence.c_str());
		else
				mh.buildSketchByNoSeedAAHash(sequence.c_str());

		auto& sketch = mh.getSektch();
		
	
		{
				std::lock_guard<std::mutex> lock(mtx2);
				seq_id = num_seqs.fetch_add(1);
				hashes.emplace_back(sketch.hashes);
				seq_ids.emplace_back(seq_id);
		}
	}
}

void consumer_cluster(int tid, gzFile fp, kseq_t* ks, int k, int m, bool xxhash_flag, int min_len) {
    while (true) {
        std::string sequence;
        int seq_id;
		{
				std::lock_guard<std::mutex> lock(mtx1);
				int length = kseq_read(ks);
				if (length < 0) break;
				if (length < min_len) continue;
				sequence = ks->seq.s;//direct copy?
				
	//			cout << ks->name.s << " " << seq_id << endl;
		}

		// FIXME: store char* for simply using cdhit to cluster
		const char* cstr = sequence.c_str();
//		char* buffer = (char*)malloc(sequence.size()+1);
//		char* buffer = new char[sequence.size() + 1];
//		strcpy(buffer, cstr);
		KHFMinHash mh;
		mh.setK(k);
		mh.setM(m);

		if (xxhash_flag)
				mh.buildSketch(sequence.c_str());
		else
				mh.buildSketchByNoSeedAAHash(sequence.c_str());

		auto& sketch = mh.getSektch();
		
	
		{
				std::lock_guard<std::mutex> lock(mtx2);
				seq_id = num_seqs.fetch_add(1);
				hashes.emplace_back(sketch.hashes);
				seq_ids.emplace_back(seq_id);
				fa_map.emplace(seq_id, sequence);
		}
	}
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

	option_input->required();

	bool xxhash_flag = false;
	auto option_xxhash = app.add_flag("-x, --xxhash", xxhash_flag, "Default hash is aahash, if this flag is enabled, use xxhash");

	bool block_on = false;
	auto option_block = app.add_flag("-b, --block-on", block_on, "If this flat is enabled, sort in block mode. eg. m cols unites m / r times");
	option_block->needs("-r");

	cluster_on = false;
	auto option_cluster = app.add_flag("-c, --cluster", cluster_on, "If this flat is enabled, clustering during the grouping");

	second_group = false;
	auto option_second = app.add_flag("-S, --second", second_group, "If this flat is enabled, second grouping during the grouping");


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
	cerr << "==========End Paramters==========" << endl;


	gzFile fp1;
	kseq_t *ks1;

	fp1 = gzopen(filename.c_str(),"r");

	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return 0;
	}

	
	ks1 = kseq_init(fp1);

//  FIXME:临时的验证方式 为了把其他聚类软件的结果和seqid对应起来
//	ostringstream seq_id_name;
//	seq_id_name << "k" << k << "m" << m << "seq-id.txt";
//	ofstream seq_id(seq_id_name.str());
	ofstream seq_id("sequence-id.txt");
	streambuf* origin_cout = cout.rdbuf();
	cout.rdbuf(seq_id.rdbuf());


	cerr << "Start Building sketches!" << endl;
	auto generation_start = chrono::high_resolution_clock::now();
	
    std::vector<std::thread> threads;
	if(cluster_on) {
    	for (int i = 0; i < num_threads; ++i) {
    	    threads.emplace_back(consumer_cluster, i, fp1, ks1, k, m, xxhash_flag, min_len);
    	}
	}else {

    	for (int i = 0; i < num_threads; ++i) {
    	    threads.emplace_back(consumer, i, fp1, ks1, k, m, xxhash_flag, min_len);
    	}
	}
    for (auto& t : threads) {
        t.join();
    }
    

	auto generation_end = chrono::high_resolution_clock::now();
	auto generation_duration = chrono::duration_cast<chrono::seconds>(generation_end - generation_start).count();

	cerr << "Sketching time: " << generation_duration << endl;
	cerr << "Total number of seqs: " << num_seqs.load() << endl;

	cerr << "Generated Map Size: " << fa_map.size() << endl;

    gzclose(fp1);
    kseq_destroy(ks1);

	//grouping
	cerr << "Start grouping!" << endl;
	GroupStream gs(num_seqs.load(), m, r, 1);
	gs.setIDs(seq_ids);
	gs.setNumThreads(num_threads);
	if(cluster_on) {
		gs.setClusterOn();
	}
	if(block_on) 
		gs.setSlideOff();
	if (second_group) {
		gs.setSecondGroup();
	}
	unordered_map<int, vector<int>> group_map;
	gs.Group(hashes, group_map);

	//输出每个seq和他的root
	cout.rdbuf(origin_cout);
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
        if (minHeap.size() > 10) {
            minHeap.pop(); // 保持堆的大小为 10
        }
//		打印id-rootid
//		for(const auto& node : pair.second)
//			cout << node << " " << pair.first << endl;
	}

//	输出代表序列
//	std::copy(rep_ids.begin(), rep_ids.end(), std::ostream_iterator<int>(rep_file, "\n"));

	while(!minHeap.empty()){
//		string file_name = "nr/" + to_string(minHeap.top().second) + ".fa";
//		cout << "results output to: " << file_name << endl;
//		ofstream ofile(file_name);
//		int rootid = minHeap.top().second;
//		for(int i : group_map[rootid]) {
//			ofile << ">" << i << endl;
//			ofile << fa_map[i] << endl;
//		}
//		ofile.close();
		cerr << minHeap.top().first << endl;
		minHeap.pop();
	}

	cerr << "Total Group Nums: " << group_map.size() << endl;
	return 0;
	
}
