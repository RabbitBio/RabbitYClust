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

#include "../kernels/GroupStream.h"

#define CLI11

using namespace std;
using namespace Sketch;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[])
{
	gzFile fp1;
	kseq_t *ks1;

#ifdef CLI11
	CLI::App app{"yclust v.0.0.1, extremely fast and scalable protein clustering"};
	int threads = 1;
	int min_len = 50;
	int k = 8;
	int m = 15;
	int r = 2;
	float similarity = 0.9;
	string filename = "";
	string res_file = "";

	auto option_threads = app.add_option("-t, --threads", threads,  "set the thread number, default 1 thread");
	auto option_min_len = app.add_option("--min-length", min_len, "set the filter minimum length (minLen), protein length less than minLen will be ignore, default 50");
	auto option_min_similarity = app.add_option("-s, --min-similarity", similarity, "set the minimum similarity for clustering, default 0.9");
	auto option_kmer_size = app.add_option("-k, --kmer-size", k, "set the kmer size, default 8");
	auto option_m_size = app.add_option("-m, --m-size", m, "set the number of hash functions will be used, default 15");
	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	auto option_r = app.add_option("-r, --r-size", r, "set the number of block");
#ifdef OUTPUT
	auto option_output = app.add_option("-o, --output", res_file, "output file, 64bit binary hashes");
	option_output->required();
#endif
#ifdef FAI
	string fa_file_name = "";
	bool comment_ignore = false;
	auto option_output_fa = app.add_option("-f, --new-fa", fa_file_name, "output new formated fasta file");
	auto option_comment = app.add_flag("-c, --comment", comment_ignore, "If this flat is enabled, ignore the comment.");
	option_output_fa->required();
#endif
	option_input->required();

	bool xxhash_flag = false;
	auto option_xxhash = app.add_flag("-x, --xxhash", xxhash_flag, "Default hash is aahash, if this flag is enabled, use xxhash");

	CLI11_PARSE(app, argc, argv);

	if(threads < 1)
	{
		cerr << "Invalid thread number: " << threads << endl;
		return 1;
	}
    const unsigned h = 1 ;  // 哈希数量（每个种子）
    const unsigned hash_num_per_seed = m / h;  // 每个种子生成的哈希数
    
    // 随机生成种子
    std::vector<std::string> seed_strings;
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_int_distribution<> dis(0, 3);
//    
//    for (int i = 0; i < h; ++i) {
//        std::string seed = "";
//        for (int j = 0; j < k; ++j) {
//            seed += std::to_string(dis(gen));  // 随机生成 0-3 的数字
//            // seed+=std::to_string(i);
//        }
//        seed_strings.push_back(seed);
//    }
//    gyj新写的代码，不用随机种子
//	std::string seed = "11111111";
//	seed_strings.push_back(seed);


	cerr << "==========Paramters==========" << endl;
	cerr << "Threads: " << threads << endl;
	cerr << "K: " << k << endl;
	cerr << "M: " << m << endl;
	cerr << "R: " << r << endl;
	cerr << "Min_len: " << min_len << endl;
	cerr << "Similarity Threshold:" << similarity << endl;
	cerr << "Input: " << filename << endl;
	#ifdef OUTPUT
	cerr << "Ouput: " << res_file << endl;
	#endif
	#ifdef FAI
	cerr << "New Fasta File: " << fa_file_name << endl;
	#endif
	cerr << "seed: " << h << endl;
	cerr << "hashNum per seed: " << hash_num_per_seed << endl;
	cerr << "==========End Paramters==========" << endl;
#elif defined(OUTPUT) && defined(FAI)
	if(argc != 4)
	{
		cerr << "Usage: ./yclust input.fa output.hash new-fasta-name" << endl;
	}
#elif defined(OUTPUT)
	if(argc != 3)
	{
		cerr << "Usage: ./yclust input.fa output.hash" << endl;
		cerr << "Note: output.hash is in binary format for uint64_t hash values!" << endl;
		cerr << "Note: the names of sequences are printed to stdout!" << endl;
		return 1;
	}
#elif defined(FAI)
	if(argc != 3)
	{
		cerr << "Usage: ./yclust input.fa new-fasta-name" << endl;
		return 1;
	}
#else
	if(argc != 2)
	{
		cerr << "Usage: ./yclust input.fa" << endl;
		return 1;
	}
#endif


#ifdef CLI11
	fp1 = gzopen(filename.c_str(),"r");

	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return 0;
	}
#else
	fp1 = gzopen(argv[1],"r");

	if(NULL == fp1){
		cerr << "Fail to open file: " << argv[1] << endl;
		return 0;
	}
#endif

#ifdef OUTPUT
#ifdef CLI11
	ofstream ofile(res_file, ios::binary);
	if( !ofile.is_open())
	{
		cerr << "Failed to open file: " << res_file << endl;
		return 1;
	}
#else
	ofstream ofile(argv[2], ios::binary);
	if( !ofile.is_open())
	{
		cerr << "Failed to open file: " << argv[2] << endl;
		return 1;
	}
#endif
	int buffer_size = 1<<20;//1MB
	int buffer_pos = 0;
	vector<char> buffer(buffer_size);

	cerr << "buffer_size: " << buffer_size << endl;
#endif

#ifdef FAI
#ifdef CLI11
	ofstream fa_output(fa_file_name, ios::binary);
#else
	ofstream fa_output(argv[3], ios::binary);
#endif
	if( !fa_output.is_open())
	{
		cerr << "Failed to open file: " << fa_file_name << endl;
		return 1;
	}

	int64_t pos = 0;
	int64_t number_seqs = 0;
	vector<uint64_t> fai;

    char *record = (char *)malloc(1<<20 * sizeof(char));  //FIXME: 1MB buffer only works for protein sequences
#endif
    
	ks1 = kseq_init(fp1);

	int count = 0;	
	vector<vector<uint64_t>> hashes;

// 输出sequence-id sequence-group groupid-sequences
//	ofstream seq_id("sequence-id.txt");
//	ostringstream group_res_name;
//	group_res_name << "k" << k << "m" << m << "group.res";
//	ofstream seqfile(group_res_name.str());
//	streambuf* origin_cout = cout.rdbuf();
//	cout.rdbuf(seq_id.rdbuf());

	auto generation_start = chrono::high_resolution_clock::now();
	while(1)
	{
		int length = kseq_read(ks1);

		if(length < 0) break;
		if (ks1->seq.l <= min_len) continue;

		char * seq1 = ks1->seq.s;
//		cout << ks1->name.s << " " << count << endl;
		////cout << ks1->comment.s << " " << ks1->seq.l << endl;
		//cout << ks1->seq.l << endl;
		//cout << ks1->seq.s << endl;
		KHFMinHash mh = KHFMinHash();			
#ifdef CLI11
		mh.setK(k);
		mh.setM(m);
#endif
		if(xxhash_flag) mh.buildSketch(seq1);
		else mh.buildSketchByNoSeedAAHash(seq1);
//		mh.buildSketch(seq1, seed_strings, h, hash_num_per_seed);
//		mh.buildSketch(seq1);

	
		auto & sketch = mh.getSektch();	
		// 构建vector<vector> hashes
		hashes.emplace_back(sketch.hashes);

		//for(int i = 0; i < sketch.hashes.size(); i++)
		//	cerr << sketch.hashes[i] << endl;	


#ifdef FAI
		fai.emplace_back(pos);

		int record_len = 0;
		int str_pos = 0;
		record[str_pos] = '>';
		str_pos++;
		for(int i = 0; i < ks1->name.l; i++) record[i+str_pos] = ks1->name.s[i];
		str_pos += ks1->name.l;
		if(!comment_ignore)
		{
			record[str_pos] = ' ';
			str_pos++;
			for(int i = 0; i < ks1->comment.l; i++) record[i+str_pos] = ks1->comment.s[i];
			str_pos += ks1->comment.l;
		}

		record[str_pos] = '\n';
		str_pos++;
		for(int i = 0; i < ks1->seq.l; i++) record[i+str_pos] = ks1->seq.s[i];
		str_pos += ks1->seq.l;
		record[str_pos] = '\n';
		str_pos++;

		//cerr << str_pos << endl;
		//cerr << string(record, ks1->name.l + ks1->comment.l + ks1->seq.l + 4);
   		if(!comment_ignore) record_len =  ks1->name.l + ks1->comment.l + ks1->seq.l + 4;
		else record_len = ks1->name.l + ks1->seq.l + 3;

		fa_output.write(record, record_len);

		pos += record_len; //note that with extra '>' ' ' and two '\n'
#endif

#ifdef OUTPUT
		const char* vec_data = reinterpret_cast<const char*>(sketch.hashes.data());
		int vec_bytes_left = sketch.hashes.size() * sizeof(uint64_t);
		int vec_pos = 0;
		while (vec_bytes_left > 0)
		{
			size_t bytes_to_copy = std::min(buffer_size - buffer_pos, vec_bytes_left);
			std::memcpy(buffer.data() + buffer_pos, vec_data + vec_pos, bytes_to_copy);
			buffer_pos += bytes_to_copy;
			vec_pos += bytes_to_copy;
			vec_bytes_left -= bytes_to_copy;
 
			if(buffer_pos == buffer_size)
			{
				cerr << "write to file" << endl;
				ofile.write(buffer.data(), buffer_size);
				buffer_pos = 0;
			}
		}
#endif
		count++;
	}
	auto generation_end = chrono::high_resolution_clock::now();
	auto generation_duration = chrono::duration_cast<chrono::seconds>(generation_end - generation_start).count();
	cerr << "generation time: " << generation_duration << endl;
	
#ifdef OUTPUT
	if(buffer_pos > 0)
		ofile.write(buffer.data(), buffer_pos);
	ofile.close();
#endif

#ifdef FAI
//	把fai写到一个文件里
//	o_file.write(reinterpret_cast<const char*>(fai.data()), fai.size() * sizeof(uint64_t));
	fai.emplace_back(pos);
	fa_output.close();
#endif
	cerr << "number of seqs: " << count << endl;

    gzclose(fp1);
    kseq_destroy(ks1);

	GroupStream gs(count, m, r);
	unordered_map<int, vector<int>> group_map;
	gs.Group(hashes, group_map);
#ifdef FAI
	int group_id = 0;
	ifstream fa_input(fa_file_name, ios::binary);
	for(const auto& pair : group_map) {
		for(const auto& seq : pair.second) {
			uint64_t start_pos = fai[seq];
			uint64_t end_pos = fai[seq+1];
			int length = (int)(end_pos - start_pos);
			fa_input.seekg(start_pos, ios::beg);
			string sequence;
			sequence.resize(length);
			fa_input.read(&sequence[0], length);
		}
	}
	fa_input.close();
#endif
// 输出每个seq和他的root
	int max_group_Size = 0;
	for(const auto& pair : group_map) {
		max_group_Size = max_group_Size > pair.second.size() ? max_group_Size : pair.second.size();
		for(const auto& node : pair.second)
			cout << node << " " << pair.first << endl;
	}

	cerr << "MAX-Group-Size: " << max_group_Size << endl;
	cerr << "Total Group Nums: " << group_map.size() << endl;
	return 0;
	
}
