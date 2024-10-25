#include "KHFMinHash.h"
#include "kseq.h"
#include "CLI11.hpp"
#include <zlib.h>

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>

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
	float similarity = 0.9;
	string filename = "";
	string res_file = "";

	auto option_threads = app.add_option("-t, --threads", threads,  "set the thread number, default 1 thread");
	auto option_min_len = app.add_option("--min-length", min_len, "set the filter minimum length (minLen), protein length less than minLen will be ignore, default 50");
	auto option_min_similarity = app.add_option("-s, --min-similarity", similarity, "set the minimum similarity for clustering, default 0.9");
	auto option_kmer_size = app.add_option("-k, --kmer-size", k, "set the kmer size, default 8");
	auto option_m_size = app.add_option("-m, --m-size", m, "set the number of hash functions will be used, default 15");
	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	auto option_output = app.add_option("-o, --output", res_file, "output file, 64bit binary hashes");

	option_input->required();
	option_output->required();

	CLI11_PARSE(app, argc, argv);

	if(threads < 1)
	{
		cerr << "Invalid thread number: " << threads << endl;
		return 1;
	}
	cerr << "==========Paramters==========" << endl;
	cerr << "Threads: " << threads << endl;
	cerr << "K: " << k << endl;
	cerr << "M: " << m << endl;
	cerr << "Min_len: " << min_len << endl;
	cerr << "Similarity Threshold:" << similarity << endl;
	cerr << "Input: " << filename << endl;
	cerr << "Ouput: " << res_file << endl;
	cerr << "==========End Paramters==========" << endl;
#else
	if(argc != 3)
	{
		cerr << "Usage: ./yclust input.fa output.hash" << endl;
		cerr << "Note: output.hash is in binary format for uint64_t hash values!" << endl;
		cerr << "Note: the names of sequences are printed to stdout!" << endl;
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

	ks1 = kseq_init(fp1);

	int buffer_size = 1<<20;//1MB
	int buffer_pos = 0;
	vector<char> buffer(buffer_size);

	int count = 0;	

	cerr << "buffer_size: " << buffer_size << endl;
	while(1)
	{
		int length = kseq_read(ks1);

		if(length < 0) break;
		if (ks1->seq.l <= min_len) continue;

		char * seq1 = ks1->seq.s;
		cout << ">" << ks1->name.s << " ";
		//cout << ks1->comment.s << " " << ks1->seq.l << endl;
		cout << ks1->seq.l << endl;
		cout << ks1->seq.s << endl;
		KHFMinHash mh = KHFMinHash();			
#ifdef CLI11
		mh.setK(k);
		mh.setM(m);
#endif
		mh.buildSketch(seq1);
	
		auto & sketch = mh.getSektch();	

		//for(int i = 0; i < sketch.hashes.size(); i++)
		//	cerr << sketch.hashes[i] << endl;	

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

		count++;
	}
	
	if(buffer_pos > 0)
		ofile.write(buffer.data(), buffer_pos);

    gzclose(fp1);
    kseq_destroy(ks1);

	ofile.close();
	cerr << "number of seqs: " << count << endl;

	return 0;
	
}
