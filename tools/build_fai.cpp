#include "kseq.h"
#include "CLI11.hpp"
#include <zlib.h>

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[])
{
	gzFile fp1;
	kseq_t *ks1;

	CLI::App app{"fai v.0.0.1, build fasta index of (gziped) protein files"};
	string filename = "";
	string res_file = "";
	string fa_file_name = "";
	bool comment_flag = false;
		
	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	auto option_output = app.add_option("-o, --output", res_file, "output fasta index file");
	auto option_output_fa = app.add_option("-f, --new-fa", fa_file_name, "output new formated fasta file");
	auto option_comment = app.add_flag("-c, --comment", comment_flag, "If this flat is enabled, ignore the comment.");

	option_input->required();
	option_output->required();
	option_output_fa->required();

	CLI11_PARSE(app, argc, argv);

	fp1 = gzopen(filename.c_str(),"r");

	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return 0;
	}

	ofstream ofile(res_file, ios::binary);
	if( !ofile.is_open())
	{
		cerr << "Failed to open file: " << res_file << endl;
		return 1;
	}

	ofstream fa_file(fa_file_name, ios::binary);
	if( !ofile.is_open())
	{
		cerr << "Failed to open file: " << fa_file_name << endl;
		return 1;
	}

	ks1 = kseq_init(fp1);

	int64_t pos = 0;
	int64_t number_seqs = 0;
	vector<uint64_t> fai;

    char *record = (char *)malloc(1<<20 * sizeof(char));  //FIXME: 1MB buffer only works for protein sequences

	while(1)
	{
		int length = kseq_read(ks1);
		int record_len = 0;

		if(length < 0) break;
		//stringstream record;
		//cerr << pos << endl;
		fai.emplace_back(pos);

		//cerr << "seq.l: " << ks1->seq.l << endl;
		//cerr << "name.l: " << ks1->name.l << endl;
		//cerr << "comment.l: " << ks1->comment.l << endl;
		//record << '>' << ks1->name.s << ' ' << ks1->comment.s << endl;
		//record << ks1->seq.s << endl;
		//cout << record.str();
		//fa_file.write(record.str().c_str(), record.str().size());
		
		int str_pos = 0;
		record[str_pos] = '>';
		str_pos++;
		for(int i = 0; i < ks1->name.l; i++) record[i+str_pos] = ks1->name.s[i];
		str_pos += ks1->name.l;
		if(!comment_flag)
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
   		if(!comment_flag) record_len =  ks1->name.l + ks1->comment.l + ks1->seq.l + 4;
		else record_len = ks1->name.l + ks1->seq.l + 3;

		fa_file.write(record, record_len);

		pos += record_len; //note that with extra '>' ' ' and two '\n'

		number_seqs++;
	}

//	for(auto &i:fai)
//		cerr << i << endl;	
	ofile.write(reinterpret_cast<const char*>(fai.data()), fai.size() * sizeof(uint64_t));

	cerr << "number of sequences: " << number_seqs << endl;
	
	ofile.close();
	fa_file.close();

	return 0;
}
