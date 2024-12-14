#include <zlib.h>
#include "CLI11.hpp"

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{
	CLI::App app{"read_fai v.0.0.1, read fasta index of protein files"};
	string filename = "";
	string fafilename = "";
	string seqIDfilename = "";
	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	auto option_fai = app.add_option("-f, --fai", fafilename, "input fa-file name");
	auto option_seqID = app.add_option("-s, --sequences", seqIDfilename, "input sequences ID file name");

	option_input->required();
	option_fai->required();
	option_seqID->required();

	CLI11_PARSE(app, argc, argv);

	ifstream faifile(fafilename, ios::binary);
	ifstream fainput(filename);
	ifstream seqIDs(seqIDfilename);

	if( !fainput )
	{
		cerr << "Fail to open file " << filename << endl;
	}
	if( !faifile )
	{
		cerr << "Fail to open file " << fafilename << endl;
	}
	if( !seqIDs )
	{
		cerr << "Fail to open file " << seqIDfilename << endl;
	}

	faifile.seekg(0, ios::end);
    size_t fileSize = faifile.tellg();
    faifile.seekg(0, ios::beg);
    size_t numElements = fileSize / sizeof(uint64_t);

    // 创建一个 vector<int> 来存储读取的数据
	std::vector<uint64_t> fa_index(numElements);
	cout << "total fa index nums: " << fa_index.size() << endl;

	// 读取数据到 vector<int>
	faifile.read(reinterpret_cast<char*>(fa_index.data()), fileSize);
	faifile.close();

	int id;
	while(seqIDs >> id)
	{
		uint64_t start_pos = fa_index[id];
		uint64_t end_pos = fa_index[id+1];
		int length = (int)(end_pos - start_pos);
		// 定位到偏移位置
		fainput.seekg(start_pos, ios::beg);
		string sequence;
		sequence.resize(length);
		fainput.read(&sequence[0], length);
		cout << id << " " << sequence << endl;
	}
	return 0;
}
