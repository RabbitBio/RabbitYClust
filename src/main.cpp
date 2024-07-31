#include "KHFMinHash.h"
#include "kseq.h"
#include <zlib.h>

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>

using namespace std;
using namespace Sketch;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char** argv)
{

	if(argc != 3)
	{
		cerr << "Usage: ./yclust input.fa output.hash" << endl;
		cerr << "Note: output.hash is in binary format for uint64_t hash values!" << endl;
		cerr << "Note: the names of sequences are printed to stdout!" << endl;
		return 1;
	}
	gzFile fp1;
	kseq_t *ks1;

	fp1 = gzopen(argv[1],"r");

	if(NULL == fp1){
		fprintf(stderr,"Fail to open file: %s\n", argv[1]);
		return 0;
	}


	ofstream ofile(string(argv[2]), ios::binary);

	if( !ofile.is_open())
	{
		cerr << "Failed to open " << argv[2] << endl;
		return 1;
	}

	ks1 = kseq_init(fp1);

	int buffer_size = 1<<20;//1MB
	int buffer_pos = 0;
	vector<char> buffer(buffer_size);

	int count = 0;	

	while(1)
	{
		int length = kseq_read(ks1);

		if(length < 0) break;

		char * seq1 = ks1->seq.s;
		cout << ks1->name.s << " ";
		cout << ks1->comment.s << endl;
		KHFMinHash mh = KHFMinHash(seq1);			
	
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
