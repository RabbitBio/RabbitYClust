
#include "KHFMinHash.h"
#include "xxhash.hpp"

#include <algorithm>
#include <queue>
#include <iostream>

#include <sys/time.h>
#include <limits.h>
#include <immintrin.h>

using namespace std;


namespace Sketch
{

KHFMinHash::KHFMinHash(char * seqNew):
	seq(seqNew)
{
	sketch();
}

void KHFMinHash::buildSketch(char * seqNew = NULL)
{
	// rebuild sketch using old data
	if(seqNew == NULL)
	{
		if(seq == NULL)
		{
			std::cerr << "WARNING: no data found" << std::endl;
			return;
		}

		sketch();
	} else {
		seq = seqNew;

		sketch();

	}
}
void KHFMinHash::sketch()
{
	sk.k = m_k;		
	sk.l = m_l;		
	sk.m = m_m;		

	sk.hashes.resize(sk.l * sk.m);
	xxhash hash;
	std::string seqStr(seq);
	
	for(int i = 0; i < sk.l * sk.m; i++) sk.hashes[i] = ULONG_MAX;
	uint64_t *ptr = sk.hashes.data();

	for(int j = 0; j < m_m; j++)
	{

		for(int i = 0; i < seqStr.size() - m_k + 1; i++)
		{
			hash.reset(j);//TODO: make it more random seed!
			hash.update(&seqStr.data()[i], m_k);
			uint64_t hash_value = hash.digest();
			ptr[j] = std::min(ptr[j], hash_value);
		}
	}


}


}// namespace Sketch
