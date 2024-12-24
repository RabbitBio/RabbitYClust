#include "int_hash.h"
#include "OrderMinHash.h"
#include "AminoEncode.h"

#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <random>
#include <cassert>
#include <cstring>

using namespace std;

namespace Sketch
{

OrderMinHash::OrderMinHash(char * seqNew):
	seq(seqNew)
{
	sketch();
}

void OrderMinHash::sketch()
{
	sk.k = m_k;		
	sk.l = m_l;		
	sk.m = m_m;		

	sk.hashes.resize(std::max(sk.l, 1) * sk.m);
	sk.positions.resize(std::max(sk.l, 1) * sk.m);

	std::string seqStr(this->seq);

	if(seqStr.size() < m_k) return;
	const uint64_t weight = m_l > 0 ? 1 : 0;
	std::vector<mer_info> mers;
	std::unordered_map<uint64_t, unsigned> occurrences;

	uint64_t * intHash = new uint64_t[seqStr.size() - m_k + 1];
	uint64_t * occ     = new uint64_t[seqStr.size() - m_k + 1];

	for(int i = 0; i < seqStr.size() - m_k + 1; i++)
	{
		//intHash[i] = hash_to_uint(&seq.data()[i], k);
		std::string kmer_seq = seqStr.substr(i, m_k);	
		intHash[i] = encodeAminoAcidsTo64Bit(kmer_seq);
	}

	//  Create list of k-mers with occurrence numbers
	for(size_t i = 0; i < seqStr.size() - m_k + 1; ++i) {
		//auto occ = occurrences[seq.substr(i, k)]++;
		auto tmpocc = occurrences[intHash[i]]++;
		//mers.emplace_back(i, occ, (uint64_t)0, hash_to_uint(&seq.data()[i], k));
		mers.emplace_back(i, tmpocc, (uint64_t)0, 0);
		occ[i] = mers[i].occ;
	}

	std::mt19937_64 gen64(mtSeed); //TODO: make 32 a parameter
	auto cmp = [](Sketch::mer_info & a, Sketch::mer_info & b){return a.hash < b.hash;};
	vector<std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> > pqueues;
	uint64_t * mseed = new uint64_t[m_m];
	for(int i = 0; i < m_m; i++)
	{
		pqueues.emplace_back(cmp);
		mseed[i] = gen64();
	}

	std::vector<mer_info> lmers;
	lmers.reserve(m_l);

	for( int id = 0; id < mers.size(); id++)
	{
		for(unsigned i = 0; i < m_m; ++i) 
		{
			std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> & pqueue = pqueues[i];
			uint64_t kmer_int = intHash[id];
		    kmer_int += occ[id] * weight;

			uint64_t kmer_hash = murmur3_fmix(kmer_int, mseed[i]);
			if( pqueue.size() < m_l || kmer_hash < pqueue.top().hash)
			{
				//pqueue.push(one_mer);
				pqueue.emplace(id, occ[id], kmer_hash, 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}

		}
	}

	for(unsigned i = 0; i < m_m; ++i) 
	{
		std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> & pqueue = pqueues[i];
		lmers.clear();
		while(!pqueue.empty())
		{
			lmers.push_back(pqueue.top());
			pqueue.pop();
		}
		std::sort(lmers.begin(), lmers.end(), [&](const mer_info& x, const mer_info& y) { return x.pos < y.pos; });
		assert(lmers.size() == m_l);

		//	block(i, j, lmers[j].pos);
		//char *ptr = sk.data.data();
		for(int j = 0; j < m_l; ++j)
		{
			sk.hashes[i*m_l + j] = lmers[j].hash;
			sk.positions[i*m_l + j] = lmers[j].pos;
			//memcpy(ptr, &seqStr.data()[lmers[j].pos], m_k);
			//ptr += m_k;
		}
	}

	delete[] intHash;
	delete[] occ;
	delete[] mseed;
	return;

}


}
