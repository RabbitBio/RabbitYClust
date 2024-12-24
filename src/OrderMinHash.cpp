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
#include <iostream>
#include <immintrin.h>

using namespace std;
#if defined __AVX512F__ && defined __AVX512DQ__
static void inspect(__m512i va)
{
	uint64_t a[8];
	_mm512_storeu_si512(a, va);
	cout << "[ ";
	for(int i = 0; i<8;i++)
		cout << a[i] << " ";
	cout << "]" << endl;
}
#else
#ifdef __AVX2__
// implement by avx2
inline __m256i avx2_mullo_epi64(__m256i a1, __m256i b1)
{
    __m256i albl = _mm256_mul_epu32(a1, b1);

    const int shuffle= 0xb1;
    const int blendmask = 0x55;
    __m256i a1shuffle = _mm256_shuffle_epi32(a1, shuffle);
    __m256i b1shuffle = _mm256_shuffle_epi32(b1, shuffle);

    __m256i ahbl = _mm256_mul_epu32(a1shuffle, b1);
    __m256i albh = _mm256_mul_epu32(a1, b1shuffle);

    ahbl = _mm256_add_epi64(ahbl, albh);//sum of albh ahbl
    ahbl = _mm256_shuffle_epi32(ahbl, shuffle);

    albh = _mm256_blend_epi32(ahbl, albl, blendmask);//res without add hi32 of albl
    __m256i zero = _mm256_set1_epi32(0);
    albl = _mm256_blend_epi32(albl, zero, blendmask);

    albh = _mm256_add_epi32(albh, albl);
    return albh;

}
inline void inspect_avx2(__m256i va)
{
	cout << "[ ";
	for(int i = 0; i < 4; i++)
		cout << ((uint64_t*)&va)[i] << " ";
	cout << "]" << endl;
}
#endif

#endif



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
		std::string kmer_seq = seqStr.substr(i, m_k);	
		intHash[i] = encodeAminoAcidsTo64Bit(kmer_seq);
	}

	//  Create list of k-mers with occurrence numbers
	for(size_t i = 0; i < seqStr.size() - m_k + 1; ++i) {
		auto tmpocc = occurrences[intHash[i]]++;
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

#define lanes 8
	int pend_size = (mers.size() / lanes) * lanes;
	uint64_t hashBuffer[lanes];


	//AVX512
	//__m512i vweight = _mm512_set1_epi64(weight);

#if defined __AVX512F__  && defined __AVX512DQ__
	__m512i vzero   = _mm512_set1_epi64(0);
	__m512i vconst0 = _mm512_set1_epi64(0xff51afd7ed558ccd);
	__m512i vconst1 = _mm512_set1_epi64(0xc4ceb9fe1a85ec53);
	__mmask8 weight_msk = m_l > 0 ? 0xFF : 0x00;
#else
#if defined __AVX2__
	//AVX2
	__m256i vweight = _mm256_set1_epi64x(weight);
	__m256i vconst0 = _mm256_set1_epi64x(0xff51afd7ed558ccd);
	__m256i vconst1 = _mm256_set1_epi64x(0xc4ceb9fe1a85ec53);
#else
#endif
#endif

	//main sketch loop
	//for(unsigned i = 0; i < m; ++i) 
		for( int id = 0; id < pend_size; id+=lanes)
	{
		//std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> & pqueue = pqueues[i];
		//body
	for(unsigned i = 0; i < m_m; ++i) 
		//for( int id = 0; id < pend_size; id+=lanes)
		{
			std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> & pqueue = pqueues[i];


#if defined __AVX512F__ && __AVX512DQ__
			//using AVX512
			__m512i va = _mm512_loadu_si512((void *)&intHash[id + 0 * 8]);
			__m512i vocc = _mm512_loadu_si512((void *)&occ[id + 0 * 8]);
			__m512i vb = _mm512_mask_add_epi64(va, weight_msk, va, vocc);
			__m512i vseed = _mm512_set1_epi64(mseed[i]);
			va = _mm512_xor_epi64(vb, vseed);
			__m512i vtmp = _mm512_srli_epi64(va, 33);
			vb = _mm512_xor_epi64(va, vtmp);
			va = _mm512_mullo_epi64(vb, vconst0);

			vtmp = _mm512_srli_epi64(va, 33);
			vb = _mm512_xor_epi64(va, vtmp);
			va = _mm512_mullo_epi64(vb, vconst1);
			vtmp = _mm512_srli_epi64(va, 33);
			vb = _mm512_xor_epi64(va, vtmp);
			_mm512_storeu_si512(hashBuffer + 0 * 8, vb);
#else
#if defined __AVX2__
			//using intrinsics using avx2 and unroll2
			//TODO: not heavily tested
			__m256i va0, va1;
			__m256i vb0, vb1;
			__m256i vseed = _mm256_set1_epi64x(mseed[i]);

			__m256i vinthash0 = _mm256_loadu_si256((__m256i*)(intHash + id));
			__m256i vinthash1 = _mm256_loadu_si256((__m256i*)(intHash + id + 4));
			__m256i vocc0     = _mm256_loadu_si256((__m256i*)(occ + id));
			__m256i vocc1     = _mm256_loadu_si256((__m256i*)(occ + id + 4));
			va0 = avx2_mullo_epi64(vocc0, vweight);
			va1 = avx2_mullo_epi64(vocc1, vweight);
			vb0 = _mm256_add_epi64(vinthash0, va0); 
			vb1 = _mm256_add_epi64(vinthash1, va1);
			va0 = _mm256_xor_si256(vb0, vseed); 
			va1 = _mm256_xor_si256(vb1, vseed); 

			vinthash0 = _mm256_srli_epi64(va0, 33);
			vinthash1 = _mm256_srli_epi64(va1, 33);
			vb0 = _mm256_xor_si256(va0, vinthash0);
			vb1 = _mm256_xor_si256(va1, vinthash1);
			va0 = avx2_mullo_epi64(vconst0, vb0);
			va1 = avx2_mullo_epi64(vconst0, vb1);

			vinthash0 = _mm256_srli_epi64(va0, 33);
			vinthash1 = _mm256_srli_epi64(va1, 33);
			vb0 = _mm256_xor_si256(va0, vinthash0);
			vb1 = _mm256_xor_si256(va1, vinthash1);
			va0 = avx2_mullo_epi64(vconst1, vb0);
			va1 = avx2_mullo_epi64(vconst1, vb1);

			vinthash0 = _mm256_srli_epi64(va0, 33);
			vinthash1 = _mm256_srli_epi64(va1, 33);
			vb0 = _mm256_xor_si256(va0, vinthash0);
			vb1 = _mm256_xor_si256(va1, vinthash1);
			//inspect_avx2(vb0);
			//inspect_avx2(vb1);
			_mm256_storeu_si256((__m256i*) hashBuffer    , vb0);
			_mm256_storeu_si256((__m256i*) (hashBuffer + 4), vb1);
			//cout << "avx2 results:" << endl;
			//for(int vid = 0; vid < lanes; vid++)
			//	cout << hashBuffer[vid] << " ";
			//cout << endl;

#else
			for(int vid = 0; vid < lanes; ++vid)
			{
			uint64_t kmer_int = intHash[id + vid];
			kmer_int += occ[id + vid] * weight;
			hashBuffer[vid] = mc::murmur3_fmix(kmer_int, mseed[i]);
			}
			//cout << "correct results:" << endl;
			//for(int vid = 0; vid < lanes; vid++)
			//	cout << hashBuffer[vid] << " ";
			//cout << endl;
#endif 
#endif
			//exit(0);
			//for(int vid = 0; vid < lanes; ++vid)
			//{
			//if( pqueue.size() < l || hashBuffer[vid] < pqueue.top().hash)
			//{
			//	pqueue.emplace(id+vid, occ[id + vid], hashBuffer[vid], 0);
			//	//pqueue.push(one_mer);
			//	if(pqueue.size() > l) pqueue.pop();
			//}
			//}

			if( pqueue.size() < m_l || hashBuffer[0] < pqueue.top().hash)
			{
				pqueue.emplace(id+0, occ[id + 0], hashBuffer[0], 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}
			if( pqueue.size() < m_l || hashBuffer[1] < pqueue.top().hash)
			{
				pqueue.emplace(id+1, occ[id + 1], hashBuffer[1], 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}
			if( pqueue.size() < m_l || hashBuffer[2] < pqueue.top().hash)
			{
				pqueue.emplace(id+2, occ[id + 2], hashBuffer[2], 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}
			if( pqueue.size() < m_l || hashBuffer[3] < pqueue.top().hash)
			{
				pqueue.emplace(id+3, occ[id + 3], hashBuffer[3], 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}
			if( pqueue.size() < m_l || hashBuffer[4] < pqueue.top().hash)
			{
				pqueue.emplace(id+4, occ[id + 4], hashBuffer[4], 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}
			if( pqueue.size() < m_l || hashBuffer[5] < pqueue.top().hash)
			{
				pqueue.emplace(id+5, occ[id + 5], hashBuffer[5], 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}
			if( pqueue.size() < m_l || hashBuffer[6] < pqueue.top().hash)
			{
				pqueue.emplace(id+6, occ[id + 6], hashBuffer[6], 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}
			if( pqueue.size() < m_l || hashBuffer[7] < pqueue.top().hash)
			{
				pqueue.emplace(id+7, occ[id + 7], hashBuffer[7], 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}

		}
	}
	//for(unsigned i = 0; i < m; ++i) 
		for( int id = pend_size; id < mers.size(); id++)
	{
		//std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> & pqueue = pqueues[i];
		//tail
	for(unsigned i = 0; i < m_m; ++i) 
		//for( int id = pend_size; id < mers.size(); id++)
		{
			std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> & pqueue = pqueues[i];
			uint64_t kmer_int = intHash[id];
		    kmer_int += occ[id] * weight;
			//one_mer.pos = id;
			//one_mer.occ = occ[id];
			uint64_t kmer_hash = murmur3_fmix(kmer_int, mseed[i]);
			if( pqueue.size() < m_l || kmer_hash < pqueue.top().hash)
			{
				//pqueue.push(one_mer);
				pqueue.emplace(id, occ[id], kmer_hash, 0);
				if(pqueue.size() > m_l) pqueue.pop();
			}

		}
	}
//	for( int id = 0; id < mers.size(); id++)
//	{
//		for(unsigned i = 0; i < m_m; ++i) 
//		{
//			std::priority_queue<mer_info, std::vector<mer_info>, decltype(cmp)> & pqueue = pqueues[i];
//			uint64_t kmer_int = intHash[id];
//		    kmer_int += occ[id] * weight;
//
//			uint64_t kmer_hash = murmur3_fmix(kmer_int, mseed[i]);
//			if( pqueue.size() < m_l || kmer_hash < pqueue.top().hash)
//			{
//				pqueue.emplace(id, occ[id], kmer_hash, 0);
//				if(pqueue.size() > m_l) pqueue.pop();
//			}
//
//		}
//	}

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

		for(int j = 0; j < m_l; ++j)
		{
			sk.hashes[i*m_l + j] = lmers[j].hash;
			sk.positions[i*m_l + j] = lmers[j].pos;
		}
	}

	delete[] intHash;
	delete[] occ;
	delete[] mseed;
	return;

}


}
