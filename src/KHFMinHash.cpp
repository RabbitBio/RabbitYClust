#include "KHFMinHash.h"
#include "xxhash.hpp"
#include "int_hash.h"
#include "aahash/aahash.hpp"
#include "aahash/helpers.hpp"
#include "AminoEncode.h"

#include <algorithm>
#include <queue>
#include <iostream>
#include <random>
#include <map>

#include <sys/time.h>
#include <limits.h>
#include <immintrin.h>

using namespace std;
extern std::map<std::string, int> top_kmer_map;
namespace Sketch
{

KHFMinHash::KHFMinHash(const char * seqNew):
	seq(seqNew)
{
	sketch();
}

void KHFMinHash::buildSketch(const char * seqNew = NULL)
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
void KHFMinHash::buildSketch(const char * seqNew, std::vector<std::string>& seed_strings, unsigned h, unsigned hash_num_per_seed)
{
	if(seqNew == NULL)
	{
		if(seq == NULL)
		{
			std::cerr << "WARNING: no data found" << std::endl;
			return;
		}

		sketchByAAHash(seed_strings, h, hash_num_per_seed);
	} else {
		seq = seqNew;

		sketchByAAHash(seed_strings, h, hash_num_per_seed);

	}
}
void KHFMinHash::buildSketchByNoSeedAAHash(const char * seqNew) {	
	if(seqNew == NULL)
	{
		if(seq == NULL)
		{
			std::cerr << "WARNING: no data found" << std::endl;
			return;
		}
		sketchByNoSeedAAHash();
	} else {
		seq = seqNew;

		sketchByNoSeedAAHash();

	}
}
void KHFMinHash::sketchByAAHash(std::vector<std::string>& seed_strings, unsigned h, unsigned hash_num_per_seed)
{
	sk.k = m_k;		
	sk.l = m_l;		
	sk.m = m_m;		


	sk.hashes.resize(sk.l * sk.m);
	std::string seqStr(seq);

	for(int i = 0; i < sk.l * sk.m; i++) sk.hashes[i] = ULONG_MAX;
	uint64_t *ptr = sk.hashes.data();

 	// 使用随机生成的种子初始化 SeedAAHash
    btllib::SeedAAHash seedaahash(seqStr.data(), seed_strings, hash_num_per_seed, sk.k);
    size_t pos = seedaahash.get_pos();
    // 迭代通过所有 k-mer，直到不能继续滚动
    bool success = true;
    while (success) {
        success = seedaahash.roll();  // 计算下一个 k-mer 的哈希值
        if (success) {
            // 获取当前的哈希值
            const uint64_t* hashes = seedaahash.hashes();
            // Hashes for position : seedaahash.get_pos()
            for(unsigned i = 0; i < h * hash_num_per_seed; i++){
				ptr[i] = std::min(ptr[i], hashes[i]);
			}
        }
    }
}
void KHFMinHash::sketchByNoSeedAAHash() 
{

	sk.k = m_k;		
	sk.l = m_l;		
	sk.m = m_m;		
	int count=0;

	sk.hashes.resize(sk.l * sk.m);
	std::string seqStr(seq);

	for(int i = 0; i < sk.l * sk.m; i++) sk.hashes[i] = ULONG_MAX;
	std::vector<std::vector<uint64_t>> kmers;

	uint64_t *ptr = sk.hashes.data();
    btllib::AAHash aahash(seqStr.data(), m_m, m_k, 1, 0); // level 1 2 3

	size_t pos = aahash.get_pos();

    bool success = true;

    while (success) {
        success = aahash.roll();  // 计算下一个 k-mer 的哈希值
		 std::vector<uint64_t> hashes_vec;
        if (success) { 
            // 获取当前的哈希值
            const uint64_t* hashes = aahash.hashes();
            // Hashes for position : aahash.get_pos()
			 std::string current_kmer = seqStr.substr(aahash.get_pos(), m_k);  // 取当前 k-mer 序列
			 if (top_kmer_map.find(current_kmer) != top_kmer_map.end()&&count<(m_l+m_k-1)) {
				count++;
                continue; // 直接跳过高频 k-mer
            }
            for(unsigned i = 0; i < m_m; i++){
				hashes_vec.push_back( hashes[i]);
			}
			kmers.emplace_back(hashes_vec);
        }
    }

	for (int i = 0; i <  (int)kmers.size(); i++) {
		for(unsigned j = 0; j < m_m; j++){
			ptr[j] = std::min(ptr[j], kmers[i][j]);
		}
    }
}
void KHFMinHash::sketch()
{
	sk.k = m_k;		
	sk.l = m_l;		
	sk.m = m_m;		

	sk.hashes.resize(sk.l * sk.m);
	std::string seqStr(seq);

	for(int i = 0; i < sk.l * sk.m; i++) sk.hashes[i] = ULONG_MAX;
	uint64_t *ptr = sk.hashes.data();

	mt19937_64 rng;
	rng.seed(mtSeed);
	uint64_t random_seed;
    if(inthash)
	{
		if(seqStr.size() < m_k)
		{
			for(int j = 0; j < m_m; j++)
			{
				random_seed = rng();
				uint64_t int_value = encodeAminoAcidsTo64Bit(seqStr);
				//FIXME: using random seed!
				uint64_t hash_value = jenkins_hash_seed(int_value, random_seed);
				ptr[j] = std::min(ptr[j], hash_value);
			}
			return;
		}

		for(int j = 0; j < m_m; j++)
		{
			random_seed = rng();
			for(int i = 0; i < seqStr.size() - m_k + 1; i++)
			{
				std::string kmer_seq = seqStr.substr(i, m_k);	
				uint64_t int_value = encodeAminoAcidsTo64Bit(kmer_seq);
				uint64_t hash_value = jenkins_hash_seed(int_value, random_seed);
				ptr[j] = std::min(ptr[j], hash_value);
		
			}
		}
	}else{
		xxhash hash;
		if(seqStr.size() < m_k)
		{
			for(int j = 0; j < m_m; j++)
			{
				random_seed = rng();
				hash.reset(random_seed);//TODO: make it more random seed!
				hash.update(seqStr.data(), seqStr.size());
				uint64_t hash_value = hash.digest();
				ptr[j] = std::min(ptr[j], hash_value);
			}
			return;
		}		
	
		for(int j = 0; j < m_m; j++)
		{
	
			random_seed = rng();
			for(int i = 0; i < seqStr.size() - m_k + 1; i++)
			{
				hash.reset(random_seed);//TODO: make it more random seed!
				hash.update(&seqStr.data()[i], m_k);
				uint64_t hash_value = hash.digest();
				ptr[j] = std::min(ptr[j], hash_value);
			}
		}
	}

	return;
}


}// namespace Sketch

