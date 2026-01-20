#include "KHFMinHash.h"
#include "xxhash.hpp"
#include "int_hash.h"
#include "aahash/aahash.hpp"
#include "aahash/helpers.hpp"
#include "AminoEncode.h"

#include <algorithm>
#include <cstdint>
#include <queue>
#include <iostream>
#include <random>

#include <sys/time.h>
#include <limits.h>
#include <immintrin.h>

using namespace std;


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
void KHFMinHash::buildSketchByNoSeedAAHash(const char * seqNew)
{
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

	sk.hashes.resize(sk.l * sk.m);
	std::string seqStr(seq);

	for(int i = 0; i < sk.l * sk.m; i++) sk.hashes[i] = ULONG_MAX;
	uint64_t *ptr = sk.hashes.data();

	// 预计算每个 k-mer 位置是否含有 'X' 或其他不确定氨基酸
	size_t numKmers = (seqStr.size() >= m_k) ? (seqStr.size() - m_k + 1) : 0;
	std::vector<bool> kmerHasX(numKmers, false);
	if (numKmers > 0) {
		int xCount = 0;
		// 初始化第一个窗口
		for (size_t i = 0; i < m_k; i++) {
			char c = seqStr[i];
			if (c == 'X' || c == 'B' || c == 'Z' || c == 'J' || c == 'U' || c == 'O') xCount++;
		}
		kmerHasX[0] = (xCount > 0);
		// 滑动窗口计算后续位置
		for (size_t i = 1; i < numKmers; i++) {
			char c_out = seqStr[i - 1];
			if (c_out == 'X' || c_out == 'B' || c_out == 'Z' || c_out == 'J' || c_out == 'U' || c_out == 'O') xCount--;        // 离开窗口的字符
			char c_in = seqStr[i + m_k - 1];
			if (c_in == 'X' || c_in == 'B' || c_in == 'Z' || c_in == 'J' || c_in == 'U' || c_in == 'O') xCount++;  // 进入窗口的字符
			kmerHasX[i] = (xCount > 0);
		}
	}

	btllib::AAHash aahash(seqStr.data(), m_m, m_k, 1, 0); // level 1 2 3

	bool success = true;
	while (success) {
		success = aahash.roll();  // 计算下一个 k-mer 的哈希值
		if (success) {
			size_t pos = aahash.get_pos();

			// 跳过含有 'X' 的 k-mer（使用预计算结果）
			if (pos < numKmers && kmerHasX[pos]) continue;

			// 检查是否在黑名单中
			std::string kmer_seq = seqStr.substr(pos, m_k);
			uint64_t int_value = encodeAminoAcidsTo64Bit(kmer_seq);
			if (isInBlacklist(int_value)) continue;

			// 获取当前的哈希值
			const uint64_t* hashes = aahash.hashes();
			for(unsigned i = 0; i < m_m; i++){
				//using second round hash function to avoid the shortcoming of AAHash
				uint64_t hash_value = murmur3_fmix(hashes[i], (uint64_t)i);
				ptr[i] = std::min(ptr[i], hash_value);
				//ptr[i] = std::min(ptr[i], hashes[i]);
			}
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

	// 预计算每个 k-mer 位置是否含有 'X' 或其他不确定氨基酸（滑动窗口优化）
	size_t numKmers = (seqStr.size() >= m_k) ? (seqStr.size() - m_k + 1) : 0;
	std::vector<bool> kmerHasX(numKmers, false);
	if (numKmers > 0) {
		int xCount = 0;
		// 初始化第一个窗口
		for (size_t i = 0; i < m_k; i++) {
			char c = seqStr[i];
			if (c == 'X' || c == 'B' || c == 'Z' || c == 'J' || c == 'U' || c == 'O') xCount++;
		}
		kmerHasX[0] = (xCount > 0);
		// 滑动窗口计算后续位置
		for (size_t i = 1; i < numKmers; i++) {
			char c_out = seqStr[i - 1];
			if (c_out == 'X' || c_out == 'B' || c_out == 'Z' || c_out == 'J' || c_out == 'U' || c_out == 'O') xCount--;
			char c_in = seqStr[i + m_k - 1];
			if (c_in == 'X' || c_in == 'B' || c_in == 'Z' || c_in == 'J' || c_in == 'U' || c_in == 'O') xCount++;
			kmerHasX[i] = (xCount > 0);
		}
	}

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
			for(size_t i = 0; i < numKmers; i++)
			{
				if (kmerHasX[i]) continue;  // 跳过含有 'X' 的 k-mer
				std::string kmer_seq = seqStr.substr(i, m_k);	
				uint64_t int_value = encodeAminoAcidsTo64Bit(kmer_seq);
				// 检查是否在黑名单中
				if (isInBlacklist(int_value)) continue;
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
			for(size_t i = 0; i < numKmers; i++)
			{
				if (kmerHasX[i]) continue;  // 跳过含有 'X' 的 k-mer
				// 检查是否在黑名单中
				std::string kmer_seq = seqStr.substr(i, m_k);
				uint64_t int_value = encodeAminoAcidsTo64Bit(kmer_seq);
				if (isInBlacklist(int_value)) continue;
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

