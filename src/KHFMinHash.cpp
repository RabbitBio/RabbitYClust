
#include "KHFMinHash.h"
#include "xxhash.hpp"
#include "int_hash.h"

#include <algorithm>
#include <queue>
#include <iostream>
#include <random>

#include <sys/time.h>
#include <limits.h>
#include <immintrin.h>

// Precomputed lookup table for 5-bit binary codes
// UJBZ X are all used!
uint8_t aminoAcidTable[26] = {
		0b00000,  // A - Alanine
		0b10111,  // B - Placeholder for non-standard amino acid (not used here)
		0b00001,  // C - Cysteine
		0b00010,  // D - Aspartic Acid
		0b00011,  // E - Glutamic Acid
		0b00100,  // F - Phenylalanine
		0b00101,  // G - Glycine
		0b00110,  // H - Histidine
		0b00111,  // I - Isoleucine
		0b10110,  // J - Placeholder for non-standard amino acid (not used here)
		0b01000,  // K - Lysine
		0b01001,  // L - Leucine
		0b01010,  // M - Methionine
		0b01011,  // N - Asparagine
		0b11111,  // O - Placeholder for non-standard amino acid (not used here)
		0b01100,  // P - Proline
		0b01101,  // Q - Glutamine
		0b01110,  // R - Arginine
		0b01111,  // S - Serine
		0b10000,  // T - Threonine
		0b10101,  // U - Placeholder for non-standard amino acid (not used here)
		0b10001,  // V - Valine
		0b10010,  // W - Tryptophan
		0b10100,  // X - For any kind of amino acid
		0b10011,  // Y - Tyrosine
		0b11000   // Z - Placeholder for non-standard amino acid (not used here)
};


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
				hash.reset(j);//TODO: make it more random seed!
				hash.update(seqStr.data(), seqStr.size());
				uint64_t hash_value = hash.digest();
				ptr[j] = std::min(ptr[j], hash_value);
			}
			return;
		}		
	
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

	return;
}


}// namespace Sketch

// Function to encode a string of amino acids to a 64-bit integer
uint64_t encodeAminoAcidsTo64Bit(const std::string& sequence) {
	if (sequence.length() > 12) {
		std::cerr << "Sequence length exceeds 12 amino acids." << std::endl;
		return 0;
	}

	uint64_t encoded = 0;

	for (size_t i = 0; i < sequence.length(); ++i) {
		char aminoAcid = sequence[i];
		if (aminoAcid < 'A' || aminoAcid > 'Z' || aminoAcidTable[aminoAcid - 'A'] == 0b11111) {
				std::cerr << "Invalid amino acid: " << aminoAcid << std::endl;
				return 0;
		}

		uint8_t code = aminoAcidTable[aminoAcid - 'A'];
		// Shift the current code into its correct position in the 64-bit integer
		encoded |= (static_cast<uint64_t>(code) << (5 * (11 - i)));
	}

	return encoded;
}
