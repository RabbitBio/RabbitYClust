/* This is an efficient implementation of KHF MinHash sketching algorithms. */
/* Attention: This is only used for protein sequences! */
#pragma once

#include <string>
#include <stdint.h>
#include <vector>


// Function to encode a string of amino acids to a 64-bit integer
uint64_t encodeAminoAcidsTo64Bit(const std::string& sequence);


namespace Sketch{
	//OMinHash
	struct KHFSketch {
		std::string       name;
		//k: kmer size && m: m hash functions
		//l: number of hashes is used for comparing
		int               k, l, m;
		//std::vector<uint32_t> data32;
		//TODO: using 64bit for test
		std::vector<uint64_t> hashes;
	};

	/// Sketching and compare sequences or strings using KHF MinHash algorithm.
	class KHFMinHash{

		public:
			/// KHFMinHash constructor
			KHFMinHash() : seq(nullptr), m_k(12), m_l(1), m_m(15), mtSeed(32)  {};
			/// KHFMinHash constructor for sketching sequences using default parameters
			KHFMinHash(char * seqNew);
			~KHFMinHash() {};

			/// return sketch result in `OSketch` type
			KHFSketch & getSektch(){ return sk;}

			/** \rst
			  Build a `KHFMinHash` sketch.
			  `seqNew` is NULL pointer in default.
			  If seqNew is NULL pointer, buildSketch() will rebuild sketh using old data.
			  This is useful when chaning parameters and build a new sketch.
			 \endrst
			*/
			void buildSketch(char * seqNew);

			/// Set parameter `kmerSize`: default 21.
			void setK(int k){ m_k = k; }

			/// Set parameter `l`: default 2 (normally 2 - 5).
			void setL(int l){ m_l = l; }

			/// Set parameter `m`: default 500.
			void setM(int m){ m_m = m; }

			/// Set seed value for random generator: default 32.
			void setSeed(uint64_t seedNew) { mtSeed = seedNew; }

			/// Return parameter `kmerSize`.
			int getK(){return m_k;}	

			/// Return parameter `l`.
			int getL(){return m_l;}	

			/// Return parameter `m`.
			int getM(){return m_m;}		

			/// Return random generator seed value.
			uint64_t getSeed() { return mtSeed; }

		private:

			char * seq = NULL;
			//kmer: 12 for protein in default, if using int hash, kmer max size is 12
			//m: 15 
			int m_k = 12, m_m = 15;
			//For the classical KHF implementation l should be 1
			int m_l = 1; 
            //choose whether using int hash
            bool inthash = true;

			KHFSketch sk;

			uint64_t mtSeed = 32; //default value

			void sketch();

	};
}
