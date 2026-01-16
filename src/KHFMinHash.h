/* This is an efficient implementation of KHF MinHash sketching algorithms. */
/* Attention: This is only used for protein sequences! */
#ifndef __KHFMINHASH_H__
#define __KHFMINHASH_H__

#include <string>
#include <stdint.h>
#include <vector>
#include <algorithm>

namespace Sketch{

	// Simple flat hash set for cache-friendly O(1) lookup
	// Uses open addressing with linear probing
	class FlatHashSet64 {
	public:
		static constexpr uint64_t EMPTY = UINT64_MAX;
		
		FlatHashSet64() : data_(nullptr), capacity_(0), size_(0) {}
		
		~FlatHashSet64() { 
			delete[] data_; 
		}
		
		// Build from a vector of values (call once after loading all kmers)
		void build(const std::vector<uint64_t>& values) {
			size_ = values.size();
			if (size_ == 0) {
				capacity_ = 0;
				delete[] data_;
				data_ = nullptr;
				return;
			}
			
			// Use ~2x capacity for good load factor (~0.5)
			capacity_ = 1;
			while (capacity_ < size_ * 2) capacity_ <<= 1;
			mask_ = capacity_ - 1;
			
			delete[] data_;
			data_ = new uint64_t[capacity_];
			std::fill(data_, data_ + capacity_, EMPTY);
			
			for (uint64_t val : values) {
				if (val == EMPTY) continue;  // Skip UINT64_MAX values
				uint64_t idx = hash(val) & mask_;
				while (data_[idx] != EMPTY) {
					idx = (idx + 1) & mask_;
				}
				data_[idx] = val;
			}
		}
		
		// O(1) average lookup, cache-friendly
		inline bool contains(uint64_t val) const {
			if (capacity_ == 0 || val == EMPTY) return false;
			uint64_t idx = hash(val) & mask_;
			while (data_[idx] != EMPTY) {
				if (data_[idx] == val) return true;
				idx = (idx + 1) & mask_;
			}
			return false;
		}
		
		void clear() {
			delete[] data_;
			data_ = nullptr;
			capacity_ = 0;
			size_ = 0;
		}
		
		size_t size() const { return size_; }
		bool empty() const { return size_ == 0; }
		
	private:
		// Fast hash mixing for uint64_t (splitmix64-style)
		static inline uint64_t hash(uint64_t x) {
			x ^= x >> 33;
			x *= 0xff51afd7ed558ccdULL;
			x ^= x >> 33;
			x *= 0xc4ceb9fe1a85ec53ULL;
			x ^= x >> 33;
			return x;
		}
		
		uint64_t* data_;
		size_t capacity_;
		size_t size_;
		size_t mask_;
	};


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
			KHFMinHash() : seq(NULL), m_k(8), m_l(1), m_m(15), mtSeed(32)  {};
			/// KHFMinHash constructor for sketching sequences using default parameters
			KHFMinHash(const char * seqNew);
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
			void buildSketch(const char * seqNew);
			void buildSketch(const char * seqNew, std::vector<std::string>& seed_strings, unsigned h, unsigned hash_num_per_seed);
			void buildSketchByNoSeedAAHash(const char * seqNew);

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

			/// Check if a kmer (encoded as uint64_t) is in the blacklist (O(1) flat hash lookup)
			inline bool isInBlacklist(uint64_t kmerHash) const { 
				if (external_blacklist_) return external_blacklist_->contains(kmerHash);
				return blacklist.contains(kmerHash); 
			}

			/// Clear the blacklist
			void clearBlacklist() { blacklist.clear(); external_blacklist_ = nullptr; }

			/// Get blacklist size
			size_t getBlacklistSize() const { 
				if (external_blacklist_) return external_blacklist_->size();
				return blacklist.size(); 
			}

			/// Set external blacklist (shared among multiple instances, avoids reloading)
			void setExternalBlacklist(const FlatHashSet64* bl) { external_blacklist_ = bl; }

		private:

			const char * seq = NULL;
			//kmer: 8 for protein in default, if using int hash, kmer max size is 12
			//m: 15 
			int m_k = 8, m_m = 15;
			//For the classical KHF implementation l should be 1
			int m_l = 1; 
            //choose whether using int hash
            bool inthash = false;

			KHFSketch sk;

			uint64_t mtSeed = 32; //default value

			// Kmer blacklist (flat hash set for O(1) cache-friendly lookup)
			FlatHashSet64 blacklist;
			// External blacklist pointer (shared, not owned)
			const FlatHashSet64* external_blacklist_ = nullptr;

			void sketch();

			void sketchByAAHash(std::vector<std::string>& seed_strings, unsigned h, unsigned hash_num_per_seed);
			void sketchByNoSeedAAHash();
	};
}

#endif
