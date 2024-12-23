#ifndef __ORDERMINHASH_H__
#define __ORDERMINHASH_H__

#include <cstdint>
#include <vector>

namespace Sketch{

struct OSketch {
	//std::string       name;
	int               k, l, m;
	std::vector<char> data;

};

class OrderMinHash{

	public:

		/// OrderMinHash constructor for sketching sequences using default parameters
		OrderMinHash() : seq(nullptr), m_k(9), m_l(2), m_m(15), mtSeed(32)  {};
		OrderMinHash(char * seqNew);
		~OrderMinHash() {};

		/// return sketch result in `OSketch` type
		OSketch getSektch(){ return sk;}

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

		//Default parameters;
		int m_k = 9, m_l = 2, m_m = 15;

		OSketch sk;
		uint64_t mtSeed = 32; //default value

		void sketch();

		inline void compute_sketch(char * ptr, const char * seq);

};

struct mer_info {
	size_t pos;
	uint64_t hash;
	uint64_t int_hash;
	unsigned occ;
	mer_info(){}
	mer_info(size_t p, unsigned o, uint64_t h, uint64_t oh)
		: pos(p)
		  , hash(h)
		  , occ(o)
		  , int_hash(oh)
	{ }
};

}

#endif
