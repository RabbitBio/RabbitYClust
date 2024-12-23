#ifndef __AMINOENCODE_H__
#define __AMINOENCODE_H__

#include <string>

extern uint8_t aminoAcidTable[];
// Function to encode a string of amino acids to a 64-bit integer
uint64_t encodeAminoAcidsTo64Bit(const std::string& sequence);

#endif
