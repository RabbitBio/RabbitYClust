#include "AminoEncode.h"
#include <iostream>

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
