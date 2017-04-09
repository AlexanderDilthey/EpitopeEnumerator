/*
 * readFiles.h
 *
 *  Created on: Mar 31, 2017
 *      Author: diltheyat
 */

#ifndef READFILES_H_
#define READFILES_H_

#include <string>
#include <vector>
#include <map>

#include "Utilities.h"

class variantFromVCF {
public:
	std::string chromosomeID;
	unsigned int position;
	std::string referenceString;
	std::vector<std::string> sampleAlleles;
};

class transcriptExon {
public:
	bool valid;
	unsigned int firstPos;
	unsigned int lastPos;
	transcriptExon()
	{
		valid = false;
	}
};

class transcript {
public:
	std::string transcriptID;
	std::string geneName;
	unsigned char strand;
	std::string chromosomeID;
	std::vector<transcriptExon> exons;
};

std::map<std::string, std::map<int, variantFromVCF>> readVariants(std::string VCF);
std::vector<transcript> readTranscripts(std::string transcriptsFile);
unsigned int countCharacters_noGaps(const std::string& S);

std::string translateCodon2AA(const std::string& codon);

void extendAsNecessary(std::string& S, unsigned int desiredLength);

std::vector<transcript> getPlusStrandTranscripts(const std::vector<transcript>& transcripts);
std::vector<transcript> getMinusStrandTranscripts(const std::vector<transcript>& transcripts, const std::map<std::string, std::string>& referenceGenome_minus);
std::map<std::string, std::string> getMinusStrandReferenceGenome(const std::map<std::string, std::string>& referenceGenome);
std::map<std::string, std::map<int, variantFromVCF>> getMinusStrandVariants(const std::map<std::string, std::map<int, variantFromVCF>>& variants, const std::map<std::string, std::string>& referenceGenome_minus);
void checkVariantsConsistentWithReferenceGenome(const std::map<std::string, std::map<int, variantFromVCF>>& variants, const std::map<std::string, std::string>& referenceGenome);

std::string seq_reverse_complement(const std::string& sequence);
char reverse_char_nucleotide(char c);

#endif /* READFILES_H_ */
