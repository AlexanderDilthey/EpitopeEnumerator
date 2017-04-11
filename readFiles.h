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
#include <set>

class variantFromVCF {
public:
	std::string chromosomeID;
	unsigned int position;
	std::string referenceString;
	std::vector<std::string> sampleAlleles;
	std::vector<bool> sampleAlleles_interesting;
	bool isRef() const;
	void allAllelesNotInteresting();

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
	void print() const;
};

std::map<std::string, std::map<int, variantFromVCF>> readVariants(std::string VCF, const std::map<std::string, std::string>& referenceGenome);
std::vector<transcript> readTranscripts(std::string transcriptsFile);

std::vector<transcript> getPlusStrandTranscripts(const std::vector<transcript>& transcripts);
std::vector<transcript> getMinusStrandTranscripts(const std::vector<transcript>& transcripts, const std::map<std::string, std::string>& referenceGenome_minus);
std::map<std::string, std::string> getMinusStrandReferenceGenome(const std::map<std::string, std::string>& referenceGenome);
std::map<std::string, std::map<int, variantFromVCF>> getMinusStrandVariants(const std::map<std::string, std::map<int, variantFromVCF>>& variants, const std::map<std::string, std::string>& referenceGenome_minus);
void extendAsNecessary(std::string& S, unsigned int desiredLength);

void checkVariantsConsistentWithReferenceGenome(const std::map<std::string, std::map<int, variantFromVCF>>& variants, const std::map<std::string, std::string>& referenceGenome);
void checkTranscriptsTranslate(const std::vector<transcript>& transcripts, const std::map<std::string, std::string>& referenceGenome);

std::map<std::string, std::map<int, variantFromVCF>> combineVariants(const std::map<std::string, std::map<int, variantFromVCF>>& variants_normalGenome, const std::map<std::string, std::map<int, variantFromVCF>>& additionalVariants_tumourGenome, const std::map<std::string, std::string>& referenceGenome);

#endif /* READFILES_H_ */
