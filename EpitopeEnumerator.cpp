//============================================================================
// Name        : EpitopeEnumerator.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include <utility>
#include <algorithm>

#include "Utilities.h"
#include "readFiles.h"

using namespace std;

std::vector<transcript> getPlusStrandTranscripts(const std::vector<transcript>& transcripts);
void enumeratePeptides(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus);
std::pair<std::string, std::vector<std::string>> get_reference_and_variantAlleles(const variantFromVCF& v, unsigned int lastReferencePos);

int main(int argc, char *argv[]) {

	std::vector<std::string> ARG (argv + 1, argv + argc + !argc);
	std::map<std::string, std::string> arguments;

	arguments["referenceGenome"] = "data/GRCh38_full_analysis_set_plus_decoy_hla.fa.chr20";
	arguments["normalVCF"] = "data/NA12878.vcf.chr21";
	arguments["transcripts"] = "data/gencode.v26.annotation.gff3";

	for(unsigned int i = 0; i < ARG.size(); i++)
	{
		if((ARG.at(i).length() > 2) && (ARG.at(i).substr(0, 2) == "--"))
		{
			std::string argname = ARG.at(i).substr(2);
			std::string argvalue = ARG.at(i+1);
			arguments[argname] = argvalue;
		}
	}

	std::vector<transcript> transcripts = readTranscripts(arguments.at("transcripts"));

	assert( 4 == 5 );

	std::map<std::string, std::map<int, variantFromVCF>> variants = readVariants(arguments.at("normalVCF"));

	assert( 2 == 4 );

	std::map<std::string, std::string> referenceGenome = Utilities::readFASTA(arguments.at("referenceGenome"));

	std::vector<transcript> transcripts_plus = getPlusStrandTranscripts(transcripts);

	enumeratePeptides(referenceGenome, transcripts_plus, variants);

	assert("Are the GFF stop coordinates inclusive?" == "");

	return 0;
}

void enumeratePeptides(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus)
{
	class runningHaplotype {
	public:
		int fraction;

		std::string AApart;
		std::vector<bool> AApart_interesting;

		std::string nucleotidePart;
		std::vector<bool> nucleotidePart_interesting;
		int nucleotidePart_nonGap;
		bool canExtendFurther;

		runningHaplotype()
		{
			fraction = 1;
			nucleotidePart_nonGap = 0;
			canExtendFurther = true;
		}

		void extendWithNucleotides(const std::string& nucleotides, const std::vector<bool> interesting)
		{
			assert(nucleotides.size() == interesting.size());
			assert(nucleotidePart.size() == nucleotidePart_interesting.size());
			assert(countCharacters_noGaps(nucleotides) == nucleotidePart_nonGap);

			nucleotidePart.insert(nucleotidePart.end(), nucleotides.begin(), nucleotides.end());
			nucleotidePart_interesting.insert(nucleotidePart_interesting.end(), interesting.begin(), interesting.end());

			int added_nonGaps = countCharacters_noGaps(nucleotides);
			nucleotidePart_nonGap += added_nonGaps;

			while(nucleotidePart_nonGap >= 3)
			{
				std::string codon; codon.reserve(3);
				int consumedIndex_inNucleotidePart = 0;
				while(codon.length() < 3)
				{
					unsigned char charForConsumption = nucleotidePart.at(consumedIndex_inNucleotidePart);
					consumedIndex_inNucleotidePart++;
					if((charForConsumption != '-') && (charForConsumption != '_'))
					{
						codon.push_back(charForConsumption);
					}
				}
				nucleotidePart = nucleotidePart.substr(consumedIndex_inNucleotidePart);
				nucleotidePart_interesting = std::vector<bool>(nucleotidePart_interesting.begin()+consumedIndex_inNucleotidePart, nucleotidePart_interesting.end());
				nucleotidePart_nonGap -= 3;

				assert(nucleotidePart.size() == nucleotidePart_interesting.size());
				assert(countCharacters_noGaps(nucleotides) == nucleotidePart_nonGap);
			}

		}
	};

	for(unsigned int transcriptI = 0; transcriptI < transcripts_plus.size(); transcriptI++)
	{
		const transcript& transcript = transcripts_plus.at(transcriptI);
		assert(transcript.strand == '+');
		std::string chromosomeID = transcript.chromosomeID;
		if(referenceGenome_plus.count(chromosomeID) == 0)
			continue;

		int exons = transcript.exons.size();

		// to make sure that the exons go from left to right
		for(int exonI = 1; exonI < exons; exonI++)
		{
			assert(transcript.exons.at(exonI-1).lastPos < transcript.exons.at(exonI-1).firstPos);
		}

		runningHaplotype referenceHaplotype;
		runningHaplotype firstHaplotype;
		std::vector<runningHaplotype> sampleHaplotypes = {firstHaplotype};


		auto doExtension = [&](std::string referenceSequence, std::vector<std::string> sampleAlleles) -> void {
			referenceHaplotype.extendWithNucleotides(referenceSequence);
			size_t existingSampleHaplotypes_maxI = sampleHaplotypes.size();
			if(sampleAlleles.size() == 1)
			{
				for(unsigned int existingHaplotypeI = 0; existingHaplotypeI < existingSampleHaplotypes_maxI; existingHaplotypeI++)
				{
					sampleHaplotypes.at(existingHaplotypeI).extendWithNucleotides(sampleAlleles.at(0));
				}
			}
			else
			{
				for(unsigned int existingHaplotypeI = 0; existingHaplotypeI < existingSampleHaplotypes_maxI; existingHaplotypeI++)
				{
					runningHaplotype preExtension = sampleHaplotypes.at(existingHaplotypeI);
					sampleHaplotypes.at(existingHaplotypeI).extendWithNucleotides(sampleAlleles.at(0));
					for(unsigned int sampleAlleleI = 1; sampleAlleleI < sampleAlleles.size(); sampleAlleleI++)
					{
						sampleHaplotypes.push_back(preExtension);
						sampleHaplotypes.back().extendWithNucleotides(sampleAlleles.at(sampleAlleleI));
					}
				}
			}
		};

		//std::string collectedReferenceSequence;
		//std::vector<> runningSampleHaplotypes;

		for(int exonI = 0; exonI < exons; exonI++)
		{
			const transcriptExon& exon = transcript.exons.at(exonI);
			for(unsigned int referencePos = exon.firstPos; referencePos <= exon.lastPos; referencePos++)
			{
				if(variants_plus.count(chromosomeID) && variants_plus.at(chromosomeID).count(referencePos))
				{
					std::pair<std::string, std::vector<std::string>> reference_and_variantAlleles = get_reference_and_variantAlleles(variants_plus.at(chromosomeID).at(referencePos), exon.lastPos);

					doExtension(reference_and_variantAlleles.first, reference_and_variantAlleles.second);

					int reference_extension_length_noGaps = countCharacters_noGaps(reference_and_variantAlleles.first);
					referencePos += (reference_extension_length_noGaps - 1);
				}
				else
				{
					std::string referenceCharacter = referenceGenome_plus.at(chromosomeID).substr(referencePos, 1);
					std::vector<std::string> extendWith = {referenceCharacter};

					doExtension(referenceCharacter, extendWith);
				}
			}
		}

	}
}

std::pair<std::string, std::vector<std::string>> get_reference_and_variantAlleles(const variantFromVCF& v, unsigned int lastReferencePos)
{
	assert(v.position <= lastReferencePos);

	// check that we have variant alleles and that they are "aligned"
	assert(v.sampleAlleles.size());
	for(const std::string& variantAllele : v.sampleAlleles)
	{
		assert(variantAllele.length() == v.referenceString.length());
	}

	// if everything within last-position constraints, do nothing ...
	unsigned int v_referenceString_length_noGaps = countCharacters_noGaps(v.referenceString);
	if((v.position + v_referenceString_length_noGaps - 1) <= lastReferencePos)
	{
		return make_pair(v.referenceString, v.sampleAlleles);
	}
	// else: cut variant alleles accordingly
	else
	{
		std::string forReturn_ref;
		std::vector<std::string> forReturn_alleles;
		forReturn_alleles.resize(v.sampleAlleles.size());

		int consumePositions = lastReferencePos - v.position + 1;
		assert(consumePositions >= 1);
		int consumedPositions = 0;
		int consumeI = 0;
		while(consumedPositions < consumePositions)
		{
			unsigned char refC = v.referenceString.at(consumeI);

			forReturn_ref.push_back(refC);
			for(unsigned int alternativeAlleleI = 0; alternativeAlleleI < v.sampleAlleles.size(); alternativeAlleleI++)
			{
				forReturn_alleles.at(alternativeAlleleI).push_back(v.sampleAlleles.at(alternativeAlleleI).at(consumeI));
			}
			if((refC != '-') && (refC != '_'))
			{
				consumedPositions++;
			}
			consumeI++;
		}

		for(const std::string& alternativeAllele : forReturn_alleles)
		{
			assert(alternativeAllele.size() == forReturn_ref.length());
		}

		return make_pair(forReturn_ref, forReturn_alleles);
	}
}

std::vector<transcript> getPlusStrandTranscripts(const std::vector<transcript>& transcripts)
{
	std::vector<transcript> forReturn;
	forReturn.reserve(transcripts.size());
	for(transcript t : transcripts)
	{
		if(t.strand == '+')
		{
			forReturn.push_back(t);
		}
	}

	std::cout << "getPlusStrandTranscripts(..): " << forReturn.size() << " plus-strand transcripts.\n" << std::flush;

	return forReturn;
}



