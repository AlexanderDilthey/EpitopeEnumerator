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
void enumeratePeptides(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool isTumour);
std::pair<std::string, std::vector<std::string>> get_reference_and_variantAlleles(const variantFromVCF& v, unsigned int lastReferencePos);

int main(int argc, char *argv[]) {

	std::vector<std::string> ARG (argv + 1, argv + argc + !argc);
	std::map<std::string, std::string> arguments;

	arguments["referenceGenome"] = "data/GRCh38_full_analysis_set_plus_decoy_hla.fa.chr20";
	arguments["normalVCF"] = "data/NA12878.vcf.chr20";
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

	std::map<std::string, std::map<int, variantFromVCF>> variants = readVariants(arguments.at("normalVCF"));

	std::map<std::string, std::string> referenceGenome = Utilities::readFASTA(arguments.at("referenceGenome"));

	std::vector<transcript> transcripts_plus = getPlusStrandTranscripts(transcripts);

	for(const transcript& t : transcripts_plus)
	{
		std::string t_referenceSequence;
		for(const transcriptExon& e : t.exons)
		{
			if(e.valid)
			{
				if(referenceGenome.count(t.chromosomeID) == 0)
				{
					std::cerr << "Unknown chromosome ID: "+t.chromosomeID << "\n" << std::flush;
					throw std::runtime_error("Unknown chromosome ID: "+t.chromosomeID);
				}
				const std::string& chromosomeSequence = referenceGenome.at(t.chromosomeID);
				assert(e.firstPos <= e.lastPos);
				assert(e.firstPos >= 0);
				assert(e.lastPos < chromosomeSequence.length());

				t_referenceSequence += referenceGenome.at(t.chromosomeID).substr(e.firstPos, e.lastPos - e.firstPos + 1);
			}
		}
		std::string translation;
		assert((t_referenceSequence.length() % 3) == 0);
		for(unsigned int i = 0; i < t_referenceSequence.length(); i += 3)
		{
			std::string codon = t_referenceSequence.substr(i, 3);
			assert(codon.length() == 3);
			translation += translateCodon2AA(codon);
		}
		if(translation.back() != '!')
		{
			// std::cout << t.transcriptID << " " << translation << "\n" << std::flush;
		}
		// assert(translation.back() == '!');
	}

	// assert("Are the GFF stop coordinates inclusive?" == "");


	enumeratePeptides(referenceGenome, transcripts_plus, variants, true);


	return 0;
}

void enumeratePeptides(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool isTumour)
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

		void extendWithNucleotides(const std::string& nucleotides, const std::vector<bool>& interesting)
		{
			assert(nucleotides.size() == interesting.size());
			assert(nucleotidePart.size() == nucleotidePart_interesting.size());
			assert(countCharacters_noGaps(nucleotidePart) == nucleotidePart_nonGap); // paranoid

			nucleotidePart.insert(nucleotidePart.end(), nucleotides.begin(), nucleotides.end());
			nucleotidePart_interesting.insert(nucleotidePart_interesting.end(), interesting.begin(), interesting.end());

			int added_nonGaps = countCharacters_noGaps(nucleotides);
			nucleotidePart_nonGap += added_nonGaps;

			while(nucleotidePart_nonGap >= 3)
			{
				std::string codon; codon.reserve(3);
				bool codon_interesting = false;
				int consumedIndex_inNucleotidePart = 0;
				while(codon.length() < 3)
				{
					unsigned char charForConsumption = nucleotidePart.at(consumedIndex_inNucleotidePart);
					codon_interesting = (codon_interesting || nucleotidePart_interesting.at(consumedIndex_inNucleotidePart));
					consumedIndex_inNucleotidePart++;
					if((charForConsumption != '-') && (charForConsumption != '_'))
					{
						codon.push_back(charForConsumption);
					}
					assert(consumedIndex_inNucleotidePart <= nucleotidePart.size());
				}

				// translate codon -> AA part
				AApart.append(translateCodon2AA(codon));
				AApart_interesting.push_back(codon_interesting);

				// translate nucleotide interesting -> AA
				bool haveRemainder = (consumedIndex_inNucleotidePart < (int)nucleotidePart.size());
				nucleotidePart = haveRemainder ? nucleotidePart.substr(consumedIndex_inNucleotidePart) : "";
				nucleotidePart_interesting = haveRemainder ? std::vector<bool>(nucleotidePart_interesting.begin()+consumedIndex_inNucleotidePart, nucleotidePart_interesting.end()) : std::vector<bool>();
				nucleotidePart_nonGap -= 3;

				assert(nucleotidePart.size() == nucleotidePart_interesting.size());
				//std::cout << nucleotides << "\n";
				//std::cout << nucleotidePart_nonGap << "\n\n" << std::flush;
				assert((int)countCharacters_noGaps(nucleotidePart) == nucleotidePart_nonGap); // paranoid
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
			if(transcript.exons.at(exonI-1).valid)
				assert(transcript.exons.at(exonI-1).firstPos <= transcript.exons.at(exonI-1).lastPos);
		}

		runningHaplotype referenceHaplotype;
		runningHaplotype firstHaplotype;
		std::vector<runningHaplotype> sampleHaplotypes = {firstHaplotype};


		auto doExtension = [&](std::string referenceSequence, std::vector<std::string> sampleAlleles, std::vector<bool> sampleAlleles_interesting) -> void {
			assert(sampleAlleles.size() == sampleAlleles_interesting.size());
			std::vector<bool> referenceSequence_interesting;
			referenceSequence_interesting.resize(referenceSequence.length(), false);
			referenceHaplotype.extendWithNucleotides(referenceSequence, referenceSequence_interesting);
			size_t existingSampleHaplotypes_maxI = sampleHaplotypes.size();

			std::vector<std::vector<bool>> sampleAlleles_perCharacter_interesting;
			sampleAlleles_perCharacter_interesting.reserve(sampleAlleles.size());
			for(const std::string& oneSampleAllele : sampleAlleles)
			{
				std::vector<bool> oneSampleAllele_interesting;
				oneSampleAllele_interesting.resize(oneSampleAllele.length(), false);
				sampleAlleles_perCharacter_interesting.push_back(oneSampleAllele_interesting);
			}
			if(sampleAlleles.size() == 1)
			{
				for(unsigned int existingHaplotypeI = 0; existingHaplotypeI < existingSampleHaplotypes_maxI; existingHaplotypeI++)
				{
					sampleHaplotypes.at(existingHaplotypeI).extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_perCharacter_interesting.at(0));
				}
			}
			else
			{
				std::cout << "!" << std::flush;
				for(unsigned int existingHaplotypeI = 0; existingHaplotypeI < existingSampleHaplotypes_maxI; existingHaplotypeI++)
				{
					runningHaplotype preExtension = sampleHaplotypes.at(existingHaplotypeI);
					sampleHaplotypes.at(existingHaplotypeI).extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_perCharacter_interesting.at(0));
					for(unsigned int sampleAlleleI = 1; sampleAlleleI < sampleAlleles.size(); sampleAlleleI++)
					{
						sampleHaplotypes.push_back(preExtension);
						sampleHaplotypes.back().extendWithNucleotides(sampleAlleles.at(sampleAlleleI), sampleAlleles_perCharacter_interesting.at(sampleAlleleI));
					}
				}
			}

			if(sampleHaplotypes.size() > 1)
				std::cout << "sampleHaplotypes.size()" << ": " << sampleHaplotypes.size() << "\n" << std::flush;
		};

		//std::string collectedReferenceSequence;
		//std::vector<> runningSampleHaplotypes;

		for(int exonI = 0; exonI < exons; exonI++)
		{
			const transcriptExon& exon = transcript.exons.at(exonI);
			if(exon.valid == false)
				continue;

			for(unsigned int referencePos = exon.firstPos; referencePos <= exon.lastPos; referencePos++)
			{
				if(variants_plus.count(chromosomeID) && variants_plus.at(chromosomeID).count(referencePos))
				{
					std::pair<std::string, std::vector<std::string>> reference_and_variantAlleles = get_reference_and_variantAlleles(variants_plus.at(chromosomeID).at(referencePos), exon.lastPos);

					std::vector<bool> extendWith_interesting;
					if(isTumour)
					{
						extendWith_interesting.reserve(reference_and_variantAlleles.second.size());
						for(const std::string& sampleAllele : reference_and_variantAlleles.second)
						{
							extendWith_interesting.push_back(!(sampleAllele == reference_and_variantAlleles.first));
						}
					}
					else
					{
						extendWith_interesting.resize(reference_and_variantAlleles.second.size(), false);
					}

					std::cout << "V: " << reference_and_variantAlleles.second.size() << "\n" << std::flush;

					doExtension(reference_and_variantAlleles.first, reference_and_variantAlleles.second, extendWith_interesting);

					int reference_extension_length_noGaps = countCharacters_noGaps(reference_and_variantAlleles.first);
					referencePos += (reference_extension_length_noGaps - 1);
				}
				else
				{
					std::string referenceCharacter = referenceGenome_plus.at(chromosomeID).substr(referencePos, 1);
					std::vector<std::string> extendWith = {referenceCharacter};
					std::vector<bool> extendWith_interesting = {false};
					doExtension(referenceCharacter, extendWith, extendWith_interesting);
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



