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
#include "EpitopeEnumerator.h"

#include "tests.h"

using namespace std;

int main(int argc, char *argv[]) {

	std::vector<std::string> ARG (argv + 1, argv + argc + !argc);
	std::map<std::string, std::string> arguments;

	int testing = 2;
	if(testing == 1)
	{
		some_simple_tests();
	}
	else if(testing == 2)
	{
		randomTests_withVariants_2();
		randomTests();
		randomTests_withVariants();
	}
	else
	{
		arguments["referenceGenome"] = "data/GRCh38_full_analysis_set_plus_decoy_hla.fa.chr20";
		arguments["normalVCF"] = "data/NA12878.vcf.chr20";
		arguments["transcripts"] = "data/gencode.v26.annotation.gff3";

		/*
		std::map<std::pair<int, int>, int> testP;
		testP[make_pair(1,2)] = 4;
		assert(testP.count(make_pair(1,2)));
		assert(! testP.count(make_pair(1,3)));
		std::cout << testP[make_pair(1,2)] << "\n" << std::flush;
		assert(2 == 10);
		 */

		for(unsigned int i = 0; i < ARG.size(); i++)
		{
			if((ARG.at(i).length() > 2) && (ARG.at(i).substr(0, 2) == "--"))
			{
				std::string argname = ARG.at(i).substr(2);
				std::string argvalue = ARG.at(i+1);
				arguments[argname] = argvalue;
			}
		}

		std::map<std::string, std::string> referenceGenome = Utilities::readFASTA(arguments.at("referenceGenome"));
		std::vector<transcript> transcripts = readTranscripts(arguments.at("transcripts"));
		std::map<std::string, std::map<int, variantFromVCF>> variants = readVariants(arguments.at("normalVCF"), referenceGenome);

		std::map<std::string, std::map<int, variantFromVCF>> variants_tumour;
		std::map<std::string, std::map<int, variantFromVCF>> variants_combined = combineVariants(variants, variants_tumour, referenceGenome);

		std::map<int, std::map<std::string, std::map<double, std::set<std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>> haplotypeStore;
		haplotypeStore[6] = std::map<std::string, std::map<double, std::set<std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>();
		//enumeratePeptideHaplotypes(referenceGenome, transcripts_plus, variants, true, haplotypeStore);

		std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
		std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;
		enumeratePeptideHaplotypes(referenceGenome, transcripts, variants, {6}, p_per_epitope, p_per_epitope_locations);

		for(auto haplotypesOneLength : p_per_epitope)
		{
			std::cout << "Length " << haplotypesOneLength.first << ": " << haplotypesOneLength.second.size() << "\n" << std::flush;
			for(auto fragment : haplotypesOneLength.second)
			{
				const std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>& probability_and_positions = fragment.second;
				double p = probability_and_positions.first;
				std::cout << "\t" << fragment.first << " " << p << "\n" << std::flush;
			}
		}
		assert(1 == 2);

		assert("minus-strand transcripts!" == "");
		assert("add stop codons!" == "");
		assert("add filter for PASS" == "");
	}

	return 0;
}

void enumeratePeptideHaplotypes(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet)
{
	std::map<std::string, std::string> referenceGenome_minus = getMinusStrandReferenceGenome(referenceGenome);
	std::map<std::string, std::map<int, variantFromVCF>> variants_minus = getMinusStrandVariants(variants, referenceGenome_minus);

	checkVariantsConsistentWithReferenceGenome(variants, referenceGenome);
	checkVariantsConsistentWithReferenceGenome(variants_minus, referenceGenome_minus);

	std::vector<transcript> transcripts_plus = getPlusStrandTranscripts(transcripts);
	std::vector<transcript> transcripts_minus = getMinusStrandTranscripts(transcripts, referenceGenome_minus);

	checkTranscriptsTranslate(transcripts_plus, referenceGenome);
	checkTranscriptsTranslate(transcripts_minus, referenceGenome_minus);

	p_per_epitope_forRet.clear();
	p_per_epitope_locations_forRet.clear();
	for(int k : haplotypeLengths)
	{
		p_per_epitope_forRet[k].count("");
		p_per_epitope_locations_forRet[k].count("");
	}

	enumeratePeptideHaplotypes_plus(referenceGenome, transcripts_plus, variants, false, p_per_epitope_forRet, p_per_epitope_locations_forRet);
	enumeratePeptideHaplotypes_plus(referenceGenome_minus, transcripts_minus, variants_minus, true, p_per_epitope_forRet, p_per_epitope_locations_forRet);
}

// this function returns:
// p_per_epitope_forRet: for each occurring epitope, the maximum probability of occurrence over all transcripts and ALL possible locations (without probabilities)
// p_per_epitope_locations_forRet: for each occurring epitope, possible locations and their marginal maximum probabilities over all transcripts
void enumeratePeptideHaplotypes_plus(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool invertPositions, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet)
{
	checkVariantsConsistentWithReferenceGenome(variants_plus, referenceGenome_plus); // paranoid
	assert(p_per_epitope_forRet.size() == p_per_epitope_locations_forRet.size());

	std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> oneTranscript_p_per_epitope_forRet;
	std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> oneTranscript_p_per_epitope_locations_forRet;
	for(auto k : p_per_epitope_forRet)
	{
		assert(p_per_epitope_locations_forRet.count(k.first));
		oneTranscript_p_per_epitope_forRet[k.first].count("");
		oneTranscript_p_per_epitope_locations_forRet[k.first].count("");
	}
	for(unsigned int transcriptI = 0; transcriptI < transcripts_plus.size(); transcriptI++)
	{
		const transcript& transcript = transcripts_plus.at(transcriptI);
		assert(transcript.strand == '+');
		std::string chromosomeID = transcript.chromosomeID;
		if(referenceGenome_plus.count(chromosomeID) == 0)
			continue;

		// the oneTranscript_* maps are clear'ed within enumeratePeptideHaplotypes_oneTranscript
		enumeratePeptideHaplotypes_oneTranscript(transcript, referenceGenome_plus, variants_plus, oneTranscript_p_per_epitope_forRet, oneTranscript_p_per_epitope_locations_forRet);

		for(auto k : p_per_epitope_forRet)
		{
			for(auto epitope : oneTranscript_p_per_epitope_forRet.at(k.first))
			{
				std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>> positions = epitope.second.second;

				if(invertPositions)
				{
					std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>> inverted_positions;
					for(auto p : positions)
					{
						std::pair<std::vector<std::pair<int, int>>, std::vector<bool>> inverted_p = p;
						for(std::pair<int, int>& pP : inverted_p.first)
						{
							size_t referenceContig_length = referenceGenome_plus.at(chromosomeID).length();
							pP.first = referenceContig_length - pP.first - 1;
							pP.second = referenceContig_length - pP.second - 1;
						}
						inverted_positions.insert(p);
					}
					positions = inverted_positions;
				}

				if(p_per_epitope_forRet.at(k.first).count(epitope.first))
				{
					if(p_per_epitope_forRet.at(k.first).at(epitope.first).first < epitope.second.first)
					{
						p_per_epitope_forRet.at(k.first).at(epitope.first).first = epitope.second.first;
					}
					p_per_epitope_forRet.at(k.first).at(epitope.first).second.insert(positions.begin(), positions.end());
				}
				else
				{
					p_per_epitope_forRet.at(k.first)[epitope.first].first = epitope.second.first;
					p_per_epitope_forRet.at(k.first)[epitope.first].second = positions;
				}
			}
			for(auto epitope : oneTranscript_p_per_epitope_locations_forRet.at(k.first))
			{
				for(auto location : epitope.second)
				{
					if(p_per_epitope_locations_forRet.at(k.first)[epitope.first].count(location.first))
					{
						if(p_per_epitope_locations_forRet.at(k.first).at(epitope.first).at(location.first) < location.second)
						{
							p_per_epitope_locations_forRet.at(k.first).at(epitope.first).at(location.first) = location.second;
						}
					}
					else
					{
						p_per_epitope_locations_forRet.at(k.first)[epitope.first][location.first] = location.second;
					}
				}
			}
		}
	}
}

// this function returns:
// p_per_epitope_forRet: for each occurring epitope, the probability of occurrence and possible locations (without probabilities)
// p_per_epitope_locations_forRet: for each occurring epitope, possible locations and their marginal probabilities
// It would also be possible to return, for each epitope, the probabilities over sets of locations, but this doesn't currently seem very useful
void enumeratePeptideHaplotypes_oneTranscript(const transcript& transcript, const std::map<std::string, std::string> referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet)
{
	// validate exons - left to right, non-overlapping
	int n_exons = transcript.exons.size();
	// make sure that the exons go from left to right in a non-overlapping fashion
	int lastValidExonI = -1;
	for(int exonI = 0; exonI < n_exons; exonI++)
	{
		if(transcript.exons.at(exonI).valid)
		{
			assert(transcript.exons.at(exonI).firstPos <= transcript.exons.at(exonI).lastPos);
			if(!(transcript.exons.at(exonI).lastPos < referenceGenome_plus.at(transcript.chromosomeID).length()))
			{
				std::cerr << "transcript.exons.at(exonI).lastPos" << ": " << transcript.exons.at(exonI).lastPos << "\n";
				std::cerr << "referenceGenome_plus.at(transcript.chromosomeID).length()" << ": " << referenceGenome_plus.at(transcript.chromosomeID).length() << "\n";
				std::cerr << std::flush;
			}
			assert(transcript.exons.at(exonI).lastPos < referenceGenome_plus.at(transcript.chromosomeID).length());
			if(lastValidExonI != -1)
			{
				assert(transcript.exons.at(lastValidExonI).lastPos < transcript.exons.at(exonI).firstPos);
			}
			lastValidExonI = exonI;
		}
	}

	// build sequenceFragments array
	std::string referenceSequence_forControl;
	std::vector<int> referenceSequence_forControl_referencePositions;
	std::vector<std::vector<std::string>> sequenceFragments;
	std::vector<std::vector<std::vector<int>>> sequenceFragments_referencePositions;
	std::vector<std::vector<std::vector<bool>>> sequenceFragments_interesting;

	auto addHomozygousStretch = [&](const std::string& sequence, const std::vector<int> referencePositions, const std::vector<bool>& interesting) -> void {
		assert(sequence.size() == referencePositions.size());
		assert(sequence.size() == interesting.size());

		if(sequenceFragments.size() && (sequenceFragments.back().size() == 1))
		{
			sequenceFragments.back().at(0).append(sequence);
			sequenceFragments_referencePositions.back().at(0).insert(sequenceFragments_referencePositions.back().at(0).end(), referencePositions.begin(), referencePositions.end());
			sequenceFragments_interesting.back().at(0).insert(sequenceFragments_interesting.back().at(0).end(), interesting.begin(), interesting.end());
		}
		else
		{
			std::vector<std::string> homStretch = {sequence};
			std::vector<std::vector<int>> homStretch_coordinates({referencePositions});
			std::vector<std::vector<bool>> homStretch_interesting({interesting});

			sequenceFragments.push_back(homStretch);
			sequenceFragments_referencePositions.push_back(homStretch_coordinates);
			sequenceFragments_interesting.push_back(homStretch_interesting);
		}
	};

	auto addHeterozygousStretch = [&](const std::vector<std::string>& sequences, const std::vector<std::vector<int>> referencePositions, const std::vector<std::vector<bool>>& interesting) -> void {

		assert(sequences.size() == 2);
		assert(sequences.size() == referencePositions.size());
		assert(sequences.size() == interesting.size());
		for(unsigned int sI = 0; sI < sequences.size(); sI++)
		{
			assert(sequences.at(sI).size() == sequences.at(0).size());
			assert(sequences.at(sI).size() == referencePositions.at(sI).size());
			assert(sequences.at(sI).size() == interesting.at(sI).size());
		}

		sequenceFragments.push_back(sequences);
		sequenceFragments_referencePositions.push_back(referencePositions);
		sequenceFragments_interesting.push_back(interesting);
	};

	const std::string& chromosomeID = transcript.chromosomeID;
	for(int exonI = 0; exonI < n_exons; exonI++)
	{
		const transcriptExon& exon = transcript.exons.at(exonI);
		if(exon.valid == false)
			continue;

		for(unsigned int referencePos = exon.firstPos; referencePos <= exon.lastPos; referencePos++)
		{
			if(variants_plus.count(chromosomeID) && variants_plus.at(chromosomeID).count(referencePos))
			{
				std::tuple<std::string, std::vector<int>, std::vector<bool>, std::vector<std::string>, std::vector<std::vector<int>>, std::vector<std::vector<bool>>> reference_and_variantAlleles = get_reference_and_variantAlleles(variants_plus.at(chromosomeID).at(referencePos), referencePos, exon.lastPos);

				assert(std::get<3>(reference_and_variantAlleles).size() == 2);

				bool homozygous =
						(std::get<3>(reference_and_variantAlleles).at(0) == std::get<3>(reference_and_variantAlleles).at(1)) &&
						(std::get<4>(reference_and_variantAlleles).at(0) == std::get<4>(reference_and_variantAlleles).at(1)) &&
						(std::get<5>(reference_and_variantAlleles).at(0) == std::get<5>(reference_and_variantAlleles).at(1));

				// this is not really very interesting
				/*
				std::vector<std::vector<bool>> extendWith_interesting;
				for(const std::string& sampleAllele : std::get<3>(reference_and_variantAlleles))
				{
					if(isTumour)
					{
						std::vector<bool> potentiallyInteresting;
						potentiallyInteresting.resize(sampleAllele.size(), sampleAllele != std::get<0>(reference_and_variantAlleles));
						extendWith_interesting.push_back(potentiallyInteresting);
					}
					else
					{
						std::vector<bool> notInteresting;
						notInteresting.resize(sampleAllele.size(), false);
						extendWith_interesting.push_back(notInteresting);
					}
				}
				*/

				if(homozygous)
				{
					addHomozygousStretch(std::get<3>(reference_and_variantAlleles).at(0), std::get<4>(reference_and_variantAlleles).at(0), std::get<5>(reference_and_variantAlleles).at(0));
				}
				else
				{
					addHeterozygousStretch(std::get<3>(reference_and_variantAlleles), std::get<4>(reference_and_variantAlleles), std::get<5>(reference_and_variantAlleles));
				}

				referenceSequence_forControl.append(std::get<0>(reference_and_variantAlleles));
				referenceSequence_forControl_referencePositions.insert(referenceSequence_forControl_referencePositions.end(), std::get<1>(reference_and_variantAlleles).begin(), std::get<1>(reference_and_variantAlleles).end());
				int reference_extension_length_noGaps = countCharacters_noGaps(std::get<0>(reference_and_variantAlleles));
				referencePos += (reference_extension_length_noGaps - 1);
			}
			else
			{
				std::string referenceCharacter = referenceGenome_plus.at(chromosomeID).substr(referencePos, 1);
				addHomozygousStretch(referenceCharacter, std::vector<int>({(int)referencePos}), std::vector<bool>({false}));
				referenceSequence_forControl.append(referenceCharacter);
				referenceSequence_forControl_referencePositions.push_back(referencePos);

			}
		}
	}

	// sanity checks
	size_t totalLength = 0;
	int n_fragments_diploid = 0;
	int lastStretchStop = -1;
	std::vector<std::pair<unsigned int, unsigned int>> stretchBoundaries;
	int first_heterozygous_position = -2;
	for(unsigned int stretchI = 0; stretchI < sequenceFragments.size(); stretchI++)
	{
		if(stretchI > 0)
		{
			if(sequenceFragments.at(stretchI).size() == 1)
			{
				assert(sequenceFragments.at(stretchI-1).size() != 1);
			}
		}

		assert((sequenceFragments.at(stretchI).size() == 2) || (sequenceFragments.at(stretchI).size() == 1));
		assert(sequenceFragments.at(stretchI).size() == sequenceFragments_referencePositions.at(stretchI).size());
		assert(sequenceFragments.at(stretchI).size() == sequenceFragments_interesting.at(stretchI).size());
		for(unsigned int sI = 0; sI < sequenceFragments.at(stretchI).size(); sI++)
		{
			assert(sequenceFragments.at(stretchI).at(sI).size() == sequenceFragments_referencePositions.at(stretchI).at(sI).size());
			assert(sequenceFragments.at(stretchI).at(sI).size() == sequenceFragments_interesting.at(stretchI).at(sI).size());
		}
		totalLength += sequenceFragments.at(stretchI).at(0).size();

		if(sequenceFragments.at(stretchI).size() == 2)
		{
			n_fragments_diploid++;
			if(first_heterozygous_position == -2)
			{
				first_heterozygous_position = stretchI;
			}
		}

		int thisStretchStop = lastStretchStop + 1 + sequenceFragments.at(stretchI).at(0).size() - 1;
		stretchBoundaries.push_back(make_pair(lastStretchStop+1,thisStretchStop));
		lastStretchStop = thisStretchStop;
	}
	assert(referenceSequence_forControl.length() == totalLength);
	assert(referenceSequence_forControl_referencePositions.size() == totalLength);
	assert(stretchBoundaries.back().second == (totalLength-1));

	size_t n_haplotype_pairs = 1;
	if(n_fragments_diploid > 1)
	{
		n_haplotype_pairs = std::pow(2, n_fragments_diploid-1);
	}

	std::vector<int> utilizingIndex;
	utilizingIndex.resize(sequenceFragments.size(), 0);
	std::string nucleotideHaplotype_1;
	std::string nucleotideHaplotype_2;
	nucleotideHaplotype_1.resize(totalLength, '#');
	nucleotideHaplotype_2.resize(totalLength, '#');
	std::vector<int> nucleotideHaplotype_1_positions;
	std::vector<int> nucleotideHaplotype_2_positions;
	std::vector<bool> nucleotideHaplotype_1_interesting;
	std::vector<bool> nucleotideHaplotype_2_interesting;
	nucleotideHaplotype_1_positions.resize(totalLength);
	nucleotideHaplotype_2_positions.resize(totalLength);
	nucleotideHaplotype_1_interesting.resize(totalLength);
	nucleotideHaplotype_2_interesting.resize(totalLength);

	auto populateNucleotideHaplotypes = [&](std::string& nucleotideHaplotype_1, std::string& nucleotideHaplotype_2, std::vector<int>& nucleotideHaplotype_1_positions, std::vector<int>& nucleotideHaplotype_2_positions, std::vector<bool>& nucleotideHaplotype_1_interesting, std::vector<bool>& nucleotideHaplotype_2_interesting, const std::vector<int>& utilizingIndex) -> void {
		if(first_heterozygous_position >= 0)
			assert(utilizingIndex.at(first_heterozygous_position) == 0);

		for(unsigned int stretchI = 0; stretchI < utilizingIndex.size(); stretchI++)
		{
			const std::pair<unsigned int, unsigned int>& boundaries = stretchBoundaries.at(stretchI);
			unsigned int L =  boundaries.second - boundaries.first + 1;

			int index_for_1;
			int index_for_2;
			if(utilizingIndex.at(stretchI) == 0)
			{
				if(sequenceFragments.at(stretchI).size() == 1)
				{
					index_for_1 = 0;
					index_for_2 = 0;
				}
				else
				{
					index_for_1 = 0;
					index_for_2 = 1;
				}
			}
			else
			{
				index_for_1 = 1;
				index_for_2 = 0;
			}

			const std::string& sequence_for_1 = sequenceFragments.at(stretchI).at(index_for_1);
			const std::string& sequence_for_2 = sequenceFragments.at(stretchI).at(index_for_2);
			const std::vector<int>& positions_for_1 = sequenceFragments_referencePositions.at(stretchI).at(index_for_1);
			const std::vector<int>& positions_for_2 = sequenceFragments_referencePositions.at(stretchI).at(index_for_2);
			const std::vector<bool>& interesting_for_1 = sequenceFragments_interesting.at(stretchI).at(index_for_1);
			const std::vector<bool>& interesting_for_2 = sequenceFragments_interesting.at(stretchI).at(index_for_2);

			assert(sequence_for_1.length() == L); assert(sequence_for_2.length() == L);
			assert(positions_for_1.size() == L); assert(positions_for_2.size() == L);
			assert(interesting_for_1.size() == L); assert(interesting_for_2.size() == L);

			nucleotideHaplotype_1.replace(boundaries.first, L, sequence_for_1);
			nucleotideHaplotype_2.replace(boundaries.first, L, sequence_for_2);
			for(unsigned int pI = boundaries.first; pI <= boundaries.second; pI++)
			{
				nucleotideHaplotype_1_positions.at(pI) = positions_for_1.at(pI - boundaries.first);
				nucleotideHaplotype_2_positions.at(pI) = positions_for_2.at(pI - boundaries.first);
				nucleotideHaplotype_1_interesting.at(pI) = interesting_for_1.at(pI - boundaries.first);
				nucleotideHaplotype_2_interesting.at(pI) = interesting_for_2.at(pI - boundaries.first);
			}
		}
	};

	populateNucleotideHaplotypes(nucleotideHaplotype_1, nucleotideHaplotype_2, nucleotideHaplotype_1_positions, nucleotideHaplotype_2_positions, nucleotideHaplotype_1_interesting, nucleotideHaplotype_2_interesting, utilizingIndex);
	assert(nucleotideHaplotype_1.find("#") == std::string::npos);
	assert(nucleotideHaplotype_2.find("#") == std::string::npos);

	assert(p_per_epitope_forRet.size() == p_per_epitope_locations_forRet.size());
	for(auto k : p_per_epitope_forRet)
	{
		p_per_epitope_forRet.at(k.first).clear();
		p_per_epitope_locations_forRet.at(k.first).clear();
	}

	bool done = false;
	int consideredHaplotypes = 0;
	while(!done)
	{
		// do something with the haplotype

		populateFragmentStorageFromNucleotideHaplotypePair_stopAware_additive(
			nucleotideHaplotype_1,
			nucleotideHaplotype_1_positions,
			nucleotideHaplotype_1_interesting,
			nucleotideHaplotype_2,
			nucleotideHaplotype_2_positions,
			nucleotideHaplotype_2_interesting,
			1.0/(double)n_haplotype_pairs,
			p_per_epitope_forRet,
			p_per_epitope_locations_forRet
		);
		consideredHaplotypes++;

		// increase by 1
		int currentDigit = utilizingIndex.size() - 1;
		while((currentDigit >= 0) && ((utilizingIndex.at(currentDigit) >= ((int)sequenceFragments.at(currentDigit).size()-1)) || (currentDigit == first_heterozygous_position)))
		{
			currentDigit--;
		}
		if(currentDigit < 0)
		{
			done = true;
			break;
		}
		else
		{
			utilizingIndex.at(currentDigit)++;
			assert(utilizingIndex.at(currentDigit) < (int)sequenceFragments.at(currentDigit).size());
			for(int digitI = currentDigit+1; digitI < (int)utilizingIndex.size(); digitI++)
			{
				utilizingIndex.at(digitI) = 0;
			}
			populateNucleotideHaplotypes(nucleotideHaplotype_1, nucleotideHaplotype_2, nucleotideHaplotype_1_positions, nucleotideHaplotype_2_positions, nucleotideHaplotype_1_interesting, nucleotideHaplotype_2_interesting, utilizingIndex);
		}
	}
	assert(done);

	std::cout << "Expected haplotype pairs: " << n_haplotype_pairs << " // processed: " << consideredHaplotypes << "\n" << std::flush;

	assert((int)n_haplotype_pairs == consideredHaplotypes);
}

void populateFragmentStorageFromNucleotideHaplotypePair_stopAware_additive(const std::string& sequence_1, const std::vector<int>& positions_1, const std::vector<bool>& interesting_1, const std::string& sequence_2, const std::vector<int>& positions_2, const std::vector<bool>& interesting_2, double p, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet)
{

	assert(p >= 0);
	assert(p <= 1);

	assert(sequence_1.length() == positions_1.size());
	assert(sequence_1.length() == interesting_1.size());
	assert(sequence_2.length() == positions_2.size());
	assert(sequence_2.length() == interesting_2.size());

	std::set<int> fragmentSizes;
	assert(p_per_epitope_forRet.size() == p_per_epitope_locations_forRet.size());
	for(auto f : p_per_epitope_forRet)
	{
		assert(f.first > 0);
		fragmentSizes.insert(f.first);
		assert(p_per_epitope_locations_forRet.count(f.first));
	}

	for(int fragmentSize : fragmentSizes)
	{
		fragmentT h1_fragment = AAHaplotypeFromSequence_stopAware(sequence_1, positions_1, interesting_1);
		fragmentT h2_fragment = AAHaplotypeFromSequence_stopAware(sequence_2, positions_2, interesting_2);

		//std::cout << "One pair:\n";
		//printFragment(h1_fragment);
		//printFragment(h2_fragment);
		//std::cout << "--" << "\n" << std::flush;

		assert(std::get<0>(h1_fragment).size() == std::get<1>(h1_fragment).size());
		assert(std::get<0>(h2_fragment).size() == std::get<2>(h2_fragment).size());
		std::vector<fragmentT> fragments_1 = AAHaplotypeIntoFragments(fragmentSize, h1_fragment);
		std::vector<fragmentT> fragments_2 = AAHaplotypeIntoFragments(fragmentSize, h2_fragment);

		std::set<fragmentT> combined_fragments;
		combined_fragments.insert(fragments_1.begin(), fragments_1.end());
		combined_fragments.insert(fragments_2.begin(), fragments_2.end());

		std::map<std::string, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>> peptide_2_position;
		for(const fragmentT& fragment : combined_fragments)
		{
			peptide_2_position[std::get<0>(fragment)].insert(make_pair(std::get<1>(fragment), std::get<2>(fragment)));
		}

		for(auto peptide_and_positions : peptide_2_position)
		{
			const std::string& peptide = peptide_and_positions.first;

			if(peptide == "SSSSSS")
			{
				//std::cout << "p: " << p << "\n" << std::flush;
				//printFragment(fragment);
			}

			if(p_per_epitope_forRet.at(fragmentSize).count(peptide))
			{
				p_per_epitope_forRet.at(fragmentSize).at(peptide).first += p;
				p_per_epitope_forRet.at(fragmentSize).at(peptide).second.insert(peptide_and_positions.second.begin(), peptide_and_positions.second.end());
			}
			else
			{
				p_per_epitope_forRet.at(fragmentSize)[peptide] = make_pair(p, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>({peptide_and_positions.second}));
			}

			for(auto thisLocation : peptide_and_positions.second)
			{
				if(p_per_epitope_locations_forRet.at(fragmentSize)[peptide].count(thisLocation))
				{
					p_per_epitope_locations_forRet.at(fragmentSize)[peptide].at(thisLocation) += p;
				}
				else
				{
					p_per_epitope_locations_forRet.at(fragmentSize)[peptide][thisLocation] = p;
				}
			}
		}
	}
}

std::vector<fragmentT> AAHaplotypeIntoFragments(int k, const fragmentT& haplotype)
{
	assert(std::get<0>(haplotype).size() == std::get<1>(haplotype).size());
	assert(std::get<0>(haplotype).size() == std::get<2>(haplotype).size());
	int n_fragments = std::get<0>(haplotype).length() - k + 1;
	std::vector<fragmentT> forReturn;
	if(n_fragments > 0)
	{
		forReturn.reserve(n_fragments);
		for(int pI = 0; pI < n_fragments; pI++)
		{
			std::string S = std::get<0>(haplotype).substr(pI, k);
			std::vector<std::pair<int, int>> S_p = std::vector<std::pair<int, int>>(std::get<1>(haplotype).begin()+pI, std::get<1>(haplotype).begin()+pI+k);
			std::vector<bool> S_i = std::vector<bool>(std::get<2>(haplotype).begin()+pI, std::get<2>(haplotype).begin()+pI+k);
			assert((int)S.length() == k); assert((int)S_p.size() == k); assert((int)S_i.size() == k);
			forReturn.push_back(make_tuple(S, S_p, S_i));
		}
	}
	return forReturn;
}

fragmentT AAHaplotypeFromSequence_stopAware(const std::string& sequence, const std::vector<int>& positions, const std::vector<bool>& interesting)
{
	assert(sequence.length() == positions.size()); assert(sequence.length() == interesting.size());
	std::vector<fragmentT> forReturn;

	std::string AAs; AAs.reserve(sequence.length()/3);
	std::vector<std::pair<int, int>> AAs_firstLast; AAs_firstLast.reserve(AAs.capacity());
	std::vector<bool> AAs_interesting; AAs_interesting.reserve(AAs.capacity());


	std::string runningCodon;
	std::vector<int> runningPositions;
	std::vector<bool> runningInteresting;

	int currentPos = 0;
	while(currentPos < (int)sequence.length())
	{
		unsigned char thisC = sequence.at(currentPos);
		if((thisC != '-') && (thisC != '_'))
		{
			runningCodon.push_back(thisC);
			runningPositions.push_back(positions.at(currentPos));
			runningInteresting.push_back(interesting.at(currentPos));
		}

		if(runningCodon.length() == 3)
		{
			std::string translation = translateCodon2AA(runningCodon);
			bool interesting = (runningInteresting.at(0) || runningInteresting.at(1) || runningInteresting.at(2));
			int lastPos = -2;
			int minPos = -2;
			int maxPos = -2;
			for(auto p : runningPositions)
			{
				if(p != -1)
				{
					if(minPos == -2)
						minPos = p;

					maxPos = p;

					if(lastPos != -2)
					{
						assert(p > lastPos);
					}
					lastPos = p;
				}
			}
			assert(maxPos >= minPos);
			assert(((minPos == -2) && (maxPos == -2)) || ((minPos >= 0) && (maxPos >= 0)));
			if((minPos == -2) && (maxPos == -2))
			{
				minPos = -1;
				maxPos = -1;
			}

			if(translation == "!")
			{
				break;
			}
			else
			{
				AAs.append(translation);
				AAs_firstLast.push_back(make_pair(minPos, maxPos));
				AAs_interesting.push_back(interesting);
			}

			runningCodon.clear();
			runningPositions.clear();
			runningInteresting.clear();
		}
		currentPos++;
	}

	assert(AAs.size() == AAs_firstLast.size());
	assert(AAs.size() == AAs_interesting.size());
	return make_tuple(AAs, AAs_firstLast, AAs_interesting);
}


void enumeratePeptideHaplotypes_improperFrequencies_oneTranscript(const transcript& transcript, const std::map<std::string, std::string> referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet)
{
	for(auto k_and_stored_fragments : p_per_epitope_forRet)
	{
		k_and_stored_fragments.second.clear();

		int peptideHaplotypeLength = k_and_stored_fragments.first;

		auto addToFragmentsFromOneExtension = [](std::map<fragmentT, double>& fragmentsStore, const std::set<fragmentT>& fragmentsThisExtension, double p) -> void {
			for(const fragmentT& fragmentFromExtension : fragmentsThisExtension)
			{
				if(fragmentsStore.count(fragmentFromExtension) == 0)
				{
					fragmentsStore[fragmentFromExtension] = 0;
				}
				fragmentsStore.at(fragmentFromExtension) += p;
			}
		};

		using runningHaplotypeKey = std::tuple<std::string, std::vector<std::pair<int, int>>, std::vector<bool>, std::string, std::vector<bool>, std::vector<int>>;
		using runningHaplotypePairKey = std::pair<runningHaplotypeKey, runningHaplotypeKey>;

		class runningHaplotype {
		public:
			std::string AApart;
			std::vector<std::pair<int, int>> AApart_firstLast;
			std::vector<bool> AApart_interesting;

			std::string nucleotidePart;
			std::vector<bool> nucleotidePart_interesting;
			std::vector<int> nucleotidePart_refCoordinates;
			int nucleotidePart_nonGap;

			bool have_been_stopped;

			runningHaplotypeKey getKey()
			{
				return make_tuple(
						AApart,
						AApart_firstLast,
						AApart_interesting,
						nucleotidePart,
						nucleotidePart_interesting,
						nucleotidePart_refCoordinates
				);
			}

			runningHaplotype()
			{
				nucleotidePart_nonGap = 0;
				have_been_stopped = false;
			}

			void extendWithNucleotides(const std::string& nucleotides, const std::vector<int>& nucleotides_refCoordinates, const std::vector<bool>& nucleotides_interesting)
			{
				assert(! have_been_stopped);
				assert(nucleotides.size() == nucleotides_refCoordinates.size());
				assert(nucleotides.size() == nucleotides_interesting.size());

				assert(nucleotidePart.size() == nucleotidePart_refCoordinates.size());
				assert(nucleotidePart.size() == nucleotidePart_interesting.size());
				assert((int)countCharacters_noGaps(nucleotidePart) == nucleotidePart_nonGap); // paranoid

				nucleotidePart.insert(nucleotidePart.end(), nucleotides.begin(), nucleotides.end());
				nucleotidePart_refCoordinates.insert(nucleotidePart_refCoordinates.end(), nucleotides_refCoordinates.begin(), nucleotides_refCoordinates.end());
				nucleotidePart_interesting.insert(nucleotidePart_interesting.end(), nucleotides_interesting.begin(), nucleotides_interesting.end());

				bool OK = true;
				int lastRefPos = -2;
				for(int i = 0; i < (int)nucleotidePart_refCoordinates.size(); i++)
				{
					if((lastRefPos == -2) || (nucleotidePart_refCoordinates.at(i) == -1) || (nucleotidePart_refCoordinates.at(i) > lastRefPos))
					{
						if(nucleotidePart_refCoordinates.at(i) != -1)
						{
							lastRefPos = nucleotidePart_refCoordinates.at(i);
						}
					}
					else
					{
						OK = false;
					}
				}

				if(! OK)
				{
					std::cerr << "Problem!\n";
					for(auto p : nucleotidePart_refCoordinates)
					{
						std::cerr << p << " ";
					}
					std::cerr << "\n";
					for(auto p : nucleotides_refCoordinates)
					{
						std::cerr << p << " ";
					}
					std::cerr << "\n" << std::flush;
					throw std::runtime_error("Position problem");
				}

				int added_nonGaps = countCharacters_noGaps(nucleotides);
				nucleotidePart_nonGap += added_nonGaps;

				while(nucleotidePart_nonGap >= 3)
				{
					// find the next codon, its coordinates, and whether it's interesting
					std::string codon; codon.reserve(3);
					int codon_firstPos = -2;
					int codon_lastPos = -2;
					bool codon_interesting = false;
					int consumedIndex_inNucleotidePart = 0;
					while(codon.length() < 3)
					{
						unsigned char charForConsumption = nucleotidePart.at(consumedIndex_inNucleotidePart);
						int refCoordinateForConsumption = nucleotidePart_refCoordinates.at(consumedIndex_inNucleotidePart);

						codon_interesting = (codon_interesting || nucleotidePart_interesting.at(consumedIndex_inNucleotidePart));
						consumedIndex_inNucleotidePart++;
						if((charForConsumption != '-') && (charForConsumption != '_'))
						{
							codon.push_back(charForConsumption);
							if(refCoordinateForConsumption != -1)
							{
								if((codon_firstPos == -2))
									codon_firstPos = refCoordinateForConsumption;
								codon_lastPos = refCoordinateForConsumption;
							}
						}

						assert(consumedIndex_inNucleotidePart <= (int)nucleotidePart.size());
					}

					if(codon_firstPos == -2)
					{
						assert(codon_lastPos == -2);
						codon_firstPos = -1;
						codon_lastPos = -1;
					}
					else
					{
						assert(codon_lastPos >= 0);
						assert(codon_firstPos <= codon_lastPos);
					}

					// translate codon -> AA part
					AApart.append(translateCodon2AA(codon));
					AApart_interesting.push_back(codon_interesting);
					if(!(codon_firstPos <= codon_lastPos))
					{
						std::cerr << "codon_firstPos" << ": " << codon_firstPos << "\n";
						std::cerr << "codon_lastPos" << ": " << codon_lastPos << "\n";
						for(auto nP : nucleotidePart_refCoordinates)
						{
							std::cerr << nP << " ";
						}
						std::cerr << "\n" << std::flush;
					}
					assert(codon_firstPos <= codon_lastPos);
					AApart_firstLast.push_back(make_pair(codon_firstPos, codon_lastPos));

					// remove the consumed components
					bool haveRemainder = (consumedIndex_inNucleotidePart < (int)nucleotidePart.size());
					nucleotidePart = haveRemainder ? nucleotidePart.substr(consumedIndex_inNucleotidePart) : "";
					nucleotidePart_interesting = haveRemainder ? std::vector<bool>(nucleotidePart_interesting.begin()+consumedIndex_inNucleotidePart, nucleotidePart_interesting.end()) : std::vector<bool>();
					nucleotidePart_refCoordinates = haveRemainder ? std::vector<int>(nucleotidePart_refCoordinates.begin()+consumedIndex_inNucleotidePart, nucleotidePart_refCoordinates.end()) : std::vector<int>();
					nucleotidePart_nonGap -= 3;

					assert(nucleotidePart.size() == nucleotidePart_refCoordinates.size());
					assert(nucleotidePart.size() == nucleotidePart_interesting.size());
					assert((int)countCharacters_noGaps(nucleotidePart) == nucleotidePart_nonGap); // paranoid
				}
			}

			void shortenAA(int peptideHaplotypeLength, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& fragmentsStore, unsigned int sampleHaplotypes_size, int openedHaplotypes, bool have_deleted_haplotype)
			{
				size_t position_stop = AApart.find("!");
				if(position_stop != std::string::npos)
				{
					AApart = AApart.substr(0, position_stop);
					AApart_firstLast = std::vector<std::pair<int, int>>(AApart_firstLast.begin(), AApart_firstLast.begin()+position_stop);
					AApart_interesting = std::vector<bool>(AApart_interesting.begin(), AApart_interesting.begin()+position_stop);
					assert(AApart.size() == AApart_firstLast.size()); assert(AApart.size() == AApart_interesting.size());
					assert(AApart.find("!") == std::string::npos);
				}

				while((int)AApart.length() >= peptideHaplotypeLength)
				{
					std::string extract_AA = AApart.substr(0, peptideHaplotypeLength);
					std::vector<std::pair<int, int>> extract_positions = std::vector<std::pair<int, int>>(AApart_firstLast.begin(), AApart_firstLast.begin()+peptideHaplotypeLength);
					std::vector<bool> extract_interesting = std::vector<bool>(AApart_interesting.begin(), AApart_interesting.begin()+peptideHaplotypeLength);
					assert((int)extract_AA.size() == peptideHaplotypeLength); assert((int)extract_positions.size() == peptideHaplotypeLength); assert((int)extract_interesting.size() == peptideHaplotypeLength);

					AApart = AApart.substr(1);
					AApart_firstLast = std::vector<std::pair<int, int>>(AApart_firstLast.begin()+1, AApart_firstLast.end());
					AApart_interesting = std::vector<bool>(AApart_interesting.begin()+1, AApart_interesting.end());
					assert(AApart.size() == AApart_firstLast.size()); assert(AApart.size() == AApart_interesting.size());

					double p;
					if((sampleHaplotypes_size <= 2) && (! have_deleted_haplotype))
					{
						p = 1;
					}
					else
					{
						if(have_deleted_haplotype && (openedHaplotypes <= 2))
						{
							assert(sampleHaplotypes_size <= 2);
							p = 1;
						}
						else
						{
							p = -1;
						}
					}
					assert(extract_AA.length() == peptideHaplotypeLength);
					if(fragmentsStore.at(peptideHaplotypeLength).count(extract_AA))
					{
						if(p > fragmentsStore.at(peptideHaplotypeLength)[extract_AA].first)
						{
							fragmentsStore.at(peptideHaplotypeLength)[extract_AA].first = p;
						}
					}
					else
					{
						fragmentsStore.at(peptideHaplotypeLength)[extract_AA].first = p;
					}
					fragmentsStore.at(peptideHaplotypeLength).at(extract_AA).second.insert(make_pair(extract_positions, extract_interesting));
				}

				if(position_stop != std::string::npos)
				{
					have_been_stopped = true;
				}

			}
		};

		assert(transcript.strand == '+');
		std::string chromosomeID = transcript.chromosomeID;
		assert(referenceGenome_plus.count(chromosomeID));

		int exons = transcript.exons.size();

		// make sure that the exons go from left to right in a non-overlapping fashion
		int lastValidExonI = -1;
		for(int exonI = 0; exonI < exons; exonI++)
		{
			if(transcript.exons.at(exonI).valid)
			{
				assert(transcript.exons.at(exonI).firstPos <= transcript.exons.at(exonI).lastPos);
				if(lastValidExonI != -1)
				{
					assert(transcript.exons.at(lastValidExonI).lastPos < transcript.exons.at(exonI).firstPos);
				}
				lastValidExonI = exonI;
			}
		}

		runningHaplotype referenceHaplotype;
		runningHaplotype firstHaplotype;
		std::vector<runningHaplotype> sampleHaplotypes = {firstHaplotype};
		int openedHaplotypes = 1;
		bool have_deleted_haplotype = false;

		auto doExtension = [&](std::string referenceSequence, std::vector<int> referenceSequence_refCoordinates, std::vector<std::string> sampleAlleles, std::vector<std::vector<int>> sampleAlleles_refCoordinates, std::vector<std::vector<bool>> sampleAlleles_interesting) -> void {
			assert(sampleAlleles.size() == 2);
			assert(sampleAlleles.size() == sampleAlleles_refCoordinates.size());
			assert(sampleAlleles.size() == sampleAlleles_interesting.size());
			assert(sampleAlleles.at(0).length() == referenceSequence.length());
			assert(sampleAlleles.at(1).length() == referenceSequence.length());

			std::vector<bool> referenceSequence_interesting;
			referenceSequence_interesting.resize(referenceSequence.length(), false);


			referenceHaplotype.extendWithNucleotides(referenceSequence, referenceSequence_refCoordinates, referenceSequence_interesting);

			std::vector<std::vector<bool>> sampleAlleles_perCharacter_interesting;
			sampleAlleles_perCharacter_interesting.reserve(sampleAlleles.size());
			for(const std::string& oneSampleAllele : sampleAlleles)
			{
				std::vector<bool> oneSampleAllele_interesting;
				oneSampleAllele_interesting.resize(oneSampleAllele.length(), false);
				sampleAlleles_perCharacter_interesting.push_back(oneSampleAllele_interesting);
			}

			bool heterozygous = (sampleAlleles.at(0) != sampleAlleles.at(1)) || (sampleAlleles_interesting.at(0) != sampleAlleles_interesting.at(1)) || (sampleAlleles_refCoordinates.at(0) != sampleAlleles_refCoordinates.at(1));
			size_t existingSampleHaplotypes_maxI = sampleHaplotypes.size();
			if(! heterozygous)
			{
				assert(sampleAlleles.at(0) == sampleAlleles.at(1));
				for(unsigned int existingHaplotypeI = 0; existingHaplotypeI < existingSampleHaplotypes_maxI; existingHaplotypeI++)
				{
					sampleHaplotypes.at(existingHaplotypeI).extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_perCharacter_interesting.at(0));
				}
			}
			else
			{
				for(unsigned int existingHaplotypeI = 0; existingHaplotypeI < existingSampleHaplotypes_maxI; existingHaplotypeI++)
				{
					runningHaplotype preExtensionHaplotype = sampleHaplotypes.at(existingHaplotypeI);

					sampleHaplotypes.at(existingHaplotypeI).extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_perCharacter_interesting.at(0));

					sampleHaplotypes.push_back(preExtensionHaplotype);
					sampleHaplotypes.back().extendWithNucleotides(sampleAlleles.at(1), sampleAlleles_refCoordinates.at(1), sampleAlleles_perCharacter_interesting.at(1));

					openedHaplotypes *= 2;
				}

			}


			//if(sampleHaplotypes.size() > 1)
			//	std::cout << "sampleHaplotypes.size()" << ": " << sampleHaplotypes.size() << "\n" << std::flush;


			std::map<fragmentT, double> newFragments_thisExtension;
			for(runningHaplotype& h : sampleHaplotypes)
			{
				h.shortenAA(peptideHaplotypeLength, p_per_epitope_forRet, sampleHaplotypes.size(), openedHaplotypes, have_deleted_haplotype);
			}

			std::set<runningHaplotypeKey> haveHaplotypes;
			std::set<unsigned int> delete_indices;
			for(unsigned int runningHaplotypeI = 0; runningHaplotypeI < sampleHaplotypes.size(); runningHaplotypeI++)
			{
				runningHaplotype& h = sampleHaplotypes.at(runningHaplotypeI);
				if(h.have_been_stopped)
				{
					delete_indices.insert(runningHaplotypeI);
					have_deleted_haplotype = true;
				}
				else
				{
					runningHaplotypeKey hKey = h.getKey();
					if(haveHaplotypes.count(hKey))
					{
						delete_indices.insert(runningHaplotypeI);
					}
					else
					{
						haveHaplotypes.insert(hKey);
					}
				}

			}

			if(sampleHaplotypes.size() > 1)
			{
				/*
				std::cout << "sampleHaplotypes.size()" << ": " << sampleHaplotypes.size() << "\n" << std::flush;
				std::cout << "haveHaplotypes.size()" << ": " <<  haveHaplotypes.size() << "\n";
				std::cout << std::flush;
				*/
			}
			unsigned int sampleHaplotypes_size_beforeDeletion = sampleHaplotypes.size();

			if(delete_indices.size())
			{
				//std::cout << "Delete " << delete_indices.size() << "\n" << std::flush;
			}
			if(delete_indices.size() >= 1)
			{
				//td::cout << "Size before: " << sampleHaplotypes.size() << "\n" << std::flush;

				unsigned int first = *(delete_indices.begin());
				unsigned int last = *(delete_indices.rbegin());
				assert(first <= last);
				int lastIndexDeleted = -1;
				for(std::set<unsigned int>::reverse_iterator deleteIndexIt = delete_indices.rbegin(); deleteIndexIt != delete_indices.rend(); deleteIndexIt++)
				{
					unsigned int indexToDelete = *deleteIndexIt;
					//std::cout << "\t\tNow delete " << indexToDelete << ", size " << sampleHaplotypes.size() << "\n" << std::flush;
					assert((lastIndexDeleted == -1) || ((int)indexToDelete < lastIndexDeleted));
					sampleHaplotypes.erase(sampleHaplotypes.begin() + indexToDelete);
					lastIndexDeleted = (int)indexToDelete;
				}

				//std::cout << "Size after: " << sampleHaplotypes.size() << "\n" << std::flush;

			}

			if(!(sampleHaplotypes.size() == (sampleHaplotypes_size_beforeDeletion - delete_indices.size())))
			{
				std::cerr << "sampleHaplotypes.size()" << ": " << sampleHaplotypes.size() << "\n";
				std::cerr << "delete_indices.size()" << ": " << delete_indices.size() << "\n";
				std::cerr << "sampleHaplotypes_size_beforeDeletion" << ": " << sampleHaplotypes_size_beforeDeletion << "\n";
				std::cerr << std::flush;
			}

			assert(sampleHaplotypes.size() == (sampleHaplotypes_size_beforeDeletion - delete_indices.size()));

			if(sampleHaplotypes.size() == 0)
			{
				return;
			}
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
				assert(sampleHaplotypes.size() >= 1);
				if(variants_plus.count(chromosomeID) && variants_plus.at(chromosomeID).count(referencePos))
				{
					std::tuple<std::string, std::vector<int>, std::vector<bool>, std::vector<std::string>, std::vector<std::vector<int>>, std::vector<std::vector<bool>>> reference_and_variantAlleles = get_reference_and_variantAlleles(variants_plus.at(chromosomeID).at(referencePos), referencePos, exon.lastPos);

					doExtension(std::get<0>(reference_and_variantAlleles), std::get<1>(reference_and_variantAlleles), std::get<3>(reference_and_variantAlleles), std::get<4>(reference_and_variantAlleles), std::get<5>(reference_and_variantAlleles));

					int reference_extension_length_noGaps = countCharacters_noGaps(std::get<0>(reference_and_variantAlleles));
					referencePos += (reference_extension_length_noGaps - 1);
				}
				else
				{
					//std::cerr << "\tnV\n" << std::flush;
					std::string referenceCharacter = referenceGenome_plus.at(chromosomeID).substr(referencePos, 1);
					std::vector<std::string> extendWith = {referenceCharacter, referenceCharacter};
					std::vector<int> referenceAllele_coordinates = {(int)referencePos};
					std::vector<std::vector<int>> sampleAlleles_refCoordinates = {{(int)referencePos},{(int)referencePos}};
					std::vector<std::vector<bool>> extendWith_interesting = {std::vector<bool>({false}), std::vector<bool>({false})};
					doExtension(referenceCharacter, referenceAllele_coordinates, extendWith, sampleAlleles_refCoordinates, extendWith_interesting);
				}
			}
		}


	}

}

void enumeratePeptideHaplotypes_improperFrequencies_plus(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool invertPositions, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet)
{
	checkVariantsConsistentWithReferenceGenome(variants_plus, referenceGenome_plus); // paranoid

	std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> oneTranscript_p_per_epitope_forRet;
	std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> oneTranscript_p_per_epitope_locations_forRet;
	for(auto k : p_per_epitope_forRet)
	{
		oneTranscript_p_per_epitope_forRet[k.first].count("");
		oneTranscript_p_per_epitope_locations_forRet[k.first].count("");
	}
	for(unsigned int transcriptI = 0; transcriptI < transcripts_plus.size(); transcriptI++)
	{
		const transcript& transcript = transcripts_plus.at(transcriptI);
		assert(transcript.strand == '+');
		std::string chromosomeID = transcript.chromosomeID;
		if(referenceGenome_plus.count(chromosomeID) == 0)
			continue;

		// the oneTranscript_* maps are clear'ed within enumeratePeptideHaplotypes_oneTranscript
		enumeratePeptideHaplotypes_improperFrequencies_oneTranscript(transcript, referenceGenome_plus, variants_plus, oneTranscript_p_per_epitope_forRet);

		for(auto k : p_per_epitope_forRet)
		{
			for(auto epitope : oneTranscript_p_per_epitope_forRet.at(k.first))
			{
				std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>> positions = epitope.second.second;

				if(invertPositions)
				{
					std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>> inverted_positions;
					for(auto p : positions)
					{
						std::pair<std::vector<std::pair<int, int>>, std::vector<bool>> inverted_p = p;
						for(std::pair<int, int>& pP : inverted_p.first)
						{
							size_t referenceContig_length = referenceGenome_plus.at(chromosomeID).length();
							pP.first = referenceContig_length - pP.first - 1;
							pP.second = referenceContig_length - pP.second - 1;
						}
						inverted_positions.insert(p);
					}
					positions = inverted_positions;
				}

				if(p_per_epitope_forRet.at(k.first).count(epitope.first))
				{
					if(p_per_epitope_forRet.at(k.first).at(epitope.first).first < epitope.second.first)
					{
						p_per_epitope_forRet.at(k.first).at(epitope.first).first = epitope.second.first;
					}
					p_per_epitope_forRet.at(k.first).at(epitope.first).second.insert(positions.begin(), positions.end());
				}
				else
				{
					p_per_epitope_forRet.at(k.first)[epitope.first].first = epitope.second.first;
					p_per_epitope_forRet.at(k.first)[epitope.first].second = positions;
				}
			}
		}
	}
}

void enumeratePeptideHaplotypes_improperFrequencies(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet)
{
	std::map<std::string, std::string> referenceGenome_minus = getMinusStrandReferenceGenome(referenceGenome);
	std::map<std::string, std::map<int, variantFromVCF>> variants_minus = getMinusStrandVariants(variants, referenceGenome_minus);

	checkVariantsConsistentWithReferenceGenome(variants, referenceGenome);
	checkVariantsConsistentWithReferenceGenome(variants_minus, referenceGenome_minus);

	std::vector<transcript> transcripts_plus = getPlusStrandTranscripts(transcripts);
	std::vector<transcript> transcripts_minus = getMinusStrandTranscripts(transcripts, referenceGenome_minus);

	checkTranscriptsTranslate(transcripts_plus, referenceGenome);
	checkTranscriptsTranslate(transcripts_minus, referenceGenome_minus);

	p_per_epitope_forRet.clear();
	for(int k : haplotypeLengths)
	{
		p_per_epitope_forRet[k].count("");
	}

	enumeratePeptideHaplotypes_improperFrequencies_plus(referenceGenome, transcripts_plus, variants, false, p_per_epitope_forRet);
	enumeratePeptideHaplotypes_improperFrequencies_plus(referenceGenome_minus, transcripts_minus, variants_minus, true, p_per_epitope_forRet);
}

/*
void enumeratePeptides(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool isTumour)
{
	int peptideHaplotypeLength = *(peptideLengths.rbegin()) * 2 - 1;
	peptideHaplotypeLength = 4;

	using fragmentT = std::tuple<std::string, std::vector<std::pair<int, int>>, std::vector<bool>>;
	std::map<fragmentT, double> fragments;

	auto printOneFragment = [](const fragmentT& f, double p) -> void {
		std::cout << std::get<0>(f) << "\n";
		for(bool b : std::get<2>(f))
		{
			std::cout << (int)b;
		}
		std::cout << "\n";
		int minPos = -1;
		int maxPos = -1;
		for(auto p : std::get<1>(f))
		{
			assert(p.first <= p.second);
			assert(p.first != -1);
			assert(p.second != -1);

			if(minPos == -1)
			{
				minPos = p.first;
			}
			assert(minPos <= p.first);

			maxPos = p.second;
		}
		std::cout << minPos << " - " << maxPos << "\n\n" << std::flush;
	};

	auto addToTotalFragmentStore = [&fragments](std::map<fragmentT, double>& fragmentsStore_oneExtension) -> void {
		for(auto oneFragment : fragmentsStore_oneExtension)
		{
			// QC
			int lastCodonEnd = -1;
			for(unsigned int pI = 0; pI < std::get<0>(oneFragment.first).size(); pI++)
			{
				assert(
						((std::get<1>(oneFragment.first).at(pI).first == -1) && (std::get<1>(oneFragment.first).at(pI).second == - 1)) ||
						(std::get<1>(oneFragment.first).at(pI).first <= std::get<1>(oneFragment.first).at(pI).second)
				);
				if(std::get<1>(oneFragment.first).at(pI).first != -1)
				{
					if(lastCodonEnd != -1)
					{
						assert(lastCodonEnd < std::get<1>(oneFragment.first).at(pI).first);
					}
					lastCodonEnd = std::get<1>(oneFragment.first).at(pI).second;
				}
			}
			if(! fragments.count(oneFragment.first))
			{
				fragments[oneFragment.first] = 0;
			}
			if(fragments.at(oneFragment.first) < oneFragment.second)
			{
				fragments.at(oneFragment.first) = oneFragment.second;
			}
		}
	};


	auto addToFragmentsFromOneExtension = [](std::map<fragmentT, double>& fragmentsStore, const std::set<fragmentT>& fragmentsThisExtension, double p) -> void {
		for(const fragmentT& fragmentFromExtension : fragmentsThisExtension)
		{
			if(fragmentsStore.count(fragmentFromExtension) == 0)
			{
				fragmentsStore[fragmentFromExtension] = 0;
			}
			fragmentsStore.at(fragmentFromExtension) += p;
		}
	};

	using runningHaplotypeKey = std::tuple<std::string, std::vector<std::pair<int, int>>, std::vector<bool>, std::string, std::vector<bool>, std::vector<int>>;
	using runningHaplotypePairKey = std::pair<runningHaplotypeKey, runningHaplotypeKey>;

	class runningHaplotype {
	public:
		std::string AApart;
		std::vector<std::pair<int, int>> AApart_firstLast;
		std::vector<bool> AApart_interesting;

		std::string nucleotidePart;
		std::vector<bool> nucleotidePart_interesting;
		std::vector<int> nucleotidePart_refCoordinates;
		int nucleotidePart_nonGap;

		runningHaplotypeKey getKey()
		{
			return make_tuple(
					AApart,
					AApart_firstLast,
					AApart_interesting,
					nucleotidePart,
					nucleotidePart_interesting,
					nucleotidePart_refCoordinates
			);
		}

		runningHaplotype()
		{
			nucleotidePart_nonGap = 0;
		}

		void extendWithNucleotides(const std::string& nucleotides, const std::vector<int>& nucleotides_refCoordinates, const std::vector<bool>& nucleotides_interesting)
		{
			assert(nucleotides.size() == nucleotides_refCoordinates.size());
			assert(nucleotides.size() == nucleotides_interesting.size());

			assert(nucleotidePart.size() == nucleotidePart_refCoordinates.size());
			assert(nucleotidePart.size() == nucleotidePart_interesting.size());
			assert((int)countCharacters_noGaps(nucleotidePart) == nucleotidePart_nonGap); // paranoid

			nucleotidePart.insert(nucleotidePart.end(), nucleotides.begin(), nucleotides.end());
			nucleotidePart_refCoordinates.insert(nucleotidePart_refCoordinates.end(), nucleotides_refCoordinates.begin(), nucleotides_refCoordinates.end());
			nucleotidePart_interesting.insert(nucleotidePart_interesting.end(), nucleotides_interesting.begin(), nucleotides_interesting.end());

			bool OK = true;
			int lastRefPos = -2;
			for(int i = 0; i < (int)nucleotidePart_refCoordinates.size(); i++)
			{
				if((lastRefPos == -2) || (nucleotidePart_refCoordinates.at(i) == -1) || (nucleotidePart_refCoordinates.at(i) > lastRefPos))
				{
					if(nucleotidePart_refCoordinates.at(i) != -1)
					{
						lastRefPos = nucleotidePart_refCoordinates.at(i);
					}
				}
				else
				{
					OK = false;
				}
			}

			if(! OK)
			{
				std::cerr << "Problem!\n";
				for(auto p : nucleotidePart_refCoordinates)
				{
					std::cerr << p << " ";
				}
				std::cerr << "\n";
				for(auto p : nucleotides_refCoordinates)
				{
					std::cerr << p << " ";
				}
				std::cerr << "\n" << std::flush;
				throw std::runtime_error("Position problem");
			}

			int added_nonGaps = countCharacters_noGaps(nucleotides);
			nucleotidePart_nonGap += added_nonGaps;

			while(nucleotidePart_nonGap >= 3)
			{
				// find the next codon, its coordinates, and whether it's interesting
				std::string codon; codon.reserve(3);
				int codon_firstPos = -2;
				int codon_lastPos = -2;
				bool codon_interesting = false;
				int consumedIndex_inNucleotidePart = 0;
				while(codon.length() < 3)
				{
					unsigned char charForConsumption = nucleotidePart.at(consumedIndex_inNucleotidePart);
					int refCoordinateForConsumption = nucleotidePart_refCoordinates.at(consumedIndex_inNucleotidePart);

					codon_interesting = (codon_interesting || nucleotidePart_interesting.at(consumedIndex_inNucleotidePart));
					consumedIndex_inNucleotidePart++;
					if((charForConsumption != '-') && (charForConsumption != '_'))
					{
						codon.push_back(charForConsumption);
						if(refCoordinateForConsumption != -1)
						{
							if((codon_firstPos == -2))
								codon_firstPos = refCoordinateForConsumption;
							codon_lastPos = refCoordinateForConsumption;
						}
					}

					assert(consumedIndex_inNucleotidePart <= (int)nucleotidePart.size());
				}

				if(codon_firstPos == -2)
				{
					assert(codon_lastPos == -2);
					codon_firstPos = -1;
					codon_lastPos = -1;
				}
				else
				{
					assert(codon_lastPos >= 0);
				}

				// translate codon -> AA part
				AApart.append(translateCodon2AA(codon));
				AApart_interesting.push_back(codon_interesting);
				if(!(codon_firstPos <= codon_lastPos))
				{
					std::cerr << "codon_firstPos" << ": " << codon_firstPos << "\n";
					std::cerr << "codon_lastPos" << ": " << codon_lastPos << "\n";
					for(auto nP : nucleotidePart_refCoordinates)
					{
						std::cerr << nP << " ";
					}
					std::cerr << "\n" << std::flush;
				}
				assert(codon_firstPos <= codon_lastPos);
				AApart_firstLast.push_back(make_pair(codon_firstPos, codon_lastPos));

				// remove the consumed components
				bool haveRemainder = (consumedIndex_inNucleotidePart < (int)nucleotidePart.size());
				nucleotidePart = haveRemainder ? nucleotidePart.substr(consumedIndex_inNucleotidePart) : "";
				nucleotidePart_interesting = haveRemainder ? std::vector<bool>(nucleotidePart_interesting.begin()+consumedIndex_inNucleotidePart, nucleotidePart_interesting.end()) : std::vector<bool>();
				nucleotidePart_refCoordinates = haveRemainder ? std::vector<int>(nucleotidePart_refCoordinates.begin()+consumedIndex_inNucleotidePart, nucleotidePart_refCoordinates.end()) : std::vector<int>();
				nucleotidePart_nonGap -= 3;

				assert(nucleotidePart.size() == nucleotidePart_refCoordinates.size());
				assert(nucleotidePart.size() == nucleotidePart_interesting.size());
				assert((int)countCharacters_noGaps(nucleotidePart) == nucleotidePart_nonGap); // paranoid
			}
		}

		void shortenAA(int peptideHaplotypeLength, std::set<fragmentT>& fragmentsStore)
		{
			while((int)AApart.length() >= peptideHaplotypeLength)
			{
				std::string extract_AA = AApart.substr(0, peptideHaplotypeLength);
				std::vector<std::pair<int, int>> extract_positions = std::vector<std::pair<int, int>>(AApart_firstLast.begin(), AApart_firstLast.begin()+peptideHaplotypeLength);
				std::vector<bool> extract_interesting = std::vector<bool>(AApart_interesting.begin(), AApart_interesting.begin()+peptideHaplotypeLength);
				assert((int)extract_AA.size() == peptideHaplotypeLength); assert((int)extract_positions.size() == peptideHaplotypeLength); assert((int)extract_interesting.size() == peptideHaplotypeLength);

				AApart = AApart.substr(1);
				AApart_firstLast = std::vector<std::pair<int, int>>(AApart_firstLast.begin()+1, AApart_firstLast.end());
				AApart_interesting = std::vector<bool>(AApart_interesting.begin()+1, AApart_interesting.end());
				assert(AApart.size() == AApart_firstLast.size()); assert(AApart.size() == AApart_interesting.size());

				fragmentT key = make_tuple(extract_AA, extract_positions, extract_interesting);
				fragmentsStore.insert(key);
			}
		}

		void storeRemainingAA(int peptideHaplotypeLength, std::set<fragmentT>& fragmentsStore)
		{
			assert((int)AApart.length() < peptideHaplotypeLength);
			assert(nucleotidePart.length() < 3);
			fragmentT key = make_tuple(AApart, AApart_firstLast, AApart_interesting);
			fragmentsStore.insert(key);
		}
	};

	class runningPairOfHaplotypes {
	public:
		double p;
		runningHaplotype h1;
		runningHaplotype h2;
		runningPairOfHaplotypes()
		{
			p = 1;
		}
		void extendWithNucleotides(const std::string& h1_nucleotides, const std::vector<int>& h1_nucleotides_refCoordinates, const std::vector<bool>& h1_nucleotides_interesting, const std::string& h2_nucleotides, const std::vector<int>& h2_nucleotides_refCoordinates, const std::vector<bool>& h2_nucleotides_interesting)
		{
			h1.extendWithNucleotides(h1_nucleotides, h1_nucleotides_refCoordinates, h1_nucleotides_interesting);
			h2.extendWithNucleotides(h2_nucleotides, h2_nucleotides_refCoordinates, h2_nucleotides_interesting);
		}

		runningHaplotypePairKey getKey()
		{
			return make_pair(h1.getKey(), h2.getKey());
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

		// make sure that the exons go from left to right in a non-overlapping fashion
		int lastValidExonI = -1;
		for(int exonI = 0; exonI < exons; exonI++)
		{
			if(transcript.exons.at(exonI).valid)
			{
				assert(transcript.exons.at(exonI).firstPos <= transcript.exons.at(exonI).lastPos);
				if(lastValidExonI != -1)
				{
					assert(transcript.exons.at(lastValidExonI).lastPos < transcript.exons.at(exonI).firstPos);
				}
				lastValidExonI = exonI;
			}
		}

		runningHaplotype referenceHaplotype;
		runningPairOfHaplotypes firstSampleHaplotype;
		std::vector<runningPairOfHaplotypes> sampleHaplotypePairs = {firstSampleHaplotype};

		auto doExtension = [&](std::string referenceSequence, std::vector<int> referenceSequence_refCoordinates, std::vector<std::string> sampleAlleles, std::vector<std::vector<int>> sampleAlleles_refCoordinates, std::vector<bool> sampleAlleles_interesting) -> void {
			assert(sampleAlleles.size() == 2);
			assert(sampleAlleles.size() == sampleAlleles_refCoordinates.size());
			assert(sampleAlleles.size() == sampleAlleles_interesting.size());
			assert(sampleAlleles.at(0).length() == referenceSequence.length());
			assert(sampleAlleles.at(1).length() == referenceSequence.length());

			std::vector<bool> referenceSequence_interesting;
			referenceSequence_interesting.resize(referenceSequence.length(), false);


			referenceHaplotype.extendWithNucleotides(referenceSequence, referenceSequence_refCoordinates, referenceSequence_interesting);

			std::vector<std::vector<bool>> sampleAlleles_perCharacter_interesting;
			sampleAlleles_perCharacter_interesting.reserve(sampleAlleles.size());
			for(const std::string& oneSampleAllele : sampleAlleles)
			{
				std::vector<bool> oneSampleAllele_interesting;
				oneSampleAllele_interesting.resize(oneSampleAllele.length(), false);
				sampleAlleles_perCharacter_interesting.push_back(oneSampleAllele_interesting);
			}


			bool nonIdenticalRefCoordinates = (sampleAlleles_refCoordinates.at(0).size() != sampleAlleles_refCoordinates.at(1).size());
			if(!nonIdenticalRefCoordinates)
			{
				for(unsigned int i = 0; i < sampleAlleles_refCoordinates.at(0).size(); i++)
				{
					if(sampleAlleles_refCoordinates.at(0).at(i) != sampleAlleles_refCoordinates.at(1).at(i))
					{
						nonIdenticalRefCoordinates = true;
						break;
					}
				}
			}
			bool heterozygous = (sampleAlleles.at(0) != sampleAlleles.at(1)) || (sampleAlleles_interesting.at(0) != sampleAlleles_interesting.at(1)) || (nonIdenticalRefCoordinates);
			size_t existingSampleHaplotypes_maxI = sampleHaplotypePairs.size();
			if(! heterozygous)
			{
				assert(sampleAlleles.at(0) == sampleAlleles.at(1));
				for(unsigned int existingHaplotypePairI = 0; existingHaplotypePairI < existingSampleHaplotypes_maxI; existingHaplotypePairI++)
				{
					sampleHaplotypePairs.at(existingHaplotypePairI).h1.extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_perCharacter_interesting.at(0));
					sampleHaplotypePairs.at(existingHaplotypePairI).h2.extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_perCharacter_interesting.at(0));
				}
			}
			else
			{
				assert(2 == 4);
				if(existingSampleHaplotypes_maxI == 1)
				{
					sampleHaplotypePairs.at(0).h1.extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_perCharacter_interesting.at(0));
					sampleHaplotypePairs.at(0).h2.extendWithNucleotides(sampleAlleles.at(1), sampleAlleles_refCoordinates.at(1), sampleAlleles_perCharacter_interesting.at(1));
				}
				else
				{
					for(unsigned int existingHaplotypePairI = 0; existingHaplotypePairI < existingSampleHaplotypes_maxI; existingHaplotypePairI++)
					{
						runningPairOfHaplotypes preExtensionPair = sampleHaplotypePairs.at(existingHaplotypePairI);

						sampleHaplotypePairs.at(existingHaplotypePairI).h1.extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_perCharacter_interesting.at(0));
						sampleHaplotypePairs.at(existingHaplotypePairI).h2.extendWithNucleotides(sampleAlleles.at(1), sampleAlleles_refCoordinates.at(1), sampleAlleles_perCharacter_interesting.at(1));
						sampleHaplotypePairs.at(existingHaplotypePairI).p *= 0.5;

						for(unsigned int sampleAlleleI = 1; sampleAlleleI < sampleAlleles.size(); sampleAlleleI++)
						{
							sampleHaplotypePairs.push_back(preExtensionPair);

							sampleHaplotypePairs.back().h1.extendWithNucleotides(sampleAlleles.at(1), sampleAlleles_refCoordinates.at(1), sampleAlleles_perCharacter_interesting.at(1));
							sampleHaplotypePairs.back().h2.extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_perCharacter_interesting.at(0));
							sampleHaplotypePairs.back().p *= 0.5;
						}
					}
				}
			}

			if(sampleHaplotypePairs.size() > 1)
				std::cout << "sampleHaplotypes.size()" << ": " << sampleHaplotypePairs.size() << "\n" << std::flush;

			double p_sum = 0;
			for(auto h : sampleHaplotypePairs)
			{
				p_sum += h.p;
			}
			assert(abs(p_sum - 1) <= 1e-4);

			std::map<fragmentT, double> newFragments_thisExtension;
			for(runningPairOfHaplotypes& hP : sampleHaplotypePairs)
			{
				std::set<fragmentT> thisPair_newFragments;
				hP.h1.shortenAA(peptideHaplotypeLength, thisPair_newFragments);
				hP.h2.shortenAA(peptideHaplotypeLength, thisPair_newFragments);
				addToFragmentsFromOneExtension(newFragments_thisExtension, thisPair_newFragments, hP.p);
			}
			addToTotalFragmentStore(newFragments_thisExtension);

			std::map<runningHaplotypePairKey, unsigned int> havePair;
			std::set<unsigned int> delete_indices;
			for(unsigned int runningHaplotypePairI = 0; runningHaplotypePairI < sampleHaplotypePairs.size(); runningHaplotypePairI++)
			{
				runningPairOfHaplotypes& hP = sampleHaplotypePairs.at(runningHaplotypePairI);
				runningHaplotypePairKey pairKey = hP.getKey();
				if(havePair.count(pairKey))
				{
					sampleHaplotypePairs.at(havePair.at(pairKey)).p += hP.p;
					delete_indices.insert(runningHaplotypePairI);
				}
				else
				{
					havePair[pairKey] = runningHaplotypePairI;
				}
			}

			if(sampleHaplotypePairs.size() > 1)
			{
				std::cout << "sampleHaplotypePairs.size()" << ": " <<  sampleHaplotypePairs.size() << "\n";
				std::cout << "havePair.size()" << ": " <<  havePair.size() << "\n";
				std::cout << std::flush;
			}
			unsigned int sampleHaplotypePairs_size_beforeDeletion = sampleHaplotypePairs.size();

			if(delete_indices.size())
			{
				assert(2 == 4);
				std::cout << "Delete " << delete_indices.size() << "\n" << std::flush;
			}
			if(delete_indices.size() >= 2)
			{
				unsigned int first = *(delete_indices.begin());
				unsigned int last = *(delete_indices.rbegin());
				assert(first < last);
				int lastIndexDeleted = -1;
				for(std::set<unsigned int>::reverse_iterator deleteIndexIt = delete_indices.rbegin(); deleteIndexIt != delete_indices.rend(); deleteIndexIt++)
				{
					unsigned int indexToDelete = *deleteIndexIt;
					assert((lastIndexDeleted == -1) || ((int)indexToDelete < lastIndexDeleted));
					sampleHaplotypePairs.erase(sampleHaplotypePairs.begin() + indexToDelete);
					lastIndexDeleted = (int)indexToDelete;
				}
			}
			assert(sampleHaplotypePairs.size() == (sampleHaplotypePairs_size_beforeDeletion - delete_indices.size()));

			p_sum = 0;
			for(auto h : sampleHaplotypePairs)
			{
				p_sum += h.p;
			}
			assert(abs(p_sum - 1) <= 1e-4);
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
				//assert(2 == 10);
				//std::cerr << "exon " << exonI << " from " << exon.firstPos << " to " << exon.lastPos << " pos: " << referencePos << "\n" << std::flush;
				if(variants_plus.count(chromosomeID) && variants_plus.at(chromosomeID).count(referencePos))
				{
					//std::cerr << "\tv\n" << std::flush;
					std::tuple<std::string, std::vector<int>, std::vector<std::string>, std::vector<std::vector<int>>> reference_and_variantAlleles = get_reference_and_variantAlleles(variants_plus.at(chromosomeID).at(referencePos), referencePos, exon.lastPos);

					std::vector<bool> extendWith_interesting;
					if(isTumour)
					{
						extendWith_interesting.reserve(std::get<2>(reference_and_variantAlleles).size());
						for(const std::string& sampleAllele : std::get<2>(reference_and_variantAlleles))
						{
							extendWith_interesting.push_back(!(sampleAllele == std::get<0>(reference_and_variantAlleles)));
						}
					}
					else
					{
						extendWith_interesting.resize(std::get<2>(reference_and_variantAlleles).size(), false);
					}

					//std::cout << "V: " << std::get<2>(reference_and_variantAlleles).size() << "\n" << std::flush;

					std::vector<int> referenceAllele_coordinates;
					std::vector<std::vector<int>> sampleAlleles_refCoordinates;
					//std::cerr << "\t" << std::get<0>(reference_and_variantAlleles) << "\n" << std::flush;
					/*
					for(unsigned int vI = 0; vI < std::get<2>(reference_and_variantAlleles).size(); vI++)
					{
						std::cerr << "\t" << std::get<2>(reference_and_variantAlleles).at(vI) << "\n" << std::flush;
						std::cerr << "\t";
						for(auto pI : std::get<3>(reference_and_variantAlleles).at(vI))
						{
							std::cerr << pI << " ";
						}
						std::cerr << "\n" << std::flush;
					}

					doExtension(std::get<0>(reference_and_variantAlleles), std::get<1>(reference_and_variantAlleles), std::get<2>(reference_and_variantAlleles), std::get<3>(reference_and_variantAlleles), extendWith_interesting);

					int reference_extension_length_noGaps = countCharacters_noGaps(std::get<0>(reference_and_variantAlleles));

					referencePos += (reference_extension_length_noGaps - 1);
				}
				else
				{
					//std::cerr << "\tnV\n" << std::flush;
					std::string referenceCharacter = referenceGenome_plus.at(chromosomeID).substr(referencePos, 1);
					std::vector<std::string> extendWith = {referenceCharacter, referenceCharacter};
					std::vector<int> referenceAllele_coordinates = {(int)referencePos};
					std::vector<std::vector<int>> sampleAlleles_refCoordinates = {{(int)referencePos},{(int)referencePos}};
					std::vector<bool> extendWith_interesting = {false, false};
					doExtension(referenceCharacter, referenceAllele_coordinates, extendWith, sampleAlleles_refCoordinates, extendWith_interesting);
				}
			}
		}


		std::map<fragmentT, double> newFragments_last;
		for(runningPairOfHaplotypes& hP : sampleHaplotypePairs)
		{
			std::set<fragmentT> thisPair_newFragments;
			hP.h1.storeRemainingAA(peptideHaplotypeLength, thisPair_newFragments);
			hP.h2.storeRemainingAA(peptideHaplotypeLength, thisPair_newFragments);
			addToFragmentsFromOneExtension(newFragments_last, thisPair_newFragments, hP.p);
		}
		addToTotalFragmentStore(newFragments_last);

	}

	std::cout << "Have " << fragments.size() << " fragments.\n" << std::flush;

	/*
	unsigned int printFragments = 10;
	unsigned int printedFragments = 0;
	for(auto f : fragments)
	{
		std::cout << "." << std::flush;
		printOneFragment(f.first, f.second);
		printedFragments++;
		std::cout << printedFragments << " " << printFragments << " " << (printedFragments == printFragments) << "\n" << std::flush;
		if(printedFragments == printFragments)
			break;
	}



	assert(1 == 11);
}

*/

std::tuple<std::string, std::vector<int>, std::vector<bool>, std::vector<std::string>, std::vector<std::vector<int>>, std::vector<std::vector<bool>>> get_reference_and_variantAlleles(const variantFromVCF& v, unsigned int startReferencePos, unsigned int lastReferencePos)
{
	assert(v.position <= lastReferencePos);

	// check that we have variant alleles and that they are "aligned"
	assert(v.sampleAlleles.size() == 2);
	assert(v.sampleAlleles.size() == v.sampleAlleles_interesting.size());
	for(const std::string& variantAllele : v.sampleAlleles)
	{
		assert(variantAllele.length() == v.referenceString.length());
	}

	auto getReferencePositionVector = [](const std::string& S, unsigned int firstPos) -> std::vector<int> {
		std::vector<int> forReturn;
		forReturn.resize(S.length(), -1);
		unsigned int runningRefPosition = firstPos;
		for(unsigned int i = 0; i < S.length(); i++)
		{
			unsigned char thisC = S.at(i);
			if((thisC != '-') && (thisC != '_'))
			{
				forReturn.at(i) = runningRefPosition;
				runningRefPosition++;
			}
		}

		int lastRefPos = -2;
		for(int i = 0; i < (int)forReturn.size(); i++)
		{
			assert((lastRefPos == -2) || (forReturn.at(i) == -1) || (forReturn.at(i) > lastRefPos));
			if(forReturn.at(i) != -1)
			{
				lastRefPos = forReturn.at(i);
			}
			// assert(forReturn.at(i-1) <= forReturn.at(i));
		}
		return forReturn;
	};

	auto getInterestingVector = [](int length, bool interesting) -> std::vector<bool> {
		std::vector<bool> forReturn;
		forReturn.resize(length, interesting);
		return forReturn;
	};


	assert(v.position == startReferencePos);
	// if everything within last-position constraints, do nothing ...
	unsigned int v_referenceString_length_noGaps = countCharacters_noGaps(v.referenceString);
	if((v.position + v_referenceString_length_noGaps - 1) <= lastReferencePos)
	{
		std::vector<std::vector<int>> referencePositions_per_sampleAllele;
		referencePositions_per_sampleAllele.reserve(v.sampleAlleles.size());
		std::vector<int> referencePositionVector = getReferencePositionVector(v.referenceString, startReferencePos);
		for(const std::string& sampleAllele : v.sampleAlleles)
		{
			referencePositions_per_sampleAllele.push_back(referencePositionVector);
		}

		std::vector<std::vector<bool>> forReturn_sampleAlleles_interesting({getInterestingVector(v.referenceString.length(), v.sampleAlleles_interesting.at(0)), getInterestingVector(v.referenceString.length(), v.sampleAlleles_interesting.at(1))});
		return make_tuple(v.referenceString, referencePositionVector, getInterestingVector(v.referenceString.length(), false), v.sampleAlleles, referencePositions_per_sampleAllele, forReturn_sampleAlleles_interesting);
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

		std::vector<std::vector<int>> referencePositions_per_returnAllele;
		referencePositions_per_returnAllele.reserve(v.sampleAlleles.size());

		std::vector<int> referencePositionVector = getReferencePositionVector(forReturn_ref, startReferencePos);

		for(const std::string& sampleAllele : forReturn_alleles)
		{
			referencePositions_per_returnAllele.push_back(referencePositionVector);
		}

		std::vector<std::vector<bool>> forReturn_sampleAlleles_interesting({getInterestingVector(forReturn_ref.length(), v.sampleAlleles_interesting.at(0)), getInterestingVector(forReturn_ref.length(), v.sampleAlleles_interesting.at(1))});
		return make_tuple(forReturn_ref, referencePositionVector, getInterestingVector(forReturn_ref.length(), false), forReturn_alleles, referencePositions_per_returnAllele, forReturn_sampleAlleles_interesting);
	}
}

void printFragment(const fragmentT& f)
{
	std::cout << "Fragment print:\n";
	std::cout << "\t" << std::get<0>(f) << "\n";
	const std::vector<std::pair<int, int>>& positions =  std::get<1>(f);
	const std::vector<bool>& interesting =  std::get<2>(f);
	std::cout << "\t";
	for(auto p : positions)
	{
		std::cout << p.first << "-" << p.second << " ";
	}
	std::cout << "\n";
	std::cout << "\t";
	for(auto i : interesting)
	{
		std::cout << (int)i;
	}
	std::cout << "\n" << std::flush;
}



