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
//void enumeratePeptides(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool isTumour);
std::tuple<std::string, std::vector<int>, std::vector<std::string>, std::vector<std::vector<int>>> get_reference_and_variantAlleles(const variantFromVCF& v, unsigned int startReferencePos, unsigned int lastReferencePos);

void enumeratePeptideHaplotypes_oneTranscript(const transcript& transcript, const std::map<std::string, std::string> referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool isTumour);
void enumeratePeptideHaplotypes(int haplotypeLength, const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool isTumour);

std::set<int> peptideLengths = {8, 9, 10, 11, 15, 16};
using fragmentT = std::tuple<std::string, std::vector<std::pair<int, int>>, std::vector<bool>>;

std::vector<fragmentT> AAHaplotypeIntoFragments(int k, const fragmentT& haplotype);
fragmentT AAHaplotypeFromSequence_stopAware(const std::string& sequence, const std::vector<int>& positions, const std::vector<bool>& interesting);

int main(int argc, char *argv[]) {

	std::vector<std::string> ARG (argv + 1, argv + argc + !argc);
	std::map<std::string, std::string> arguments;

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


	enumeratePeptideHaplotypes(10, referenceGenome, transcripts_plus, variants, true);

	assert(1 == 2);

	assert("minus-strand transcripts!" == "");
	assert("add stop codons!" == "");

	return 0;
}

void enumeratePeptideHaplotypes(int haplotypeLength, const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool isTumour)
{
	for(unsigned int transcriptI = 0; transcriptI < transcripts_plus.size(); transcriptI++)
	{
		const transcript& transcript = transcripts_plus.at(transcriptI);
		assert(transcript.strand == '+');
		std::string chromosomeID = transcript.chromosomeID;
		if(referenceGenome_plus.count(chromosomeID) == 0)
			continue;

		enumeratePeptideHaplotypes_oneTranscript(transcript, referenceGenome_plus, variants_plus, isTumour);
	}

}

void enumeratePeptideHaplotypes_oneTranscript(const transcript& transcript, const std::map<std::string, std::string> referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool isTumour, std::map<int, std::map<fragmentT, double>>& haplotypeStore)
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
				std::tuple<std::string, std::vector<int>, std::vector<std::string>, std::vector<std::vector<int>>> reference_and_variantAlleles = get_reference_and_variantAlleles(variants_plus.at(chromosomeID).at(referencePos), referencePos, exon.lastPos);

				assert(std::get<2>(reference_and_variantAlleles).size() == 2);
				bool homozygous = (std::get<2>(reference_and_variantAlleles).at(0) == std::get<2>(reference_and_variantAlleles).at(1)) && (std::get<3>(reference_and_variantAlleles).at(0) == std::get<3>(reference_and_variantAlleles).at(1));

				// this is not really very interesting
				std::vector<std::vector<bool>> extendWith_interesting;
				for(const std::string& sampleAllele : std::get<2>(reference_and_variantAlleles))
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

				if(homozygous)
				{
					addHomozygousStretch(std::get<2>(reference_and_variantAlleles).at(0),  std::get<3>(reference_and_variantAlleles).at(0), extendWith_interesting.at(0));
				}
				else
				{
					addHeterozygousStretch(std::get<2>(reference_and_variantAlleles), std::get<3>(reference_and_variantAlleles), extendWith_interesting);
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


	bool done = false;
	int consideredHaplotypes = 0;
	while(!done)
	{
		// do something with the haplotype
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

	assert(n_haplotype_pairs == consideredHaplotypes);
}

void populateFragmentStorageFromNucleotideHaplotypePair_stopAware_additive(const std::string& sequence_1, const std::vector<int>& positions_1, const std::vector<bool>& interesting_1, const std::string& sequence_2, const std::vector<int>& positions_2, const std::vector<bool>& interesting_2, double p, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& fragmentStore)
{
	std::set<int> fragmentSizes;
	for(auto f : fragmentStore)
	{
		assert(f.first > 0);
		fragmentSizes.insert(f.first);
	}

	for(int fragmentSize : fragmentSizes)
	{
		fragmentT h1_fragment = make_tuple(sequence_1, positions_1, interesting_1);
		fragmentT h2_fragment = make_tuple(sequence_2, positions_2, interesting_2);
		std::vector<fragmentT> fragments_1 = AAHaplotypeIntoFragments(fragmentSize, h1_fragment);
		std::vector<fragmentT> fragments_2 = AAHaplotypeIntoFragments(fragmentSize, h2_fragment);

		std::set<fragmentT> combined_fragments;
		combined_fragments.insert(fragments_1.begin(), fragments_1.end());
		combined_fragments.insert(fragments_2.begin(), fragments_2.end());

		for(const fragmentT& fragment : combined_fragments)
		{
			if(fragmentStore.at(fragmentSize).count(std::get<0>(fragment)))
			{
				fragmentStore.at(fragmentSize).at(std::get<0>(fragment)).first += p;
				fragmentStore.at(fragmentSize).at(std::get<0>(fragment)).second.insert(make_pair(std::get<1>(fragment), std::get<2>(fragment)));
			}
			else
			{
				fragmentStore.at(fragmentSize)[std::get<0>(fragment)] = make_pair(p, {make_pair(std::get<1>(fragment), std::get<2>(fragment))});
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
	forReturn.reserve(n_fragments);
	for(unsigned int pI = 0; pI < n_fragments; pI++)
	{
		std::string S = std::get<0>(haplotype).substr(pI, k);
		std::vector<std::pair<int, int>> S_p = std::vector<std::pair<int, int>>(std::get<1>(haplotype).begin()+pI, std::get<1>(haplotype).begin()+pI+k);
		std::vector<bool> S_i = std::vector<bool>(std::get<2>(haplotype).begin()+pI, std::get<2>(haplotype).begin()+pI+k);
		assert(S.length() == k); assert(S_p.size() == k); assert(S_i.size() == k);
		forReturn.push_back(make_tuple(S, S_p, S_i));
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
	while(currentPos < sequence.length())
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
			bool interesting = runningInteresting.at(0) || runningPositions.at(1) || runningPositions.at(2);
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
		}
		currentPos++;
	}

	return make_tuple(AAs, AAs_firstLast, AAs_interesting);
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

std::tuple<std::string, std::vector<int>, std::vector<std::string>, std::vector<std::vector<int>>> get_reference_and_variantAlleles(const variantFromVCF& v, unsigned int startReferencePos, unsigned int lastReferencePos)
{
	assert(v.position <= lastReferencePos);

	// check that we have variant alleles and that they are "aligned"
	assert(v.sampleAlleles.size());
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
		return make_tuple(v.referenceString, referencePositionVector, v.sampleAlleles, referencePositions_per_sampleAllele);
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

		return make_tuple(forReturn_ref, referencePositionVector, forReturn_alleles, referencePositions_per_returnAllele);
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



