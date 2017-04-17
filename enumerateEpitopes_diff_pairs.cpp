/*
 * enumerateEpitopesdiffproper.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: diltheyat
 */

#include "enumerateEpitopes_diff_pairs.h"

/*
 * enumerateEpitopeshaplotypePairs.cpp
 *
 *  Created on: Apr 11, 2017
 *      Author: diltheyat
 */

#include <iostream>
#include <assert.h>
#include <cmath>
#include <omp.h>

#include "Util.h"
#include "enumerateEpitopes_haplotypePairs.h"

void forBaseline_populateFragmentStorageFromNucleotideHaplotypePair_stopAware_additive(int k, const std::string& sequence_1, const std::vector<int>& positions_1, const std::vector<bool>& interesting_1, const std::string& sequence_2, const std::vector<int>& positions_2, const std::vector<bool>& interesting_2, double p, std::map<std::string, double>& epitopes_and_p)
{
	assert(p >= 0);
	assert(p <= 1);

	assert(sequence_1.length() == positions_1.size());
	assert(sequence_1.length() == interesting_1.size());
	assert(sequence_2.length() == positions_2.size());
	assert(sequence_2.length() == interesting_2.size());

	fragmentT h1_fragment = AAHaplotypeFromSequence_stopAware(sequence_1, positions_1, interesting_1);
	fragmentT h2_fragment = AAHaplotypeFromSequence_stopAware(sequence_2, positions_2, interesting_2);

	/*
	std::cout << "One pair:\n";
	std::cout << sequence_1 << "\n";
	printFragment(h1_fragment);
	std::cout << sequence_2 << "\n";
	printFragment(h2_fragment);
	std::cout << "--" << "\n" << std::flush;
	*/
	assert(std::get<0>(h1_fragment).size() == std::get<1>(h1_fragment).size());
	assert(std::get<0>(h2_fragment).size() == std::get<2>(h2_fragment).size());
	std::vector<fragmentT> fragments_1 = AAHaplotypeIntoFragments(k, h1_fragment);
	std::vector<fragmentT> fragments_2 = AAHaplotypeIntoFragments(k, h2_fragment);

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

		if(epitopes_and_p.count(peptide))
		{
			epitopes_and_p.at(peptide) += p;
		}
		else
		{
			epitopes_and_p[peptide] = p;
		}
	}
}

std::map<std::string, double> enumeratePeptideHaplotypes_baseLine_oneTranscript(int k, const transcript& transcript, const std::map<std::string, std::string>& referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus)
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
		stretchBoundaries.push_back(std::make_pair(lastStretchStop+1,thisStretchStop));
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

	std::cout << n_haplotype_pairs << "  " << totalLength << "       " << "\n" << std::flush;
	
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

	bool done = false;
	int consideredHaplotypes = 0;
	std::map<std::string, double> forReturn_thisTranscript;
	while(!done)
	{
		// do something with the haplotype

		forBaseline_populateFragmentStorageFromNucleotideHaplotypePair_stopAware_additive(
			k,
			nucleotideHaplotype_1,
			nucleotideHaplotype_1_positions,
			nucleotideHaplotype_1_interesting,
			nucleotideHaplotype_2,
			nucleotideHaplotype_2_positions,
			nucleotideHaplotype_2_interesting,
			1.0/(double)n_haplotype_pairs,
			forReturn_thisTranscript
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

	// std::cout << "Expected haplotype pairs: " << n_haplotype_pairs << " // processed: " << consideredHaplotypes << "\n" << std::flush;

	assert((int)n_haplotype_pairs == consideredHaplotypes);

	return forReturn_thisTranscript;
}

std::map<std::string, double> enumeratePeptideHaplotypes_baseLine_plus(int k, const std::map<std::string, std::string>& referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus)
{
	std::map<std::string, double> forReturn;
	checkVariantsConsistentWithReferenceGenome(variants_plus, referenceGenome_plus);

	for(unsigned int transcriptI = 0; transcriptI < transcripts_plus.size(); transcriptI++)
	{
		std::cout << "\r" << transcriptI << " / " << transcripts_plus.size() << "         " << std::flush;
		const transcript& transcript = transcripts_plus.at(transcriptI);
		assert(transcript.strand == '+');
		std::string chromosomeID = transcript.chromosomeID;
		if(referenceGenome_plus.count(chromosomeID) == 0)
			continue;

		// the oneTranscript_* maps are clear'ed within enumeratePeptideHaplotypes_oneTranscript
		std::map<std::string, double> epitopes_this_transcript = enumeratePeptideHaplotypes_baseLine_oneTranscript(k, transcript, referenceGenome_plus, variants_plus);

		for(auto epitope : epitopes_this_transcript)
		{
			if(forReturn.count(epitope.first) == 0)
			{
				forReturn[epitope.first] = epitope.second;
			}
			else
			{
				if(epitope.second > forReturn.at(epitope.first))
				{
					forReturn.at(epitope.first) = epitope.second;
				}
			}
		}
	}
	
	std::cout << "\n";

	return forReturn;
}

std::map<std::string, double> enumeratePeptideHaplotypes_baseLine(int threads, int k, const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants)
{
	assert(threads == 1);

	std::cout << timestamp() << "\t\t Set up minus-strand data..\n" << std::flush;
	std::map<std::string, std::string> referenceGenome_minus = getMinusStrandReferenceGenome(referenceGenome);
	std::map<std::string, std::map<int, variantFromVCF>> variants_minus = getMinusStrandVariants(variants, referenceGenome_minus);
	std::cout << timestamp() << "\t\t ... done.\n" << std::flush;

	checkVariantsConsistentWithReferenceGenome(variants, referenceGenome);
	checkVariantsConsistentWithReferenceGenome(variants_minus, referenceGenome_minus);

	std::vector<transcript> transcripts_plus = getPlusStrandTranscripts(transcripts);
	std::vector<transcript> transcripts_minus = getMinusStrandTranscripts(transcripts, referenceGenome_minus);

	checkTranscriptsTranslate(transcripts_plus, referenceGenome);

	std::map<std::string, double> forReturn = enumeratePeptideHaplotypes_baseLine_plus(k, referenceGenome, transcripts_plus, variants);
	std::map<std::string, double> minus_epitopes = enumeratePeptideHaplotypes_baseLine_plus(k, referenceGenome_minus, transcripts_minus, variants_minus);

	for(auto epitope : minus_epitopes)
	{
		if(forReturn.count(epitope.first) == 0)
		{
			forReturn[epitope.first] = epitope.second;
		}
		else
		{
			if(epitope.second > forReturn.at(epitope.first))
			{
				forReturn.at(epitope.first) = epitope.second;
			}
		}
	}

	return forReturn;
}

