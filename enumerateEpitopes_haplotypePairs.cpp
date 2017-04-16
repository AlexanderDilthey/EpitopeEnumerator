/*
 * enumerateEpitopeshaplotypePairs.cpp
 *
 *  Created on: Apr 11, 2017
 *      Author: diltheyat
 */

#include "enumerateEpitopes_haplotypePairs.h"

#include <iostream>
#include <assert.h>

#include "Util.h"
#include <cmath>
#include <omp.h>

std::map<int, std::set<std::string>> enumeratePeptideHaplotypes_properFrequencies_easy(int threads, const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, bool limitToCertainEpitopes)
{
	std::map<int, std::set<std::string>> forReturn;

	std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
	std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;
	enumeratePeptideHaplotypes(threads, referenceGenome, transcripts, variants, haplotypeLengths, p_per_epitope, p_per_epitope_locations);

	for(auto haplotypesOneLength : p_per_epitope)
	{
		int k = haplotypesOneLength.first;
		for(auto fragment : haplotypesOneLength.second)
		{
			std::string peptide = fragment.first;
			double p = fragment.second.first;

			if(limitToCertainEpitopes)
			{
				if(abs(p - 1) <= 1e-5)
				{
					forReturn[k].insert(peptide);
				}
			}
			else
			{
				forReturn[k].insert(peptide);
			}
		}
	}

	return forReturn;
}

void enumeratePeptideHaplotypes(int threads, const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet)
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

	if(threads == 1)
	{
		enumeratePeptideHaplotypes_plus(referenceGenome, transcripts_plus, variants, false, p_per_epitope_forRet, p_per_epitope_locations_forRet);
		enumeratePeptideHaplotypes_plus(referenceGenome_minus, transcripts_minus, variants_minus, true, p_per_epitope_forRet, p_per_epitope_locations_forRet);
	}
	else
	{
		enumeratePeptideHaplotypes_plus_mt(threads, referenceGenome, transcripts_plus, variants, false, p_per_epitope_forRet, p_per_epitope_locations_forRet);
		enumeratePeptideHaplotypes_plus_mt(threads, referenceGenome_minus, transcripts_minus, variants_minus, true, p_per_epitope_forRet, p_per_epitope_locations_forRet);
	}
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

void enumeratePeptideHaplotypes_plus_mt(int threads, const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool invertPositions, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet)
{
	assert(threads > 0);

	checkVariantsConsistentWithReferenceGenome(variants_plus, referenceGenome_plus); // paranoid
	assert(p_per_epitope_forRet.size() == p_per_epitope_locations_forRet.size());

	std::vector<std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>> p_per_epitope_perThread;
	std::vector<std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>> p_per_epitope_locations_perThread;

	p_per_epitope_perThread.resize(threads);
	p_per_epitope_locations_perThread.resize(threads);
	for(auto k : p_per_epitope_forRet)
	{
		assert(p_per_epitope_locations_forRet.count(k.first));
	}

	size_t transcripts_n = transcripts_plus.size();
	omp_set_num_threads(threads);
	#pragma omp parallel for
	for(unsigned int transcriptI = 0; transcriptI < transcripts_plus.size(); transcriptI++)
	{
		int tI = omp_get_thread_num();
		if(tI == 0)
		{
			std::cout << "Thread " << tI << ", transcript " << transcriptI << " / " << transcripts_n << " (this strand)." << std::flush;
		}
		std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> oneTranscript_p_per_epitope;
		std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> oneTranscript_p_per_epitope_locations;

		for(const auto& k : p_per_epitope_forRet)
		{
			oneTranscript_p_per_epitope[k.first].count("");
			oneTranscript_p_per_epitope_locations[k.first].count("");
		}

		const transcript& transcript = transcripts_plus.at(transcriptI);
		assert(transcript.strand == '+');
		std::string chromosomeID = transcript.chromosomeID;
		if(referenceGenome_plus.count(chromosomeID) == 0)
			continue;

		// the oneTranscript_* maps are clear'ed within enumeratePeptideHaplotypes_oneTranscript
		enumeratePeptideHaplotypes_oneTranscript(transcript, referenceGenome_plus, variants_plus, oneTranscript_p_per_epitope, oneTranscript_p_per_epitope_locations);

		for(const auto& k : p_per_epitope_forRet)
		{
			for(const auto& epitope : oneTranscript_p_per_epitope.at(k.first))
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

				if(p_per_epitope_perThread.at(tI)[k.first].count(epitope.first))
				{
					if(p_per_epitope_perThread.at(tI).at(k.first).at(epitope.first).first < epitope.second.first)
					{
						p_per_epitope_perThread.at(tI).at(k.first).at(epitope.first).first = epitope.second.first;
					}
					p_per_epitope_perThread.at(tI).at(k.first).at(epitope.first).second.insert(positions.begin(), positions.end());
				}
				else
				{
					p_per_epitope_perThread.at(tI)[k.first][epitope.first].first = epitope.second.first;
					p_per_epitope_perThread.at(tI).at(k.first)[epitope.first].second = positions;
				}
			}
			for(const auto& epitope : oneTranscript_p_per_epitope_locations.at(k.first))
			{
				for(const auto& location : epitope.second)
				{
					if(p_per_epitope_locations_perThread.at(tI)[k.first][epitope.first].count(location.first))
					{
						if(p_per_epitope_locations_perThread.at(tI).at(k.first).at(epitope.first).at(location.first) < location.second)
						{
							p_per_epitope_locations_perThread.at(tI).at(k.first).at(epitope.first).at(location.first) = location.second;
						}
					}
					else
					{
						p_per_epitope_locations_perThread.at(tI)[k.first][epitope.first][location.first] = location.second;
					}
				}
			}
		}
	}

	for(const auto& k : p_per_epitope_forRet)
	{
		for(int tI = 0; tI < threads; tI++)
		{
			if(p_per_epitope_perThread.at(tI).count(k.first))
			{
				for(const auto& epitope : p_per_epitope_perThread.at(tI).at(k.first))
				{
					const std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>& positions = epitope.second.second;

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

			if(p_per_epitope_locations_perThread.at(tI).count(k.first))
			{
				for(const auto& epitope : p_per_epitope_locations_perThread.at(tI).at(k.first))
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
	bool runningCodon_interesting = false;
	while(currentPos < (int)sequence.length())
	{
		unsigned char thisC = sequence.at(currentPos);
		if((thisC != '-') && (thisC != '_'))
		{
			runningCodon.push_back(thisC);
			runningPositions.push_back(positions.at(currentPos));
			runningInteresting.push_back(interesting.at(currentPos));
		}
		runningCodon_interesting = (runningCodon_interesting || interesting.at(currentPos));

		if(runningCodon.length() == 3)
		{
			std::string translation = translateCodon2AA(runningCodon);
			bool interesting = (runningInteresting.at(0) || runningInteresting.at(1) || runningInteresting.at(2));
			interesting = runningCodon_interesting;
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
				AAs_firstLast.push_back(std::make_pair(minPos, maxPos));
				AAs_interesting.push_back(interesting);
			}

			runningCodon.clear();
			runningPositions.clear();
			runningInteresting.clear();

			runningCodon_interesting = false;
		}
		currentPos++;
	}

	assert(AAs.size() == AAs_firstLast.size());
	assert(AAs.size() == AAs_interesting.size());
	return make_tuple(AAs, AAs_firstLast, AAs_interesting);
}


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
		for(unsigned int i = 0; i < v.sampleAlleles.size(); i++)
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

		for(unsigned int i = 0; i < v.sampleAlleles.size(); i++)
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


