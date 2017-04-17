/*
 * enumerateEpitopesnoHaplotypePairs.cpp
 *
 *  Created on: Apr 11, 2017
 *      Author: diltheyat
 */

#include <iostream>
#include <assert.h>

#include "enumerateEpitopes_noHaplotypePairs.h"
#include "Util.h"


std::map<int, std::set<std::string>> enumeratePeptideHaplotypes_improperFrequencies_easy(const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, bool limitToCertainEpitopes, const std::set<std::string>* ignoreCoreEpitopes, int corePadding)
{
	std::map<int, std::set<std::string>> forReturn;

	std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
	enumeratePeptideHaplotypes_improperFrequencies(referenceGenome, transcripts, variants, haplotypeLengths, p_per_epitope, ignoreCoreEpitopes, corePadding);

	for(auto haplotypesOneLength : p_per_epitope)
	{
		int k = haplotypesOneLength.first;
		forReturn[k].count("");

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

void enumeratePeptideHaplotypes_improperFrequencies_oneTranscript(const transcript& transcript, const std::map<std::string, std::string>& referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, const std::set<std::string>* ignoreCoreEpitopes, int corePadding)
{
	for(auto k_and_stored_fragments : p_per_epitope_forRet)
	{
		k_and_stored_fragments.second.clear();

		int peptideHaplotypeLength = k_and_stored_fragments.first;

		/*
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
		using runningHaplotypePairKey = std::pair<runningHaplotypeKey, runningHaplotypeKey>;

		 */
		using runningHaplotypeKey = std::tuple<std::string, std::vector<std::pair<int, int>>, std::vector<bool>, std::string, std::vector<bool>, std::vector<int>>;

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
					AApart_firstLast.push_back(std::make_pair(codon_firstPos, codon_lastPos));

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

			void shortenAA(int peptideHaplotypeLength, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& fragmentsStore, size_t sampleHaplotypes_size, size_t openedHaplotypes, bool have_deleted_haplotype, const std::set<std::string>* ignoreCoreEpitopes, int corePadding)
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
					assert((int)extract_AA.length() == peptideHaplotypeLength);

					bool doStore = true;

					if(ignoreCoreEpitopes != 0)
					{
						assert(corePadding >= 0);
						assert(extract_AA.length() > 2 * corePadding);
						std::string coreEpitope = extract_AA.substr(corePadding, extract_AA.length() - 2 * corePadding);
						assert(coreEpitope.length() > 0);
						if(ignoreCoreEpitopes->count(coreEpitope))
						{
							doStore = false;
						}
					}

					if(doStore)
					{
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
						fragmentsStore.at(peptideHaplotypeLength).at(extract_AA).second.insert(std::make_pair(extract_positions, extract_interesting));
					}
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
		size_t openedHaplotypes = 1;
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

			/*
			std::vector<std::vector<bool>> sampleAlleles_perCharacter_interesting;
			sampleAlleles_perCharacter_interesting.reserve(sampleAlleles.size());
			for(bool oneSampleAllele_isInteresting : sampleAlleles_interesting)
			{
				std::vector<bool> oneSampleAllele_interesting;
				oneSampleAllele_interesting.resize(referenceSequence.length(), oneSampleAllele_isInteresting);
				sampleAlleles_perCharacter_interesting.push_back(oneSampleAllele_interesting);
			}
			*/

			bool heterozygous = (sampleAlleles.at(0) != sampleAlleles.at(1)) || (sampleAlleles_interesting.at(0) != sampleAlleles_interesting.at(1)) || (sampleAlleles_refCoordinates.at(0) != sampleAlleles_refCoordinates.at(1));
			size_t existingSampleHaplotypes_maxI = sampleHaplotypes.size();
			if(! heterozygous)
			{
				assert(sampleAlleles.at(0) == sampleAlleles.at(1));
				for(unsigned int existingHaplotypeI = 0; existingHaplotypeI < existingSampleHaplotypes_maxI; existingHaplotypeI++)
				{
					sampleHaplotypes.at(existingHaplotypeI).extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_interesting.at(0));
				}
			}
			else
			{
				for(unsigned int existingHaplotypeI = 0; existingHaplotypeI < existingSampleHaplotypes_maxI; existingHaplotypeI++)
				{
					runningHaplotype preExtensionHaplotype = sampleHaplotypes.at(existingHaplotypeI);

					sampleHaplotypes.at(existingHaplotypeI).extendWithNucleotides(sampleAlleles.at(0), sampleAlleles_refCoordinates.at(0), sampleAlleles_interesting.at(0));

					sampleHaplotypes.push_back(preExtensionHaplotype);
					sampleHaplotypes.back().extendWithNucleotides(sampleAlleles.at(1), sampleAlleles_refCoordinates.at(1), sampleAlleles_interesting.at(1));
				}
				
				if(openedHaplotypes < 1e6)
				{	
					openedHaplotypes *= 2;
				}
			}

			/*
			if(!(sampleHaplotypes.size() <= openedHaplotypes))
			{
				std::cerr << "sampleHaplotypes.size()" << ": " << sampleHaplotypes.size() << "\n";
				std::cerr << "openedHaplotypes" << ": " << openedHaplotypes << "\n";
				std::cerr << std::flush;
			}
			*/
			assert(sampleHaplotypes.size() <= openedHaplotypes);


			//if(sampleHaplotypes.size() > 1)
			//	std::cout << "sampleHaplotypes.size()" << ": " << sampleHaplotypes.size() << "\n" << std::flush;


			std::map<fragmentT, double> newFragments_thisExtension;
			for(runningHaplotype& h : sampleHaplotypes)
			{
				h.shortenAA(peptideHaplotypeLength, p_per_epitope_forRet, sampleHaplotypes.size(), openedHaplotypes, have_deleted_haplotype, ignoreCoreEpitopes, corePadding);
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
				if(sampleHaplotypes.size() >= 1)
				{
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
}

void enumeratePeptideHaplotypes_improperFrequencies_plus(const std::map<std::string, std::string>& referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool invertPositions, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, const std::set<std::string>* ignoreCoreEpitopes, int corePadding)
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

		std::cout << "\r" << transcriptI << " / " << transcripts_plus.size() << "         " << std::flush;

		// the oneTranscript_* maps are clear'ed within enumeratePeptideHaplotypes_oneTranscript
		enumeratePeptideHaplotypes_improperFrequencies_oneTranscript(transcript, referenceGenome_plus, variants_plus, oneTranscript_p_per_epitope_forRet, ignoreCoreEpitopes, corePadding);

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
	
	std::cout << "\n";	
}

void enumeratePeptideHaplotypes_improperFrequencies(const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, const std::set<std::string>* ignoreCoreEpitopes, int corePadding)
{
	std::cout << timestamp() << "\t\t Set up minus-strand data..\n" << std::flush;

	std::map<std::string, std::string> referenceGenome_minus = getMinusStrandReferenceGenome(referenceGenome);
	std::map<std::string, std::map<int, variantFromVCF>> variants_minus = getMinusStrandVariants(variants, referenceGenome_minus);
	
	std::cout << timestamp() << "\t\t ... done.\n" << std::flush;

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

	enumeratePeptideHaplotypes_improperFrequencies_plus(referenceGenome, transcripts_plus, variants, false, p_per_epitope_forRet, ignoreCoreEpitopes, corePadding);
	enumeratePeptideHaplotypes_improperFrequencies_plus(referenceGenome_minus, transcripts_minus, variants_minus, true, p_per_epitope_forRet, ignoreCoreEpitopes, corePadding);
}
