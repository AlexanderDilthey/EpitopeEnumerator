/*
 * findDifferentEpitopes.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: diltheyat
 */

#include "findEpitopeDifferences.h"
#include "enumerateEpitopes_haplotypePairs.h"
#include "enumerateEpitopes_diff_pairs.h"
#include "enumerateEpitopes_noHaplotypePairs.h"
#include "Util.h"

#include <assert.h>
#include <iostream>


void produceDifferencesFile (const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants_normal,  const std::map<std::string, std::map<int, variantFromVCF>>& variants_tumour, int coreEpitopeLength, int additionalBuffer, std::ofstream* outputStream)
{
	std::cout << timestamp() << "\t Combine variants.\n" << std::flush;

	std::map<std::string, std::map<int, variantFromVCF>> variants_combined = combineVariants(variants_normal, variants_tumour, referenceGenome);

	std::cout << timestamp() << "\t ... done.\n" << std::flush;
		
	std::cout << timestamp() << "\t Start germline epitope enumeration...\n" << std::flush;

	std::map<std::string, double> epitopes_normal = enumeratePeptideHaplotypes_baseLine(1, coreEpitopeLength, referenceGenome, transcripts, variants_normal, false);
	std::set<std::string> ignorePeptides;
	for(auto epitope : epitopes_normal)
	{
		ignorePeptides.insert(epitope.first);
	}

	std::set<std::string> chromosomes;
	for(const transcript& t : transcripts)
	{
		chromosomes.insert(t.chromosomeID);
	}


	std::map<std::string, double> max_p_per_epitope;
	std::map<std::string, std::map<std::string, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>> locations_per_epitope;

	for(std::string chromosomeID : chromosomes)
	{
		std::vector<transcript> transcripts_thisChromosome;
		transcripts_thisChromosome.reserve(transcripts.size());
		for(const transcript& t : transcripts)
		{
			if(t.chromosomeID == chromosomeID)
			{
				transcripts_thisChromosome.push_back(t);
			}
		}

		int tumour_epitope_length = coreEpitopeLength + 2 * additionalBuffer;
		std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
		enumeratePeptideHaplotypes_improperFrequencies(referenceGenome, transcripts_thisChromosome, variants_combined, {tumour_epitope_length}, p_per_epitope, &ignorePeptides, additionalBuffer);
		assert(p_per_epitope.count(tumour_epitope_length));

		for(auto foundEpitope : p_per_epitope.at(tumour_epitope_length))
		{
			std::string epitope = foundEpitope.first;
			double maxP = foundEpitope.second.first;
			const std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>& locations = foundEpitope.second.second;

			if(max_p_per_epitope.count(epitope) == 0)
			{
				max_p_per_epitope[epitope] = maxP;
			}
			else
			{
				if(max_p_per_epitope.at(epitope) < maxP)
				{
					max_p_per_epitope.at(epitope) = maxP;
				}
			}

			for(auto location : locations)
			{
				locations_per_epitope[epitope][chromosomeID].insert(location);
			}
		}
	}

	for(auto epitope_and_p : max_p_per_epitope)
	{
		const std::string& epitope = epitope_and_p.first;
		int locationIndex = 0;
		for(auto locationsPerChromosome : locations_per_epitope.at(epitope))
		{
			const std::string& chromosomeID = locationsPerChromosome.first;

			for(const std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>& location : locations_per_epitope.at(epitope).at(chromosomeID))
			{
				*outputStream <<epitope << "\t" <<
								coreEpitopeLength << "\t" <<
								additionalBuffer << "\t" <<
								locationIndex << "\t" <<
								chromosomeID << "\t";

								for(auto interesting : location.second)
								{
									*outputStream << (int)interesting;
								}
								*outputStream << "\t";

								bool firstPosition = true;
								for(auto position : location.first)
								{
									if(not firstPosition)
									{
										*outputStream << ";";
									}
									*outputStream << position.first << "-" << position.second;
									firstPosition = false;
								}
								*outputStream << "\t";

				*outputStream << max_p_per_epitope.at(epitope) << "\n";

				locationIndex++;
			}
		}
	}
}


std::set<std::string> identifyDifferences_faster(const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants_normal,  const std::map<std::string, std::map<int, variantFromVCF>>& variants_tumour, int coreEpitopeLength, int additionalBuffer)
{

	std::cout << timestamp() << "\t Combine variants.\n" << std::flush;


	std::map<std::string, std::map<int, variantFromVCF>> variants_combined = combineVariants(variants_normal, variants_tumour, referenceGenome);

	std::cout << timestamp() << "\t ... done.\n" << std::flush;


	std::set<int> epitopeLengths_tumour;
	std::set<int> epitopeLengths_normal;

	int combinedTumourEpitopeLength = coreEpitopeLength + 2 * additionalBuffer;
	epitopeLengths_tumour.insert(coreEpitopeLength + 2 * additionalBuffer);
	epitopeLengths_normal.insert(coreEpitopeLength);

	std::cout << timestamp() << "\t Start germline epitope enumeration...\n" << std::flush;

	std::map<std::string, double> epitopes_normal = enumeratePeptideHaplotypes_baseLine(1, coreEpitopeLength, referenceGenome, transcripts, variants_normal, false);

	std::set<std::string> ignorePeptides;
	for(auto epitope : epitopes_normal)
	{
		ignorePeptides.insert(epitope.first);
	}

	std::cout << timestamp() << "\t Start tumour epitope enumeration...\n" << std::flush;

	std::map<int, std::set<std::string>> epitopes_tumour = enumeratePeptideHaplotypes_improperFrequencies_easy(referenceGenome, transcripts, variants_combined, epitopeLengths_tumour, true, &ignorePeptides, additionalBuffer);

	std::cout << timestamp() << "\t identifyDifferences_faster(..): Searched for length " << coreEpitopeLength << " + 2x " << additionalBuffer << "\n" << std::flush;
	assert(epitopes_tumour.count(combinedTumourEpitopeLength));

	std::cout << "\t" << epitopes_normal.size() << " normal epitopes.\n";
	std::cout << "\t" << epitopes_tumour.at(combinedTumourEpitopeLength).size() << " tumour epitopes of length " << combinedTumourEpitopeLength << " with a tumour-exclusive core of " << coreEpitopeLength << " AAs.\n" << std::flush;

	std::set<std::string> forReturn;

	for(auto epitope : epitopes_tumour.at(combinedTumourEpitopeLength))
	{
		forReturn.insert(epitope);
	}

	return forReturn;
}

std::set<std::string> identifyDifferences_naive(const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants_normal,  const std::map<std::string, std::map<int, variantFromVCF>>& variants_tumour, int coreEpitopeLength, int additionalBuffer)
{
	std::map<std::string, std::map<int, variantFromVCF>> variants_combined = combineVariants(variants_normal, variants_tumour, referenceGenome);

	std::set<int> epitopeLengths_tumour;
	std::set<int> epitopeLengths_normal;

	epitopeLengths_tumour.insert(coreEpitopeLength + 2 * additionalBuffer);
	epitopeLengths_normal.insert(coreEpitopeLength);

	std::map<int, std::set<std::string>> epitopes_normal = enumeratePeptideHaplotypes_properFrequencies_easy(1, referenceGenome, transcripts, variants_normal, epitopeLengths_normal, false);
	std::map<int, std::set<std::string>> epitopes_tumour = enumeratePeptideHaplotypes_properFrequencies_easy(1, referenceGenome, transcripts, variants_combined, epitopeLengths_tumour, true);

	int combinedTumourEpitopeLength = coreEpitopeLength + 2 * additionalBuffer;

	std::cout << "Search for length " << coreEpitopeLength << " + 2x" << additionalBuffer << "\n" << std::flush;
	assert(epitopes_normal.count(coreEpitopeLength));
	assert(epitopes_tumour.count(combinedTumourEpitopeLength));

	std::cout << "\t" << epitopes_normal.at(coreEpitopeLength).size() << " normal epitopes.\n";
	std::cout << "\t" << epitopes_tumour.at(combinedTumourEpitopeLength).size() << " tumour epitopes.\n" << std::flush;

	std::set<std::string> forReturn;

	for(const std::string& extendedTumourEpitope : epitopes_tumour.at(combinedTumourEpitopeLength))
	{
		std::string coreEpitope = extendedTumourEpitope.substr(additionalBuffer, coreEpitopeLength);
		assert((int)coreEpitope.length() == coreEpitopeLength);

		if(epitopes_normal.at(coreEpitopeLength).count(coreEpitope) == 0)
		{
			forReturn.insert(extendedTumourEpitope);
		}
	}

	std::cout << "identifyDifferences_naive(..): " << forReturn.size() << " tumour-only epitopes.\n" << std::flush;

	return forReturn;
}
