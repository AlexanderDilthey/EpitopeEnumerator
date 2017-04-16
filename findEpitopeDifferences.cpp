/*
 * findDifferentEpitopes.cpp
 *
 *  Created on: Apr 16, 2017
 *      Author: diltheyat
 */

#include "findEpitopeDifferences.h"
#include "enumerateEpitopes_haplotypePairs.h"
#include "enumerateEpitopes_diff_pairs.h"

#include <assert.h>
#include <iostream>


std::set<std::string> identifyDifferences_faster(const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants_normal,  const std::map<std::string, std::map<int, variantFromVCF>>& variants_tumour, int coreEpitopeLength, int additionalBuffer)
{
	std::map<std::string, std::map<int, variantFromVCF>> variants_combined = combineVariants(variants_normal, variants_tumour, referenceGenome);

	std::set<int> epitopeLengths_tumour;
	std::set<int> epitopeLengths_normal;

	int combinedTumourEpitopeLength = coreEpitopeLength + 2 * additionalBuffer;
	epitopeLengths_tumour.insert(coreEpitopeLength + 2 * additionalBuffer);
	epitopeLengths_normal.insert(coreEpitopeLength);

	std::map<std::string, double> epitopes_normal = enumeratePeptideHaplotypes_baseLine(1, coreEpitopeLength, referenceGenome, transcripts, variants_normal);

	std::set<std::string> ignorePeptides;
	for(auto epitope : epitopes_normal)
	{
		ignorePeptides.insert(epitope.first);
	}

	std::map<int, std::set<std::string>> epitopes_tumour = enumeratePeptideHaplotypes_properFrequencies_easy(1, referenceGenome, transcripts, variants_combined, epitopeLengths_tumour, true, &ignorePeptides, additionalBuffer);


	std::cout << "identifyDifferences_faster: Search for length " << coreEpitopeLength << " + 2x" << additionalBuffer << "\n" << std::flush;
	assert(epitopes_tumour.count(combinedTumourEpitopeLength));

	std::cout << "\t" << epitopes_normal.size() << " normal epitopes.\n";
	std::cout << "\t" << epitopes_tumour.at(combinedTumourEpitopeLength).size() << " tumour epitopes.\n" << std::flush;

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
