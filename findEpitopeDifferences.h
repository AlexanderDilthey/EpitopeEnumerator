/*
 * findDifferentEpitopes.h
 *
 *  Created on: Apr 16, 2017
 *      Author: diltheyat
 */

#ifndef FINDEPITOPEDIFFERENCES_H_
#define FINDEPITOPEDIFFERENCES_H_

#include <set>
#include <map>
#include <vector>
#include <string>
#include <fstream>

#include "EpitopeEnumerator.h"

std::set<std::string> identifyDifferences_naive(const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants_normal,  const std::map<std::string, std::map<int, variantFromVCF>>& variants_tumour, int coreEpitopeLength, int additionalBuffer);
std::set<std::string> identifyDifferences_faster(const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants_normal,  const std::map<std::string, std::map<int, variantFromVCF>>& variants_tumour, int coreEpitopeLength, int additionalBuffer);

void produceDifferencesFile (const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants_normal,  const std::map<std::string, std::map<int, variantFromVCF>>& variants_tumour, int coreEpitopeLength, int additionalBuffer, std::ofstream* outputStream);

#endif /* FINDEPITOPEDIFFERENCES_H_ */
