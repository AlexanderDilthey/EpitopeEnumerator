/*
 * enumerateEpitopesnoHaplotypePairs.h
 *
 *  Created on: Apr 11, 2017
 *      Author: diltheyat
 */

#ifndef ENUMERATEEPITOPES_NOHAPLOTYPEPAIRS_H_
#define ENUMERATEEPITOPES_NOHAPLOTYPEPAIRS_H_

#include <vector>
#include <map>
#include <set>

#include "enumerateEpitopes_haplotypePairs.h"
#include "EpitopeEnumerator.h"

void enumeratePeptideHaplotypes_improperFrequencies_oneTranscript(const transcript& transcript, const std::map<std::string, std::string> referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet);
void enumeratePeptideHaplotypes_improperFrequencies_plus(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool invertPositions, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet);
void enumeratePeptideHaplotypes_improperFrequencies(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet);

std::map<int, std::set<std::string>> enumeratePeptideHaplotypes_improperFrequencies_easy(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, bool limitToCertainEpitopes);


#endif /* ENUMERATEEPITOPES_NOHAPLOTYPEPAIRS_H_ */
