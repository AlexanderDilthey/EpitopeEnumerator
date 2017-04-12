/*
 * enumerateEpitopeshaplotypePairs.h
 *
 *  Created on: Apr 11, 2017
 *      Author: diltheyat
 */

#ifndef ENUMERATEEPITOPES_HAPLOTYPEPAIRS_H_
#define ENUMERATEEPITOPES_HAPLOTYPEPAIRS_H_

#include <vector>
#include <map>
#include <set>
#include <utility>
#include <algorithm>

#include "readFiles.h"

void populateFragmentStorageFromNucleotideHaplotypePair_stopAware_additive(const std::string& sequence_1, const std::vector<int>& positions_1, const std::vector<bool>& interesting_1, const std::string& sequence_2, const std::vector<int>& positions_2, const std::vector<bool>& interesting_2, double p, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet);
void enumeratePeptideHaplotypes_oneTranscript(const transcript& transcript, const std::map<std::string, std::string> referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet);
void enumeratePeptideHaplotypes_plus(const std::map<std::string, std::string> referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus, bool invertPositions, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet);
void enumeratePeptideHaplotypes(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>& p_per_epitope_forRet, std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>>& p_per_epitope_locations_forRet);


using fragmentT = std::tuple<std::string, std::vector<std::pair<int, int>>, std::vector<bool>>;
std::vector<fragmentT> AAHaplotypeIntoFragments(int k, const fragmentT& haplotype);
fragmentT AAHaplotypeFromSequence_stopAware(const std::string& sequence, const std::vector<int>& positions, const std::vector<bool>& interesting);
void printFragment(const fragmentT& f);


std::tuple<std::string, std::vector<int>, std::vector<bool>, std::vector<std::string>, std::vector<std::vector<int>>, std::vector<std::vector<bool>>> get_reference_and_variantAlleles(const variantFromVCF& v, unsigned int startReferencePos, unsigned int lastReferencePos);

std::map<int, std::set<std::string>> enumeratePeptideHaplotypes_properFrequencies_easy(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths, bool limitToCertainEpitopes);

#endif /* ENUMERATEEPITOPES_HAPLOTYPEPAIRS_H_ */
