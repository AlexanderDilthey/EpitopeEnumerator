/*
 * enumerateEpitopesdiffproper.h
 *
 *  Created on: Apr 16, 2017
 *      Author: diltheyat
 */

#ifndef ENUMERATEEPITOPES_DIFF_PAIRS_H_
#define ENUMERATEEPITOPES_DIFF_PAIRS_H_

#include <map>
#include <vector>
#include <set>
#include <string>

#include "readFiles.h"

void forBaseline_populateFragmentStorageFromNucleotideHaplotypePair_stopAware_additive(int k, const std::string& sequence_1, const std::vector<int>& positions_1, const std::vector<bool>& interesting_1, const std::string& sequence_2, const std::vector<int>& positions_2, const std::vector<bool>& interesting_2, double p, std::map<std::string, double>& epitopes_and_p);

std::map<std::string, double> enumeratePeptideHaplotypes_baseLine_oneTranscript(int k, const transcript& transcript, const std::map<std::string, std::string>& referenceGenome_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus);
std::map<std::string, double> enumeratePeptideHaplotypes_baseLine_plus(int k, const std::map<std::string, std::string>& referenceGenome_plus, const std::vector<transcript>& transcripts_plus, const std::map<std::string, std::map<int, variantFromVCF>>& variants_plus);
std::map<std::string, double> enumeratePeptideHaplotypes_baseLine(int threads, int k, const std::map<std::string, std::string>& referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants);

#endif /* ENUMERATEEPITOPES_DIFF_PAIRS_H_ */
