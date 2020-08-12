/*
 * Util.h
 *
 *  Created on: Apr 11, 2017
 *      Author: diltheyat
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <map>
#include <string>
#include <vector>
#include <set>

std::map<std::string, std::string> readFASTA(std::string file, bool fullIdentifier = false);
void eraseNL(std::string& s);

std::vector<std::string> split(const std::string& input, const std::string& delimiter);
int StrtoI(const std::string& s);

double randomDouble();

std::string join(std::vector<std::string> parts, std::string delim);

unsigned int randomNumber(unsigned int max_inclusive);
std::string generateRandomNucleotideSequence(int length);


std::string timestamp();

std::vector<std::string> partitionStringIntokMers(const std::string& str, int k);


std::string seq_reverse_complement(const std::string& sequence);
char reverse_char_nucleotide(char c);

std::string removeGaps(const std::string& in);

unsigned int countCharacters_noGaps(const std::string& S);

void fillTranslationTables();
std::string nucleotide2AA(const std::string& nucleotides);
std::set<std::string> translateAA2Codon(const std::string& AA);
std::string translateCodon2AA(const std::string& codon);
char randomNucleotide();
std::string randomAA();
std::string translateAASequence2Codons(const std::string& AAs);
std::string generateRandomAASequence(int length);



#endif /* UTIL_H_ */
