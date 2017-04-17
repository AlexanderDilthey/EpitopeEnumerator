/*
 * Util.cpp
 *
 *  Created on: Apr 11, 2017
 *      Author: diltheyat
 */

#include "Util.h"

#include <fstream>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <assert.h>
#include <sstream>
#include <time.h>

int StrtoI(const std::string& s)
{
	  std::stringstream ss(s);
	  int i;
	  ss >> i;
	  return i;
}

void eraseNL(std::string& s)
{
	while (!s.empty() && ((s[s.length()-1] == '\r') || (s[s.length()-1] == '\n'))) {
		    s.erase(s.length()-1);
	}
}

std::string timestamp()
{
    time_t ltime; /* calendar time */
    ltime=time(NULL); /* get current cal time */
    char* timeCString = asctime( localtime(&ltime) );
    std::string forReturn(timeCString);
    return " [ "+ forReturn.substr(0, forReturn.length() - 1)+" ] ";
}


std::map<std::string, std::string> readFASTA(std::string file, bool fullIdentifier)
{
	std::map<std::string, std::string> forReturn;

	std::ifstream FASTAstream;

	FASTAstream.open(file.c_str());
	if(! FASTAstream.is_open())
	{
		throw std::runtime_error("readFASTA(): Cannot open file "+file);
	}

	while(FASTAstream.good())
	{
		std::string line;
		size_t lineCounter = 0;

		std::string currentSequenceIdentifier;
		while(FASTAstream.good())
		{
			std::getline(FASTAstream, line);
			eraseNL(line);

			lineCounter++;

			if(line.substr(0, 1) == ">")
			{
				std::string ident = line.substr(1);

				if(! fullIdentifier)
				{
					for(size_t i = 0; i < ident.size(); i++)
					{
						if((ident.at(i) == ' ') || (ident.at(i) == '\t'))
						{
							ident = ident.substr(0, i);
							break;
						}
					}
				}
				
				currentSequenceIdentifier = ident;
				assert(forReturn.count(ident) == 0);
			}
			else
			{
				forReturn[currentSequenceIdentifier] += line;
			}
		}
	}

	return forReturn;
}

std::string removeGaps(const std::string& in)
{
	std::string out;
	out.reserve(in.size());
	for(size_t i = 0; i < in.size(); i++)
	{
		if((in.at(i) != '_') && (in.at(i) != '-'))
		{
			out.push_back(in.at(i));
		}
	}
	return out;
}

std::vector<std::string> split(const std::string& input, const std::string& delimiter)
{
	std::vector<std::string> output;
	if(input.length() == 0)
	{
		return output;
	}

	if(delimiter == "")
	{
		output.reserve(input.size());
		for(unsigned int i = 0; i < input.length(); i++)
		{
			output.push_back(input.substr(i, 1));
		}
	}
	else
	{
		if(input.find(delimiter) == std::string::npos)
		{
			output.push_back(input);
		}
		else
		{
			int s = 0;
			size_t p = input.find(delimiter);

			do {
				output.push_back(input.substr(s, p - s));
				s = p + delimiter.size();
				p = input.find(delimiter, s);
			} while (p != std::string::npos);
			output.push_back(input.substr(s));
		}
	}

	return output;
}


double randomDouble()
{
	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);
	return f;
}



unsigned int countCharacters_noGaps(const std::string& S)
{
	unsigned int forReturn = 0;
	for(unsigned int i = 0; i < S.length(); i++)
	{
		char c = S.at(i);
		if((c != '-') && (c != '_'))
			forReturn++;
	}
	return forReturn;
}

std::map<std::string, std::string> codon2AA;
std::map<std::string, std::set<std::string>> AA2codon;

void fillTranslationTables()
{
	if(codon2AA.size() == 0)
	{
		codon2AA["CCT"] = "P";
		codon2AA["CAC"] = "H";
		codon2AA["CTG"] = "L";
		codon2AA["CAG"] = "Q";
		codon2AA["GGA"] = "G";
		codon2AA["CGG"] = "R";
		codon2AA["TAT"] = "Y";
		codon2AA["GAT"] = "D";
		codon2AA["ATT"] = "I";
		codon2AA["AAC"] = "N";
		codon2AA["CCG"] = "P";
		codon2AA["TCC"] = "S";
		codon2AA["CGA"] = "R";
		codon2AA["GTG"] = "V";
		codon2AA["GTC"] = "V";
		codon2AA["CTA"] = "L";
		codon2AA["AAG"] = "K";
		codon2AA["CGT"] = "R";
		codon2AA["TTA"] = "L";
		codon2AA["AAT"] = "N";
		codon2AA["ACA"] = "T";
		codon2AA["GGT"] = "G";
		codon2AA["GGC"] = "G";
		codon2AA["GCC"] = "A";
		codon2AA["GCA"] = "A";
		codon2AA["GAG"] = "E";
		codon2AA["CAT"] = "H";
		codon2AA["TGT"] = "C";
		codon2AA["ATG"] = "M";
		codon2AA["ATC"] = "I";
		codon2AA["TTC"] = "F";
		codon2AA["TTT"] = "F";
		codon2AA["CAA"] = "Q";
		codon2AA["AGC"] = "S";
		codon2AA["TGG"] = "W";
		codon2AA["GCT"] = "A";
		codon2AA["GAC"] = "D";
		codon2AA["CGC"] = "R";
		codon2AA["CCC"] = "P";
		codon2AA["TTG"] = "L";
		codon2AA["ACT"] = "T";
		codon2AA["ATA"] = "I";
		codon2AA["AGA"] = "R";
		codon2AA["AGT"] = "S";
		codon2AA["CTT"] = "L";
		codon2AA["GCG"] = "A";
		codon2AA["AGG"] = "R";
		codon2AA["AAA"] = "K";
		codon2AA["ACG"] = "T";
		codon2AA["TGA"] = "!";
		codon2AA["CCA"] = "P";
		codon2AA["GTT"] = "V";
		codon2AA["GGG"] = "G";
		codon2AA["TCG"] = "S";
		codon2AA["GTA"] = "V";
		codon2AA["TCA"] = "S";
		codon2AA["CTC"] = "L";
		codon2AA["TGC"] = "C";
		codon2AA["TAC"] = "Y";
		codon2AA["GAA"] = "E";
		codon2AA["TAG"] = "!";
		codon2AA["ACC"] = "T";
		codon2AA["TAA"] = "!";
		codon2AA["TCT"] = "S";

		for(auto c2A : codon2AA)
		{
			AA2codon[c2A.second].insert(c2A.first);
		}
	}
}
std::set<std::string> translateAA2Codon(const std::string& AA)
{
	fillTranslationTables();

	if(AA2codon.count(AA) == 0)
	{
		throw std::runtime_error("AA "+AA+" undefined.");
	}

	return AA2codon.at(AA);
}

std::string translateAASequence2Codons(const std::string& AAs)
{
	std::string forReturn;
	forReturn.reserve(AAs.size() * 3);
	for(unsigned int AAi = 0; AAi < AAs.size(); AAi++)
	{
		std::string AA = AAs.substr(AAi, 1);
		std::set<std::string> codons = translateAA2Codon(AA);
		if(codons.size() > 1)
		{
			std::vector<std::string> codons_vec(codons.begin(), codons.end());
			int n = rand() % codons_vec.size();
			assert((n >= 0) && (n < (int)codons_vec.size()));
			forReturn.append(codons_vec.at(n));
		}
		else
		{
			std::string codon = *(codons.begin());
			forReturn.append(codon);
		}
	}

	return forReturn;
}

std::string translateCodon2AA(const std::string& codon)
{
	fillTranslationTables();
	assert(codon.length() == 3);

	if(codon.find("N") != std::string::npos)
	{
		return "?";
	}
	
	if(codon2AA.count(codon) == 0)
	{
		throw std::runtime_error("Codon "+codon+" undefined.");
	}

	return codon2AA.at(codon);
}



std::string seq_reverse_complement(const std::string& sequence)
{
	int length = sequence.size();
	std::string forReturn;
	forReturn.resize(length);
    for(int k=0; k < length; k++)
    {
        forReturn[k] = reverse_char_nucleotide(sequence.at(length-k-1));
    }
    return forReturn;
}

char reverse_char_nucleotide(char c)
{
    switch (c)
    {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'N':
			return 'N';
		case 'a':
			return 't';
		case 'c':
			return 'g';
		case 'g':
			return 'c';
		case 't':
			return 'a';
		case 'n':
			return 'n';
		case '_':
			return '_';
		case '-':
			return '-';
		case '*':
			return '*';
		default:
			std::string errorString = "reverse_char_nucleotide: nucleotide not existing!";
			errorString.push_back(c);
			throw std::runtime_error(errorString);
    }
}


char randomNucleotide()
{
	char nucleotides[4] = {'A', 'C', 'G', 'T'};
	int n = rand() % 4;
	assert((n >= 0) && (n <= 3));
	return nucleotides[n];
}

std::vector<std::string> AAs_for_randomAA;
std::string randomAA()
{
	if(AAs_for_randomAA.size() == 0)
	{
		fillTranslationTables();
		for(auto AA : AA2codon)
		{
			if(AA.first != "!")
				AAs_for_randomAA.push_back(AA.first);
		}
	}

	int n = rand() % AAs_for_randomAA.size();
	assert((n >= 0) && (n < (int)AAs_for_randomAA.size()));
	return AAs_for_randomAA.at(n);
}

std::string generateRandomAASequence(int length)
{
	std::string forReturn;
	forReturn.reserve(length);
	for(int i = 0; i < length; i++)
	{
		forReturn.append(randomAA());
	}
	assert(forReturn.length() == length);
	return forReturn;
}





std::vector<std::string> partitionStringIntokMers(const std::string& str, int k)
{
	std::vector<std::string> forReturn;
	if((int)str.length() >= k)
	{
		forReturn.reserve(str.length() - k + 1);
		for(int i = 0; i <= (str.length() - k); i++)
		{
			std::string kMer = str.substr(i, k);
			assert((int)kMer.length() == k);
			forReturn.push_back(kMer);
		}
	}
	return forReturn;
}


unsigned int randomNumber(unsigned int max_inclusive)
{
	int n = rand() % (max_inclusive+1);
	assert((n >= 0) && (n <= max_inclusive));
	return n;
}

std::string generateRandomNucleotideSequence(int length)
{
	std::string forReturn;
	forReturn.reserve(length);
	for(int i = 0; i < length; i++)
	{
		forReturn.push_back(randomNucleotide());
	}
	assert(forReturn.length() == length);
	return forReturn;
}





