/*
 * readFiles.cpp
 *
 *  Created on: Mar 31, 2017
 *      Author: diltheyat
 */

#include "readFiles.h"

#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdexcept>
#include <exception>
#include <set>

std::map<std::string, std::map<int, variantFromVCF>> readVariants(std::string VCF)
{
	std::map<std::string, std::map<int, variantFromVCF>> forReturn;

	std::ifstream inputStream;
	inputStream.open(VCF.c_str());
	assert(inputStream.is_open());
	bool sawHeader = false;
	std::string line;
	size_t read_variants = 0;
	while(inputStream.good())
	{
		std::getline(inputStream, line);
		Utilities::eraseNL(line);
		if(line.length())
		{
			if((line.length() >= 2) && (line.substr(0, 2) == "##"))
			{
				continue;
			}

			std::vector<std::string> line_fields = Utilities::split(line, "\t");

			if(line.substr(0, 1) == "#")
			{
				assert(line_fields.at(0) == "#CHROM");
				assert(line_fields.at(1) == "POS");
				assert(line_fields.at(3) == "REF");
				assert(line_fields.at(4) == "ALT");
				assert(line_fields.size() >= 10);
				sawHeader = true;
			}
			else
			{
				assert(sawHeader);

				variantFromVCF thisVariant;
				thisVariant.chromosomeID = line_fields.at(0);
				thisVariant.position = Utilities::StrtoI(line_fields.at(1)) - 1;
				thisVariant.referenceString = line_fields.at(3);

				std::vector<std::string> alternativeAlleles = Utilities::split(line_fields.at(4), ",");
				std::map<int, std::string> number_2_allele;
				number_2_allele[0] = thisVariant.referenceString;
				for(unsigned int aI = 0; aI < alternativeAlleles.size(); aI++)
				{
					std::string allele = alternativeAlleles.at(aI);
					number_2_allele[aI+1] = allele;
				}

				std::vector<std::string> thisSample_alleles_unphased = Utilities::split(line_fields.at(9), "/");
				std::vector<std::string> thisSample_alleles_phased = Utilities::split(line_fields.at(9), "|");
				assert(((thisSample_alleles_unphased.size() == 2) && (thisSample_alleles_phased.size() == 1)) || ((thisSample_alleles_unphased.size() == 1) && (thisSample_alleles_phased.size() == 2)));
				std::vector<std::string> thisSample_alleles = ((thisSample_alleles_unphased.size() == 2) && (thisSample_alleles_phased.size() == 1)) ? thisSample_alleles_unphased : thisSample_alleles_phased;
				assert(thisSample_alleles.size() == 2);

				for(std::string a : thisSample_alleles)
				{
					int a_numeric = Utilities::StrtoI(a);
					assert(number_2_allele.count(a_numeric));
					thisVariant.sampleAlleles.push_back(number_2_allele.at(a_numeric));
				}

				unsigned int maxLength = thisVariant.referenceString.length();
				for(auto a : thisVariant.sampleAlleles)
				{
					if(maxLength < a.length())
						maxLength = a.length();
				}

				extendAsNecessary(thisVariant.referenceString, maxLength);
				for(std::string& a : thisVariant.sampleAlleles)
				{
					extendAsNecessary(a, maxLength);
				}
				assert(forReturn[thisVariant.chromosomeID].count(thisVariant.position) == 0);
				assert(thisVariant.sampleAlleles.size() == 2);
				forReturn[thisVariant.chromosomeID][thisVariant.position] = thisVariant;
				read_variants++;
			}
		}
	}

	std::cout << "readVariants(..): Have " << read_variants << " variants.\n" << std::flush;

	return forReturn;
}

std::vector<transcript> readTranscripts(std::string transcriptsFile)
{
	std::ifstream inputStream;
	inputStream.open(transcriptsFile.c_str());
	assert(inputStream.is_open());
	std::string line;
	size_t read_exons = 0;

	std::map<std::string, transcript> transcript_storage;
	std::map<std::string, std::map<int, transcriptExon>> exons_per_transcript;
	std::set<std::string> ignored_transcriptIDs;

	long long lineI = -1;
	while(inputStream.good())
	{
		lineI++;

		std::getline(inputStream, line);
		Utilities::eraseNL(line);
		if(line.length())
		{
			std::vector<std::string> line_fields = Utilities::split(line, "\t");
			if(line_fields.at(0) != "chr20") // todo remove
				continue;

			if(line_fields.at(2) != "CDS")
				continue;

			int startPos = Utilities::StrtoI(line_fields.at(3)) - 1;
			int stopPos = Utilities::StrtoI(line_fields.at(4)) - 1;
			assert(startPos <= stopPos);

			std::string dataFields_string = line_fields.at(8);
			std::vector<std::string> dataFields_vector = Utilities::split(dataFields_string, ";");
			std::map<std::string, std::string> dataFields;
			for(std::string oneDataField : dataFields_vector)
			{
				std::vector<std::string> thisDataField_vector = Utilities::split(oneDataField, "=");
				assert(thisDataField_vector.size() == 2);
				dataFields[thisDataField_vector.at(0)] = thisDataField_vector.at(1);
			}

			assert(dataFields.count("ID"));
			std::string ID = dataFields.at("ID");

			std::string geneName;
			if(dataFields.count("gene_name"))
			{
				geneName = dataFields.at("gene_name");
			}

			if((line.find("cds_end_NF") != std::string::npos) || (line.find("cds_start_NF") != std::string::npos))
			{
				ignored_transcriptIDs.insert(ID);
				continue;
			}

			if(dataFields.count("transcript_type"))
			{
				if(dataFields.at("transcript_type") != "protein_coding")
				{
					ignored_transcriptIDs.insert(ID);
					continue;
				}
			}

			if(dataFields.count("gene_type"))
			{
				if(dataFields.at("gene_type") != "protein_coding")
				{
					std::cerr << "Warning - gene_type for CDS is " << dataFields.at("gene_type") << " (in file " << transcriptsFile << ", line " << lineI << ")\n" << std::flush;
					ignored_transcriptIDs.insert(ID);
					continue;
				}
				assert(dataFields.at("gene_type") == "protein_coding");
			}
			assert(line_fields.at(6).size() == 1);

			if(transcript_storage.count(ID) == 0)
			{
				transcript_storage[ID].transcriptID = ID;
				transcript_storage[ID].chromosomeID = line_fields.at(0);
				transcript_storage[ID].geneName = geneName;
				transcript_storage[ID].strand = line_fields.at(6).at(0);
			}
			else
			{
				assert(transcript_storage[ID].transcriptID == ID);
				assert(transcript_storage[ID].chromosomeID == line_fields.at(0));
				assert(transcript_storage[ID].geneName == geneName);
				assert(transcript_storage[ID].strand == line_fields.at(6).at(0));
			}

			assert(dataFields.count("exon_number"));
			int exonI = Utilities::StrtoI(dataFields.at("exon_number"));

			assert(exons_per_transcript[ID].count(exonI) == 0);
			transcriptExon exon;
			exon.firstPos = startPos;
			exon.lastPos = stopPos;

			exons_per_transcript[ID][exonI] = exon;

			read_exons++;
		}
	}

	size_t ignored_missingExons = 0;
	size_t ignored_non3Dividable = 0;
	std::vector<transcript> forReturn;
	for(auto transcript_and_id : transcript_storage)
	{
		std::string transcriptID = transcript_and_id.first;
		transcript transcriptCopy = transcript_and_id.second;
		assert(exons_per_transcript.count(transcriptID));
		unsigned int maxCollectedExon = 0;
		for(auto eI : exons_per_transcript.at(transcriptID))
		{
			if(eI.first > maxCollectedExon)
				maxCollectedExon = eI.first;
		}
		transcriptCopy.exons.resize(maxCollectedExon);
		bool allExonsOK = true;
		size_t allExons_length = 0;
		for(unsigned int i = 1; i <= maxCollectedExon; i++)
		{
			/*
			if(!exons_per_transcript.at(transcriptID).count(i))
			{
				allExonsOK = false;
				//std::cerr << "TranscriptID " << transcriptID << " exon " << i << " is missing. Have:\n";
				for(auto e : exons_per_transcript.at(transcriptID))
				{
					//std::cerr << "\t" << e.first << "\n";
				}
				std::cerr << std::flush;
			}
			//assert(exons_per_transcript.at(transcriptID).count(i));
			 *
			 */
			if(exons_per_transcript.at(transcriptID).count(i))
			{
				transcriptCopy.exons.at(i-1) = exons_per_transcript.at(transcriptID).at(i);
				transcriptCopy.exons.at(i-1).valid = true;
				allExons_length += (transcriptCopy.exons.at(i-1).lastPos - transcriptCopy.exons.at(i-1).firstPos + 1);
			}
			else
			{

			}
		}

		if(allExonsOK)
		{
			if((allExons_length % 3) == 0)
			{
				assert((allExons_length % 3) == 0);
				forReturn.push_back(transcriptCopy);
			}
			else
			{
				ignored_non3Dividable++;
			}
		}
		else
		{
			//std::cerr << "Missing exons: " << transcriptID << "\n";
			ignored_missingExons++;
		}

		/*

		if(transcriptCopy.strand == '+')
		{
			int lastPos = -1;
			for(transcriptExon e : transcriptCopy.exons)
			{
				if(e.valid)
				{
					if(lastPos != -1)
					{
						assert(e.firstPos > lastPos);
					}
					lastPos = e.lastPos;
				}
			}
		}

		*/
	}

	std::cout << "readTranscripts(..): Have " << forReturn.size() << " transcripts with " << read_exons << " exons; ignored because of missing exons " << ignored_missingExons << "; ignored because length not multiple of 3: " << ignored_non3Dividable << "; ignored because of other criteria: " << ignored_transcriptIDs.size() << ".\n" << std::flush;

	return forReturn;
}

unsigned int countCharacters_noGaps(const std::string& S)
{
	unsigned int forReturn = 0;
	for(unsigned int i = 0; i < S.length(); i++)
	{
		unsigned char c = S.at(i);
		if((c != '-') && (c != '_'))
			forReturn++;
	}
	return forReturn;
}

std::map<std::string, std::string> codon2AA;
std::string translateCodon2AA(const std::string& codon)
{
	assert(codon.length() == 3);
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
	}

	if(codon2AA.count(codon) == 0)
	{
		throw std::runtime_error("Codon "+codon+" undefined.");
	}

	return codon2AA.at(codon);
}

void extendAsNecessary(std::string& S, unsigned int desiredLength)
{
	assert(S.length() <= desiredLength);
	int missing = desiredLength - S.length();
	if(missing > 0)
	{
		std::string append;
		append.resize(missing, '-');
		S.append(append);
	}

}
