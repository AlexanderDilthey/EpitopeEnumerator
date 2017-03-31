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

				for(std::string a : thisSample_alleles_unphased)
				{
					int a_numeric = Utilities::StrtoI(a);
					assert(number_2_allele.count(a_numeric));
					thisVariant.sampleAlleles.push_back(number_2_allele.at(a_numeric));
				}

				assert(forReturn[thisVariant.chromosomeID].count(thisVariant.position) == 0);
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

	long long lineI = -1;
	while(inputStream.good())
	{
		lineI++;

		std::getline(inputStream, line);
		Utilities::eraseNL(line);
		if(line.length())
		{
			std::vector<std::string> line_fields = Utilities::split(line, "\t");
			if(line_fields.at(0) != "chr21") // todo remove
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

			if(dataFields.count("gene_type"))
			{
				if(dataFields.at("gene_type") != "protein_coding")
				{
					std::cerr << "Warning - gene_type for CDS is " << dataFields.at("gene_type") << " (in file " << transcriptsFile << ", line " << lineI << ")\n" << std::flush;
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

	std::vector<transcript> forReturn;
	for(auto transcript_and_id : transcript_storage)
	{
		std::string transcriptID = transcript_and_id.first;
		transcript transcriptCopy = transcript_and_id.second;
		assert(exons_per_transcript.count(transcriptID));
		unsigned int collectedExons = exons_per_transcript.at(transcriptID).size();
		transcriptCopy.exons.resize(collectedExons);
		for(unsigned int i = 0; i < collectedExons; i++)
		{
			assert(exons_per_transcript.at(transcriptID).count(i));
			transcriptCopy.exons.at(i) = exons_per_transcript.at(transcriptID).at(i);
		}

		forReturn.push_back(transcriptCopy);
	}

	std::cout << "readTranscripts(..): Have " << forReturn.size() << " transcripts with " << read_exons << " exons.\n" << std::flush;

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
