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

std::map<std::string, std::map<int, variantFromVCF>> readVariants(std::string VCF, const std::map<std::string, std::string>& referenceGenome)
{
	std::map<std::string, std::map<int, variantFromVCF>> forReturn;
	std::map<std::string, std::vector<bool>> positionsCoveredByVariant;

	for(auto r : referenceGenome)
	{
		positionsCoveredByVariant[r.first].resize(r.second.length(), false);
	}

	std::ifstream inputStream;
	inputStream.open(VCF.c_str());
	assert(inputStream.is_open());
	bool sawHeader = false;
	std::string line;
	size_t read_variants = 0;
	size_t skipped_variants = 0;
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

				if(positionsCoveredByVariant.count(thisVariant.chromosomeID) == 0)
				{
					std::cerr << "Problem with variant: assumedly wrong reference allele. Did you generate your VCF from the specified reference genome?" << "\n";
					std::cerr << "\t" << "VCF: " << VCF << "\n";
					std::cerr << "\t" << "Chromosome: " << thisVariant.chromosomeID << "\n";
					std::cerr << "\t" << "Position (0-based): " << thisVariant.position << "\n";
				}

				unsigned int thisVariant_lastPosition = thisVariant.position + thisVariant.referenceString.length() - 1;
				bool allReferencePositionsFree = true;
				for(unsigned int pI = thisVariant.position; pI <= thisVariant_lastPosition; pI++)
				{
					if(positionsCoveredByVariant.at(thisVariant.chromosomeID).at(pI))
					{
						allReferencePositionsFree = false;
					}
				}
				//

				if(allReferencePositionsFree)
				{
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

					for(unsigned int pI = thisVariant.position; pI <= thisVariant_lastPosition; pI++)
					{
						positionsCoveredByVariant.at(thisVariant.chromosomeID).at(pI) = true;
					}
				}
				else
				{
					skipped_variants++;
					std::cerr << "Warning: variant " << thisVariant.chromosomeID << ":" << thisVariant.position << " was ignored because it overlapped with other variants in the same VCF.\n" << std::flush;
				}
				read_variants++;

			}
		}
	}

	std::cout << "readVariants(..): Have " << read_variants << " variants; of which " << skipped_variants << " were skipped because they overlapped with existing variants.\n" << std::flush;

	checkVariantsConsistentWithReferenceGenome(forReturn, referenceGenome);
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


std::vector<transcript> getPlusStrandTranscripts(const std::vector<transcript>& transcripts)
{
	std::vector<transcript> forReturn;
	forReturn.reserve(transcripts.size());
	for(transcript t : transcripts)
	{
		if(t.strand == '+')
		{
			forReturn.push_back(t);
		}
	}

	// std::cout << "getPlusStrandTranscripts(..): " << forReturn.size() << " plus-strand transcripts.\n" << std::flush;

	return forReturn;
}

std::vector<transcript> getMinusStrandTranscripts(const std::vector<transcript>& transcripts, const std::map<std::string, std::string>& referenceGenome_minus)
{
	std::vector<transcript> forReturn;
	for(const transcript& t : transcripts)
	{
		if(t.strand == '-')
		{
			transcript t_minus;
			t_minus.chromosomeID = t.chromosomeID;
			t_minus.geneName = t.geneName;
			t_minus.strand = '+';
			t_minus.exons.resize(t.exons.size());

			size_t referenceContigLength = referenceGenome_minus.at(t.chromosomeID).length();

			// make sure exons go from right to left in non-overlapping fashion
			int lastValidExonI = -1;
			for(int exonI = 0; exonI < (int)t.exons.size(); exonI++)
			{
				if(t.exons.at(exonI).valid)
				{
					assert(t.exons.at(exonI).firstPos <= t.exons.at(exonI).lastPos);
					if(lastValidExonI != -1)
					{
						assert(t.exons.at(exonI).lastPos < t.exons.at(lastValidExonI).firstPos);
					}
					lastValidExonI = exonI;

					t_minus.exons.at(exonI).valid = true;

					long long minus_firstPos = referenceContigLength - t.exons.at(exonI).firstPos - 1;
					long long minus_lastPos = referenceContigLength - t.exons.at(exonI).lastPos - 1;
					if(!(minus_firstPos >= minus_lastPos))
					{
						t.print();
					}
					assert(minus_firstPos >= minus_lastPos);
					long long third = minus_firstPos;
					minus_firstPos = minus_lastPos;
					minus_lastPos = third;
					assert(minus_firstPos <= minus_lastPos);

					t_minus.exons.at(exonI).firstPos = minus_firstPos;
					t_minus.exons.at(exonI).lastPos = minus_lastPos;
				}
			}

			forReturn.push_back(t_minus);
		}
	}

	return forReturn;
}



std::map<std::string, std::map<int, variantFromVCF>> getMinusStrandVariants(const std::map<std::string, std::map<int, variantFromVCF>>& variants, const std::map<std::string, std::string>& referenceGenome_minus)
{
	std::map<std::string, std::map<int, variantFromVCF>> forReturn;
	for(auto chromosomeData : variants)
	{
		for(auto positionAndVariant : chromosomeData.second)
		{
			const variantFromVCF& existingVariant = positionAndVariant.second;

			size_t referenceContigLength = referenceGenome_minus.at(existingVariant.chromosomeID).length();
			unsigned int existingVariant_lastVariantPosition = existingVariant.position + countCharacters_noGaps(existingVariant.referenceString) - 1;
			assert(existingVariant_lastVariantPosition < referenceContigLength);

			variantFromVCF minusVariant;
			minusVariant.chromosomeID = existingVariant.chromosomeID;
			minusVariant.position = referenceContigLength - existingVariant_lastVariantPosition - 1;
			minusVariant.referenceString = seq_reverse_complement(existingVariant.referenceString);
			assert(referenceGenome_minus.at(minusVariant.chromosomeID).substr(minusVariant.position, countCharacters_noGaps(minusVariant.referenceString)) == removeGaps(minusVariant.referenceString));
			for(auto sA : existingVariant.sampleAlleles)
			{
				minusVariant.sampleAlleles.push_back(seq_reverse_complement(sA));
			}
			forReturn[chromosomeData.first][minusVariant.position] = minusVariant;
		}
	}

	return forReturn;
}

std::map<std::string, std::string> getMinusStrandReferenceGenome(const std::map<std::string, std::string>& referenceGenome)
{
	std::map<std::string, std::string> forReturn;
	for(auto referenceGenomeEntry : referenceGenome)
	{
		forReturn[referenceGenomeEntry.first] = seq_reverse_complement(referenceGenomeEntry.second);
	}
	return forReturn;
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
			std::string errorString = "Utilities::reverse_char_nucleotide: nucleotide not existing!";
			errorString.push_back(c);
			throw std::runtime_error(errorString);
    }
}

void checkVariantsConsistentWithReferenceGenome(const std::map<std::string, std::map<int, variantFromVCF>>& variants, const std::map<std::string, std::string>& referenceGenome)
{
	for(auto chromosomeData : variants)
	{
		for(auto positionAndVariant : chromosomeData.second)
		{
			const variantFromVCF& variant = positionAndVariant.second;
			assert(variant.chromosomeID == chromosomeData.first);
			assert((int)variant.position == positionAndVariant.first);
			for(auto sA : variant.sampleAlleles)
			{
				assert(sA.length() == variant.referenceString.length());
			}

			std::string supposedReferenceString = referenceGenome.at(variant.chromosomeID).substr(variant.position, countCharacters_noGaps(variant.referenceString));
			std::string variant_referenceString_noGaps = removeGaps(variant.referenceString);
			if(supposedReferenceString != variant_referenceString_noGaps)
			{
				std::cerr << "Problem with variant: assumedly wrong reference allele. Did you generate your VCF from the specified reference genome?" << "\n";
				std::cerr << "\t" << "Chromosome: " << variant.chromosomeID << "\n";
				std::cerr << "\t" << "Position (0-based): " << variant.position << "\n";
				std::cerr << "\t" << "Variant reference string" << ": " << variant.referenceString << "\n";
				std::cerr << "\t" << "Variant reference string, no gaps" << ": " << variant_referenceString_noGaps << "\n";
				std::cerr << "\t" << "Expected reference string" << ": " << supposedReferenceString << "\n";
				std::cerr << std::flush;
			}
			assert(variant_referenceString_noGaps == supposedReferenceString);

		}
	}
}

std::map<std::string, std::map<int, variantFromVCF>> combineVariants(const std::map<std::string, std::map<int, variantFromVCF>>& variants_normalGenome, const std::map<std::string, std::map<int, variantFromVCF>>& additionalVariants_tumourGenome, const std::map<std::string, std::string>& referenceGenome)
{
	std::map<std::string, std::map<int, variantFromVCF>> forReturn;

	checkVariantsConsistentWithReferenceGenome(variants_normalGenome, referenceGenome);
	checkVariantsConsistentWithReferenceGenome(additionalVariants_tumourGenome, referenceGenome);

	std::map<std::string, std::vector<bool>> positionsCoveredByVariant;
	for(auto r : referenceGenome)
	{
		positionsCoveredByVariant[r.first].resize(r.second.length(), false);
	}

	for(auto variantsPerChromosome : additionalVariants_tumourGenome)
	{
		std::string chromosomeID = variantsPerChromosome.first;
		assert(positionsCoveredByVariant.count(chromosomeID));
		for(auto positionAndVariant : variantsPerChromosome.second)
		{
			const variantFromVCF& tumourVariant = positionAndVariant.second;
			unsigned int lastPosition = tumourVariant.position + countCharacters_noGaps(tumourVariant.referenceString) - 1;

			forReturn[tumourVariant.chromosomeID][tumourVariant.position] = tumourVariant;
			for(unsigned int pI = tumourVariant.position; pI <= lastPosition; pI++)
			{
				positionsCoveredByVariant.at(tumourVariant.chromosomeID).at(pI) = true;
			}
		}
	}

	for(auto variantsPerChromosome : variants_normalGenome)
	{
		std::string chromosomeID = variantsPerChromosome.first;
		assert(positionsCoveredByVariant.count(chromosomeID));
		for(auto positionAndVariant : variantsPerChromosome.second)
		{
			const variantFromVCF& normalGenomeVariant = positionAndVariant.second;
			unsigned int lastPosition = normalGenomeVariant.position + countCharacters_noGaps(normalGenomeVariant.referenceString) - 1;

			bool noExistingVariantCover = true;

			for(unsigned int pI = normalGenomeVariant.position; pI <= lastPosition; pI++)
			{
				if(positionsCoveredByVariant.at(normalGenomeVariant.chromosomeID).at(pI))
				{
					noExistingVariantCover = false;
				}
			}

			if(noExistingVariantCover)
			{
				forReturn[normalGenomeVariant.chromosomeID][normalGenomeVariant.position] = normalGenomeVariant;
				for(unsigned int pI = normalGenomeVariant.position; pI <= lastPosition; pI++)
				{
					positionsCoveredByVariant.at(normalGenomeVariant.chromosomeID).at(pI) = true;
				}
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

void transcript::print() const
{
	std::cout << "Transcript " << transcriptID << " / " << geneName << "\n";
	std::cout << "==================================================================\n";
	std::cout << "Chromosome: " << chromosomeID << "\n";
	std::cout << "Strand: " << strand << "\n";
	std::cout << "Exons:\n";
	for(transcriptExon e : exons)
	{
		std::cout << "\t" << e.firstPos << " " << e.lastPos << "\n";
	}
	std::cout << "\n" << std::flush;
}

void checkTranscriptsTranslate(const std::vector<transcript>& transcripts, const std::map<std::string, std::string>& referenceGenome)
{
	for(const transcript& t : transcripts)
	{
		assert(t.strand == '+');
		std::string t_referenceSequence;
		for(const transcriptExon& e : t.exons)
		{
			if(e.valid)
			{
				if(referenceGenome.count(t.chromosomeID) == 0)
				{
					std::cerr << "Unknown chromosome ID: "+t.chromosomeID << "\n" << std::flush;
					throw std::runtime_error("Unknown chromosome ID: "+t.chromosomeID);
				}
				const std::string& chromosomeSequence = referenceGenome.at(t.chromosomeID);
				assert(e.firstPos <= e.lastPos);
				assert(e.firstPos >= 0);
				assert(e.lastPos < chromosomeSequence.length());

				t_referenceSequence += referenceGenome.at(t.chromosomeID).substr(e.firstPos, e.lastPos - e.firstPos + 1);
			}
		}
		std::string translation;
		assert((t_referenceSequence.length() % 3) == 0);
		for(unsigned int i = 0; i < t_referenceSequence.length(); i += 3)
		{
			std::string codon = t_referenceSequence.substr(i, 3);
			assert(codon.length() == 3);
			translation += translateCodon2AA(codon);
		}
		if(translation.back() != '!')
		{
			// std::cout << t.transcriptID << " " << translation << "\n" << std::flush;
		}
		// assert(translation.back() == '!');
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


