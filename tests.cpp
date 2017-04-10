/*
 * tests.cpp
 *
 *  Created on: Apr 9, 2017
 *      Author: diltheyat
 */

#include "tests.h"

void randomTests()
{
	// recovery of expected AA-mers
	for(unsigned int iteration = 0; iteration < 100; iteration++)
	{
		std::cout << "Testing iteration " << iteration << "\n" << std::flush;
		std::string AAsequence = generateRandomAASequence(20);
		std::string nucleotideSequnce = translateAASequence2Codons(AAsequence);

		std::cout << "\t" << AAsequence << "\n" << std::flush;

		std::set<int> AA_mers = {5, 7, 12};
		for(char strand : {'+', '-'})
		{
			// all as 1 exon, no variants
			{
				std::map<std::string, std::string> referenceGenome;
				referenceGenome["chr"] = (strand == '+') ? nucleotideSequnce : seq_reverse_complement(nucleotideSequnce);

				std::map<std::string, std::map<int, variantFromVCF>> variants;

				transcriptExon oneExon;
				oneExon.valid = true;
				oneExon.firstPos = 0;
				oneExon.lastPos = nucleotideSequnce.length() - 1;

				transcript oneTranscript;
				oneTranscript.chromosomeID = "chr";
				oneTranscript.geneName = "testGene";
				oneTranscript.strand = strand;
				oneTranscript.exons = {oneExon};

				std::vector<transcript> transcripts = {oneTranscript};

				std::map<int, std::map<std::string, std::map<double, std::set<std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>> haplotypeStore;
				for(auto k : AA_mers)
				{
					haplotypeStore[k].count("");
				}

				enumeratePeptideHaplotypes(referenceGenome, transcripts, variants, true, haplotypeStore);

				for(auto k : AA_mers)
				{
					std::vector<std::string> expected_AAs = partitionStringIntokMers(AAsequence, k);
					std::set<std::string> expected_AAs_set(expected_AAs.begin(), expected_AAs.end());
					assert(haplotypeStore.count(k));
					std::set<std::string> got_AAs;
					for(auto AAmers : haplotypeStore.at(k))
					{
						got_AAs.insert(AAmers.first);
					}

					assert_AA_sets_identical(got_AAs,  expected_AAs_set);
				}
			}

			{
				int coveredBases = 0;
				int currentExonStart = 0;
				std::vector<std::pair<int, int>> exons;
				while(coveredBases < (int)nucleotideSequnce.length())
				{
					int max_exon_length = nucleotideSequnce.length() - currentExonStart;
					int exon_length = randomNumber(max_exon_length);
					int exon_stop_pos = currentExonStart + exon_length - 1;
					assert(exon_stop_pos < (int)nucleotideSequnce.length());
					assert(exon_stop_pos >= ((int)currentExonStart - 1));

					exons.push_back(make_pair(currentExonStart, exon_stop_pos));

					coveredBases = (exon_stop_pos+1);
					currentExonStart = exon_stop_pos + 1;
				}

				std::cerr << exons.back().second << " " << ((int)nucleotideSequnce.length() - 1) << "\n" << std::flush;
				assert(exons.back().second == ((int)nucleotideSequnce.length() - 1));

				// no intermediary sequence
				{
					std::map<std::string, std::string> referenceGenome;
					referenceGenome["chr"] = (strand == '+') ? nucleotideSequnce : seq_reverse_complement(nucleotideSequnce);

					std::map<std::string, std::map<int, variantFromVCF>> variants;

					std::vector<transcriptExon> transcript_exons;
					for(auto e : exons)
					{
						transcriptExon thisExon;
						if(e.first <= e.second)
						{
							thisExon.valid = true;
							thisExon.firstPos = e.first;
							thisExon.lastPos = e.second;
						}
						transcript_exons.push_back(thisExon);
					}

					if(strand == '-')
					{
						for(transcriptExon& e : transcript_exons)
						{
							e.firstPos = referenceGenome["chr"].length() - e.firstPos - 1;
							e.lastPos = referenceGenome["chr"].length() - e.lastPos - 1;
							assert(e.firstPos >= e.lastPos);
							int third = e.firstPos;
							e.firstPos = e.lastPos;
							e.lastPos = third;
						}
					}

					transcript oneTranscript;
					oneTranscript.chromosomeID = "chr";
					oneTranscript.geneName = "testGene";
					oneTranscript.strand = strand;
					oneTranscript.exons = transcript_exons;

					std::vector<transcript> transcripts = {oneTranscript};

					std::map<int, std::map<std::string, std::map<double, std::set<std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>>> haplotypeStore;
					for(auto k : AA_mers)
					{
						haplotypeStore[k].count("");
					}

					enumeratePeptideHaplotypes(referenceGenome, transcripts, variants, true, haplotypeStore);

					for(auto k : AA_mers)
					{
						std::vector<std::string> expected_AAs = partitionStringIntokMers(AAsequence, k);
						std::set<std::string> expected_AAs_set(expected_AAs.begin(), expected_AAs.end());
						assert(haplotypeStore.count(k));
						std::set<std::string> got_AAs;
						for(auto AAmers : haplotypeStore.at(k))
						{
							got_AAs.insert(AAmers.first);
						}

						assert_AA_sets_identical(got_AAs,  expected_AAs_set);
					}
				}
			}
		}
	}
}

void assert_AA_sets_identical(const std::set<std::string>& s1, const std::set<std::string>& s2)
{
	if(s1 != s2)
	{
		std::cerr << "s1 size: " << s1.size() << "; s2: " << s2.size() << "\n";
		std::set<std::string> combined;
		combined.insert(s1.begin(), s1.end());
		combined.insert(s2.begin(), s2.end());
		for(auto e : combined)
		{
			std::cerr << "\t" << e << " " << (s1.count(e) ? 1 : 0) << " " << (s2.count(e) ? 1 : 0) << "\n";
		}
		std::cerr << std::flush;
	}
	assert(s1 == s2);
}

unsigned int randomNumber(unsigned int max_inclusive)
{
	int n = rand() % (max_inclusive+1);
	assert((n >= 0) && (n <= max_inclusive));
	return n;
}

