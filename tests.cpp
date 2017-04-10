/*
 * tests.cpp
 *
 *  Created on: Apr 9, 2017
 *      Author: diltheyat
 */

#include "tests.h"

void randomTests_withVariants()
{
	// recovery of expected AA-mers
	for(unsigned int iteration = 0; iteration < 100; iteration++)
	{
		std::cout << "Testing iteration " << iteration << "\n" << std::flush;

		std::vector<std::string> all_AA_sequences;
		//std::vector<std::string> all_nucleotide_sequences;
		std::string AAsequence = generateRandomAASequence(40);
		all_AA_sequences.push_back(AAsequence);

		std::string nucleotideSequnce = translateAASequence2Codons(AAsequence);

		std::map<int, unsigned char> have_variant_at_position;
		int substituted_codons = 0;
		int variants = 0;
		for(unsigned int variantAAi = 0; variantAAi < 4; variantAAi++)
		{
			std::string AAsequence_variant = AAsequence;
			for(unsigned int pos = 0; pos < AAsequence_variant.size(); pos++)
			{
				if(Utilities::randomDouble() <= (double)1/AAsequence_variant.size())
				{
					std::string newAA;
					do {
						newAA = randomAA();
					} while(newAA == AAsequence.substr(pos, 1));

					std::string oldCodon = nucleotideSequnce.substr(3 * pos, 3);
					assert(translateCodon2AA(oldCodon) == AAsequence.substr(pos, 1));

					std::string newCodon = translateAASequence2Codons(newAA);
					bool canSubstitute = true;
					for(unsigned int pI = 0; pI < 3; pI++)
					{
						unsigned int pI_nucleotideSequence = 3 * pos + pI;
						unsigned char oldChar = oldCodon.at(pI);
						unsigned char newChar = newCodon.at(pI);
						assert(nucleotideSequnce.at(pI_nucleotideSequence) == oldChar);
						if(oldChar != newChar)
						{
							if(have_variant_at_position.count(pI_nucleotideSequence) && (have_variant_at_position.at(pI_nucleotideSequence) != newChar))
							{
								canSubstitute = false;
							}
						}
					}

					if(canSubstitute)
					{
						substituted_codons++;
						for(unsigned int pI = 0; pI < 3; pI++)
						{
							unsigned int pI_nucleotideSequence = 3 * pos + pI;
							unsigned char oldChar = oldCodon.at(pI);
							unsigned char newChar = newCodon.at(pI);
							assert(nucleotideSequnce.at(pI_nucleotideSequence) == oldChar);
							if(oldChar != newChar)
							{
								have_variant_at_position[pI_nucleotideSequence] = newChar;
								variants++;
							}
						}
						AAsequence_variant.at(pos) = newAA.at(0);
					}
				}
			}
			all_AA_sequences.push_back(AAsequence_variant);
		}

		std::cout << "Generated " << all_AA_sequences.size() << " sequence with " << substituted_codons << " substituted codons and " << variants << " variants.\n" << std::flush;

		std::set<int> AA_mers = {5, 7, 12};
		for(char strand : {'+'})
		{
			// all as 1 exon, no variants
			{
				std::map<std::string, std::string> referenceGenome;
				referenceGenome["chr"] = (strand == '+') ? nucleotideSequnce : seq_reverse_complement(nucleotideSequnce);

				std::map<std::string, std::map<int, variantFromVCF>> variants;
				for(auto v : have_variant_at_position)
				{
					variantFromVCF vv;
					vv.chromosomeID = "chr";
					vv.position = v.first;
					vv.referenceString = nucleotideSequnce.substr(vv.position, 1);
					std::string variantAllele = {v.second};
					vv.sampleAlleles = std::vector<std::string>({vv.referenceString, variantAllele});
					variants["chr"][v.first] = vv;
				}

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
					std::set<std::string> expected_AAs;
					for(auto aaS : all_AA_sequences)
					{
						std::vector<std::string> aa = partitionStringIntokMers(aaS, k);
						expected_AAs.insert(aa.begin(), aa.end());
					}
					assert(haplotypeStore.count(k));

					std::set<std::string> got_AAs;
					for(auto AAmers : haplotypeStore.at(k))
					{
						got_AAs.insert(AAmers.first);
					}

					std::cout << "k = " << k << ", expect at least " << expected_AAs.size() << ", have " << got_AAs.size() << "\n" << std::flush;

					for(auto expectedAA : expected_AAs)
					{
						assert(got_AAs.count(expectedAA));
					}
				}
			}

			/*
			{
				for(int intermediarySequence : {0, 1})
				{
					std::map<std::string, std::string> referenceGenome;
					referenceGenome["chr"] = "";

					int coveredBases = 0;
					int currentExonStart = 0;
					std::vector<std::pair<int, int>> exons;
					int offset_added_bases = 0;
					while(coveredBases < (int)nucleotideSequnce.length())
					{
						int max_exon_length = nucleotideSequnce.length() - currentExonStart;
						int exon_length = randomNumber(max_exon_length);
						int exon_stop_pos = currentExonStart + exon_length - 1;
						assert(exon_stop_pos < (int)nucleotideSequnce.length());
						assert(exon_stop_pos >= ((int)currentExonStart - 1));


						exons.push_back(make_pair(offset_added_bases+currentExonStart, offset_added_bases+exon_stop_pos));
						// std::cerr << "Exon " << currentExonStart << " " << exon_stop_pos << "\n" << std::flush;

						if(currentExonStart <= exon_stop_pos)
						{
							referenceGenome.at("chr").append(nucleotideSequnce.substr(currentExonStart, exon_stop_pos - currentExonStart + 1));
						}

						coveredBases = (exon_stop_pos+1);
						currentExonStart = exon_stop_pos + 1;

						if(intermediarySequence)
						{
							int random_add_length = randomNumber(1000);
							if(random_add_length)
							{
								referenceGenome.at("chr").append(generateRandomNucleotideSequence(random_add_length));
								offset_added_bases += random_add_length;
							}
						}
					}

					std::cerr << exons.back().second << " " << ((int)nucleotideSequnce.length() - 1) << "\n" << std::flush;

					if(intermediarySequence == 0)
						assert(exons.back().second == ((int)nucleotideSequnce.length() - 1));

					if(strand == '-')
					{
						referenceGenome["chr"] = seq_reverse_complement(referenceGenome.at("chr"));
					}


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
			*/
		}
	}
}

void randomTests()
{
	// recovery of expected AA-mers
	for(unsigned int iteration = 0; iteration < 1000; iteration++)
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
				for(int intermediarySequence : {0, 1})
				{
					std::map<std::string, std::string> referenceGenome;
					referenceGenome["chr"] = "";

					int coveredBases = 0;
					int currentExonStart = 0;
					std::vector<std::pair<int, int>> exons;
					int offset_added_bases = 0;
					while(coveredBases < (int)nucleotideSequnce.length())
					{
						int max_exon_length = nucleotideSequnce.length() - currentExonStart;
						int exon_length = randomNumber(max_exon_length);
						int exon_stop_pos = currentExonStart + exon_length - 1;
						assert(exon_stop_pos < (int)nucleotideSequnce.length());
						assert(exon_stop_pos >= ((int)currentExonStart - 1));


						exons.push_back(make_pair(offset_added_bases+currentExonStart, offset_added_bases+exon_stop_pos));
						// std::cerr << "Exon " << currentExonStart << " " << exon_stop_pos << "\n" << std::flush;

						if(currentExonStart <= exon_stop_pos)
						{
							referenceGenome.at("chr").append(nucleotideSequnce.substr(currentExonStart, exon_stop_pos - currentExonStart + 1));
						}

						coveredBases = (exon_stop_pos+1);
						currentExonStart = exon_stop_pos + 1;

						if(intermediarySequence)
						{
							int random_add_length = randomNumber(1000);
							if(random_add_length)
							{
								referenceGenome.at("chr").append(generateRandomNucleotideSequence(random_add_length));
								offset_added_bases += random_add_length;
							}
						}
					}

					std::cerr << exons.back().second << " " << ((int)nucleotideSequnce.length() - 1) << "\n" << std::flush;

					if(intermediarySequence == 0)
						assert(exons.back().second == ((int)nucleotideSequnce.length() - 1));

					if(strand == '-')
					{
						referenceGenome["chr"] = seq_reverse_complement(referenceGenome.at("chr"));
					}


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

double randomDouble()
{
	assert(RAND_MAX != 0);
	double f = (double)rand() / RAND_MAX;
	assert(f >= 0);
	assert(f <= 1);
	return f;
}


