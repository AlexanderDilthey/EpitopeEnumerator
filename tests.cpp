/*
 * tests.cpp
 *
 *  Created on: Apr 9, 2017
 *      Author: diltheyat
 */

#include "tests.h"

void randomTests_withVariants()
{
	// this function generates a range of modified AA haplotypes and
	// generates the SNVs corresponding to the changed codons
	// The condition we impose is only whether we see all expected AAs,
	// and not whether the generated AA set is fully correct.

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
		for(char strand : {'+', '-'})
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
					vv.allAllelesNotInteresting();

					if(strand == '-')
					{
						vv.position = referenceGenome["chr"].length() - vv.position - 1;
						vv.referenceString = seq_reverse_complement(vv.referenceString);
						for(std::string& a : vv.sampleAlleles)
						{
							a = seq_reverse_complement(a);
						}
					}
					variants["chr"][vv.position] = vv;
				}

				checkVariantsConsistentWithReferenceGenome(variants, referenceGenome);

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

				std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
				std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;
				enumeratePeptideHaplotypes(referenceGenome, transcripts, variants, AA_mers, p_per_epitope, p_per_epitope_locations);

				for(auto k : AA_mers)
				{
					std::set<std::string> expected_AAs;
					for(auto aaS : all_AA_sequences)
					{
						std::vector<std::string> aa = partitionStringIntokMers(aaS, k);
						expected_AAs.insert(aa.begin(), aa.end());
					}
					assert(p_per_epitope.count(k));

					std::set<std::string> got_AAs;
					for(auto AAmers : p_per_epitope.at(k))
					{
						got_AAs.insert(AAmers.first);
					}

					std::cout << "k = " << k << ", expect at least " << expected_AAs.size() << ", have " << got_AAs.size() << "\n" << std::flush;

					for(auto expectedAA : expected_AAs)
					{
						assert(got_AAs.count(expectedAA));
					}
				}

				compare_proper_improper_peptides(referenceGenome, transcripts, variants, AA_mers);
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

void test_proper_improper_enumeration()
{
	for(unsigned int iteration = 0; iteration < 1000; iteration++)
	{
		std::cout << "Testing iteration " << iteration << "\n" << std::flush;

		std::string AAsequence = generateRandomAASequence(50);
		std::string nucleotideSequnce = translateAASequence2Codons(AAsequence);

		std::set<int> AA_mers = {5, 7, 12};
		for(char strand : {'+', '-'})
		{
			{
				std::map<std::string, std::string> referenceGenome;
				referenceGenome["chr"] = "";

				int coveredBases = 0;
				int currentExonStart = 0;
				std::vector<std::pair<int, int>> exons;
				int offset_added_bases = 0;
				std::map<int, int> translation_coordinates_beforeExons_toAfterExons;
				while(coveredBases < (int)nucleotideSequnce.length())
				{
					int max_exon_length_in_AA = (nucleotideSequnce.length() - currentExonStart)/3;
					int exon_length = randomNumber(max_exon_length_in_AA) * 3;
					int exon_stop_pos = currentExonStart + exon_length - 1;
					assert(exon_stop_pos < (int)nucleotideSequnce.length());
					assert(exon_stop_pos >= ((int)currentExonStart - 1));


					for(int pI = currentExonStart; pI <= exon_stop_pos; pI++)
					{
						translation_coordinates_beforeExons_toAfterExons[pI] = offset_added_bases+pI;
					}
					exons.push_back(make_pair(offset_added_bases+currentExonStart, offset_added_bases+exon_stop_pos));
					// std::cerr << "Exon " << currentExonStart << " " << exon_stop_pos << "\n" << std::flush;

					if(currentExonStart <= exon_stop_pos)
					{
						referenceGenome.at("chr").append(nucleotideSequnce.substr(currentExonStart, exon_stop_pos - currentExonStart + 1));
					}

					coveredBases = (exon_stop_pos+1);
					currentExonStart = exon_stop_pos + 1;

					int random_add_length = randomNumber(20);
					if(random_add_length)
					{
						referenceGenome.at("chr").append(generateRandomNucleotideSequence(random_add_length));
						offset_added_bases += random_add_length;
					}
				}

				std::string referenceGenome_plus = referenceGenome.at("chr");
				if(strand == '-')
				{
					referenceGenome["chr"] = seq_reverse_complement(referenceGenome.at("chr"));
				}


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
						if(e.valid)
						{
							e.firstPos = referenceGenome["chr"].length() - e.firstPos - 1;
							e.lastPos = referenceGenome["chr"].length() - e.lastPos - 1;
							assert(e.firstPos >= e.lastPos);
							int third = e.firstPos;
							e.firstPos = e.lastPos;
							e.lastPos = third;
						}
					}
				}

				transcript oneTranscript;
				oneTranscript.chromosomeID = "chr";
				oneTranscript.geneName = "testGene";
				oneTranscript.strand = strand;
				oneTranscript.exons = transcript_exons;

				std::vector<transcript> transcripts = {oneTranscript};

				std::map<std::string, std::map<int, variantFromVCF>> variants;

				size_t n_substitutions = 0;
				size_t n_INs = 0;
				size_t n_DELs = 0;

				for(unsigned int pI = 0; pI < referenceGenome_plus.length(); pI++)
				{
					if(Utilities::randomDouble() < 0.02)
					{
						variantFromVCF vv;

						vv.chromosomeID = "chr";
						vv.position = pI;
						vv.referenceString = referenceGenome_plus.substr(vv.position, 1);

						std::string variantAllele;

						if(Utilities::randomDouble() <= 0.5)
						{
							if(Utilities::randomDouble() <= 0.5)
							{
								variantAllele = "-";
								n_DELs++;
							}
							else
							{
								variantAllele = vv.referenceString;
								variantAllele.push_back(randomNucleotide());
								vv.referenceString.push_back('-');
								n_INs++;
							}
						}
						else
						{
							char newAllele;
							do {
								newAllele = randomNucleotide();
							} while(newAllele == vv.referenceString.at(0));
							variantAllele.push_back(newAllele);
							n_substitutions++;
						}

						assert(vv.referenceString.length() == variantAllele.length());
						assert(variantAllele != vv.referenceString);
						if(Utilities::randomDouble() < 0.2)
						{
							vv.sampleAlleles = std::vector<std::string>({variantAllele, variantAllele});
							vv.sampleAlleles_interesting = {true, true};
						}
						else
						{
							vv.sampleAlleles = std::vector<std::string>({vv.referenceString, variantAllele});
							vv.sampleAlleles_interesting = {false, true};
						}

						if(strand == '-')
						{
							vv.position = referenceGenome_plus.length() - vv.position - countCharacters_noGaps(vv.referenceString);
							vv.referenceString = seq_reverse_complement(vv.referenceString);

							for(std::string& a : vv.sampleAlleles)
							{
								a = seq_reverse_complement(a);
							}
						}


						variants["chr"][vv.position] = vv;

					}
				}
				std::cout << "Strand " << strand << ", nucleotide sequence length " << nucleotideSequnce.length() << "; generated " << n_substitutions << " substitutions " << n_INs << " insertions; " << n_DELs << " deletions.\n" << std::flush;

				checkVariantsConsistentWithReferenceGenome(variants, referenceGenome);

				compare_proper_improper_peptides(referenceGenome, transcripts, variants, AA_mers);
			}
		}
	}
}

void randomTests_withVariants_2()
{
	// similar to randomTests_withVariants(), just that we're not starting
	// with AA variants and fully check the generated AA set for correctness

	for(unsigned int iteration = 0; iteration < 100; iteration++)
	{
		std::cout << "Testing iteration " << iteration << "\n" << std::flush;

		std::string AAsequence = generateRandomAASequence(20);
		std::string nucleotideSequnce = translateAASequence2Codons(AAsequence);

		for(int doINDELs : {0, 1})
		{
			std::map<int, std::string> AAvariants;
			std::map<int, std::string> codonVariants;
			for(unsigned int i = 1; i < AAsequence.length(); i++)
			{
				if(doINDELs)
				{
					if(randomDouble() <= 0.1)
					{
						AAvariants[i] = "!";
						codonVariants[3 * i] = translateAASequence2Codons(AAvariants[i]);
					}
					else
					{
						if(randomDouble() <= 0.5)
						{
							std::string newAA = randomAA();
							std::string oldAA = AAsequence.substr(i, 1);
							AAvariants[i] = oldAA+newAA;

							std::string oldCodon = nucleotideSequnce.substr(3 * i, 3);
							assert(translateCodon2AA(oldCodon) == oldAA);
							std::string newCodon = translateAASequence2Codons(newAA);

							codonVariants[3 * i] = oldCodon + newCodon;
						}
						else
						{
							AAvariants[i] = "";
							codonVariants[3 * i] = "";
						}
					}
				}
				else
				{
					std::string newAA;
					do {
						newAA = randomAA();
					} while(newAA == AAsequence.substr(i, 1));
					AAvariants[i] = newAA;

					std::string oldCodon = nucleotideSequnce.substr(3 * i, 3);
					assert(translateCodon2AA(oldCodon) == AAsequence.substr(i, 1));
					std::string newCodon = translateAASequence2Codons(newAA);
					assert(oldCodon != newCodon);

					codonVariants[3 * i] = newCodon;
				}

				i += randomNumber(4);
			}


			std::vector<std::string> AA_haplotypes = {""};
			for(unsigned int pI = 0; pI < AAsequence.length(); pI++)
			{
				if(AAvariants.count(pI) == 0)
				{
					for(string& h : AA_haplotypes)
					{
						h.push_back(AAsequence.at(pI));
					}
				}
				else
				{
					assert(AAvariants.at(pI) != AAsequence.substr(pI, 1));
					size_t n_existing_haplotypes = AA_haplotypes.size();
					for(unsigned int hI = 0; hI < n_existing_haplotypes; hI++)
					{
						AA_haplotypes.push_back(AA_haplotypes.at(hI));
						AA_haplotypes.at(hI).append(AAsequence.substr(pI, 1));
						AA_haplotypes.back().append(AAvariants.at(pI));
					}
				}
			}

			assert(AA_haplotypes.at(0) == AAsequence);

			std::cout << "Generated " << codonVariants.size() << " substituted codons and " << AA_haplotypes.size() << " AA haplotypes.\n" << std::flush;

			std::set<int> AA_mers = {5, 7, 12};
			for(char strand : {'+', '-'})
			{
				// all as 1 exon, no variants
				{
					for(int intermediarySequence : {0, 1})
					{
						std::map<std::string, std::string> referenceGenome;
						referenceGenome["chr"] = "";

						int coveredBases = 0;
						int currentExonStart = 0;
						std::vector<std::pair<int, int>> exons;
						int offset_added_bases = 0;
						std::map<int, int> translation_coordinates_beforeExons_toAfterExons;
						while(coveredBases < (int)nucleotideSequnce.length())
						{
							int max_exon_length_in_AA = (nucleotideSequnce.length() - currentExonStart)/3;
							int exon_length = randomNumber(max_exon_length_in_AA) * 3;
							int exon_stop_pos = currentExonStart + exon_length - 1;
							assert(exon_stop_pos < (int)nucleotideSequnce.length());
							assert(exon_stop_pos >= ((int)currentExonStart - 1));


							for(int pI = currentExonStart; pI <= exon_stop_pos; pI++)
							{
								translation_coordinates_beforeExons_toAfterExons[pI] = offset_added_bases+pI;
							}
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

						if(intermediarySequence == 0)
							assert(exons.back().second == ((int)nucleotideSequnce.length() - 1));

						std::string referenceGenome_plus = referenceGenome.at("chr");
						if(strand == '-')
						{
							referenceGenome["chr"] = seq_reverse_complement(referenceGenome.at("chr"));
						}


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

						std::map<std::string, std::map<int, variantFromVCF>> variants;
						for(auto v : codonVariants)
						{
							variantFromVCF vv;
							vv.chromosomeID = "chr";
							vv.position = translation_coordinates_beforeExons_toAfterExons.at(v.first);

							vv.referenceString = nucleotideSequnce.substr(v.first, 3);
							assert(vv.referenceString == referenceGenome_plus.substr(vv.position, 3));

							std::string variantAllele = {v.second};
							if(variantAllele.length() == 0)
							{
								variantAllele = "---";
							}
							else if(variantAllele.length() == 6)
							{
								vv.referenceString.append("---");
							}
							assert(vv.referenceString.length() == variantAllele.length());

							vv.sampleAlleles = std::vector<std::string>({vv.referenceString, variantAllele});
							vv.allAllelesNotInteresting();

							if(strand == '-')
							{
								//std::cerr << "Variant reference: " << vv.referenceString << "\n";
								//std::cerr << "vv.position: " << vv.position << "\n";

								vv.position = referenceGenome["chr"].length() - vv.position - countCharacters_noGaps(vv.referenceString);
								vv.referenceString = seq_reverse_complement(vv.referenceString);

								//std::cerr << "- Variant reference: " << vv.referenceString << "\n";
								//std::cerr << "- vv.position: " << vv.position << "\n";
								//std::cerr << std::flush;

								for(std::string& a : vv.sampleAlleles)
								{
									a = seq_reverse_complement(a);
								}
							}
							variants["chr"][vv.position] = vv;
						}

						checkVariantsConsistentWithReferenceGenome(variants, referenceGenome);

						std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
						std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;
						enumeratePeptideHaplotypes(referenceGenome, transcripts, variants, AA_mers, p_per_epitope, p_per_epitope_locations);

						for(auto k : AA_mers)
						{
							std::set<std::string> expected_AAs;
							for(auto aaS : AA_haplotypes)
							{
								std::string aaS_noGaps = removeGaps(aaS);
								size_t position_stop = aaS_noGaps.find("!");
								if(position_stop != std::string::npos)
								{
									//std::cerr << "Trim " << position_stop << " / " << aaS_noGaps.size() << "\n";
									aaS_noGaps = aaS_noGaps.substr(0, position_stop);
								}
								std::vector<std::string> aa = partitionStringIntokMers(aaS_noGaps, k);
								expected_AAs.insert(aa.begin(), aa.end());
							}
							assert(p_per_epitope.count(k));

							std::set<std::string> got_AAs;
							for(auto AAmers : p_per_epitope.at(k))
							{
								got_AAs.insert(AAmers.first);
							}

							assert_AA_sets_identical(got_AAs,  expected_AAs);
						}

						compare_proper_improper_peptides(referenceGenome, transcripts, variants, AA_mers);
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

				std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
				std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;
				enumeratePeptideHaplotypes(referenceGenome, transcripts, variants, AA_mers, p_per_epitope, p_per_epitope_locations);

				for(auto k : AA_mers)
				{
					std::vector<std::string> expected_AAs = partitionStringIntokMers(AAsequence, k);
					std::set<std::string> expected_AAs_set(expected_AAs.begin(), expected_AAs.end());
					assert(p_per_epitope.count(k));
					std::set<std::string> got_AAs;
					for(auto AAmers : p_per_epitope.at(k))
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

					std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
					std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;
					enumeratePeptideHaplotypes(referenceGenome, transcripts, variants, AA_mers, p_per_epitope, p_per_epitope_locations);
					compare_proper_improper_peptides(referenceGenome, transcripts, variants, AA_mers);

					for(auto k : AA_mers)
					{
						std::vector<std::string> expected_AAs = partitionStringIntokMers(AAsequence, k);
						std::set<std::string> expected_AAs_set(expected_AAs.begin(), expected_AAs.end());
						assert(p_per_epitope.count(k));
						std::set<std::string> got_AAs;
						for(auto AAmers : p_per_epitope.at(k))
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

void some_simple_tests()
{
	{
		std::cout << "Test 1\n" << std::flush;

		std::map<std::string, std::string> referenceGenome;
		referenceGenome["chr"] = "ACGGCAGCAGCAGCAGCAGCAAAA"; // last pos: 23 AAAAAAK
		//                        0         1         2
		//                                  T  T
		referenceGenome["chr"].append("NNNNNNNNNN"); // last pos: 33
		referenceGenome["chr"].append("TCGTCGTCGTCGTCGTCGTCGTCGTCGTAA"); // last pos: 45

		std::map<std::string, std::map<int, variantFromVCF>> variants;

		variantFromVCF oneVariant;
		oneVariant.chromosomeID = "chr";
		oneVariant.position = 10;
		oneVariant.referenceString = "C";
		oneVariant.sampleAlleles = {"C", "T"};
		oneVariant.sampleAlleles_interesting = {false, true};

		variantFromVCF oneVariant_2;
		oneVariant_2.chromosomeID = "chr";
		oneVariant_2.position = 13;
		oneVariant_2.referenceString = "C";
		oneVariant_2.sampleAlleles = {"C", "T"};
		oneVariant_2.sampleAlleles_interesting = {false, true};

		variants["chr"][10] = oneVariant;
		//variants["chr"][13] = oneVariant_2;

		// std::map<int, variantFromVCF>();

		transcriptExon oneExon;
		oneExon.valid = true;
		oneExon.firstPos = 3;
		oneExon.lastPos = oneExon.firstPos + 7 * 3 - 1; // 23
		transcriptExon twoExon;
		twoExon.valid = true;
		twoExon.firstPos = 34;
		twoExon.lastPos = 63; // 23

		transcript oneTranscript;
		oneTranscript.chromosomeID = "chr";
		oneTranscript.geneName = "testGene";
		oneTranscript.strand = '+';
		oneTranscript.exons = {oneExon, twoExon};

		std::vector<transcript> transcripts_plus = {oneTranscript};

		checkVariantsConsistentWithReferenceGenome(variants, referenceGenome);

		std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;

		enumeratePeptideHaplotypes_improperFrequencies(referenceGenome, transcripts_plus, variants, {6}, p_per_epitope);

		for(auto haplotypesOneLength : p_per_epitope)
		{
			std::cout << "Length " << haplotypesOneLength.first << ": " << haplotypesOneLength.second.size() << "\n" << std::flush;
			for(auto fragment : haplotypesOneLength.second)
			{
				const std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>& probabilites_and_positions = fragment.second;
				double maxP = probabilites_and_positions.first;
				std::cout << "\t" << fragment.first << " " << maxP << "\n" << std::flush;
				for(auto positions_and_interesting : probabilites_and_positions.second)
				{
					std::vector<std::pair<int, int>> positions = positions_and_interesting.first;
					std::vector<bool> interesting = positions_and_interesting.second;
					std::cout << "\t";
					for(auto i : interesting)
					{
						std::cout << (int)i;
					}
					std::cout << "\n" << std::flush;
				}
			}
		}

		compare_proper_improper_peptides(referenceGenome, transcripts_plus, variants, {6});

		//std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;
		//enumeratePeptideHaplotypes(referenceGenome, transcripts_plus, variants, {6}, p_per_epitope, p_per_epitope_locations);
	}

	if(1 == 1)
	{
		std::cout << "Test 2\n" << std::flush;

		std::map<std::string, std::string> referenceGenome;
		referenceGenome["chr"] = "ACGGCAGCAGCAGCAGCAGCATAATTT";
		//                        0         1         2
		std::map<std::string, std::map<int, variantFromVCF>> variants;

		variantFromVCF oneVariant;
		oneVariant.chromosomeID = "chr";
		oneVariant.position = 10;
		oneVariant.referenceString = "C";
		oneVariant.sampleAlleles = {"C", "T"};
		oneVariant.sampleAlleles_interesting = {false, true};

		variantFromVCF oneVariant_2;
		oneVariant_2.chromosomeID = "chr";
		oneVariant_2.position = 13;
		oneVariant_2.referenceString = "C";
		oneVariant_2.sampleAlleles = {"C", "T"};
		oneVariant_2.sampleAlleles_interesting = {false, true};

		variants["chr"][10] = oneVariant;
		variants["chr"][13] = oneVariant_2;

		// std::map<int, variantFromVCF>();

		transcriptExon oneExon;
		oneExon.valid = true;
		oneExon.firstPos = 3;
		oneExon.lastPos = oneExon.firstPos + 7 * 3 - 1;
		transcript oneTranscript;
		oneTranscript.chromosomeID = "chr";
		oneTranscript.geneName = "testGene";
		oneTranscript.strand = '+';
		oneTranscript.exons = {oneExon};

		std::vector<transcript> transcripts_plus = {oneTranscript};
		std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;



		enumeratePeptideHaplotypes_improperFrequencies(referenceGenome, transcripts_plus, variants, {6}, p_per_epitope);

		for(auto haplotypesOneLength : p_per_epitope)
		{
			std::cout << "Length " << haplotypesOneLength.first << ": " << haplotypesOneLength.second.size() << "\n" << std::flush;
			for(auto fragment : haplotypesOneLength.second)
			{
				const std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>& probabilites_and_positions = fragment.second;
				double maxP = probabilites_and_positions.first;
				std::cout << "\t" << fragment.first << " " << maxP << "\n" << std::flush;
				for(auto positions_and_interesting : probabilites_and_positions.second)
				{
					std::vector<std::pair<int, int>> positions = positions_and_interesting.first;
					std::vector<bool> interesting = positions_and_interesting.second;
					std::cout << "\t";
					for(auto i : interesting)
					{
						std::cout << (int)i;
					}
					std::cout << "\n" << std::flush;
				}
			}
		}

		compare_proper_improper_peptides(referenceGenome, transcripts_plus, variants, {6});

		//std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;
		//enumeratePeptideHaplotypes(referenceGenome, transcripts_plus, variants, {6}, p_per_epitope, p_per_epitope_locations);

	}
}

void compare_proper_improper_peptides(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths)
{
	std::map<int, std::map<std::string, double>> peptides_proper;
	std::map<int, std::map<std::string, double>> peptides_improper;

	std::map<int, std::map<std::string, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>> peptides_proper_positions;
	std::map<int, std::map<std::string, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>> peptides_improper_positions;
	{
		std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
		enumeratePeptideHaplotypes_improperFrequencies(referenceGenome, transcripts, variants, haplotypeLengths, p_per_epitope);

		for(auto haplotypesOneLength : p_per_epitope)
		{
			for(auto fragment : haplotypesOneLength.second)
			{
				peptides_improper[haplotypesOneLength.first][fragment.first] = fragment.second.first;
				peptides_improper_positions[haplotypesOneLength.first][fragment.first] = fragment.second.second;
			}
		}
	}

	{
		std::map<int, std::map<std::string, std::pair<double, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>>>>> p_per_epitope;
		std::map<int, std::map<std::string, std::map<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>, double>>> p_per_epitope_locations;

		enumeratePeptideHaplotypes(referenceGenome, transcripts, variants, haplotypeLengths, p_per_epitope, p_per_epitope_locations);

		for(auto haplotypesOneLength : p_per_epitope)
		{
			for(auto fragment : haplotypesOneLength.second)
			{
				peptides_proper[haplotypesOneLength.first][fragment.first] = fragment.second.first;
				peptides_proper_positions[haplotypesOneLength.first][fragment.first] = fragment.second.second;
			}
		}
	}

	assert(peptides_proper.size() == peptides_improper.size());
	for(auto l_and_sets : peptides_proper)
	{
		std::set<std::string> set_proper;
		std::set<std::string> set_improper;
		for(auto p : peptides_proper.at(l_and_sets.first))
		{
			set_proper.insert(p.first);
		}
		for(auto p : peptides_improper.at(l_and_sets.first))
		{
			set_improper.insert(p.first);
		}

		assert_AA_sets_identical(set_proper, set_improper);
		//std::cerr << "set_proper.size() for k = " << l_and_sets.first << ": " << set_proper.size() << "\n" << std::flush;
		for(std::string peptide : set_improper)
		{
			double p_proper = peptides_proper.at(l_and_sets.first).at(peptide);
			double p_improper = peptides_improper.at(l_and_sets.first).at(peptide);
			if(p_improper == 1)
			{
				assert(abs(1 - p_proper) < 1e-4);
				//std::cerr << "p_improper: " << p_improper << "\n" << std::flush;
			}

			if(!(peptides_proper_positions.at(l_and_sets.first).at(peptide) == peptides_improper_positions.at(l_and_sets.first).at(peptide)))
			{
				std::cout << "Print for " << peptide << ", proper:\n";
				print_positions_and_interesting_sets(peptide, peptides_proper_positions.at(l_and_sets.first).at(peptide));
				std::cout << "Print for " << peptide << ", improper:\n";
				print_positions_and_interesting_sets(peptide, peptides_improper_positions.at(l_and_sets.first).at(peptide));
			}
			assert(peptides_proper_positions.at(l_and_sets.first).at(peptide) == peptides_improper_positions.at(l_and_sets.first).at(peptide));
		}
	}


	//enumeratePeptideHaplotypes(referenceGenome, transcripts_plus, variants, {6}, p_per_epitope, p_per_epitope_locations);
}

void print_positions_and_interesting_sets(const std::string& peptide, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>> s)
{
	std::cout << "\tPositions and interesting for " << peptide << "\n";
	for(auto p_and_i : s)
	{
		std::cout << "\t\t";
		for(auto p : p_and_i.first)
		{
			std::cout << "[" << p.first << "-" << p.second << "] ";
		}
		std::cout << "\n";

		std::cout << "\t\t";
		for(auto i : p_and_i.second)
		{
			std::cout << i;
		}
		std::cout << "\n\n";
	}
	std::cout << std::flush;
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


