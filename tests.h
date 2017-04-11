/*
 * tests.h
 *
 *  Created on: Apr 9, 2017
 *      Author: diltheyat
 */

#ifndef TESTS_H_
#define TESTS_H_

#include "EpitopeEnumerator.h"
#include <set>
#include <string>

void randomTests();
void randomTests_withVariants_2();
void randomTests_withVariants();
void some_simple_tests();

void print_positions_and_interesting_sets(const std::string& peptide, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>> s);

void test_proper_improper_enumeration();

void assert_AA_sets_identical(const std::set<std::string>& s1, const std::set<std::string>& s2);

unsigned int randomNumber(unsigned int max_inclusive);

std::string generateRandomNucleotideSequence(int length);
double randomDouble();

void compare_proper_improper_peptides(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths);

#endif /* TESTS_H_ */
