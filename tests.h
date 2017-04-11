/*
 * tests.h
 *
 *  Created on: Apr 9, 2017
 *      Author: diltheyat
 */

#ifndef TESTS_H_
#define TESTS_H_

#include <set>
#include <string>

#include "EpitopeEnumerator.h"

void test_proper_improper_enumeration();
void randomTests();
void randomTests_withVariants_2();
void randomTests_withVariants();
void some_simple_tests();

void assert_AA_sets_identical(const std::set<std::string>& s1, const std::set<std::string>& s2);
void compare_proper_improper_peptides(const std::map<std::string, std::string> referenceGenome, const std::vector<transcript>& transcripts, const std::map<std::string, std::map<int, variantFromVCF>>& variants, std::set<int> haplotypeLengths);

void print_positions_and_interesting_sets(const std::string& peptide, std::set<std::pair<std::vector<std::pair<int, int>>, std::vector<bool>>> s);

#endif /* TESTS_H_ */
