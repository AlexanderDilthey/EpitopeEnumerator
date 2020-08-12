#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Long; 
use List::Util qw/all/;
use List::MoreUtils qw/mesh/;
use Storable qw/store retrieve dclone/;
use FindBin;
use lib "$FindBin::Bin";
use ini;

$| = 1;

my $paths_href = getIni('');

my $netMHCpan_bin = $paths_href->{'netMHCpan_bin'};
my $netMHCIIpan_bin = $paths_href->{'netMHCIIpan_bin'};

die "NetMHCpan ($netMHCpan_bin) not executable - please check that path is specified correctly in Perl file" unless(-x $netMHCpan_bin);
die "NetMHCIIpan ($netMHCIIpan_bin) not executable - please check that path is specified correctly in Perl file" unless(-x $netMHCIIpan_bin);

unless(-e '_codon_2_aaShort')
{
	system('perl getCodonTables.pl') and die "Could not get codon tables";
	die unless(-e '_codon_2_aaShort');
}
my $codon_2_aaShort_href = retrieve '_codon_2_aaShort';
my $codon_2_aaLong_href = retrieve '_codon_2_aaLong';

my $aaShort_2_codon_href = {};
foreach my $codon (keys %$codon_2_aaShort_href)
{
	my $AAshort = $codon_2_aaShort_href->{$codon};
	push(@{$aaShort_2_codon_href->{$AAshort}}, $codon);
}

my $AD5_href = readFASTA('AD5.fa');
my $prepost_href = readFASTA('prepost.fa');
my $E3_begin_1based = 27858;
my $E3_end_1based = 30839;
my $targetLength_nucleotides = ($E3_end_1based - $E3_begin_1based + 1) - (length($prepost_href->{pre}) + length($prepost_href->{post}));
my $targetLength = $targetLength_nucleotides / 3;

my $prefix = '';
my $useAll = 0;
my $HLAtypes = '';
GetOptions (
	'prefix:s' => \$prefix, 
	'HLAtypes:s' => \$HLAtypes, 
	'useAll:s' => \$useAll, 
);

my $peptides_file = $prefix . 'peptides.txt';

die "Output from C++ component - file name $peptides_file - cannot be found. This file name is constructed as the string given by --prefix PLUS the string 'peptides.txt'" unless(-e $peptides_file);
die "Please specify output from HLA*PRG:LA for normal genome via --HLAtypes" unless(-e $HLAtypes);

my %lengths = (
	'classI' => [8, 9, 10, 11],
	'classII' => [15],
);

my %hla = %{translate_HLAPRG_to_netMHC($HLAtypes)};

print "Sample HLA types, translated:\n";
print "\tclass I : ", join(", ", @{$hla{classI}}), "\n";
print "\tclass II: ", join(", ", @{$hla{classII}}), "\n";

my $tmpDir = 'tmp';
mkdir($tmpDir) unless(-d $tmpDir);
die unless(-d $tmpDir);

my %peptides_minimum_rank;
my %peptides_minimum_rank_sourceAllele;
foreach my $class (qw/classI classII/)
{
	die unless(($class eq 'classI') or ($class eq 'classII'));
	
	foreach my $length (@{$lengths{$class}})
	{
		my $fn = $tmpDir . '/' . $prefix . $class . '_' . $length . '.fa';
		open(OUT, '>', $fn) or die "Cannot open $fn";
		open(PEPTIDES, '<', $peptides_file) or die "Cannot open file $peptides_file";
		my $headerLine = <PEPTIDES>;
		chomp($headerLine);
		my @header_fields = split(/\t/, $headerLine);
		my $took_peptides = 0;
		while(<PEPTIDES>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line);
			die unless($#line_fields == $#header_fields);
			my %line = (mesh @header_fields, @line_fields);
			if(not $useAll)
			{
				 next unless($line{epitopeMaxP_allLocations} == 1);
			}
			if($line{coreLength} == $length)
			{
				my $additionalPadding = $line{additionalPadding};
				my $peptide = $line{peptide};
				$peptide = substr($peptide, $additionalPadding, length($peptide) - 2 * $additionalPadding);
				die unless(length($peptide) == $length);
				print OUT '>', $peptide, "\n", $peptide, "\n";
				$took_peptides++;
			}
		}
		close(PEPTIDES);
		close(OUT);
		
		print "Produced file $fn with $took_peptides peptides -- now apply prediction algorithm!\n";
		
		foreach my $HLAtype (@{$hla{$class}})
		{
			print "\t", $HLAtype, "\n";
			
			my $HLAtype_noExtraCharacters = $HLAtype;
			$HLAtype_noExtraCharacters =~ s/[^\w]//g;
			
			my $fn_output = $fn . '.' . $HLAtype_noExtraCharacters;			

			my $cmd;
			if($class eq 'classI')
			{
				if(-e $HLAtype)
				{
					$cmd = qq($netMHCpan_bin -f $fn -l $length -hlaseq $HLAtype > $fn_output);				
				}
				else
				{
					$cmd = qq($netMHCpan_bin -f $fn -l $length -a $HLAtype > $fn_output);
				}
			} 
			else
			{
				if(-e $HLAtype)
				{
					$cmd = qq($netMHCIIpan_bin -f $fn -length $length -hlaseq $HLAtype > $fn_output);				
				}
				else
				{
					$cmd = qq($netMHCIIpan_bin -f $fn -length $length -a $HLAtype > $fn_output);
				}
			}

			unless(-e $fn_output)
			{
				system($cmd) and die "Command $cmd failed";
			}
			
			open(OUTPUT, '<', $fn_output) or die "Cannot open $fn_output";
			my @header_fields;
			my $n_peptides_this_protein = 0;
			my $n_peptides_total = 0;			
			my $inResults = 0;
			while(<OUTPUT>)
			{
				my $line = $_;
				chomp($line);
				next unless($line);
				if((not $inResults) and ($line =~ /^-+$/))
				{
					my $header_fields = <OUTPUT>;
					chomp($header_fields);
					$header_fields =~ s/^\s+//;
					@header_fields = split(/\s+/, $header_fields);
					$n_peptides_this_protein = 0;
					die "Weird number header files, maximum index $#header_fields, line $. of $fn_output:\n$line" unless((scalar(@header_fields) == 15) || (scalar(@header_fields) == 12));
					my $nextLine = <OUTPUT>;
					chomp($nextLine);
					die "Unexpected non-dash line in $fn_output line $.:\n$nextLine" unless($nextLine =~ /^-+$/);
					$inResults = 1;
				}
				else
				{
					if($inResults)
					{
						if($line =~ /^-+$/)
						{
							warn "More/less (${n_peptides_this_protein}) than expected binding predictions, ending line $. of $fn_output" unless($n_peptides_this_protein == 1);
							
							if($class eq 'classI')
							{
								<OUTPUT>;
								my $proteinLine = <OUTPUT>;
								die unless($proteinLine =~ /^Protein/);
								<OUTPUT>;
								my $spacerLine = <OUTPUT>;
								chomp($spacerLine);
								die unless($spacerLine =~ /^-+$/);
							}
							else
							{
								my $binderLine = <OUTPUT>;
								die unless($binderLine =~ /^Number/);
								my $spacerLine = <OUTPUT>;
								chomp($spacerLine);
								die unless($spacerLine =~ /^-+$/);							
							}
							
							$n_peptides_this_protein = 0;
							$inResults = 0;							
						}
						else
						{
							$line =~ s/^\s+//;
							$line =~ s/(<=) (.B)/$1$2/;
							
							my @line_fields = split(/\s+/, $line);
							die "Weird field count in line $. of $fn_output" unless(($#line_fields == $#header_fields) || ($#line_fields == ($#header_fields-1)));
							if($#line_fields == ($#header_fields-1))
							{
								push(@line_fields, '');
							}
							my %line = (mesh @header_fields, @line_fields);
							
							my $peptide = $line{Peptide};
							my $rank = $line{'%Rank'};
							
							$n_peptides_this_protein++;
							$n_peptides_total++;
							
							unless($rank =~ /^[\d\.]+$/)
							{
								die "Weird rank $rank from line in $fn_output:\n$line" 
							}
							if((not exists $peptides_minimum_rank{$peptide}) or ($peptides_minimum_rank{$peptide} > $rank))
							{
								$peptides_minimum_rank{$peptide} = $rank;
								$peptides_minimum_rank_sourceAllele{$peptide} = $HLAtype;
							}

						}
					}
				}
			}
			close(OUTPUT);
			
			print "\t\tProcessed $n_peptides_total peptides\n";
		}
	}
}


my %peptide_lengths;
foreach my $peptide (keys %peptides_minimum_rank)
{
	$peptide_lengths{length($peptide)}++;
}
my @all_lengths = sort keys %peptide_lengths;

# Thosea asigna virus 2A (T2A of TaV)
# EGRGSLLTCGDVEENPGP
# 


my $spacerSequence = 'EGRGSLLTCGDVEENPGP';
my $spacer_between_individual_targets = length($spacerSequence);
my $spacerSequence_translated = join('', map {randomChoice($aaShort_2_codon_href->{$_})} split(//, $spacerSequence));
my $currentLength_noSpacer = 0;
my %selected_explicitly;


my @peptides_by_rank = sort {$peptides_minimum_rank{$a} <=> $peptides_minimum_rank{$b}} keys %peptides_minimum_rank;
for(my $peptideI = 0; $peptideI <= $#peptides_by_rank; $peptideI++)
{
	my $peptide = $peptides_by_rank[$peptideI];
	my $n_current_selection = scalar(keys %selected_explicitly);
	# next if(exists $selected_explicitly{$peptide});  # for multiple rounds
	
	print "Alternative $peptideI / $#peptides_by_rank, rank $peptides_minimum_rank{$peptide}, sequence $peptide -- currently in: $n_current_selection peptides \n";
	
	my $current_spacer_requirement = ($n_current_selection <= 1) ? 0 : ((scalar(keys %selected_explicitly)-1) * $spacer_between_individual_targets);
	my $current_total_space = $currentLength_noSpacer + $current_spacer_requirement;
	
	print "\tRunning sequence length $current_total_space\n";
	
	my %thisOption_contains_explicitSelection;
	foreach my $k (@all_lengths)
	{
		my @kmers = kMers($peptide, $k);
		foreach my $kmer (@kmers)
		{
			if(exists $selected_explicitly{$kmer})
			{
				$thisOption_contains_explicitSelection{$kmer}++;
			}
		}
	}
	
	print "\tImplicitly contained in this alternative: ", scalar(keys %thisOption_contains_explicitSelection), " existing peptides\n";
	
	my $ifAdd_length_noSpacer = $currentLength_noSpacer;
	foreach my $existingPeptide (keys %thisOption_contains_explicitSelection)
	{
		$ifAdd_length_noSpacer -= length($existingPeptide);
	}
	$ifAdd_length_noSpacer += length($peptide);
	my $ifAdd_n_selection = $n_current_selection - scalar(keys %thisOption_contains_explicitSelection) + 1; die unless($ifAdd_n_selection > 0);
	my $ifAdd_spacer_requirement = ($ifAdd_n_selection-1) * $spacer_between_individual_targets;
	my $ifAdd_total_space = $ifAdd_length_noSpacer + $ifAdd_spacer_requirement;
	
	if($ifAdd_total_space <= $targetLength)
	{
		die if($selected_explicitly{$peptide});
		$selected_explicitly{$peptide} = 1;
		$currentLength_noSpacer += length($peptide);
		print "\tAdding $peptide and kicking out ", scalar(keys %thisOption_contains_explicitSelection), " existing peptides.\n";
		foreach my $existingPeptide (keys %thisOption_contains_explicitSelection)
		{
			die unless($selected_explicitly{$existingPeptide});
			delete $selected_explicitly{$existingPeptide};
			$currentLength_noSpacer -= length($existingPeptide);
		}		
	}	
}

my $targetSequence = join($spacerSequence, keys %selected_explicitly);

print "\n\nDone - target sequence length ", length($targetSequence), "\n";

my @bindingranks_explicitly = map {$peptides_minimum_rank{$_}} keys %selected_explicitly;
print "\tSelected explicitly: ", scalar(@bindingranks_explicitly), " peptides\n";
print "\t\tMin-median-max binding rank: ", join(", ", min_max_median(\@bindingranks_explicitly)), "\n";

my %selected_implicitly;
foreach my $peptide (keys %selected_explicitly)
{
	foreach my $k (@all_lengths)
	{
		my @kmers = kMers($peptide, $k);
		foreach my $kmer (@kmers)
		{
			if(exists $peptides_minimum_rank{$kmer})
			{
				$selected_implicitly{$kmer}++;
			}
		}
	}	
}

my @bindingranks_implicitly = map {$peptides_minimum_rank{$_}} keys %selected_implicitly;

print "\tSelected implicitly: ", scalar(@bindingranks_implicitly), " peptides\n";
print "\t\tMin-median-max binding rank: ", join(", ", min_max_median(\@bindingranks_implicitly)), "\n";

# generate virus genome

die unless(substr($prepost_href->{pre}, length($prepost_href->{pre})-4, 4) eq 'ATGG');
$prepost_href->{pre} = substr($prepost_href->{pre}, 0, length($prepost_href->{pre})-1);
die unless(substr($prepost_href->{pre}, length($prepost_href->{pre})-3, 3) eq 'ATG');

die unless($codon_2_aaLong_href->{substr($prepost_href->{post}, 0, 3)} eq 'Stop');

my %peptides_translated;
my $peptide_starts_withT;
foreach my $selectedPeptide (keys %selected_explicitly)
{
	my $remaingPeptideForTranslation = $selectedPeptide;
	my $runningTranslation;
	if(not $peptide_starts_withT)
	{
		my $firstAA = substr($selectedPeptide, 0, 1);
		die unless(exists $aaShort_2_codon_href->{$firstAA});
		my @possibleTranslations = @{$aaShort_2_codon_href->{$firstAA}};
		my @possibleTranslations_startWithG = grep {substr($_, 0, 1) eq 'G'} @possibleTranslations;
		if(scalar(@possibleTranslations_startWithG))
		{
			$runningTranslation .= randomChoice(\@possibleTranslations_startWithG);
			$remaingPeptideForTranslation = substr($remaingPeptideForTranslation, 1);
			$peptide_starts_withT = $selectedPeptide;
		}
	}
	my @remainingAAs = split(//, $remaingPeptideForTranslation);
	foreach my $AA (@remainingAAs)
	{
		die unless(exists $aaShort_2_codon_href->{$AA});
		my @possibleTranslations = @{$aaShort_2_codon_href->{$AA}};
		$runningTranslation .= randomChoice(\@possibleTranslations);
	}
	
	die unless(length($runningTranslation) == (3 * length($selectedPeptide)));
	$peptides_translated{$selectedPeptide} = $runningTranslation;
}

my @peptides_in_order = $peptide_starts_withT ?
	($peptide_starts_withT, (grep {$_ ne $peptide_starts_withT} keys %selected_explicitly)) :
	(keys %selected_explicitly);
die unless(scalar(@peptides_in_order) == scalar(keys %selected_explicitly));


my $combined_nucleotide_sequence = join($spacerSequence_translated, map {my $t = $peptides_translated{$_}; die unless($t); $t} @peptides_in_order);
if($peptide_starts_withT)
{
	die unless(substr($combined_nucleotide_sequence, 0, 1) eq 'G');
}
else
{
	warn "No proper Kozak consensus sequence";
}
die unless($combined_nucleotide_sequence =~ /^[ACGT]+$/);

my $fn_output_peptides = $prefix . '.peptides';
my $fn_output_peptides_2A = $prefix . '.peptidesAsDNAWith2A';
my $fn_output_withPromotorAndEnd = $prefix . '.encodedPeptides.withPromotorAndTail';
my $fn_output_completeGenome = $prefix . '.completeGenome';

my $string_encoded_peptides = $combined_nucleotide_sequence;
my $string_with_promotor = $prepost_href->{pre} . $combined_nucleotide_sequence . $prepost_href->{post};


die unless(scalar(keys %$AD5_href) == 1);
my $AD5_genome = (values %$AD5_href)[0];

my $string_genome = $AD5_genome;
substr($string_genome, $E3_begin_1based - 1, $E3_end_1based - $E3_begin_1based + 1) = $string_with_promotor;

open(OUTP, '>', $fn_output_peptides) or die "Cannot open $fn_output_peptides";
print OUTP join("\n", keys %selected_explicitly), "\n";
close(OUTP);

open(OUTPEP, '>', $fn_output_peptides_2A) or die "Cannot open $fn_output_peptides_2A";
print OUTPEP $string_encoded_peptides, "\n";
close(OUTPEP);

open(OUTPEPPROM, '>', $fn_output_withPromotorAndEnd) or die "Cannot open $fn_output_withPromotorAndEnd";
print OUTPEPPROM $string_with_promotor, "\n";
close(OUTPEPPROM);

open(OUTGENOME, '>', $fn_output_completeGenome) or die "Cannot open $fn_output_completeGenome";
print OUTGENOME $string_genome, "\n";
close(OUTGENOME);

print "\nGenerated files:\n";
print "\t - $fn_output_peptides\n";
print "\t - $fn_output_withPromotorAndEnd\n";
print "\t - $fn_output_completeGenome\n";

sub randomChoice
{
	my $aref = shift;
	die unless(scalar(@{$aref}) >= 1);
	if(scalar(@{$aref}) == 1)
	{
		return $aref->[0];
	}	
	my $i = int(rand(scalar(@$aref)));
	die unless($i >= 0);
	die unless($i <= $#{$aref});
	return $aref->[$i];
}

sub min_max_median
{
	my $values_aref = shift;
	my @v = @$values_aref;
	die unless(all {defined($_)} @v);
	die unless(scalar(@v));
	@v = sort @v;
	my $middle_index = int($#v/2);
	return ($v[0], $v[$middle_index], $v[$#v]);
}


sub kMers
{
	my $string = shift;
	my $k = shift;
	
	die unless($k > 0);
	
	my $num_kMers = length($string) - $k + 1;
	$num_kMers = 0 if($num_kMers < 0);
	
	my @kMers;
	$#kMers = ($num_kMers - 1);
	
	for(my $i = 0; $i < $num_kMers; $i++)
	{
		my $kMer = substr($string, $i, $k);
		die unless(length($kMer) == $k);
		$kMers[$i] = $kMer;
	}
	
	die unless(scalar(@kMers) == $num_kMers);
	
	return @kMers;
}

sub readFASTA
{
	my $file = shift;	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{
		if(($. % 1000000) == 0)
		{
		# 	print "\r", $.;
		}
		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			$currentSequence = substr($line, 1);
		}
		else
		{
			$R{$currentSequence} .= $line;
		}
	}	
	close(F);
		
	return \%R;
}


sub translate_HLAPRG_to_netMHC
{
	my $fn_HLAPRG = shift;
	die "File $fn_HLAPRG not there" unless(-e $fn_HLAPRG);
	

	my $types_classI_str = `$netMHCpan_bin -listMHC`;
	my @types_classI = split(/\n/, $types_classI_str);
	my %have_types_classI = map {$_ => 1} @types_classI;
	
	my $types_classII_str = `$netMHCIIpan_bin -list`;
	my @types_classII = split(/\n/, $types_classII_str);
	my %have_types_classII = map {$_ => 1} @types_classII;
	
	my @classI_alleles;
	my @classII_alleles;
	my %types = ();
	open(F, '<', $fn_HLAPRG) or die "Cannot open $fn_HLAPRG";
	my $hL = <F>; chomp($hL);
	my @hF = split(/\t/, $hL);
	while(<F>)
	{
		my $l = $_;
		chomp($l);
		next unless($l);
		my @l = split(/\t/, $l);
		my %l = (mesh @hF, @l);
		my $allele = $l{Allele};
		$types{$l{Locus}}{$allele}++;
	}
	close(F);
	
	foreach my $classILocus (qw/A B C/)
	{
		unless(exists $types{$classILocus})
		{
			warn "No HLA types for locus $classILocus?";
			next;
		}	
		
		foreach my $allele (keys %{$types{$classILocus}})
		{
			die unless($allele =~ /\w\*(\d+):(\d+)/);
			my $format_like_netMHC = 'HLA-' . $classILocus . $1 . ':' . $2;
			if(exists $have_types_classI{$format_like_netMHC})
			{
				push(@classI_alleles, $format_like_netMHC);
			}	
			else
			{
				warn "Can't translate allele $allele to NetMCHpan format";
			}
		}
	}
	
	foreach my $classIILocus (qw/DRB1/)
	{
		unless(exists $types{$classIILocus})
		{
			warn "No HLA types for locus $classIILocus?";
			next;
		}	
		
		foreach my $allele (keys %{$types{$classIILocus}})
		{
			die unless($allele =~ /\w\*(\d+):(\d+)/);
			my $format_like_netMHC = $classIILocus . '_' . $1 . '' . $2;
			if(exists $have_types_classII{$format_like_netMHC})
			{
				push(@classII_alleles, $format_like_netMHC);
			}	
			else
			{
				warn "Can't translate allele $allele to NetMCHIIpan format (no $format_like_netMHC?)";
			}
		}
	}	
	
	foreach my $classII (qw/DQ/)
	{
		unless((exists $types{$classII.'A1'}) and (exists $types{$classII.'B1'}))
		{
			warn "No HLA types for locus $classII A/B?";
			next;
		}	
		
		foreach my $alleleA (keys %{$types{$classII.'A1'}})
		{
			die unless($alleleA =~ /\w\*(\d+):(\d+)/);
			my $alleleA_format_like_netMHC = $classII.'A1' . '' . $1 . '' . $2;
			
			foreach my $alleleB (keys %{$types{$classII.'B1'}})
			{
				die unless($alleleB =~ /\w\*(\d+):(\d+)/);
				my $alleleB_format_like_netMHC = $classII.'B1' . '' . $1 . '' . $2;

				my $AB_like_netMHC = 'HLA-' . $alleleA_format_like_netMHC . '-' . $alleleB_format_like_netMHC;
				if(exists $have_types_classII{$AB_like_netMHC})
				{
					push(@classII_alleles, $AB_like_netMHC);
				}	
				else
				{
					warn "Can't translate allele pair $alleleA / $alleleB to NetMCHIIpan format (no $AB_like_netMHC?)";
				}
			}		
		}
	}		

	return {
		'classI' => \@classI_alleles,
		'classII' => \@classII_alleles,
	};		
}