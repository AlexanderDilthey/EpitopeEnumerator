#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Long; 
use List::Util qw/all/;
use List::MoreUtils qw/mesh/;

$| = 1;

my $netMHCpan_bin = qq(/data/projects/phillippy/projects/Vaccination/prediction/netMHCpan-3.0/netMHCpan);
my $netMHCIIpan_bin = qq(/data/projects/phillippy/projects/Vaccination/prediction/netMHCIIpan-3.1/netMHCIIpan);

my $peptides_file = 'peptides_intermediate.txt';
GetOptions (
	'peptides_file:s' => \$peptides_file, 
);

my %lengths = (
	'classI' => [8, 9, 10, 11],
	'classII' => [15],
);

my %hla = (
	'classI' => [qw/HLA-A03:01 HLA-A11:01 HLA-B27:05 HLA-B39:01 HLA-C02:02 HLA-C07:02/],
	# 'classII' => [qw/HLA-DQA10101-DQB10501 HLA-DQA10101-DQB10502 HLA-DQA10102-DQB10501 HLA-DQA10102-DQB10502 DRB1_0101/, '/data/projects/phillippy/projects/Vaccination/prediction/netMHCIIpan-3.1/DRB11596.fsa'],
	'classII' => [qw/HLA-DQA10101-DQB10501 HLA-DQA10101-DQB10502 HLA-DQA10102-DQB10501 HLA-DQA10102-DQB10502 DRB1_0101/],
);

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
		my $fn = $tmpDir . '/' . $class . '_' . $length . '.fa';
		open(OUT, '>', $fn) or die "Cannot open $fn";
		open(PEPTIDES, '<', $peptides_file) or die "Cannot open file $peptides_file";
		my $headerLine = <PEPTIDES>;
		chomp($headerLine);
		my @header_fields = split(/\t/, $headerLine);
		while(<PEPTIDES>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @line_fields = split(/\t/, $line);
			die unless($#line_fields == $#header_fields);
			my %line = (mesh @header_fields, @line_fields);
			if($line{coreLength} == $length)
			{
				my $additionalPadding = $line{additionalPadding};
				my $peptide = $line{peptide};
				$peptide = substr($peptide, $additionalPadding, length($peptide) - 2 * $additionalPadding);
				die unless(length($peptide) == $length);
				print OUT '>', $peptide, "\n", $peptide, "\n";
			}
		}
		close(PEPTIDES);
		close(OUT);
		print "Produced file $fn -- now apply prediction algorithm!\n";
		
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

my $targetLength = 200;
my $spacerSequence = '!!';
my $spacer_between_individual_targets = length($spacerSequence);

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

