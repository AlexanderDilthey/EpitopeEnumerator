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

my $peptides_file = '../covid_peptides.txt';
my $peptides_file_out = '../covid_peptides.candidates.txt';

my $netMHCpan_bin = '/home/dilthey/netMHCpan/netMHCpan-4.1/netMHCpan';
die "NetMHCpan ($netMHCpan_bin) not executable - please check that path is specified correctly in Perl file" unless(-x $netMHCpan_bin);

my $tmpDir = 'tmp';
mkdir($tmpDir) unless(-d $tmpDir);
die unless(-d $tmpDir);

my @HLAtypes = (['H-2-Db', 'classI'], ['H-2-Kb', 'classI']);

my %peptide_prediction_data;
foreach my $HLAtype_tupel (@HLAtypes)
{
	my $HLAtype = $HLAtype_tupel->[0];
	my $class = $HLAtype_tupel->[1];
	
	my $HLAtype_noExtraCharacters = $HLAtype;
	$HLAtype_noExtraCharacters =~ s/[^\w]//g;
	
	my %peptide_lengths;
	{
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

			die unless($#line_fields == $#header_fields);		
			my %line = (mesh @header_fields, @line_fields);

			die unless($line{coreLength});
			$peptide_lengths{$line{coreLength}}++;
		}
		close(PEPTIDES);	
	}
	
	foreach my $length (sort {$a <=> $b} keys %peptide_lengths)
	{
		my $fn = $tmpDir . '/' . 'COVID_' . $class . '_' . $length . '_' . $HLAtype_noExtraCharacters . '.fa';
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
			my $additionalPadding = $line{additionalPadding};
			
			my $peptide = $line{peptide};
			$peptide = substr($peptide, $additionalPadding, length($peptide) - 2 * $additionalPadding);
			die unless(length($peptide) == $line{coreLength});
			if(length($peptide) == $length)
			{
				print OUT '>', $peptide, "\n", $peptide, "\n";
			}
		}
		close(PEPTIDES);
		close(OUT);
		
		my $fn_output = $fn . '.predictedAffinities';
		my $cmd;
		if($class eq 'classI')
		{
			$cmd = qq($netMHCpan_bin -f $fn -l $length -a $HLAtype -BA > $fn_output);
		} 
		else
		{
			die unless($class eq 'classII');
			die;
			#$cmd = qq($netMHCIIpan_bin -f $fn -length $length -a $HLAtype > $fn_output);
		}

		if(not -e $fn_output)
		{
			system($cmd) and die "Command $cmd failed";
		}		
		
		print "Generated $fn_output\n";
		
		my %peptides_minimum_rank;
		
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
				die "Weird number header files, maximum index $#header_fields, line $. of $fn_output:\n$line" unless((scalar(@header_fields) == 17) || (scalar(@header_fields) == 12));
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
						
						
						$n_peptides_this_protein++;
						$n_peptides_total++;
						
						my $Score_EL = $line{'Score_EL'};
						my $rank_EL = $line{'%Rank_EL'};

						my $Score_BA = $line{'Score_BA'};
						my $rank_BA = $line{'%Rank_BA'};
						
						my $Aff = $line{'Aff(nM)'};
						my $BindingStatus = $line{'BindLevel'};
						
						if($BindingStatus)
						{
							$BindingStatus =~ s/<=//;

						}
												
						unless($rank_EL =~ /^[\d\.]+$/)
						{
							die "Weird rank $rank_EL from line in $fn_output:\n$line" 
						}
						
						$peptide_prediction_data{$peptide}{$HLAtype}{Score_EL} = $Score_EL;
						$peptide_prediction_data{$peptide}{$HLAtype}{rank_EL} = $rank_EL;
						$peptide_prediction_data{$peptide}{$HLAtype}{Score_BA} = $Score_BA;
						$peptide_prediction_data{$peptide}{$HLAtype}{rank_BA} = $rank_BA;
						$peptide_prediction_data{$peptide}{$HLAtype}{Aff} = $Aff;
						$peptide_prediction_data{$peptide}{$HLAtype}{BindingStatus} = $BindingStatus;
					}
				}
			}
		}
		close(OUTPUT);
			
		print "\t\tProcessed $n_peptides_total peptides\n";		
	}
}

open(OUT, '>', $peptides_file_out) or die "Cannot open $peptides_file_out";
open(PEPTIDES, '<', $peptides_file) or die "Cannot open file $peptides_file";
my $headerLine = <PEPTIDES>;
chomp($headerLine);
my @header_fields = split(/\t/, $headerLine);
my $took_peptides = 0;
my @header_fields_out = @header_fields;
my @output_fields_predictedAffinities = qw/Score_EL rank_EL Score_BA rank_BA Aff BindingStatus/;
foreach my $HLAtype_tupel (@HLAtypes)
{
	push(@header_fields_out, map {$HLAtype_tupel->[0] . '_' . $_} @output_fields_predictedAffinities);
}
print OUT join("\t", @header_fields_out), "\n";
while(<PEPTIDES>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @line_fields = split(/\t/, $line);
	
	die unless($#line_fields == $#header_fields);		
	my %line = (mesh @header_fields, @line_fields);
	my $additionalPadding = $line{additionalPadding};
	
	my $peptide = $line{peptide};
	$peptide = substr($peptide, $additionalPadding, length($peptide) - 2 * $additionalPadding);
	die unless(length($peptide) == $line{coreLength});
	
	my $strongBinders_HLA = 0;
	foreach my $HLAtype_tupel (@HLAtypes)
	{
		unless(exists $peptide_prediction_data{$peptide}{$HLAtype_tupel->[0]})
		{
			warn "Missing binding predictions for $peptide // $HLAtype_tupel->[0]";
		}
		
		if($peptide_prediction_data{$peptide}{$HLAtype_tupel->[0]}{BindingStatus} eq 'SB')
		{	
			$strongBinders_HLA++;
		}
	}
	
	if($strongBinders_HLA)
	{
		my @output_fields = @line_fields;		
		foreach my $HLAtype_tupel (@HLAtypes)
		{		
			push(@output_fields, map {$peptide_prediction_data{$peptide}{$HLAtype_tupel->[0]}{$_}} @output_fields_predictedAffinities);
		}
		
		print OUT join("\t", @output_fields), "\n";
	}
}
close(OUT);#

print "\nDone. Output all strong binders to $peptides_file_out\n";



