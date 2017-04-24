#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Long; 
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
	'classII' => [qw/HLA-DQA10101-DQB10501 HLA-DQA10101-DQB10502 HLA-DQA10102-DQB10501 HLA-DQA10102-DQB10502 DRB1_0101/, '/data/projects/phillippy/projects/Vaccination/prediction/netMHCIIpan-3.1/DRB11596.fsa'],
);

my $tmpDir = 'tmp';
mkdir($tmpDir) unless(-d $tmpDir);
die unless(-d $tmpDir);

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

			system($cmd) and die "Command $cmd failed";
		}
	}
}