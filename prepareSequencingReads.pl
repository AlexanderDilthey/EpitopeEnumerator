#!/usr/bin/perl

# example command: perl prepareSequencingReads.pl --normalBAM /data/projects/phillippy/projects/Vaccination/tumourData/case_007/TCRBOA7-N-WEX.bam --cancerBAM /data/projects/phillippy/projects/Vaccination/tumourData/case_007/TCRBOA7-T-WEX.bam --outputDirectory reproduction
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";
use ini;
$| = 1;

my $paths_href = getIni('all');

my $normalBAM;
my $cancerBAM;
my $outputDirectory;
my $threads = 16;

GetOptions (
	'normalBAM:s' => \$normalBAM,
	'cancerBAM:s' => \$cancerBAM,
	'outputDirectory:s' => \$outputDirectory,
	'threads:s' => \$threads,
);

unless(($normalBAM) and (-e $normalBAM))
{
	die "Input file --normalBAM $normalBAM not found";
}
unless(($cancerBAM) and (-e $cancerBAM))
{
	die "Input file --cancerBAM $cancerBAM not found";
}

die "Please specify two different files for --normalBAM and --cancerBAM" if($normalBAM eq $cancerBAM);


unless($outputDirectory)
{
	die "Please specify --outputDirectory";
}
unless(-d $outputDirectory)
{
	mkdir($outputDirectory) or die "Cannot mkdir $outputDirectory";
}
my $normalOutputDir = $outputDirectory . '/normal';
my $cancerOutputDir = $outputDirectory . '/cancer';
unless(-d $normalOutputDir)
{
	mkdir($normalOutputDir) or die "Cannot mkdir $normalOutputDir";
}
unless(-d $cancerOutputDir)
{
	mkdir($cancerOutputDir) or die "Cannot mkdir $cancerOutputDir";
}
make_sure_ref_indexed($paths_href->{HumanGenomeReference}, $paths_href->{'bwa_bin'});

# generic preparation files

my $finalNormalBAM;
my $finalCancerBAM;

foreach my $cfg ([$normalBAM, $normalOutputDir, \$finalNormalBAM, 'N'], [$cancerBAM, $cancerOutputDir, \$finalCancerBAM, 'T'])
{
	my $inputBAM = $cfg->[0];
	my $thisOutputDirectory = $cfg->[1];
	my $newReadGroup = $cfg->[3];
	die unless(defined $newReadGroup);
	
	my $r1 = $thisOutputDirectory . '/R1.fa';
	my $r2 = $thisOutputDirectory . '/R2.fa';
	my $rBAM = $thisOutputDirectory . '/rawreads.bam';
	my $rSortedBAM = $thisOutputDirectory . '/rawreads.sorted.bam';
	my $rSortedBAMWithRG = $thisOutputDirectory . '/rawreads.sorted.withRG.bam';
	my $recalibrationTable = $thisOutputDirectory . '/recalibrationTable';
	my $rSortedBAMRecalibrated = $thisOutputDirectory . '/rawreads.sorted.withRG.recalibrated.bam';

	# generic preparation

	my $cmd_extract = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{Picard_jar} SamToFastq I=$inputBAM F=$r1 F2=$r2 VALIDATION_STRINGENCY=LENIENT);
	exc($cmd_extract);

	my $cmd_align = qq($paths_href->{bwa_bin} mem -t $threads $paths_href->{HumanGenomeReference} $r1 $r2 | $paths_href->{samtools_bin} view -Sb - > $rBAM );
	exc($cmd_align);

	my $cmd_sort = qq($paths_href->{samtools_bin} sort $rBAM > $rSortedBAM);
	exc($cmd_sort);

	my $cmd_index = qq($paths_href->{samtools_bin} index $rSortedBAM);
	exc($cmd_index);

	my $cmd_RG = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{Picard_jar} AddOrReplaceReadGroups I=$rSortedBAM O=$rSortedBAMWithRG RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$newReadGroup);
	exc($cmd_RG);

	my $cmd_index_2 = qq($paths_href->{samtools_bin} index $rSortedBAMWithRG);
	exc($cmd_index_2);

	my $cmd_recalibrate_1 = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T BaseRecalibrator -R $paths_href->{HumanGenomeReference} -I $rSortedBAMWithRG -knownSites $paths_href->{knownPositions} -o $recalibrationTable);
	exc($cmd_recalibrate_1);

	my $cmd_recalibrate_2 = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T PrintReads -R $paths_href->{HumanGenomeReference} -I $rSortedBAMWithRG -BQSR $recalibrationTable -o $rSortedBAMRecalibrated);
	exc($cmd_recalibrate_2);

	my $cmd_index_3 = qq($paths_href->{samtools_bin} index $rSortedBAMRecalibrated);
	exc($cmd_index_3);
	
	${$cfg->[2]} = $rSortedBAMRecalibrated;
} 

my $finalNormalVCF = $normalOutputDir . '/calls_filtered.pass.vcf';
my $finalCancerVCF = $cancerOutputDir . '/calls_somatic_filtered.vcf';

{
	my $VCFcalls = $normalOutputDir . '/calls.vcf';
	my $SNPcalls = $normalOutputDir . '/calls.SNPs.vcf';
	my $INDELcalls = $normalOutputDir . '/calls.INDELs.vcf';
	my $SNPcallsFiltered = $normalOutputDir . '/calls.SNPs.filtered.vcf';
	my $INDELcallsFiltered = $normalOutputDir . '/calls.INDELs.filtered.vcf';
	my $VCFcallsFiltered = $normalOutputDir . '/calls_filted.vcf';
	my $VCFpass = $finalNormalVCF;
	
	my $cmd_haplotypecaller = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T HaplotypeCaller -R $paths_href->{HumanGenomeReference} -I $finalNormalBAM -o $VCFcalls);	
	exc($cmd_haplotypecaller);

	my $cmd_select1 = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T SelectVariants -R $paths_href->{HumanGenomeReference} -V $VCFcalls -selectType SNP -o $SNPcalls);
	exc($cmd_select1);

	my $cmd_select2 = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T SelectVariants -R $paths_href->{HumanGenomeReference} -V $VCFcalls -selectType INDEL -o $INDELcalls);
	exc($cmd_select2);

	my $cmd_filtration1 = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T VariantFiltration -R $paths_href->{HumanGenomeReference} -V $SNPcalls --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o $SNPcallsFiltered);
	exc($cmd_filtration1);

	my $cmd_filtration2 = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T VariantFiltration -R $paths_href->{HumanGenomeReference} -V $INDELcalls --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter" -o $INDELcallsFiltered);
	exc($cmd_filtration2);

	my $cmd_combine = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T CombineVariants -R $paths_href->{HumanGenomeReference} -V $SNPcallsFiltered -V $INDELcallsFiltered -o $VCFcallsFiltered --genotypemergeoption UNIQUIFY);
	exc($cmd_combine);
	
	my $cmd_PASS = qq(grep '\\(PASS\\)\\|\\(^#\\)' $VCFcallsFiltered > $VCFpass);
	exc($cmd_PASS);
	
	# die join("\n\n", $cmd_haplotypecaller, $cmd_select1, $cmd_select2, $cmd_filtration1, $cmd_filtration2, $cmd_combine, $cmd_PASS);
}

{
	my $VCFcallsSomatic = $cancerOutputDir . '/calls_somatic.vcf';
	my $VCFcallsSomaticFiltered = $finalCancerVCF;

	my $cmd_call = qq($paths_href->{Java_bin_and_arguments} -jar $paths_href->{GATK_jar} -T MuTect2 -R $paths_href->{HumanGenomeReference} -I:tumor $finalCancerBAM -I:normal $finalNormalBAM -o $VCFcallsSomatic);	
	exc($cmd_call);
	
	my $cmd_PASS = qq(grep '\\(PASS\\)\\|\\(^#\\)' $VCFcallsSomatic > $VCFcallsSomaticFiltered);
	exc($cmd_PASS);

	# die join("\n\n", $cmd_call, $cmd_PASS);
}

die "Expected output file missing" unless(-e $finalNormalVCF);
die "Expected output file missing" unless(-e $finalCancerVCF);

print "\n\nProcessing done. Produced VCFs:\n";
print "\t - $finalNormalVCF\n";
print "\t - $finalCancerVCF\n";
print "\n\n";

sub exc
{
	my $cmd = shift;
	print "Now executing:\n\t$cmd\n\n";
	system($cmd) and die "Command $cmd failed";
}

sub make_sure_ref_indexed
{
	my $ref_fa = shift;
	my $bwa_bin = shift;
	unless(-e $ref_fa . '.bwt')
	{
		my $cmd = qq($bwa_bin index $ref_fa);
		system($cmd) and die "bwa index command $cmd failed";
	}
}