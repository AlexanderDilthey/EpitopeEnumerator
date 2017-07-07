use strict;
use File::Which; 
use FindBin;

sub getIni
{
	my $what = shift;
	
	my %ini;
	my $paths_file = $FindBin::Bin . '/paths.ini';
	open(INI, '<', $paths_file) or die "Cannot open $paths_file - please execute me from the directory in which I am located";
	while(<INI>)
	{
		chomp;
		$_ =~ s/[\n\r]//g;
		next unless($_);
		die "Invalid format line $. paths.ini" unless($_ =~ /^(.+?)=(.+)$/);
		$ini{$1} = $2;
	}
	close(INI); 
	
	# Java
	my $java_version_command = $ini{'Java_bin_and_arguments'} . ' -version';
	# system($java_version_command) and "Command $java_version_command failed - please specify correct path to java in $paths_file";
		
	$ini{'bwa_bin'} = existing_or_which(\%ini, 'bwa_bin', $paths_file) if($what eq 'all');
	$ini{'samtools_bin'} = existing_or_which(\%ini, 'samtools_bin', $paths_file) if($what eq 'all');
	$ini{'netMHCpan_bin'} = existing_or_which(\%ini, 'netMHCpan_bin', $paths_file);
	$ini{'netMHCIIpan_bin'} = existing_or_which(\%ini, 'netMHCIIpan_bin', $paths_file);
	
	$ini{'GATK_jar'} = existing_or_binaryCwd(\%ini, 'GATK_jar', $paths_file) if($what eq 'all');
	$ini{'Picard_jar'} = existing_or_binaryCwd(\%ini, 'Picard_jar', $paths_file) if($what eq 'all');
	$ini{'HumanGenomeReference'} = existing_or_binaryCwd(\%ini, 'HumanGenomeReference', $paths_file);
	$ini{'knownPositions'} = existing_or_binaryCwd(\%ini, 'knownPositions', $paths_file) if($what eq 'all');
	$ini{'Gencode_GFF3'} = existing_or_binaryCwd(\%ini, 'Gencode_GFF3', $paths_file);
	
	return \%ini;
}

sub existing_or_binaryCwd
{
	my $ini_href = shift;
	my $key = shift;
	my $paths_file = shift;	
	
	if(-e $ini_href->{$key})
	{
		return $ini_href->{$key};
	}
	elsif(-e $FindBin::Bin . '/' . $ini_href->{$key})
	{
		return $FindBin::Bin . '/' . $ini_href->{$key}
	}
	else
	{
		die "Please specify correct path for key $key in $paths_file";
	}
}

sub existing_or_which
{
	my $ini_href = shift;
	my $key = shift;
	my $paths_file = shift;
	
	if(-x $ini_href->{$key})
	{
		return $ini_href->{$key};
	}
	elsif(-x which($ini_href->{$key}))
	{
		return which($ini_href->{$key});
	}
	else
	{
		die "Please specify correct executable path for key $key in $paths_file";
	}

}	
1;