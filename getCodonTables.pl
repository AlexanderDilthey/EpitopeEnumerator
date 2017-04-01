use strict;

my %codon_2_aaShort;
my %codon_2_aaLong;
while(<main::DATA>)
{
	my $l = $_;
	chomp($l);
	next unless($l);
	my @fields = split(/\t/, $l);
	die unless(scalar(@fields) == 3);
	my $aaLong = $fields[0];
	my $aaShort = $fields[1];
	$fields[2] =~ s/\s//g;
	my @codons = split(/,/, $fields[2]);
	foreach my $codon (@codons)
	{
		die "Duplicate codon $codon" if($codon_2_aaLong{$codon});
		$codon_2_aaLong{$codon} = $aaLong;
		$codon_2_aaShort{$codon} = $aaShort;
	}
}	

print "std::map<std::string, std::string> codon2AA;\n";
foreach my $codon (keys %codon_2_aaShort)
{
	print qq(codon2AA["$codon"] = "$codon_2_aaShort{$codon}";), "\n";
}	
__DATA__
Ala	A	GCT, GCC, GCA, GCG
Arg	R	CGT, CGC, CGA, CGG, AGA, AGG
Asn	N	AAT, AAC
Asp	D	GAT, GAC
Cys	C	TGT, TGC
Gln	Q	CAA, CAG
Glu	E	GAA, GAG
Gly	G	GGT, GGC, GGA, GGG
His	H	CAT, CAC
Ile	I	ATT, ATC, ATA
Leu	L	TTA, TTG, CTT, CTC, CTA, CTG
Lys	K	AAA, AAG
Met	M	ATG
Phe	F	TTT, TTC
Pro	P	CCT, CCC, CCA, CCG
Ser	S	TCT, TCC, TCA, TCG, AGT, AGC
Thr	T	ACT, ACC, ACA, ACG
Trp	W	TGG
Tyr	Y	TAT, TAC
Val	V	GTT, GTC, GTA, GTG
Stop	!	TAA, TGA, TAG
