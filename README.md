# EpitopeEnumerator #

Computes (= enumerates) a list of peptide motifs (= epitopes) that could be used for creating a therapeutic cancer vaccine.

The intention here is to provide a simple proof-of-concept implementation for demonstration - there are additional data sources that could/should be considered when designing an actual therapeutic vaccine. See the accompanying blog post for further discussion.

Input:
- a GRCh38 VCF file describing the normal patient genome
- a GRCh38 VCF file (from Mutect2) specifying the tumour genome
- reference genome and annotation files

If you don't have up-to-date GRCh38 VCFs for your samples, there is also a utility script `prepareSequencingReads.pl` that will take existing cancer / normal BAMs and apply current GATK 'best practises' (from https://software.broadinstitute.org/gatk/best-practices/, June 2017) to generate VCFs from the original reads. This step includes re-alignment to the B38 reference genome and is computationally intensive (both whole-exome and whole-genome data are supported in principle, but if you're dealing with the data volumes of whole-genome data, you might consider pushing the alignment process onto your local grid -- `prepareSequencingReads.pl`, however, can't do that).

## Installation ##

### Prerequisites ###
You need a compiler with support for C++11 (e.g. g++ 4.8.5);  a modern version of Perl (>= 5.10); Java (tested on 1.8.0_66).

### Installing EpitopeEnumerator

```
mkdir ~/EpitopeEnumerator ~/EpitopeEnumerator/src ~/EpitopeEnumerator/obj ~/EpitopeEnumerator/bin
cd ~/EpitopeEnumerator/src
git clone https://github.com/AlexanderDilthey/EpitopeEnumerator .
make all
```

Check that the compilation produced a functional binary: a `../bin/EpitopeEnumerator --action testBinary` should produce the message `EpitopeEnumerator binary is functional.`.

You also need the Perl module `List::MoreUtils` - installing this should be as easy as
```
perl -MCPAN -e shell
install List::MoreUtils
```

Check that the Perl modules are functional:

```
perl -c prepareSequencingReads.pl
perl -c selectEpitopes.pl
```
### Download the data package (7.3G)

```
cd ~/EpitopeEnumerator/src/dataPackage
wget https://gembox.cbcb.umd.edu/shared/EpitopeEnumeratorData.tar.gz
tar -xvzf EpitopeEnumeratorData.tar.gz .
```

### External dependencies ###

You need locally installed versions of
- NetMHCpan (http://www.cbs.dtu.dk/services/NetMHCpan/ - tested with version 3.0)
- NetMHCIIpan (http://www.cbs.dtu.dk/services/NetMHCIIpan/ - tested with version 3.1)

For each sample that you want to apply EpitopeEnumerator to, you also need HLA types in the HLA\*PRG:LA format -- usually the easiest way to obtain these is by running HLA\*PRG:LA, but you can also encode existing HLA types in the right format (see below).
- HLA\*PRG:LA (https://github.com/AlexanderDilthey/HLA-PRG-LA)

If you want to use `prepareSequencingReads.pl` (take sequencing reads from existing cancer / normal BAMs and apply up-to-date variant calling with GATK), you also need locally installed versions of:
- bwa (https://github.com/lh3/bwa - >= 0.7.12)
- GATK (https://software.broadinstitute.org/gatk/ - tested with version 3-7)
- Picard (https://broadinstitute.github.io/picard/ - >= 2.7.1)

### Modifying paths.ini

Finally, you need to tell EpitopeEnumerator where to find the installed external software packages. This is done by modifying the file `~/EpitopeEnumerator/src/paths.ini`. If any of these programs are in your path, you can just specify the command without arguments - they will be resolved via a call to `which`. Otherwise, absolute paths are preferred. The version of `paths.ini` that comes with the source gives examples.

If you only installed NetMHCpan and NetMHCIIpan, you only need to modify the entries `netMHCpan_bin` and `netMHCIIpan_bin`. All other entries can be left unmodified.

If you installed bwa, GATK and Picard as well, you also need to modify `bwa_bin`, `GATK_jar` and `samtools_bin`. If you want to change the path of the JVM used for Picard and GATK, or if you want to specify JVM arguments, you can modify `Java_bin_and_arguments`. 

Finally, you can also use the keys `HumanGenomeReference`, `Gencode_GFF3` and `knownPositions` to specify the locations of the data package files, but usually these don't need to be modified (unless you extracted the data package into a non-standard path).

## Running EpitopeEnumerator ##

### Generate VCFs ###

If you don't have cancer and normal VCFs for your sample, run the following command:

```
cd ~/EpitopeEnumerator/src/
perl prepareSequencingReads.pl --normalBAM /path/to/NORMAL.bam --cancerBAM /path/to/CANCER.bam --outputDirectory /path/to/output/for/SAMPLE
```

Substitute `/path/to/NORMAL.bam` and `/path/to/CANCER.bam` with the paths to the normal and cancer BAM files, respectively.

Modify `--outputDirectory` as desired -- use one directory per sample.

You can pass a ``--threads`` parameter to control the number of utilized CPUs (default 16; e.g. ``--threads 4``).

When successfully finished, the program should have produced the following two files:

```
/path/to/output/for/SAMPLE/normal/calls_filtered.pass.vcf
/path/to/output/for/SAMPLE/cancer/calls_somatic.vcf.filtered
```

### Get sample HLA types ###

For each sample, you need HLA types in the format of HLA\*PRG:LA.

#### Encoding existing types ####
If you already have HLA types for your sample, you can encode these in the right format.

The format is self-explanatory and an example is given here: https://github.com/AlexanderDilthey/HLA-PRG-LA/blob/master/NA12878_example_output_G.txt.

The only required fields are `Locus`, `Chromosome` and `Allele`.

When encoding existing data, make sure to use the correct G group format for the encoded alleles (example: `A*11:01:01G`). You can use the official IMGT G groups table to translate non-G calls into G group calls: http://hla.alleles.org/wmda/hla_nom_g.txt.

#### Running HLA\*PRG:LA ####

Alternatively, you can run HLA\*PRG:LA on the normal BAM file. Example:

```
cd ~/HLA-PRG-LA/src
perl HLA-PRG-LA.pl --BAM /path/to/output/for/SAMPLE/normal/rawreads.sorted.bam --graph PRG_MHC_GRCh38_withIMGT --sampleID mySample --maxThreads 7
```

... which would (with default settings) produce a file `~/HLA-PRG-LA/working/mySample/hla/R1_bestguess_G.txt` in the right format. Full instructions can be found on the HLA\*PRG:LA repository page: https://github.com/AlexanderDilthey/HLA-PRG-LA.

As the example shows, you can run HLA\*PRG:LA on the BAM files generated during the "Generate VCFs" step.

### Enumerate candidate epitopes ###

From directory `~/EpitopeEnumerator/src/`:

```
../bin/EpitopeEnumerator --action enumerate --prefix mySample --normalVCF /path/to/output/for/SAMPLE/normal/calls_filtered.pass.vcf --tumourVCF /path/to/output/for/SAMPLE/cancer/calls_somatic.vcf.filtered --transcripts dataPackage/gencode.v25.annotation.gff3 --referenceGenome dataPackage/hs38DH.fa --mutectVCFMode 1 
```

Output: a file `mySamplepeptides.txt`, containing a list of peptide epitopes. Description of the columns:
- `peptide`: peptide string (amino acids)
- `coreLength`: core peptide length (see below)
- `additionalPadding`: non-core additional amino acids (eiter side of the core, see below)
- `epitopeLocationIndex`: index over the locations that generated the peptide in the cancer genome
- `chromosome`: chromosome of the printed locstion
- `interesting`: a binary string that specifies the positions within the peptide string that correspond to cancer-unique mutations (there are some border cases in which this string cannot be computed correctly, in which case all positions are set to 0).
- `positions`: for each amino acid listed in the peptide string, start and stop position in the cancer genome on the specified chromosome.
- `epitopeMaxP_allLocations`: over all locations of the peptide string in the cancer genome, maximum probability of occurrence. Unless the program is used to exhaustively enumerate haplotypes, there are only two possible values: `1` for "certain" and `-1` for "maybe".

Parameters:
- `--action enumerate`: required for enumeration mode
- `--prefix mySample`: prefix output files with string `mySample`. Optional.
- `--normalVCF PATH`: path to normal VCF
- `--tumourVCF PATH`: path to cancer VCF
- `--transcripts PATH`: path to transcript annotation file. Is part of data package.
- `--referenceGenome PATH`: path to reference genome. Is part of data package.
- `--mutectVCFMode 1`: Assume that the cancer VCF was generated by mutect. See below for detailed explanation.

### Select epitopes ###

From directory `~/EpitopeEnumerator/src/`:

```
perl selectEpitopes.pl --prefix mySample --HLAtypes ~/HLA-PRG-LA/working/mySample/hla/R1_bestguess_G.txt
```

Parameters:
- `--prefix mySample`: prefix for input (`mySamplepeptides.txt`) and output files.
- `--HLAtypes PATH`: path to sample HLA types in HLA\*PRG:LA format (see above).
- `--useAll 1`: use this parameter if you want to use ALL peptides from the input file (the default setting is to use only peptides with a probability of 1).

Output:
- `myPrefix.peptides`: explicitly selected peptides (amino acid sequences)
- `myPrefix.peptidesAsDNAWith2A`: peptide sequences translated into DNA (random codons to prevent homologous recombination in healthy cells) and linked with 2A sequences
- `myPrefix.encodedPeptides.withPromotorAndTail`: peptides and 2A linkers in DNA, with CMV promotor and some polyA tail sequence
- `myPrefix.completeGenome`: genome of Human Adenovirus 5, with E3 components substituted with the content of `myPrefix.encodedPeptides.withPromotorAndTail`.

#### Computational background ####

This script takes the set of cancer-exclusive peptides; for each such peptide, it uses NetMHCpan / NetMHCIIpan to predict how well the peptide would be presented by the patient's HLA proteins; finally, it tries to build a combined peptide string of a desired length, optimizing for presentability of the contained peptides.

This process of building a combined peptide string S happens in a greedy manner:
- set `S` = ''
- we sort the complete list of peptides (i.e. not distinguishing by peptide length) by by how well they're predicted to bind to any of the patient's HLA molecules (i.e., for each peptide, we take the maximum presentation "probability" across all HLA types).
- we iterate through this sorted list and consider each peptide `peptide` in turn.
- if length(`S` + `linker` + `peptide`) < desired_length, we set `S` = `S` + `linker` + `peptide`. `linker` is a sequence that we use to link multiple peptide epitopes (see below). If the just-added `peptide` contains one of the previously added components of `S` as a sub-string, we remove these previously added components of `S` from `S` before considering the next peptide in the sorted list.

#### Conceptional / biological background ####
- In the default configuration, the script searches for class I peptides of length 8, 9, 10, 11 (amino acids) and for class II peptides of length 15. This can easily configured by editing the `%lengths` hash.
- We use a 2A linker string ('EGRGSLLTCGDVEENPGP' from [Szymczak et al.](https://www.ncbi.nlm.nih.gov/pubmed/15064769)) to connect different epitopes in the same translational unit (this is referred to as 'multicistronic' constructs). Briefly, the primary amino acid sequence will be post-translationally separated at the 2A positions, the cleavage happening between the final 'G' and 'P' amino acids. As an unwanted side effect, this will lead to the formation of additional epitopes.
- We currently use random codon encoding for the translation from amino acids -> DNA. This might result in non-optimal expression, but reduces the likelihood of homologous recombination between the expression construct and the source genes in healthy cells (this is a smaller concern in non-viral delivery systems, see below). Two alternatives would be a) to use native DNA encoding for the cancer-specific epitope sequences and an optimized encoding for the 2A sequences or b) to use a human-specific optimized translation table.
- The script produces a modified [Human Adenovirus 5](https://www.ncbi.nlm.nih.gov/nuccore/AC_000008.1) genome. Specifically, the E3 components of that genome are substituted with the DNA sequence encoding the selected tumour-specific epitopes - the idea being that infection with the virus so-designed would induce a specific immune reaction against the selected epitopes. Of note, the modified virus would still be replication-competent, but the deletion of the E3 components should reduce its fitness when compared to wild-type strains. To enable expression and translation of the epitope-encoding DNA sequence, the DNA sequence is inserted into the genomic context of a CMV promotor and polyA sequence, taken from a modified Adenovirus genome published by [Gambotto et al.](https://www.ncbi.nlm.nih.gov/nucleotide/13194617)
