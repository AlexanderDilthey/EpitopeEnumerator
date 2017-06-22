# EpitopeEnumerator #

Computes (= enumerates) a list of peptide motifs (= epitopes) that could be used for creating a therapeutic cancer vaccine.

The intention here is to provide a simple proof-of-concept implementation for demonstration - there are additional data sources that could/should be considered when designing an actual therapeutic vaccine. See the accompanying blog post for further discussion.

Input:
- a GRCh38 VCF file describing the normal patient genome
- a GRCh38 VCF file (from Mutect2) specifying the tumour genome
- reference genome and annotation files

If you don't have up-to-date GRCh38 VCFs for your samples, there is also a utility script `prepareSequencingReads.pl` that will take existing cancer / normal BAMs and apply current GATK 'best-practises' (from https://software.broadinstitute.org/gatk/best-practices/, June 2017) to generate VCFs from the original reads. This step includes re-alignment to the B38 reference genome and is computationally intensive (both whole-exome and whole-genome data are supported in principle, but if you're dealing with the data volumes of whole-genome data, you might consider transferring the alignment process to your local grid -- `prepareSequencingReads.pl`, however, can't do that).

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
### Download the data package

```
cd ~/EpitopeEnumerator/src/dataPackage
wget URL_TO_DATA_PACKAGE
tar -xvzf dataPackage.tar.gz .
```

### External dependencies ###

You need locally installed versions of
- NetMHCpan (http://www.cbs.dtu.dk/services/NetMHCpan/ - tested with version 3.0)
- NetMHCIIpan (http://www.cbs.dtu.dk/services/NetMHCIIpan/ - tested with version 3.1)

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

