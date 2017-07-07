# INSTALL #

## Prerequisites ##
You need a compiler with support for C++11 (e.g. g++ 4.8.5);  a modern version of Perl (>= 5.10); Java (tested on 1.8.0_66).

## Installing EpitopeEnumerator ##

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
## Download the data package (7.3G) ##

```
cd ~/EpitopeEnumerator/src/dataPackage
wget https://gembox.cbcb.umd.edu/shared/EpitopeEnumeratorData.tar.gz
tar -xvzf EpitopeEnumeratorData.tar.gz
```

## External dependencies ##

You need locally installed versions of
- NetMHCpan (http://www.cbs.dtu.dk/services/NetMHCpan/ - tested with version 3.0)
- NetMHCIIpan (http://www.cbs.dtu.dk/services/NetMHCIIpan/ - tested with version 3.1)

For each sample that you want to apply EpitopeEnumerator to, you also need HLA types in the HLA\*PRG:LA format -- usually the easiest way to obtain these is by running HLA\*PRG:LA, but you can also encode existing HLA types in the right format (see below).
- HLA\*PRG:LA (https://github.com/AlexanderDilthey/HLA-PRG-LA)

If you want to use `prepareSequencingReads.pl` (take sequencing reads from existing cancer / normal BAMs and apply up-to-date variant calling with GATK), you also need locally installed versions of:
- bwa (https://github.com/lh3/bwa - >= 0.7.12)
- GATK (https://software.broadinstitute.org/gatk/ - tested with version 3-7)
- Picard (https://broadinstitute.github.io/picard/ - >= 2.7.1)

## Modifying paths.ini ##

Finally, you need to tell EpitopeEnumerator where to find the installed external software packages. This is done by modifying the file `~/EpitopeEnumerator/src/paths.ini`. If any of these programs are in your path, you can just specify the command without arguments - they will be resolved via a call to `which`. Otherwise, absolute paths are preferred. The version of `paths.ini` that comes with the source gives examples.

If you only installed NetMHCpan and NetMHCIIpan, you only need to modify the entries `netMHCpan_bin` and `netMHCIIpan_bin`. All other entries can be left unmodified.

If you installed bwa, GATK and Picard as well, you also need to modify `bwa_bin`, `GATK_jar` and `samtools_bin`. If you want to change the path of the JVM used for Picard and GATK, or if you want to specify JVM arguments, you can modify `Java_bin_and_arguments`. 

You can use the keys `HumanGenomeReference`, `Gencode_GFF3` and `knownPositions` to specify the locations of the data package files, but usually these don't need to be modified (unless you extracted the data package into a non-standard path).

