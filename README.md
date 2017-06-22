# EpitopeEnumerator #

Computes (= enumerates) a list of peptide motifs (= epitopes) that could be used for creating a therapeutic cancer vaccine.

The intention here is to provide a simple proof-of-concept implementation for demonstration - there are additional data sources that could/should be considered when designing an actual therapeutic vaccine. See the accompanying blog post for further discussion.

Input:
- a GRCh38 VCF file describing the normal patient genome
- a GRCh38 VCF file (from Mutect2) specifying the tumour genome
- reference genome and annotation files

If you don't have up-to-date GRCh38 VCFs for your samples, there is also a utility script prepareSequencingReads.pl that will take existing cancer / normal BAMs and apply current GATK 'best-practises' to generate VCFs from the original reads. This step includes re-alignment to the B38 reference genome and is computationally intensive (both whole-exome and whole-genome data are supported in principle, but if you're dealing with the data volumes of whole-genome data, you might consider transferring the alignment process to your local grid -- prepareSequencingReads.pl, however, can't do that).

## Installation ##

### Dependencies ###

Clone into directory of your choice. Unpack data package into dataPackage. 
### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
