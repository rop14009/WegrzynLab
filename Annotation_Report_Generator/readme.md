# Transcriptome Annotation Pipeline

## Introduction

The aim of this project is to create a software tool to improve the accuracy and speed of transcriptime (gene space) annotation of any non-model organism.

This software tool accomplishes the following:

* Determines the best possible hit from a set of search results (from 1, 2, or 3 sets of search results)
* Appends Blast2GO and Interpro information to the annotation file
* Calculates statistics on the set of search queries.
  *  N50 Statistic
  *  Median query length
  *  Average query length
  *  Shortest query length
  *  Longest query length
  *  Top 10 hits
  *  Top 10 contaminants found
  *  Number of queries with an informative/uninformative hit
  *  Number of contaminants found
  *  Number of sequences with/without domain identification
  *  Number of Components, Functions, Processes from both Blast2GO and Interpro
  *  Number of transcripts with at least 1 Component, Function, or Process
* Generates an XML version of the results that can be exported to Blast2GO
  * 1 XML file is generated per database input, and 1 XML file is generated from all of the given databases which contains the best hits from all of the given databases.
* Generates a list of all of the sequences from the query file that did not have a hit from any of the databases (saved as text file)
* Generates a list of all of the contaminant sequences from the query file (saved as text file)

## Execution Examples

Execution example for the script, after the configuration file has been filled out.

```
cd /path/to/directory/containing/script/directory_containing_script/
python report_generator.py
```

Execution example for usearch:

```
usearch -ublast query.fasta -db /path_to_udb_database/refseq_protein.udb -threads 12 -evalue 1e-9 -weak_evalue 0.0001 -blast6out results
```

Execution example for vsearch:

```
vsearch query.fasta --db /path_to_udb_database/refseq_protein.udb --threads 12 --id 1e-9 --weak_id 0.0001 --blast6out results
```

Execution example for Interpro:

Note: $input is binded to the configuration files for your interproscan run (refer to: https://code.google.com/p/interproscan/wiki/RC4HowToRun)
```
interproscan.sh -i $input -o $output -f xml -appl pfam, Panther -goterms -iprlookup
```

Execution example for Blast2GO:

Please refer to: https://wiki.hpcc.msu.edu/display/Bioinfo/Using+Blast2GO

For Blast2GO execution examples for the GUI and CLI versions of blast.


## Getting Started

### Dependencies

This script requires Python 2.X as a dependency, in order to check which version you have installed run the command:

```
python -V
```

### External Applications

The external applications that are required in order to function with this script are Blast2GO, Interpro, and usearch (or vsearch)


Usearch is a sequence analysis tool.

USearch Download Link: http://www.drive5.com/usearch/download.html

Vsearch is a sequence analysis tool that is open source, and free.

Vsearch Download Link: https://github.com/torognes/vsearch

Blast2GO is an application that allows for the functional annotation of sequences, and the analysis of annotation data.

Blast2GO Download Link: https://www.blast2go.com/blast2go-pro/download-b2g

Interpro is an application that provides functional analysis of protein sequences. It is used to generated the gene ontology of the transcripts.

Interpro Download Link: https://www.ebi.ac.uk/interpro/download.html;jsessionid=B76DDDA8BCBB1AD8AFA31F3FE0E476B5



### Configuration File

The configuration file is the only file you have to change to run this script, and is of the following form:

```
Path to query FASTA:
Source Databases: (0, 1, 2, 3)
Sequence Search Application:
Query Organism:
Path to formatted database1:
Path to FASTA-version of database1:
Path to search results from database1:
Path to formatted database2 (organism-specific):
Path to FASTA-version of database2:
Path to search results from database2:
Path to formatted database3 (full length):
Path to FASTA-version of database3:
Path to search results from database3:
Minimum Evalue: 1e-5
Full-length coverage requirement: (0 - 1)
Generate XML for Blast2GO: yes/no
Path to Interpro Results:
Path to Blast2GO Results:
Contaminant Detection Options
Insects: yes/no
Fungi: yes/no
Bacteria: yes/no
```

In the configuration file, your option is always placed after the ":". If there is no ":" on that line, then it is not an option.

Also it is important to note that the configuration file is line sensitive, it must always be of the form above, deleting any of the lines in the configuration file will lead to errors
The first line is where you will place the path to your query fasta file. 

#### Query Fasta

This option allows you to specify the multifasta file containing your search queries.
This file will be of the following format:

```
>Contig_PtHGO
ATTAAATCCAAACACTATTACTGTGACATTAACAAGGAGAGCAACGGCCTCTGTCTGTATTTGTGTTTGGCAGACGATATTTATGGATTCGCAGCACGTCAATATTTCATGAAAAAAATATCGCCCTAATATAAGATAATCATAACTTTGTTTTTCCTTTGAACTTTACTTCAACCGAAGATCTCGACTGCTGAATTCTGAATTCTGCAATTTGAATGCGCGGGGTTAAAAATGGAGAAATTGGCGAAGCTAAACATTTCAGAAGAGGACTTGGAATTCGAATTTAGCCATCTGGAATACCAGAGCGGATTTGGCAATAATTTTTCCTCGGAAGCAAAACCTGGTGCACTGCCGGTGGGGCAGAATAATCCCTTGAAGTGTCCTTATGGCCTCTACGCCGAGCAAATCTCAGGGACGGCATTTACAGTCCCCAGGAAGCTCAACCAGAGAAGTTGGCTGTACCGTATCAAGCCATCAGTAACACATGAACCTTTTCATCCTCGGGTTCCATCACATGATTATCTGATAAGTGAGTTCAATCAGTCATCAAGCTCAGCAACACCAACGCAATTGCGTTGGAATCCTTTCAACATTCCAGATACGAAAACAAATTTCATTGATGGACTGTTCACTATATGTGGAGCTGGTAGTTCATTCCTTCGGCATGGCTTTGCAGTGCACATGTATGTCGCAAATGCATCTATGGAAGGCTGTGCATTTGCTAATGCTGATGGGGACTTTCTGATAGTTCCACAACAAGGAAGATTATGGATAACGACAGAGTTAGGAAGATTGCAGGTTTCTCCTGGTGAAATAGTTGTGTTGCAACAAGGTTTCCGGTATTCCATTGATTTACCAGATGGTCCTTCTCGAGGGTACGTTGTAGAGGTTTTCAGTGGCCATTTCCAGCTTCCCGATCTTGGTCCAATCGGTGCAAATGGCCTGGCCTCACCGCCTGACTTTCTAACACCTACAGCATGGTTTGATGACAAAGCCTATCCTGGTTACACGATTGTTCATAAGTTTGGAGGTTCATTGTTTACAGCAAAACAAAATTTCTCCCCCTTCAATGTGGTTGCTTGGCATGGCAATTATGCACCTTACAAGTATGATTTGAAGAAGTTCTGTCCTTTCAATACTGTTTTATTTGATCATGGTGACCCATCAATAAACACAGTTCTGACAGTTCCATCAGAAAAACCAGGAGTAGCTGTGGTTGACTTTGTCATTTTTCCACCTCGCTGGCTGGTAGCAGAGCATACCTTCCGTCCTCCATATTACCACCGTAACTGCATGAGTGAATTTATGGGCCTGATTTATGGCAGCTATGAGGCTAAAGCAGATGGATTTCAACCAGGAGGGGCAAGCTTGCACAGCTGCATGACCCCACATGGGCCTGATACTACCACGTTTGAGAAAACTATAGCTGAAGAGGGCGATGCAAAACCTGCCAAGATAAAAGATACCATGGCTTTCATGTTTGAGTCATCATTGATTCCCAGAATCACACCATGGGTTCTCAAATCTCCCCATTTGGATAACGATTATTACAAGTGCTGGACAGGTCTCAAGTCCCATTTTCACCATGAACATCTTCCTGAAAACGGGGAATTACATAGTTAGTAAATCTGTTTCTGTAAATATGGAGAGATATGTACT

>Contig_PtHPPD
ATTTAGTGCAAGAAACAATGGGTAAGGATGGGAGTCAGAGCTTTAAACTGGTGGGATACAAGAACTTTGTTCGCCACAATCCCATGAGCGATCGCTTCCAAGTGAAGAAATTCCACCACATTGAGTTCTGGTGCGGCGATGCCACCAATACCCACAAGCGCTTCTCATGGGGTTTGGGCATGCAGGTGGTGGCCAAGTCCGATCAGAGCACCGGAAATCAGACATATTGCAGCTTCGCCATACAGTCCAACGATCTCACCTTAGTCTTCACGGCGCCTTATTCTTCGAAAATTGACCAGTCCACCACCAGAATGCCCCACCCAGCTTTCGATTCCGATGAGGCGCGGAGTTTCTTTGTCAAGCATGGGTTAGCAGTTCGAGCCGTGGGAGTGCTTGTCGAAGATGCTCGAGACGCTTTCCAGAATAGCGTGGAAAACGGTGCCGTGGCTGTGCTGGAGCCATGCGATCTTGTCGATGAGGCCACGGGCGGGAAGGTCAGTATGTCGGAGGTTAAGCTTTACGGAGATGTGGTTCTTCGTTTTGTAAGCCAGGACGGTTACAAGGGCCCCTTTTTGCCCAATTACGAGCCCGTTCAGTCGGTTCAATTGTCCTATGGAATCATTCGGGTCGACCACGCTGTGGGCAATGTGGAGAAACTGGAGGAAGCCGTCCAATaTGTGGCCAAGTTTAGCGGGTTTCACCGTTTTGCCGAGTTTACGGCAGAAGACGTCGGGACTGGAGAGAGTGGCCTGAATTCGATGGTGTTGGCGAGCAATAATGAGATGGTTCTATTGCCCATGAACGAGCCTGTGTTTGGGACCAAAAGGAAGTCACAGATACAGACCTATCTGGAGCACAATGAGGGCCCCGGATTGCAGCATTTGGCTCTTATTTGCTCTGATATATTCTCCACGCTTAAGGAGATGCGCATGCGCAGTGAAATTGGAGGGTTTGAATTCATGCCCAAGCCTCCTCCTACTTATTATAAGAATCTGAGGAATAGGGTTGCAGATGTCTTGACTGAGAAGCAGATGGAAGAGTGCGACGAGTTGGGGATTTTGGTCGATAGGGATGAtCAAGGAGTCCTCTTACAGATTTTCACAAAGCCCATTGGCGACAGGCCCACTATTTTTATTGAGATAATTCAACGTGTAGGCTGTGTGCTAACTGATGAGAATGGGAAAGCCTATCAAAAGGGTGGATGTGGAGGATTTGGCAAAGGCAATTTTTCAGAACTCTTTAAATCCATCGAGGAATACGAGAAAACACTTGAAGGAAGAGTCAATTCCTAGGCCCTTCGAAGTTCGTAGTGATCTCCTAGAGTAGATTCATCTTGGAATTAGACTGATAGAGATTTATTATATACGTGAGAGTGGCATAATGATAAAATACATGTACTACGGTTGCTTGTGTTAATGGTTTCTAGATGTCGTTTAGTGAAGATGCCTTGTAAAGAAGCAAGGTATATATGTTATACAAATATACTCTGTTTCTGGTTTGGTTCTAGGCTGCATGCCCTTAAGAGGAATTTAACTGTTATAATTTAAATATATTATGTTCAGATAGGCCCACGGCTTATCATTTTATAATTTAAATATATTCTGTTCAGATATGTCCACGGCTGGTTAAGTTCACTTGTCAAGTAACTCCATTTTTATTCACGCTGGAGCGTAATTTTGGATCTAGGGATATAGAATCAGGATGGATGCTCTCTTTAAAGCTTATAACATATTTATGGTTCGT

```

In this format the first line is prefixed with a ">" and then the sequence name is written. The line following contains the sequence.

#### Source Databases

This option contains the number of databases to be parsed by the script, and can be a number between 0 and 3.

Note, if 0 is selected then no databases will be parsed.

#### Sequence Search Application

This script currently supports 2 different sequence search applications, namely usearch and vsearch.

When running these applications, it is important to use the "-blast6out" command line option for both usearch and vsearch, to ensure that the correct format is generated.

The correct format looks like this:

```
Contig_PtHGO    gi|460412999|ref|XP_004251883.1|        78.7    432     92      0       286     1581    12      443     4.3e-212        744.2
```

#### Paths to the Databases

The current version of the script does not support formatted versions of databases (potentially will be included in future versions)
This option should simply be left blank when filling out the configuration file.


The fasta versions of the database will be in the multi-fasta file format (similar to the format of the query fasta), and the search results will be of the format shown above in the "Sequence Search Application" section.

#### Minimum E-Value and Minimum coverage requirements

The minimum E-Value can be input in either scientific notation (1e-5) or as a decimal (0.00001).

The minimum coverage requirement can only be input as a decimal number between 0 and 1 (such as 0.7)

#### Interpro and Blast2GO paths

For these two options you must enter the path to your interpro and blast2go input files.

The interpro file format is as follows:

```
c10004_g1_i1    3F8F69B8D17CEA7F        185     HMMPfam PF13976 gag_pre-integrs 37      108     1.0E-10 T       06-Nov-2014     IPR025724       GAG-pre-integrase domain
```

The blast2go file format is as follows:

```
c58101_g1_i1|m.1        uninformative   390     2       7.2E-10 0.0%    3       F:GO:0003674; P:GO:0008150; C:GO:0005575
```

#### Contaminant DB selection

This script currently has 3 databases containing species and genus information of Kingdoms known to be contaminants.

Specifically,

* Insects
* Bacteria
* Fungi

In in configuration file each of these databases can be toggled by writing either "y" or "yes" to include them when searching for contaminants, or write "no" or anything that does not include the character "y" or leave it blank to ignore those types of contaminants.


