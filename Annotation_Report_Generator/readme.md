# Transcriptome Annotation Pipeline

## Introduction

The aim of this project is to create a software tool to improve the accuracy and speed of transcriptome (gene space) annotation of any non-model organism.

This software tool accomplishes the following:

* Determines the best possible hit from a set of search results (from 1, 2, or 3 sets of search results)
* Appends Blast2GO and InterProScan information to the annotation file
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
  *  Number of Components, Functions, Processes from both Blast2GO and InterProScan
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

Execution example for InterProScan:

Note: $input is binded to the configuration files for your interproscan run (refer to: https://code.google.com/p/interproscan/wiki/RC4HowToRun)
```
interproscan.sh -i $input -o $output -f xml -appl pfam, Panther -goterms -iprlookup
```

Execution example for Blast2GO:

Please refer to: https://wiki.hpcc.msu.edu/display/Bioinfo/Using+Blast2GO

For Blast2GO execution examples for the GUI and CLI versions of blast.

### Executing combine_annotations.py independently

This script takes the Blast2GO, InterProScan, and annotation input and merges them together into one file

@PARAMETERS

--input, --blast2go, --interpro, --output


Example usages of the script that will result in correct outputs:


python combine_annotations.py --input <annotation filename> --blast2go <blast2go filename> --interpro <interpro filename> --output <output_file_name.tsv>


Note this example will generate an output with the default name combined_annotation.tsv in the directory that the script is being executed in

python combine_annotations.py --input <annotation filename> --blast2go <blast2go filename> --interpro <interpro filename>

And this example will specify a custom output name

python combine_annotations.py --input <annotation filename> --interpro <interpro filename> --output <output_file_name.tsv>


For performing merges with 2 interpro files:

python combine_annotations.py --input <annotation filename> --interpro <interpro filename> <interpro filename two> --output <output_file_name.tsv>


When the script finishes running, the output file will be availible in the output directory.



## Getting Started

### Overview of how the script functions

![Diagram of script execution](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/Annotation_Report_Generator/diagram.png)

The diagram above shows the order in which the seperate operations must be performed in order to generate a complete annotation report.

Steps:

1) You must first run your preferred search application (usearch or vsearch) to obtain a correctly formatted set of search results from your database

2) You must correctly fill out the configuration file to run the report_generator.py script (set up correct options and parameters -- see below for how to do this)

3) After completing a successful execution of the report generator script, you will have a default annotation, a log file, and if you set up your configuration file to generate XML versions of the database those will be present as well.

4) If you do not wish to generate the combined annotation file, you can stop here otherwise, input your XML files back into Blast2GO, so they can be mapped to their respective gene ontology terms.

5) You must also run InterProScan to generate your InterProScan results file (instructions below)

6) Once you have the results from InterProScan and Blast2GO, you can then run the combine_annotations.py script to generate your final, complete annotation file with all of the information from all of the databases.

Important Note:

In order for the statistics on the InterProScan results and Blast2GO results to be calculated, you must re-run the report_generator.py script with your InterProScan results and Blast2GO results.



### Dependencies

This script requires Python 2.X as a dependency, in order to check which version you have installed run the command:

```
python -V
```

Python can be downloaded for free at: https://www.python.org/downloads/

### External Applications

The external applications that are required in order to function with this script are Blast2GO, InterProScan, and usearch (or vsearch)


usearch is a sequence analysis tool.

usearch Download Link: http://www.drive5.com/usearch/download.html

vsearch is a sequence analysis tool that is open source, and free.

vsearch Download Link: https://github.com/torognes/vsearch

Blast2GO is an application that allows for the functional annotation of sequences, and the analysis of annotation data.

Blast2GO Download Link: https://www.blast2go.com/blast2go-pro/download-b2g

InterProScan is an application that provides functional analysis of protein sequences. It is used to generated the gene ontology of the transcripts.

InterProScan Download Link: https://www.ebi.ac.uk/interpro/download.html;jsessionid=B76DDDA8BCBB1AD8AFA31F3FE0E476B5



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
Path to InterProScan Results: (If no results leave blank)
Path to Blast2GO Results: (If no results leave blank)
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

#### InterProScan and Blast2GO paths

For these two options you must enter the path to your InterProScan and Blast2GO input files.

The InterProScan file format is as follows:

```
c10004_g1_i1    3F8F69B8D17CEA7F        185     HMMPfam PF13976 gag_pre-integrs 37      108     1.0E-10 T       06-Nov-2014     IPR025724       GAG-pre-integrase domain
```

The Blast2GO file format is as follows:

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

### Output from the script

The script outputs several files.

There will be a text file of the name: log_[date of execution goes here].txt 

This file will contain the log with all of the statistical data generated in respect to each database, and the collection of databases as a whole.

There will be a text file of the name: nohits_[date of execution goes here].txt

This file will contain the log of all of the nohits found while searching for matches to the search queries in the search results from the databases.

There will be a series of XML files of the name: blastxml_[db number OR combined_db]_[date of execution goes here].xml 

These are the file(s) that can be imported back into Blast2GO as an XML file.

There will be a file of the name: uninformative_[date of execution goes here].txt 

This file will contain the log of all the uninformative hits.

There will be a file of the name: default_output_annotation_2015-03-22 13:07:41.898127.tsv

This file will contain all of the search queries from the query fasta that had a hit, along with their best hit and description.

There will be a file of the name: combined_annotation_[date of execution goes here].tsv

This file will contain all of the search queries from the query fasta that had a hit, their respective best hits, descriptions, and InterProScan and Blast2GO results.
