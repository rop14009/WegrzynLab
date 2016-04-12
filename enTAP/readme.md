# Eukaryote Non-Model Transcriptome Annotation Pipeline (enTAP)

## Introduction

The software tool, enTAP, is a set of Python scripts designed to improve the accuracy and speed of functional annotations for de novo transcriptome (gene space) assemblies of non-model Eukaryotic species.

The scripts in this package enable the following:

* Integrates BLAST-style search results from up to three unique databases (and selects the most optimal annotation)
* Identifies informative and uninformative annotations consistently
* Identifies contaminants through optional filters for fungal, bacterial, and insect annotations.
* Generates BLAST2GO compatible XML files to facilitate Gene Ontology term assignments
* Integrates the results of InterProScan runs for up to two protein domain databases
* Produces a full set of spreadsheet compatible results including summary statistics on the assembly and resulting annotations
 

## Install
### Downloading the required files

enTAP can be downloaded with the following terminal command (Note: git is required to successfuly clone the repository):

```
git clone https://github.com/SamGinzburg/WegrzynLab.git
```

The files will be placed into the "WegrzynLab/enTAP" folder.

The contaminant databases can be donwloaded/updated by utilizing the update.py script (usage demonstrated below)

### Dependencies

This script requires Python 2.X.   To determine the current version, please use:

```
python -V
```

Python can be downloaded for free at: https://www.python.org/downloads/

This script requires the "filter_list.txt" file to be placed in the same directory as report_generator.py in order to function correctly. This script is included in the git repository.

### External Applications

1) Diamond is a really fast sequence-search algorithm which can compare nucleotides sequences to protein sequences for very large datasets. Alternatively blastx can also be used for small datasets. 

Diamond download:

http://github.com/bbuchfink/diamond/releases/download/v0.7.12/diamond-linux64.tar.gz

Alternatively, blastx can be used for smaller datasets where diamond is slower than blastx:

ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

2)  Blast2GO provides a full-featured GUI experience for sequence annotation.  Here, we circumvent some of the limited search options and slow network speeds by limiting its use to Gene Ontology aquisition.  Blast2GO will be run as an indepedent application during the workflow to generate a single output file.  

BLAST2GO download: https://www.blast2go.com/blast2go-pro/download-b2g

3)  InterProScan is a sequence comparison tool for the identification of protein domains and wraps around several public protein domain databases.  It is used here to identify those domains and also Ontology terms associated with them.

InterProScan Download: https://www.ebi.ac.uk/interpro/download.html;jsessionid=B76DDDA8BCBB1AD8AFA31F3FE0E476B5

## Getting Started

### Transcriptome Annotation Overview

![Diagram of script execution](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/diagram.png)

The diagram above shows the order in which the seperate operations must be performed in order to generate a complete annotation report.

Steps:

1) Run your preferred search application (diamond or blastx) to obtain a correctly formatted set of search results from your database.  A total of three database runs can be integrated in the current software version.  The Blast6out tab-delimited text output option is required as is the original FASTA database used.

Execution example for diamond:

For making diamond database file:
```
diamond makedb --in fastadatabase.fasta -d databasename
```

For running diamond, use the following syntax:

```
diamond blastx -d /path_to_diamond_database/databasename -q query.fasta -p no_of_threads -a diamond_alignment_archive.daa -t /directory_for_temp_files --sensitive --evalue evalue_score
```
For obtaining results in blast6 format, use the following command:

```
diamond view -a diamond_alignment_archive.daa -o outputfilename.txt -f tab
```

Execution example for blastx:

For making the blastx database file, use the following command:

```
makeblastdb -in fastadatabase.fasta -out fastadatabase.db -dbtype prot
```

For running blastx, use the following command:

```
blastx -query query.fasta -db fastadatabase.db -out outputfilename.txt -outfmt 6 -num_threads no_of_threads -evalue evalue_score
```


2) Run InterProScan to generate your InterProScan results file.  This can be run for up to 2 databases (specified by the appl flag).

Execution example for InterProScan:

Note: $input is binded to the configuration files for your interproscan run (refer to: https://code.google.com/p/interproscan/wiki/RC4HowToRun)
```
interproscan.sh -i $input -o $output -f xml -appl pfam, Panther -goterms -iprlookup
```

3) Prepare a configuration file with the proper paths and preferences to run the report_generator.py script.  Details on the configuration file are in the section below.

```
Path to query FASTA:
Source Databases (1, 2, 3): 
Sequence Search Application: (diamond)
Query Organism:
Database 1 score (1, 2, 3): 
Path to FASTA-version of database1:
Path to search results from database1:
Database 2 score (1, 2, 3):
Path to FASTA-version of database2:
Path to search results from database2:
Database 3 score (1, 2, 3):
Path to FASTA-version of database3:
Path to search results from database3:
Full-length coverage requirement (0 - 1):
Minimum Evalue (decimal form 0.00001 or scientific notation form 1e-5): 
Generate XML for Blast2GO (yes / no): 
Path to InterProScan Results:
Path to Blast2GO Results:
Contaminant Detection Options
Insects (yes / no): 
Fungi (yes / no): 
Bacteria (yes / no): 
Contaminant Database File Paths
Insects:
Fungi:
Bacteria:
```
4) Run the Report Generator script in the same directory with output files from the diamond/blastx runs and with the configuration file as the only parameter.

```
cd /path/to/directory/containing/script/directory_containing_script/
python report_generator.py configuration_file.txt
```

After completing a successful execution of the report generator script, you will have a tab-delimited annotation file, a log file with summary statistics, flat files for contaminants and nohit sequences, and XML files for each database (and combined) if selected

5) From the combined XML file generated, you have the option to run Blast2GO to assign ontology terms:
Execution example for Blast2GO:

GUI Version:

You begin by importing your query fasta file containing your sequence information.

![Import Sequences](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/load_sequences.png)

![Import Sequences Part 2](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/load_sequences_2.png)

After loading your sequences you must import your Blast2GO XML file containing your sequeunce descriptions.

![Import XML File](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/load_xml.png)

![Import XML File Part 2](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/load_xml_2.png)

![Import XML File Part 3](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/load_xml_3.png)

After importing the sequences and the XML files (in that order), you must run the GO mapping step

![Mapping Step](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/mapping.png)

![Mapping Step](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/mapping_2.png)

After completing these steps, you have successfully mapped the GO terms, and are ready to export the Blast2GO file containing the GO terms.

![Exporting the files](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/export_blast2go.png)

After exporting the text file, you should get a file that looks like this:

![Example File](https://raw.githubusercontent.com/SamGinzburg/WegrzynLab/master/enTAP/map_go_terms_instructions_images/blast2go_file.png)

This file can then be entered back into combine_annotations.py or report_generator.py to generate a combined annotation file.

6) Once you have the results from InterProScan and Blast2GO, you can then run the <i>combine_annotations.py</i> script to generate your final, complete annotation file and summary statistics.

This script takes the Blast2GO, InterProScan, and annotation input and merges them together into one file

@PARAMETERS

--log --input, --blast2go, --interpro, --output
```
python combine_annotations.py --input [annotation filename] --blast2go [blast2go filename] --interpro [interpro filename] --output [output_file_name.tsv]
```
Note this example will generate an output with the default name combined_annotation.tsv in the directory that the script is being executed in
```
python combine_annotations.py --input [annotation filename] --blast2go [blast2go filename] --interpro [interpro filename]
```
The following will specify a custom output name
```
python combine_annotations.py --input [annotation filename] --interpro [interpro filename] --output [output_file_name.tsv]
```
For performing merges with 2 interpro files:
```
python combine_annotations.py --input [annotation filename] --interpro [interpro filename] [interpro filename two] --output [output_file_name.tsv]
```
When the script is complete, the output file will be availible in the output directory.

APPENDING TO AN EXISTING LOG

python combine_annotations.py --log [log name] --input [annotation filename] --interpro [interpro filename] [interpro filename two] --output [output_file_name.tsv]


##Configuration File
In the configuration file, your option is always placed after the ":". If there is no ":" on that line, then it is not an option.

The configuration file is line sensitive, it must always be of the form above, deleting any of the lines in the configuration file will lead to errors or unintended behavior.

```
Path to query FASTA: This is the location in which you will place the path to your query fasta file

Source Databases: (1, 2, 3) Here you will enter the number of databases you wish to run your query fasta file against

Sequence Search Application:(diamond/blastx) Here you will specify the format of the search results from your databases.

Query Organism: This setting is optional, and in the current version of the script only affects the name of the folder in which the results are placed, ex: if "walnut" was entered here the output folder would be "walnut_output_[date]"

Database 1 score: (1, 2, 3) Here you will enter the score of database 1, higher values indicate a higher preference in best hit selection

Path to FASTA-version of database1: Here you will enter the path to the fasta version of this database

Path to search results from database1: Here you will enter the path to the search results on this database

Database 2 score: (1, 2, 3) Here you will enter the score of database 2, higher values indicate a higher preference in best hit selection

Path to FASTA-version of database2: Here you will enter the path to the fasta version of this database

Path to search results from database2: Here you will enter the path to the search results on this database

Database 3 score: (1, 2, 3) Here you will enter the score of database 3, higher values indicate a higher preference in best hit selection

Path to FASTA-version of database3: Here you will enter the path to the fasta version of this database

Path to search results from database3: Here you will enter the path to the search results on this database

Full-length coverage requirement: (0 - 1, not including 0) Here you will enter the full length coverage requirement (recommended starting value = 0.7)

Minimum Evalue: (decimal form 0.00001 or scientific notation form 1e-5) Here you will enter the minimum e-value to be used in best hit selection (recommended starting value = 1e-5)

Generate XML for Blast2GO: (yes / no)

Path to InterProScan Results: Here you will enter the path to your InterProScan results (if you have them), these results can also be appended later via the combine_annotations.py script.

Path to Blast2GO Results: Here you will enter the path to your Blast2GO results (if you have them), these results can also be appended later via the combine_annotations.py script.

Contaminant Detection Options // Select yes to include these hits when running the script, and no to not include them.

Insects: (yes / no)
Fungi: (yes / no)
Bacteria: (yes / no)

Contaminant Database File Paths // These databases can be downloaded using the update.py script

Insects: Path to the insects contaminant database
Fungi: Path to the fungi contaminant database
Bacteria: Path to the bacteria contaminant database
```

#### Query Fasta

At least one database run is necessary and the query FASTA represents the transcriptome assembly and this should be a multi-FASTA text file:

```
>Contig_PtHGO
ATTAAATCCAAACACTATTACTGTGACATTAACAAGGAGAGCAACGGCCTCTGTCTGTATTTGTGTTTGGCAGACGATATTTATGGATTCGCAGCACGTCAATATTTCATGAAAAAAATATCGCCCTAATATAAGATAATCATAACTTTGTTTTTCCTTTGAACTTTACTTCAACCGAAGATCTCGACTGCTGAATTCTGAATTCTGCAATTTGAATGCGCGGGGTTAAAAATGGAGAAATTGGCGAAGCTAAACATTTCAGAAGAGGACTTGGAATTCGAATTTAGCCATCTGGAATACCAGAGCGGATTTGGCAATAATTTTTCCTCGGAAGCAAAACCTGGTGCACTGCCGGTGGGGCAGAATAATCCCTTGAAGTGTCCTTATGGCCTCTACGCCGAGCAAATCTCAGGGACGGCATTTACAGTCCCCAGGAAGCTCAACCAGAGAAGTTGGCTGTACCGTATCAAGCCATCAGTAACACATGAACCTTTTCATCCTCGGGTTCCATCACATGATTATCTGATAAGTGAGTTCAATCAGTCATCAAGCTCAGCAACACCAACGCAATTGCGTTGGAATCCTTTCAACATTCCAGATACGAAAACAAATTTCATTGATGGACTGTTCACTATATGTGGAGCTGGTAGTTCATTCCTTCGGCATGGCTTTGCAGTGCACATGTATGTCGCAAATGCATCTATGGAAGGCTGTGCATTTGCTAATGCTGATGGGGACTTTCTGATAGTTCCACAACAAGGAAGATTATGGATAACGACAGAGTTAGGAAGATTGCAGGTTTCTCCTGGTGAAATAGTTGTGTTGCAACAAGGTTTCCGGTATTCCATTGATTTACCAGATGGTCCTTCTCGAGGGTACGTTGTAGAGGTTTTCAGTGGCCATTTCCAGCTTCCCGATCTTGGTCCAATCGGTGCAAATGGCCTGGCCTCACCGCCTGACTTTCTAACACCTACAGCATGGTTTGATGACAAAGCCTATCCTGGTTACACGATTGTTCATAAGTTTGGAGGTTCATTGTTTACAGCAAAACAAAATTTCTCCCCCTTCAATGTGGTTGCTTGGCATGGCAATTATGCACCTTACAAGTATGATTTGAAGAAGTTCTGTCCTTTCAATACTGTTTTATTTGATCATGGTGACCCATCAATAAACACAGTTCTGACAGTTCCATCAGAAAAACCAGGAGTAGCTGTGGTTGACTTTGTCATTTTTCCACCTCGCTGGCTGGTAGCAGAGCATACCTTCCGTCCTCCATATTACCACCGTAACTGCATGAGTGAATTTATGGGCCTGATTTATGGCAGCTATGAGGCTAAAGCAGATGGATTTCAACCAGGAGGGGCAAGCTTGCACAGCTGCATGACCCCACATGGGCCTGATACTACCACGTTTGAGAAAACTATAGCTGAAGAGGGCGATGCAAAACCTGCCAAGATAAAAGATACCATGGCTTTCATGTTTGAGTCATCATTGATTCCCAGAATCACACCATGGGTTCTCAAATCTCCCCATTTGGATAACGATTATTACAAGTGCTGGACAGGTCTCAAGTCCCATTTTCACCATGAACATCTTCCTGAAAACGGGGAATTACATAGTTAGTAAATCTGTTTCTGTAAATATGGAGAGATATGTACT

>Contig_PtHPPD
ATTTAGTGCAAGAAACAATGGGTAAGGATGGGAGTCAGAGCTTTAAACTGGTGGGATACAAGAACTTTGTTCGCCACAATCCCATGAGCGATCGCTTCCAAGTGAAGAAATTCCACCACATTGAGTTCTGGTGCGGCGATGCCACCAATACCCACAAGCGCTTCTCATGGGGTTTGGGCATGCAGGTGGTGGCCAAGTCCGATCAGAGCACCGGAAATCAGACATATTGCAGCTTCGCCATACAGTCCAACGATCTCACCTTAGTCTTCACGGCGCCTTATTCTTCGAAAATTGACCAGTCCACCACCAGAATGCCCCACCCAGCTTTCGATTCCGATGAGGCGCGGAGTTTCTTTGTCAAGCATGGGTTAGCAGTTCGAGCCGTGGGAGTGCTTGTCGAAGATGCTCGAGACGCTTTCCAGAATAGCGTGGAAAACGGTGCCGTGGCTGTGCTGGAGCCATGCGATCTTGTCGATGAGGCCACGGGCGGGAAGGTCAGTATGTCGGAGGTTAAGCTTTACGGAGATGTGGTTCTTCGTTTTGTAAGCCAGGACGGTTACAAGGGCCCCTTTTTGCCCAATTACGAGCCCGTTCAGTCGGTTCAATTGTCCTATGGAATCATTCGGGTCGACCACGCTGTGGGCAATGTGGAGAAACTGGAGGAAGCCGTCCAATaTGTGGCCAAGTTTAGCGGGTTTCACCGTTTTGCCGAGTTTACGGCAGAAGACGTCGGGACTGGAGAGAGTGGCCTGAATTCGATGGTGTTGGCGAGCAATAATGAGATGGTTCTATTGCCCATGAACGAGCCTGTGTTTGGGACCAAAAGGAAGTCACAGATACAGACCTATCTGGAGCACAATGAGGGCCCCGGATTGCAGCATTTGGCTCTTATTTGCTCTGATATATTCTCCACGCTTAAGGAGATGCGCATGCGCAGTGAAATTGGAGGGTTTGAATTCATGCCCAAGCCTCCTCCTACTTATTATAAGAATCTGAGGAATAGGGTTGCAGATGTCTTGACTGAGAAGCAGATGGAAGAGTGCGACGAGTTGGGGATTTTGGTCGATAGGGATGAtCAAGGAGTCCTCTTACAGATTTTCACAAAGCCCATTGGCGACAGGCCCACTATTTTTATTGAGATAATTCAACGTGTAGGCTGTGTGCTAACTGATGAGAATGGGAAAGCCTATCAAAAGGGTGGATGTGGAGGATTTGGCAAAGGCAATTTTTCAGAACTCTTTAAATCCATCGAGGAATACGAGAAAACACTTGAAGGAAGAGTCAATTCCTAGGCCCTTCGAAGTTCGTAGTGATCTCCTAGAGTAGATTCATCTTGGAATTAGACTGATAGAGATTTATTATATACGTGAGAGTGGCATAATGATAAAATACATGTACTACGGTTGCTTGTGTTAATGGTTTCTAGATGTCGTTTAGTGAAGATGCCTTGTAAAGAAGCAAGGTATATATGTTATACAAATATACTCTGTTTCTGGTTTGGTTCTAGGCTGCATGCCCTTAAGAGGAATTTAACTGTTATAATTTAAATATATTATGTTCAGATAGGCCCACGGCTTATCATTTTATAATTTAAATATATTCTGTTCAGATATGTCCACGGCTGGTTAAGTTCACTTGTCAAGTAACTCCATTTTTATTCACGCTGGAGCGTAATTTTGGATCTAGGGATATAGAATCAGGATGGATGCTCTCTTTAAAGCTTATAACATATTTATGGTTCGT

```

#### Source Databases

This option contains the number of databases to be parsed by the script, and can be a number between 0 and 3.

#### Sequence Search Application

This script currently supports 2 different sequence search applications, namely diamond and blastx.

The correct format looks like this:

```
Contig_PtHGO    gi|460412999|ref|XP_004251883.1|        78.7    432     92      0       286     1581    12      443     4.3e-212        744.2
```

#### Paths to the Databases

The current version of the script does not support indexed versions of databases so the flat file (multi-FASTA version) of the database must be supplied for each database queried.

#### Minimum E-Value and Minimum coverage requirements

The minimum E-Value can be input in either scientific notation (1e-5) or as a decimal (0.00001).

The minimum coverage requirement is a value between 0 and 1 where 70% coverage represented as 0.7 is the default value.

#### InterProScan and Blast2GO paths

For these two options you must enter the path to your InterProScan and Blast2GO output file (Sequence Table format).

The InterProScan file format is as follows:

```
c10004_g1_i1    3F8F69B8D17CEA7F        185     HMMPfam PF13976 gag_pre-integrs 37      108     1.0E-10 T       06-Nov-2014     IPR025724       GAG-pre-integrase domain
```

The Blast2GO Sequence Table format is as follows:

```
c58101_g1_i1|m.1        uninformative   390     2       7.2E-10 0.0%    3       F:GO:0003674; P:GO:0008150; C:GO:0005575
```

#### Contaminant DB selection

This script currently has three optional flat file databases containing genus and species information of organisms that could be filtered as potential contaminants in Eukaryote transcriptome studies:

* Insects
* Bacteria
* Fungi

In in configuration file, each of these can be toggled by writing either "y" or "yes" to include them when searching for contaminants, or write "no" or leave it blank to ignore these filters.

#### Downloading updated contaminant databases using update.py

The update.py script will automatically download all of the contaminant files required to run report_generator.py

The files will be deposited into a folder titled "contaminant_databases".

correct update.py usage:
```
python update.py
```

Note: An internet connection is required to run update.py.

## Output

The script outputs several files.

There will be a text file of the name: log_[date of execution goes here].txt 

This file will contain the log with all of the statistical data generated in respect to each database, and the collection of databases as a whole.
The file will be of the form:

```
Weighted Summary Statistics on the Transcriptome Assembly Input (Query)
Total Number of Query Sequences:        29785
Median Query Length:    825
Average Query Length:   1048.32
Shortest Query Length:  297
Longest Query Length:   15345
N50 Statistic:  1344
Number of queries with an informative hit:      16183
Number of queries with an uninformative hit:    9404
Number of queries with no hit:  4140
Number of contaminants:         58
The top 10 hits by species:
uninformative:  9404
PREDICTED: putative ribonuclease H protein At1g65750-like:      94
ATP binding protein, putative:  79
DNA binding protein, putative:  64
transcription factor, putative: 60
membrane protein:       59
PREDICTED: TMV resistance protein N-like:       46
protein binding protein, putative:      45
zinc finger protein, putative:  41
transcriptional regulator:      35
The top 10 contaminants:
Escherichia coli        28
Laccaria bicolor        4
Zymoseptoria tritici    3
Coprinopsis cinerea     3
Magnaporthe oryzae      2
Ustilago maydis 2
Pyrenophora teres       2
Cryptococcus gattii     2
Pyrenophora tritici-repentis    2
Mycobacterium kansasii  1

Summary Statistics on the Transcriptome Assembly Input (Query)
Total Number of Query Sequences:        29785
Median Query Length:    825
Average Query Length:   1048.32
Shortest Query Length:  297
Longest Query Length:   15345
N50 Statistic:  1344
Number of queries with an informative hit:      16183
Number of queries with an uninformative hit:    9404
Number of queries with no hit:  4140
Number of contaminants:         58
The top 10 hits by species:
uninformative:  9404
PREDICTED: putative ribonuclease H protein At1g65750-like:      94
ATP binding protein, putative:  79
DNA binding protein, putative:  64
transcription factor, putative: 60
membrane protein:       59
PREDICTED: TMV resistance protein N-like:       46
protein binding protein, putative:      45
zinc finger protein, putative:  41
transcriptional regulator:      35
The top 10 contaminants:
Escherichia coli        28
Laccaria bicolor        4
Zymoseptoria tritici    3
Coprinopsis cinerea     3
Magnaporthe oryzae      2
Ustilago maydis 2
Pyrenophora teres       2
Cryptococcus gattii     2
Pyrenophora tritici-repentis    2
Mycobacterium kansasii  1

DB:     /path/to/file/combined_genes_blast6out_ppf.out
Number of queries with an informative hit:      13723
Number of queries with an uninformative hit:    10857
Number of contaminants: 0
Number of queries with no hit:  5205
The top 10 hits by species:
uninformative:  10857
PREDICTED: putative ribonuclease H protein At1g65750-like:      85
ATP binding protein, putative:  62
DNA binding protein, putative:  48
PREDICTED: TMV resistance protein N-like:       45
transcription factor, putative: 40
protein binding protein, putative:      38
pentatricopeptide repeat-containing protein, putative:  29
PREDICTED: UPF0481 protein At3g47200-like:      29
zinc finger family protein:     28
The top 10 contaminants by species:
No contaminants present

DB:     /path/to/file/combined_genes_blast6out_rsp.out
Number of queries with an informative hit:      16103
Number of queries with an uninformative hit:    9196
Number of contaminants: 58
Number of queries with no hit:  4428
The top 10 hits by species:
uninformative:  9196
PREDICTED: putative ribonuclease H protein At1g65750-like:      85
ATP binding protein, putative:  78
DNA binding protein, putative:  64
transcription factor, putative: 60
membrane protein:       59
protein binding protein, putative:      45
PREDICTED: TMV resistance protein N-like:       44
zinc finger protein, putative:  41
transcriptional regulator:      35
The top 10 contaminants by species:
Escherichia coli:       28
Laccaria bicolor:       4
Zymoseptoria tritici:   3
Coprinopsis cinerea:    3
Magnaporthe oryzae:     2
Ustilago maydis:        2
Pyrenophora teres:      2
Cryptococcus gattii:    2
Pyrenophora tritici-repentis:   2
Mycobacterium kansasii:         1

Interpro File:  /path/to/file/interpro.raw
Blast2GO File:  /path/to/file/blast2go_walnut_transcriptome.txt
Number of sequences with Domain Identification:         25587
Number of sequences without Domain Identification:      0
Blast2GO Gene Ontology Stats
Number of Components:   4552
Number of Functions:    326
Number of Processes:    4986
Number of Transcripts with at least 1 Component:        4552
Number of Transcripts with at least 1 Function:         4986
Number of Transcripts with at least 1 Process:  326
Interpro Gene Ontology Stats (Totals)
Component:      2859
Function:       8690
Process:        13303
Number of Transcripts with at least 1 Component:        2421
Number of Transcripts with at least 1 Function:         7622
Number of Transcripts with at least 1 Process:  6826
```


There will be a text file of the name: nohits_[date of execution goes here].txt
The file will be of the form:
```
c54972_g2_i2|m.50780
c54259_g1_i5|m.2219
c37958_g2_i1|m.6596
c55455_g2_i1|m.35545
c51509_g3_i1|m.34497
c54913_g5_i1|m.11269
c29778_g1_i1|m.6050
c970_g1_i1|m.12443
c53384_g7_i1|m.30991
```
This file will contain the log of all of the nohits found while searching for matches to the search queries in the search results from the databases.

There will be a series of XML files of the name: blastxml_[db number OR combined_db]_[date of execution goes here].xml 

These are the file(s) that can be imported back into Blast2GO as an XML file.

There will be a file of the name: uninformative_[date of execution goes here].txt 

This file will contain the log of all the uninformative hits.

There will be a file of the name: default_output_annotation_2015-03-22 13:07:41.898127.tsv

This file will contain all of the search queries from the query fasta that had a hit, along with their best hit and description.

This file will be of the form:
```
Query Subject_id Identity(%) Alignment_length Mismatches Number of gap opens Query_start Query_end Subject_start Subject_end E-value Bit_score Origin Database Subject Description Species
c48227_g1_i1|m.33460 gi|255540671|ref|XP_002511400.1| 66.2 542 168 7 1 1605 1 534 1.4e-198 699.5 /group/nealedata3/Walnut_Transcriptome/large_set/combined/combined_genes_blast6out_rsp.out uninformative Ricinus communis
```

There will be a file of the name: combined_annotation_[date of execution goes here].tsv

This file will be of the form:
```
Query Subject_id Identity(%) Alignment_length Mismatches Number of gap opens Query_start Query_end Subject_start Subject_end E-value Bit_score Origin Database Subject Description Species Signature Description InterPro accession number InterPro description blast2go_process blast2go_function blast2go_component
c48227_g1_i1|m.33460 gi|255540671|ref|XP_002511400.1| 66.2 542 168 7 1 1605 1 534 1.4e-198 699.5 /group/nealedata3/Walnut_Transcriptome/large_set/combined/combined_genes_blast6out_rsp.out uninformative Ricinus communis Pex14_N IPR006785 Peroxisome membrane anchor protein Pex14p, N-terminal
c55482_g1_i1|m.33976 gi|449513131|ref|XP_004164240.1| 50.1 814 367 16 28 2406 8 803 2.5e-207 729.2 /group/nealedata3/Walnut_Transcriptome/large_set/combined/combined_genes_blast6out_rsp.out PREDICTED: probable receptor-like protein kinase At2g23200-like Cucumis sativus Pkinase_Tyr, Malectin_like IPR001245, IPR024788 Serine-threonine/tyrosine-protein kinase catalytic domain, Malectin-like carbohydrate-binding domain Molecular Function:protein kinase activity (GO:0004672), Biological Process:protein phosphorylation (GO:0006468)  F:protein serine/threonine kinase activity
```


This file will contain all of the search queries from the query fasta that had a hit, their respective best hits, descriptions, and InterProScan and Blast2GO results.
