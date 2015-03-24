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

## Execution Example


```
cd /path/to/directory/containing/script/directory_containing_script/
python report_generator.py
```

## Getting Started

### Dependencies

This script requires Python 2.0 as a dependency, in order to check which version you have installed run the command:

```
python -V
```

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




