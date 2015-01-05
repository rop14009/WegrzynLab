
Annotation Report Generator Readme

Author: Sam Ginzburg

Files Included:

report_generator.py
combine_annotations.py
bacteria_db.txt
filter_list.txt
fungi_db.txt
configuration_file.txt


*********************************************************************************************************************
Install Instructions:
*********************************************************************************************************************
In order to install the script the files above must all be copied to the desired directory on your machine.

*********************************************************************************************************************
Setup Instructions:
*********************************************************************************************************************
After installation and before running the script, the configuration_file.txt file must be filled out completely.

In configuration_file.txt you must enter the query fasta, fasta version of at least 1 database, the number of databases, the type of databse, the minimum coverage 
requirement, the minimum e-value, and the search results of a least 1 database.

If you would like to solely append interpro output you can change the number of databses to '0' in the script and then put in an interpro filepath.

You can toggle BlastXML output by typing yes/no after the option in the configuration file.

At this time the Query Organism option doesn't perform an action, as well as the paths to the formatted databases.

*NOTE* When entering values into the minimum e-value or coverage requirement, be sure to type the coverage requirement as a decimal from 0 to 1, and the minimum
e-value in either scientific notation or in decimal form.

An example of a successfull configuration:
--------------------------------------------------------------------------------------------------
Path to query FASTA: /home/dgi/question/sequences3.fasta
Source Databases: 1
Sequence Search Application:usearch
Query Organism:
Path to formatted database1:
Path to FASTA-version of database1:/share/nealedata/databases/downloadDir/refseq_protein.fasta
Path to search results from database1:/home/dgi/question/sequences3_refseq_protein
Path to formatted database2 (organism-specific):
Path to FASTA-version of database2:
Path to search results from database2:
Path to formatted database3 (full length): 
Path to FASTA-version of database3:
Path to search results from database3:
Full-length coverage requirement: 0.7
Minimum Evalue:1e-5
Generate XML for Blast2GO:yes
Path to Interpro Results: 
--------------------------------------------------------------------------------------------------



********************************************************************************************************************* 
Running the script:
*********************************************************************************************************************

After setting up the script properly it can be run by typing:
--------------------------
python report_generator.py
--------------------------
into your terminal.

*********************************************************************************************************************
Compatability
*********************************************************************************************************************

This script has been tested with

Python 2.7.6
Ubuntu 14.04 Linux kernel version 3.13.0-36-generic


and should work with other versions of python/operating systems but has not been tested on them yet.


