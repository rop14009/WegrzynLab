# Transcriptome Annotation Pipeline

Features

 * Input multiple databases
 * Determine Best hits
 * Export XML versions of the resulting databases back to Blast2GO

# Instructions on how to use the script

1) Configuration File Format

Example File


```
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
Path to Blast2GO Results:

``` 


The first line indicates the file containing the query transcripts of the following form:

```
>Contig_PtHGO
ATTAAATCCAAACACTATTACTGTGACATTAACAAGGAGAGCAACGGCCTCTGTCTGTATTTGTGTTTGGCAGACGATATTTATGGATTCGCAGCACGTCAATATTTCATGAAAAAAATATCGCCCTAATATAAGATAATCATAACTTTGTTTTTCCTTTGAACTTTACTTCAACCGAAGATCTCGACTGCTGAATTCTGAATTCTGCAATTTGAATGCGCGGGGTTAAAAATGGAGAAATTGGCGAAGCTAAACATTTCAGAAGAGGACTTGGAATTCGAATTTAGCCATCTGGAATACCAGAGCGGATTTGGCAATAATTTTTCCTCGGAAGCAAAACCTGGTGCACTGCCGGTGGGGCAGAATAATCCCTTGAAGTGTCCTTATGGCCTCTACGCCGAGCAAATCTCAGGGACGGCATTTACAGTCCCCAGGAAGCTCAACCAGAGAAGTTGGCTGTACCGTATCAAGCCATCAGTAACACATGAACCTTTTCATCCTCGGGTTCCATCACATGATTATCTGATAAGTGAGTTCAATCAGTCATCAAGCTCAGCAACACCAACGCAATTGCGTTGGAATCCTTTCAACATTCCAGATACGAAAACAAATTTCATTGATGGACTGTTCACTATATGTGGAGCTGGTAGTTCATTCCTTCGGCATGGCTTTGCAGTGCACATGTATGTCGCAAATGCATCTATGGAAGGCTGTGCATTTGCTAATGCTGATGGGGACTTTCTGATAGTTCCACAACAAGGAAGATTATGGATAACGACAGAGTTAGGAAGATTGCAGGTTTCTCCTGGTGAAATAGTTGTGTTGCAACAAGGTTTCCGGTATTCCATTGATTTACCAGATGGTCCTTCTCGAGGGTACGTTGTAGAGGTTTTCAGTGGCCATTTCCAGCTTCCCGATCTTGGTCCAATCGGTGCAAATGGCCTGGCCTCACCGCCTGACTTTCTAACACCTACAGCATGGTTTGATGACAAAGCCTATCCTGGTTACACGATTGTTCATAAGTTTGGAGGTTCATTGTTTACAGCAAAACAAAATTTCTCCCCCTTCAATGTGGTTGCTTGGCATGGCAATTATGCACCTTACAAGTATGATTTGAAGAAGTTCTGTCCTTTCAATACTGTTTTATTTGATCATGGTGACCCATCAATAAACACAGTTCTGACAGTTCCATCAGAAAAACCAGGAGTAGCTGTGGTTGACTTTGTCATTTTTCCACCTCGCTGGCTGGTAGCAGAGCATACCTTCCGTCCTCCATATTACCACCGTAACTGCATGAGTGAATTTATGGGCCTGATTTATGGCAGCTATGAGGCTAAAGCAGATGGATTTCAACCAGGAGGGGCAAGCTTGCACAGCTGCATGACCCCACATGGGCCTGATACTACCACGTTTGAGAAAACTATAGCTGAAGAGGGCGATGCAAAACCTGCCAAGATAAAAGATACCATGGCTTTCATGTTTGAGTCATCATTGATTCCCAGAATCACACCATGGGTTCTCAAATCTCCCCATTTGGATAACGATTATTACAAGTGCTGGACAGGTCTCAAGTCCCATTTTCACCATGAACATCTTCCTGAAAACGGGGAATTACATAGTTAGTAAATCTGTTTCTGTAAATATGGAGAGATATGTACT
```


The second line indicates the number of databases to be parsed, and the sequence search application can be specified to either usearch or ncbi (all lower case). This refers to the format of the file that is being input into the script.

Note: It is possible to enter 0 as an option, and if you input a number greater than 3, the script will currently not parse any of the databases.

usearch search results format example:
```
Contig_PtPAL2   gi|357150263|ref|XP_003575399.1|        72.4    58      16      0       1401    1228    15      72      1.6e-15 88.2
```

fasta version of database format example:
```
>gi|15674171|ref|NP_268346.1| 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis Il1403] >gi|116513137|ref|YP_812044.1| 30S ribosomal protein S18 [Lactococcus lactis subsp. cremo
ris SK11] >gi|125625229|ref|YP_001033712.1| 30S ribosomal protein S18 [Lactococcus lactis subsp. cremoris MG1363] >gi|281492845|ref|YP_003354825.1| 30S ribosomal protein S18 [Lactococcus
lactis subsp. lactis KF147] >gi|385831755|ref|YP_005869568.1| 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis CV56] >gi|385839508|ref|YP_005877138.1| 30S ribosomal protein S18
 [Lactococcus lactis subsp. cremoris A76] >gi|389855617|ref|YP_006357861.1| 30S ribosomal protein S18 [Lactococcus lactis subsp. cremoris NZ9000] >gi|414075194|ref|YP_007000411.1| 30S ribosomal protein S18 [Lactococcus lactis subsp. cremoris UC509.9] >gi|459286377|ref|YP_007509482.1| 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis IO-1] >gi|544395586|ref|YP_008569870.1| ribosomal protein S18 RpsR [Lactococcus lactis subsp. cremoris KW2] >gi|554464728|ref|YP_008703967.1| 30S ribosomal protein S18 [Lactococcus lactis subsp. lactis KLDS 4.0325] >gi|489223532|ref|WP_003131952.1| 30S ribosomal protein S18 [Lactococcus lactis]
MAQQRRGGFKRRKKVDFIAANKIEVVDYKDTELLKRFISERGKILPRRVTGTSAKNQRKVVNAIKRARVMALLPFVAEDQ
N
>gi|66816243|ref|XP_642131.1| hypothetical protein DDB_G0277827 [Dictyostelium discoideum AX4]
MASTQNIVEEVQKMLDTYDTNKDGEITKAEAVEYFKGKKAFNPERSAIYLFQVYDKDNDGKITIKELAGDIDFDKALKEY
KEKQAKSKQQEAEVEEDIEAFILRHNKDDNTDITKDELIQGFKETGAKDPEKSANFILTEMDTNKDGTITVKELRVYYQK
VQKLLNPDQ
>gi|66818355|ref|XP_642837.1| hypothetical protein DDB_G0276911 [Dictyostelium discoideum AX4]
MKTKSSNNIKKIYYISSILVGIYLCWQIIIQIIFLMDNSIAILEAIGMVVFISVYSLAVAINGWILVGRMKKSSKKAQYE
DFYKKMILKSKILLSTIIIVIIVVVVQDIVINFILPQNPQPYVYMIISNFIVGIADSFQMIMVIFVMGELSFKNYFKFKR
IEKQKNHIVIGGSSLNSLPVSLPTVKSNESNESNTISINSENNNSKVSTDDTINNVM
```

# Full length Coverage requirement and E-Value

These two parameters are used in the best hit detection process. The e-value can be entered either as a decimal (0.000001) or in scientific notation (1e-5). The full length coverage requirement can be entered as a decimal between 0 and 1 (inclusive).

# XML Generation

If this parameter is toggled to "yes", then XML versions of the output of the script will be generated. There will be 1 xml file per database entered into the parameters, and a combined xml file will also be generated. The files are named according to which database number they correspond to, so blastxml_db1 will refer to the database entered in the paths given for database 1.

If this parameter is toggled to "no" then no XML databases will be generated.


# Interpro / Blast2GO

These two parameters allow you to input file paths to corresponding interpro and blast2go outputs to parse the gene ontologies.

Example Interpro:
```
Wal_cds.c10004_g1_i1    3F8F69B8D17CEA7F        185     HMMPfam PF13976 gag_pre-integrs 37      108     1.0E-10 T       06-Nov-2014     IPR025724       GAG-pre-integrase domain
```
Example Blast2GO:
```
Seq. Name       Seq. Description        Seq. Length     #Hits   min. eValue     mean Similarity #GOs    GOs     Enzyme Codes    InterProScan
c26233_g2_i1|m.5        nb-arc domain-containing disease resistance protein     339     1       4.8E-10 0.0%    2       P:GO:0006952; C:GO:0005575
```

# IMPORTANT NOTE

In order for the script to function correctly, the sequence names between the usearch results, query fasta file, and interpro/blast2go files should match.

For example:

Wal_cds.c10004_g1_i1 and c10004_g1_i1 will not match with each other, so the sequence names should be: c10004_g1_i1 and c10004_g1_i1.

If isoforms are present in the databases, it is necessary that the query fasta and the databases both include the isoform identifier as c26233_g2_i1|m.5 and c26233_g2_i1 will also not match with each other.





