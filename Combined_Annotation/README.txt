This script takes the blast2go, interpro, and annotation input and merges them together into one file

@PARAMETERS

--blast2go, --interpro, --o, --h

The blast2go and interpro command line options indicate whether you will be using one, the other, or both to generate the combined annotation

If you use both of these command line arguments you should type the initial annotation file name, then the blast2go file, then the interpro file


IF the --o flag is not used the default output file name is combine_annotations.tsv


***NOTE***

In all situations if the --o command line argument is used to generate a custom file name for the output the last argument for the script MUST BE THE NEW FILE NAME, otherwise the script will terminate

In addition to this if the --h help flag is used, the the script will print out a series of examples usages but will then TERMINATE, regardless of the rest of the flags/arguments

***NOTE***

Example usages of the script that will result in correct outputs:


py combine_annotations.py --blast2go --interpro annotation.tsv blast2go.txt interpro.raw output_file_name.tsv


py combine_annotations.py --interpro annotation.tsv interpro.raw


py combine_annotations.py --blast2go annotation.tsv blast2go.txt output_file_name.tsv


When the script finishes running, the output file will be availible in the output directory.


