Author: Sam Ginzburg

This script takes the blast2go, interpro, and annotation input and merges them together into one file


@PARAMETERS

--input, --blast2go, --interpro, --output



Example usages of the script that will result in correct outputs:


py combine_annotations.py --input <annotation filename> --blast2go <blast2go filename> --interpro <interpro filename> --output <output_file_name.tsv>


#note this example will generate an output with the default name combined_annotation.tsv in the directory that the script is being executed in
py combine_annotations.py --input <annotation filename> --blast2go <blast2go filename> --interpro <interpro filename>

py combine_annotations.py --input <annotation filename> --interpro <interpro filename> --output <output_file_name.tsv>



For performing merges with 2 interpro files:

py combine_annotations.py --input <annotation filename> --interpro <interpro filename> <interpro filename two> --output <output_file_name.tsv>



When the script finishes running, the output file will be availible in the output directory.


