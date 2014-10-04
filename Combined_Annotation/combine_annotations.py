'''
Author: Sam Ginzburg

Description: An implementation of a program that combines blast2go, interpro, and previous annotation files

@param The input files should all be placed in the 'input' directory adjacent to the location of this script

'''

import time
import sys
import os
import csv

def combined_anno_id_parser(id):
    current_char = ""
    return_string = ""
    delimiter = "|"
    length = 0
    while current_char != delimiter:
        return_string += current_char
        current_char = id[length:length+1]
        length+=1
   
    return return_string

def parse_blast2go_column(row):
    return_string = "";
    split_string = row[7].split(";") # splits the column by ;
    #print(split_string)
    for string in split_string:
        string = string.strip()
        if string[:2] == "P:" or string[:2] == "F:" or string[:2] == "C:":
            return_string = return_string + string + "; " #if the string begins with the proper beginning, add it to the return_string
    return return_string[:-2] # substring taken to just remove the last ; deliminator for aesthetic reasons

def make_blast2go_walnut_combined_dict(file_name):
    blast2go_data = dict()
    with open(file_name,'r') as tsv_old:
        row = next(tsv_old) #skip the first header line
        for row in csv.reader(tsv_old, delimiter='\t'):
            id = combined_anno_id_parser(row[0]) #get the key
            blast2go_data[id] = parse_blast2go_column(row)
    print("blast2go hashtable complete")
    return blast2go_data

def make_walnut_interpro_dict(file_name):
    walnut_interpro_data = dict()
    value_list = list()
    with open(file_name,'r') as tsv_old:
        for row in csv.reader(tsv_old, delimiter='\t'):
            id = row[0][8:] #get the key
            value_list.append(row[5])
            value_list.append(row[11])
            value_list.append(row[12])
            walnut_interpro_data[id] = value_list
            value_list = list()
    print("walnut_interpro hastable complete")
    
    return walnut_interpro_data

if __name__ == '__main__':
    start_time = time.clock()
    

    file_path_combined_anno = os.path.dirname(os.path.realpath(__file__)) + '\\input\\annotation.tsv'
    file_path_blast2go_walnut = os.path.dirname(os.path.realpath(__file__)) + '\\input\\blast2go.txt'
    file_path_walnut_interpro = os.path.dirname(os.path.realpath(__file__)) + '\\input\\interpro.raw'

    #first make dicts for walnut interpro and blast2go for easier searching based off of IDs
    walnut_interpro_hashtable = make_walnut_interpro_dict(file_path_walnut_interpro)
    blast2go_hastable = make_blast2go_walnut_combined_dict(file_path_blast2go_walnut)
       
        
    #then go through old combined anno and line by line make new anno file from previous data
       
    with open(file_path_combined_anno,'r') as tsv_old, \
    open(os.path.dirname(os.path.realpath(__file__)) + '\\output\\combined_annotation.tsv', 'w') as tsv_new:
        tsv_new = csv.writer(tsv_new, delimiter='\t')
        tsv_old = csv.reader(tsv_old, delimiter='\t')
        row = next(tsv_old) # the purpose of this line is to skip the header in the csv file, so there is no need to iterate another 10K plus times through blast2go_table_walnut/walnut_interpro
        combined_row = row + ["walnut 5", "walnut 11", "walnut 12", "blast2go_process_column"]
        tsv_new.writerow(combined_row) # copy the header
        for row in tsv_old:
            id = combined_anno_id_parser(row[0])
                
            walnut_results_interpro = walnut_interpro_hashtable.get(id)
                
            if not walnut_results_interpro: #this is to create the blank spaces if there are no interpro results for corresponding IDs
                walnut_results_interpro = ["N/A","N/A","N/A"]
               
            walnut_results_blast2go = blast2go_hastable.get(id)
                
            if not walnut_results_blast2go: # if there is no match from blast2go use as filler
                walnut_results_blast2go = "N/A"
                
            combined_row = row + walnut_results_interpro + [walnut_results_blast2go]
            tsv_new.writerow(combined_row)
    print (str(time.clock() - start_time) + " seconds")
    print("complete -- file now availible in output directory")