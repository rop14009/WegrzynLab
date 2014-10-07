'''
Author: Sam Ginzburg

Description: An implementation of a program that combines blast2go, interpro, and previous annotation files

@param 	at least one of the following flags are required to run the program:  --blast2go --interpro --o , (--o : the output file name is optional, but if included MUST be present)
		
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
	process = list()
	function = list()
	component = list()
	return_string = ""
	split_string = row[7].split(";") # splits the column by ;
	#print(split_string)
	for string in split_string:
		string = string.strip()
		if string[:2] == "P:":
			process.append(string)
		elif string[:2] == "F:":
			function.append(string)
		elif string[:2] == "C:":
			component.append(string)
	return [" ".join(process), " ".join(function), " ".join(component)]

def make_blast2go_walnut_combined_dict(file_name):
	blast2go_data = dict()
	with open(file_name,'r') as tsv_old:
		row = next(tsv_old) #skip the first header line
		for row in csv.reader(tsv_old, delimiter='\t'):
			id = combined_anno_id_parser(row[0]) #get the key
			blast2go_data[id] = parse_blast2go_column(row)
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
	return walnut_interpro_data

def print_usage():
	print("To run this script use at least one of the following flags:")
	print("--blast2go --interpro")
	print("EX: py combine_annotations.py --blast2go --interpro annotation.tsv blast2go.txt interpro.raw output_file_name.tsv")
	print("EX: py combine_annotations.py --interpro annotation.tsv interpro.raw")
	print("EX: py combine_annotations.py --blast2go annotation.tsv blast2go.txt output_file_name.tsv")
	sys.exit(-1)

def parse_flags(param_list):
	blast2go = 0
	interpro = 0
	custom_output_name = 0
	for param in param_list:
		if param == "--h":
			print_usage()
		elif param == "--blast2go":
			blast2go += 1
		elif param == "--interpro":
			interpro += 1
		elif param == "--o":
			custom_output_name += 1
	
	if interpro == 1 and blast2go == 1:
		if custom_output_name == 1:
			return "all_params_custom_output"
		return "all_params"
	elif interpro == 1 and blast2go == 0:
		if custom_output_name == 1:
			return "interpro_custom_output"
		return "interpro"
	elif interpro == 0 and blast2go == 1:
		if custom_output_name == 1:
			return "blast2go_custom_output"
		return "blast2go"
	else:
		print_usage()
		return "incorrect_parameters"


		
if __name__ == '__main__':
	start_time = time.clock()
	arguments_list = sys.argv
	params = parse_flags(arguments_list)
	print (arguments_list)
	
	
	if params == "all_params" or params == "all_params_custom_output":
	
	
		if params == "all_params_custom_output":
			output = arguments_list[len(arguments_list) - 1]
			file_path_combined_anno = arguments_list[4]
			file_path_blast2go_walnut = arguments_list[5]
			file_path_walnut_interpro = arguments_list[6]
			
			if output == file_path_walnut_interpro or output == file_path_combined_anno or output == file_path_blast2go_walnut:
				print ("WARNING user entered option for custom file name and entered name of existing file (would cause overwrite) -- aborting")
				print_usage()
				sys.exit(-1)
		else:
			output = "combine_annotations.tsv"
			file_path_combined_anno = arguments_list[3]
			file_path_blast2go_walnut = arguments_list[4]
			file_path_walnut_interpro = arguments_list[5]
		

		
		#first make dicts for walnut interpro and blast2go for easier searching based off of IDs
		walnut_interpro_hashtable = make_walnut_interpro_dict(file_path_walnut_interpro)
		blast2go_hastable = make_blast2go_walnut_combined_dict(file_path_blast2go_walnut)
		
		with open(file_path_combined_anno,'r') as tsv_old, \
		open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
			tsv_new = csv.writer(tsv_new, delimiter='\t')
			tsv_old = csv.reader(tsv_old, delimiter='\t')		
			row = next(tsv_old) # the purpose of this line is to skip the header in the csv file, so there is no need to iterate another 10K plus times through blast2go_table_walnut/walnut_interpro
			combined_row = row + ["walnut 5", "walnut 11", "walnut 12", "blast2go_process","blast2go_function","blast2go_component"]
			tsv_new.writerow(combined_row) # copy the header
			for row in tsv_old:
				id = combined_anno_id_parser(row[0])
				
				walnut_results_interpro = walnut_interpro_hashtable.get(id)
				
				if not walnut_results_interpro: #this is to create the blank spaces if there are no interpro results for corresponding IDs
					walnut_results_interpro = ["N/A","N/A","N/A"]
					
				walnut_results_blast2go = blast2go_hastable.get(id)
				
				if not walnut_results_blast2go: # if there is no match from blast2go use as filler
					walnut_results_blast2go = ["N/A","N/A","N/A"]
					
				combined_row = row + walnut_results_interpro + walnut_results_blast2go
				tsv_new.writerow(combined_row)	
				
	elif params == "interpro" or params == "interpro_custom_output":
		
		print (arguments_list[len(arguments_list) - 1])
		
		if params == "interpro_custom_output":
			output = arguments_list[len(arguments_list) - 1]
			file_path_combined_anno = arguments_list[3]
			file_path_walnut_interpro = arguments_list[4]
			if output == file_path_walnut_interpro or output == file_path_combined_anno:
				print ("WARNING user entered option for custom file name and entered name of existing file (would cause overwrite) -- aborting")
				print_usage()
				sys.exit(-1)
			
		else:
			output = "combine_annotations.tsv"
			file_path_combined_anno = arguments_list[2]
			file_path_walnut_interpro = arguments_list[3]
		
		print (output)
		file_path_blast2go_walnut = "None"
		
		#first make dicts for walnut interpro and blast2go for easier searching based off of IDs
		walnut_interpro_hashtable = make_walnut_interpro_dict(file_path_walnut_interpro)
		
		with open(file_path_combined_anno,'r') as tsv_old, \
		open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
			print(os.path.realpath(__file__) + "//" + output)
			tsv_new = csv.writer(tsv_new, delimiter='\t')
			tsv_old = csv.reader(tsv_old, delimiter='\t')		
			row = next(tsv_old) # the purpose of this line is to skip the header in the csv file, so there is no need to iterate another 10K plus times through blast2go_table_walnut/walnut_interpro
			combined_row = row + ["walnut 5", "walnut 11", "walnut 12", "blast2go_process","blast2go_function","blast2go_component"]
			tsv_new.writerow(combined_row) # copy the header
			for row in tsv_old:
				id = combined_anno_id_parser(row[0])
				
				walnut_results_interpro = walnut_interpro_hashtable.get(id)
				
				if not walnut_results_interpro: #this is to create the blank spaces if there are no interpro results for corresponding IDs
					walnut_results_interpro = ["N/A","N/A","N/A"]
					
				combined_row = row + walnut_results_interpro + ["N/A","N/A","N/A"]
				tsv_new.writerow(combined_row)		
		
	elif params == "blast2go" or params == "blast2go_custom_output":
		
		if params == "blast2go_custom_output":
			output = arguments_list[len(arguments_list) - 1]
			file_path_combined_anno = arguments_list[3]
			file_path_blast2go_walnut = arguments_list[4]
			if output == file_path_blast2go_walnut or output == file_path_combined_anno:
				print ("WARNING user entered option for custom file name and entered name of existing file (would cause overwrite) -- aborting")
				print_usage()
				sys.exit(-1)
		else:
			output = "combine_annotations.tsv"
			file_path_combined_anno = arguments_list[2]
			file_path_blast2go_walnut = arguments_list[3]
		
	
		file_path_walnut_interpro = ""
		
		#first make dicts for walnut interpro and blast2go for easier searching based off of IDs
		blast2go_hastable = make_blast2go_walnut_combined_dict(file_path_blast2go_walnut)
		
		with open(file_path_combined_anno,'r') as tsv_old, \
		open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
			tsv_new = csv.writer(tsv_new, delimiter='\t')
			tsv_old = csv.reader(tsv_old, delimiter='\t')		
			row = next(tsv_old) # the purpose of this line is to skip the header in the csv file, so there is no need to iterate another 10K plus times through blast2go_table_walnut/walnut_interpro
			combined_row = row + ["walnut 5", "walnut 11", "walnut 12", "blast2go_process","blast2go_function","blast2go_component"]
			tsv_new.writerow(combined_row) # copy the header
			for row in tsv_old:
				id = combined_anno_id_parser(row[0])
				
				walnut_results_interpro = ["N/A","N/A","N/A"]
					
				walnut_results_blast2go = blast2go_hastable.get(id)
				
				if not walnut_results_blast2go: # if there is no match from blast2go use as filler
					walnut_results_blast2go = "N/A"
					
				combined_row = row + walnut_results_interpro + walnut_results_blast2go
				tsv_new.writerow(combined_row)	
	elif params == "incorrect_parameters":
		sys.exit(-1)
		

	print (str(time.clock() - start_time) + " seconds")
	print("complete -- annotation file now available")