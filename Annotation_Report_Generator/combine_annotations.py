'''
Author: Sam Ginzburg

Description: An implementation of a program that combines blast2go, interpro, and previous annotation files

@param 	at least one of the following flags are required to run the program:  --blast2go --interpro --o , (--o : the output file name is optional, but if included MUST be present)
		
'''

import time
import sys
import os
import csv


# Returns the number of components, functions, and processes (in that order)
def get_go_counts():
	global num_component
	global num_process
	global num_function

	return [num_component, num_process, num_function]

def combined_anno_id_parser(id):
    #print (id)
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
	global num_component
	global num_process
	global num_function

	process = list()
	function = list()
	component = list()
	return_string = ""

	#print (row)
	#print (len(row))
	if len(row) > 13:
		if ";" in row[13]:
			split_string = row[13].split(";") # splits the column by ;
		elif "," in row[13]:
			split_string = row[13].split(",") # splits the column by ;	
		else:
			split_string = row[13]
	else:
		split_string = ""

	#print("splitstring: " + str(split_string))
	for string in split_string:
		string = string.strip()
		if "P:" in string or "Biological Process:" in string:
			process.append(string)
			num_process += 1
		if "F:" in string or "Molecular Function:" in string:
			function.append(string)
			num_function += 1
		if "C:" in string or "Cellular Component:" in string:
			component.append(string)
			num_component += 1
	return [" ".join(process), " ".join(function), " ".join(component)]

def make_blast2go_walnut_combined_dict(file_name):
	blast2go_data = dict()
	with open(file_name,'r') as tsv_old:
		row = next(tsv_old) #skip the first header line
		for row in csv.reader(tsv_old, delimiter='\t'):
			#print (row)
			id = row[0] #oombined_anno_id_parser(row[0]) #get the key
			blast2go_data[id] = parse_blast2go_column(row)
	return blast2go_data

def make_walnut_interpro_dict(file_name):
	walnut_interpro_data = dict()
	value_list = list()
	with open(file_name,'r') as tsv_old:
		for row in csv.reader(tsv_old, delimiter='\t'):
			#print (row)
			#id = row[1][8:] #get the key
			id = str(row[0])
			#print (id)
			value_list.append(row[5])
			value_list.append(row[11])
			value_list.append(row[12])
			walnut_interpro_data[id] = value_list
			value_list = list()
	#print (walnut_interpro_data)
	return walnut_interpro_data


#legacy code from older version may/may not add back in later
def print_usage():
	print("To run this script use at least one of the following flags:")
	print("--blast2go --interpro")
	print("EX: py combine_annotations.py --blast2go --interpro annotation.tsv blast2go.txt interpro.raw output_file_name.tsv")
	print("EX: py combine_annotations.py --interpro annotation.tsv interpro.raw")
	print("EX: py combine_annotations.py --blast2go annotation.tsv blast2go.txt output_file_name.tsv")
	sys.exit(-1)
	

#helper method for parse_flags
def parse_input_params(param_list):
	input = []
	previous_param = ""
	for param in param_list:
		if param[:2] == "--" and param != "--input":
			if param == "--output":
				return input
			else:
				previous_param = param
				input.append(param)
		elif previous_param == "--blast2go" or previous_param == "--interpro":
			input.append(param)
	return input
	
#helper method for parse_flags
def parse_output_param(param_list):
	output = []
	prev_param = ""
	for param in param_list:
		if param == "--output":
			output.append(param)
			prev_param = param
		elif prev_param == "--output":
			output.append(param)
			return output
	return ""
	
#helper method for parse_flags
def determine_blast_or_interpro_input(input):
	for param in input:
		if param == "--interpro":
			return "interpro"
		if param == "--blast2go":
			return "blast2go"
	return ""
	

def parse_flags(param_list):
	#parse input paramters first
	input = parse_input_params(param_list)
	#check if --output flag is present
	output = parse_output_param(param_list)
	
	#print(input)
	#print (output)	
	if len(output) != 2: #this means use default output
		if len(input) == 4: #all params
			return "all_params"
		elif len(input) == 2: #either 1 or the other, so do a quick check
			if determine_blast_or_interpro_input(param_list) == "blast2go":
				return "blast2go"
			else:
				return "interpro"
	else: #custom output
		if len(input) == 4: #all params
			return "all_params_custom_output"
		elif len(input) == 2: #either 1 or the other, so do a quick check
			if determine_blast_or_interpro_input(param_list) == "blast2go":
				return "blast2go_custom_output"
			elif determine_blast_or_interpro_input(param_list) == "interpro":
				return "interpro_custom_output"	
	

def get_num_sequences_identification():
	global count_sequences_identification
	return count_sequences_identification

	






def main(args):
	print ("appending interpro information now...")
	global count_sequences_identification
	
	global num_component
	global num_process
	global num_function

	num_component = 0
	num_process = 0
	num_function = 0
	
	count_sequences_identification = 0
	
	start_time = time.clock()
	arguments_list = args
	print (arguments_list)
	params = parse_flags(arguments_list)
	print (params)
	if params == "all_params" or params == "all_params_custom_output":
		if params == "all_params_custom_output":
			output = arguments_list[len(arguments_list) - 1]
			file_path_combined_anno = arguments_list[2]
			file_path_blast2go_walnut = arguments_list[4]
			file_path_walnut_interpro = arguments_list[6]
			
			if output == file_path_walnut_interpro or output == file_path_combined_anno or output == file_path_blast2go_walnut:
				print ("WARNING user entered option for custom file name and entered name of existing file (would cause overwrite) -- aborting")
				#print_usage()
				sys.exit(-1)
		else:
			output = "combine_annotations.tsv"
			file_path_combined_anno = arguments_list[2]
			file_path_blast2go_walnut = arguments_list[4]
			file_path_walnut_interpro = arguments_list[6]
		

		
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
				#id = combined_anno_id_parser(row[0])
				id = combined_anno_id_parser(row[1])
				walnut_results_interpro = walnut_interpro_hashtable.get(id)
				
				if not walnut_results_interpro: #this is to create the blank spaces if there are no interpro results for corresponding IDs
					walnut_results_interpro = ["N/A","N/A","N/A"]
					
				walnut_results_blast2go = blast2go_hastable.get(id)
				
				if not walnut_results_blast2go: # if there is no match from blast2go use as filler
					walnut_results_blast2go = ["N/A","N/A","N/A"]
					
				combined_row = row + walnut_results_interpro + walnut_results_blast2go
				tsv_new.writerow(combined_row)	
	
	elif params == "interpro" or params == "interpro_custom_output":
		#print("interpro custom output")
		if params == "interpro_custom_output":
			output = arguments_list[len(arguments_list) - 1]
			#print (output)
			file_path_combined_anno = arguments_list[1]
			#print (file_path_combined_anno)
			file_path_walnut_interpro = arguments_list[3]# changed to 3 from 4
			#print (file_path_walnut_interpro)
			if output == file_path_walnut_interpro or output == file_path_combined_anno:
				print ("WARNING user entered option for custom file name and entered name of existing file (would cause overwrite) -- aborting")
				#print_usage()
				sys.exit(-1)
			
		else:
			output = "combine_annotations.tsv"
			file_path_combined_anno = arguments_list[1]
			file_path_walnut_interpro = arguments_list[4]
		
		file_path_blast2go_walnut = "None"
		#first make dicts for walnut interpro and blast2go for easier searching based off of IDs
		walnut_interpro_hashtable = make_walnut_interpro_dict(file_path_walnut_interpro)
		blast2go_hastable = make_blast2go_walnut_combined_dict(file_path_walnut_interpro)
		with open(file_path_combined_anno,'r') as tsv_old, \
		open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
			#print(os.path.realpath(__file__) + "//" + output)
			tsv_new = csv.writer(tsv_new, delimiter='\t')
			tsv_old = csv.reader(tsv_old, delimiter='\t')		
			row = next(tsv_old) # the purpose of this line is to skip the header in the csv file, so there is no need to iterate another 10K plus times through blast2go_table_walnut/walnut_interpro
			combined_row = row + ["walnut 5", "walnut 11", "walnut 12", "blast2go_process","blast2go_function","blast2go_component"]
			tsv_new.writerow(combined_row) # copy the header
			for row in tsv_old:
				#print (row)
				#id = combined_anno_id_parser(row[1])
				id = str(row[0])
				walnut_results_interpro = walnut_interpro_hashtable.get(id)
				#print (walnut_results_interpro)	
				if not walnut_results_interpro: #this is to create the blank spaces if there are no interpro results for corresponding IDs
					walnut_results_interpro = ["N/A","N/A","N/A"]
				else:
					count_sequences_identification += 1
				
				walnut_results_blast2go = blast2go_hastable.get(id)

				if not walnut_results_blast2go:
					walnut_results_blast2go = ["N/A","N/A","N/A"]
					

				combined_row = row + walnut_results_interpro + walnut_results_blast2go
				tsv_new.writerow(combined_row)		
		
	elif params == "blast2go" or params == "blast2go_custom_output":
		
		if params == "blast2go_custom_output":
			output = arguments_list[len(arguments_list) - 1]
			file_path_combined_anno = arguments_list[2]
			file_path_blast2go_walnut = arguments_list[4]
			if output == file_path_blast2go_walnut or output == file_path_combined_anno:
				print ("WARNING user entered option for custom file name and entered name of existing file (would cause overwrite) -- aborting")
				#print_usage()
				sys.exit(-1)
		else:
			output = "combine_annotations.tsv"
			file_path_combined_anno = arguments_list[2]
			file_path_blast2go_walnut = arguments_list[4]
		
	
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
					walnut_results_blast2go = ["N/A"]
					
				combined_row = row + walnut_results_interpro + walnut_results_blast2go
				tsv_new.writerow(combined_row)	

		
	
	print (str(time.clock() - start_time) + " seconds")
	print("complete -- annotation file now available")
