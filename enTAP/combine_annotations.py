'''
Author: Sam Ginzburg

Description: An implementation of a program that combines blast2go, interpro, and previous annotation files

@param 	at least one of the following flags are required to run the program:  --blast2go --interpro --output , (--output : the output file name is optional, but if included MUST be present)
		
'''

import time
import sys
import os
import csv
import argparse


def get_go_counts():
	global num_component
	global num_process
	global num_function

	return [num_component, num_process, num_function]




# Returns the number of components, functions, and processes (in that order)
def get_go_interpro_counts():
	global num_component_interpro
	global num_process_interpro
	global num_function_interpro

	return [num_component_interpro, num_process_interpro, num_function_interpro]

def combined_anno_id_parser(id):
    #print (id)
    current_char = ""
    return_string = ""
    delimiter = "|"
    length = 0
    while current_char != delimiter and current_char != "":
        return_string += current_char
        current_char = id[length:length+1]
        length+=1
    return return_string

def parse_column(row, file_type):
	global num_component_interpro
	global num_process_interpro
	global num_function_interpro


	global num_component
	global num_process
	global num_function

	#global at_least_1_c
	#global at_least_1_p
	#global at_least_1_f
	
	process = list()
	function = list()
	component = list()
	
	return_string = ""

	
	if file_type == "blast":
		row_num = 7	
	else:
		row_num = 13

	if len(row) > row_num:
		if ";" in row[row_num]:
			split_string = row[row_num].split(";") # splits the column by ;
		elif "," in row[row_num]:
			split_string = row[row_num].split(",") # splits the column by ;	
		elif " " in row[row_num]:
			split_string = row[row_num].split(" ") # splits the column by " "
		else:
			split_string = row[row_num]
	else:
		split_string = ""

	#print("splitstring: " + str(split_string))

	for string in split_string:
		string = string.strip()
		if "P:" in string or "Biological Process:" in string:
			process.append(string)
				
			if file_type == "blast":
				num_process += 1
			else:
				num_process_interpro += 1
			
		if "F:" in string or "Molecular Function:" in string:
			function.append(string)
			
			if file_type == "blast":
				num_function += 1
			else:
				num_function_interpro += 1
			
		if "C:" in string or "Cellular Component:" in string:
			component.append(string)
			
			if file_type == "blast":
				num_component += 1
			else:
				num_component_interpro += 1
		
	#print ([" ".join(process), " ".join(function), " ".join(component)])
	return [" ".join(process), " ".join(function), " ".join(component)]

def make_blast2go_walnut_combined_dict(file_name):
	blast2go_data = dict()
	with open(file_name,'r') as tsv_old:
		row = next(tsv_old) #skip the first header line
		for row in csv.reader(tsv_old, delimiter='\t'):
			#print (row)
			id = row[0].split("|")[0] #combined_anno_id_parser(row[0]) #get the key
			#print (id)
			blast2go_data[id] = parse_column(row,"blast")
			#print (blast2go_data[id])

	
	print (len(blast2go_data))
	return blast2go_data

def make_walnut_interpro_dict(file_name):
	walnut_interpro_data = dict()
	value_list = list()
	with open(file_name,'r') as tsv_old:
		csv_reader = csv.reader(tsv_old, delimiter='\t')
		csv_reader.next()
		for row in csv_reader:
			#print (row)
			#print (len(row))
			#id = row[1][8:] #get the key
			id = str(row[0])
			
			if "." in id:
				id = id[8:]
			id = id.split("|")[0]


			parse_column(row, "interpro") # this is done to calc interpro GO terms
			#print (id)
			value_list.append(row[5])
			value_list.append(row[11])
			value_list.append(row[12])
			if len(row) > 13:
				#print ("interpro go terms being added...")
				value_list.append(row[13]) #interpro go terms

			if not walnut_interpro_data.get(id) is None and not value_list[1] in walnut_interpro_data.get(id)[1]:
				walnut_interpro_data.get(id)[0] = value_list[0] + ", " + walnut_interpro_data.get(id)[0]
				walnut_interpro_data.get(id)[1] = value_list[1] + ", " + walnut_interpro_data.get(id)[1]
				walnut_interpro_data.get(id)[2] = value_list[2] + ", " + walnut_interpro_data.get(id)[2]
				#print (walnut_interpro_data.get(id))
				if len(value_list) > 3 and len(walnut_interpro_data.get(id)) > 3:
					walnut_interpro_data.get(id)[3] = value_list[3] + ", " + walnut_interpro_data.get(id)[3]
				elif len(value_list) > 3:
					walnut_interpro_data.get(id).append(value_list[3])	
			else:
				walnut_interpro_data[id] = value_list
			
			value_list = list()
	print (len(walnut_interpro_data))
	return walnut_interpro_data


#legacy code from older version may/may not add back in later
def print_usage():
	print("To run this script use at least one of the following flags:")
	print("--blast2go --interpro")
	print("EX: python combine_annotations.py --blast2go --interpro annotation.tsv blast2go.txt interpro.raw output_file_name.tsv")
	print("EX: python combine_annotations.py --interpro annotation.tsv interpro.raw")
	print("EX: python combine_annotations.py --blast2go annotation.tsv blast2go.txt output_file_name.tsv")
	sys.exit(-1)
	

#helper method for parse_flags
def parse_input_params(param_list):
	global interpro_path
	global blast2go_path
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
			if previous_param == "--blast2go":
				blast2go_path = param
			elif previous_param == "--interpro":
				interpro_path = param
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
	global append_log
	#parse input paramters first
	#nput = parse_input_params(param_list)
	#check if --output flag is present
	#output = parse_output_param(param_list)


	# check for log to append to
	
	if "--log" in param_list:
		loc = param_list.index("--log")
		append_log = param_list[loc+1]
		param_list.remove("--log")
		param_list.remove(param_list[loc+1])
	else:
		append_log = None
			

	input = parse_input_params(param_list)
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


def get_at_least_1():
	global at_least_1_c
	global at_least_1_p
	global at_least_1_f
	return [at_least_1_c, at_least_1_p, at_least_1_f]
	




def get_at_least_1_interpro():
        global at_least_1_c_interpro
        global at_least_1_p_interpro
        global at_least_1_f_interpro
        return [at_least_1_c_interpro, at_least_1_p_interpro, at_least_1_f_interpro]



def check_cpf(l, t):
	global at_least_1_c
	global at_least_1_p
	global at_least_1_f
	global at_least_1_c_interpro
	global at_least_1_p_interpro
	global at_least_1_f_interpro
	

	global num_component_interpro
	global num_process_interpro
	global num_function_interpro


	global num_component
	global num_process
	global num_function

	c = False
	p = False
	f = False
	print (type(l))
	
	if t == "interpro":	
		print ("interpro GO:")
		print (l)
	print (l)
	exit()
	for element in l:
		if "C:" in element or "Cellular Component:" in element:		
			c = True
		
		if "P:" in element or "Biological Process:" in element:
			p = True
			
		if "F:" in element or "Molecular Function:" in element:
			f = True

	if t == "blast2go":
		if c:
			at_least_1_c += 1
		if p:
			at_least_1_p += 1
		if f:
			at_least_1_f += 1
	if t == "interpro":
		if c:
			at_least_1_c_interpro += 1
		if p:
			at_least_1_p_interpro += 1
		if f:
			at_least_1_f_interpro += 1





def main(args):
	print ("appending interpro information now...")
	global count_sequences_identification
	
	global num_component_interpro
	global num_process_interpro
	global num_function_interpro
	
	global num_component
	global num_process
	global num_function
	
	
	global at_least_1_c
	global at_least_1_p
	global at_least_1_f
	global at_least_1_c_interpro
	global at_least_1_p_interpro
	global at_least_1_f_interpro

	global append_log

	global interpro_path
	global blast2go_path

	interpro_path = ""
	blast2go_path = ""

	at_least_1_c = 0
	at_least_1_p = 0
	at_least_1_f = 0


	at_least_1_c_interpro = 0
	at_least_1_p_interpro = 0
	at_least_1_f_interpro = 0


	num_component_interpro = 0
	num_process_interpro = 0
	num_function_interpro = 0
	
	num_component = 0
	num_process = 0
	num_function = 0



	count_sequences_identification = [0, 0]
	
	start_time = time.clock()
	

	parser = argparse.ArgumentParser()
	parser.add_argument('--log', nargs='?')
	parser.add_argument('--input', nargs=1)
	parser.add_argument('--interpro', nargs='*')
	parser.add_argument('--blast2go', nargs='?')
	parser.add_argument('--output', nargs='?')

	args = parser.parse_args()
	print (args)
	#print (args.log)
	#exit()

	#arguments_list = args
	#print (arguments_list)
	#params = parse_flags(arguments_list)
	#print (params)

	params = ""
	if args.blast2go and args.interpro:
		if args.output:
			params = "all_params_custom_output"
		else:
			params = "all_params"	
	elif args.blast2go:
		if args.output:
			params = "blast2go_custom_output"
		else:
			params = "blast2go"
	elif args.interpro:
		if args.output:
			params = "interpro_custom_output"
		else:
			params = "interpro"

	#print (params)
	if params == "all_params" or params == "all_params_custom_output":
		#first make dicts for walnut interpro and blast2go for easier searching based off of IDs
		if len(args.interpro) == 1:
			walnut_interpro_hashtable = make_walnut_interpro_dict(str(args.interpro[0]))
		else:
			walnut_interpro_hashtable = make_walnut_interpro_dict(str(args.interpro[0]))
			walnut_interpro_hashtable.update(make_walnut_interpro_dict(str(args.interpro[1])))

		blast2go_hastable = make_blast2go_walnut_combined_dict(str(args.blast2go))
		with open(args.input[0] ,'r') as tsv_old, \
		open(os.path.dirname(os.path.realpath(__file__)) + "//" + args.output, 'w') as tsv_new:
			tsv_new = csv.writer(tsv_new, delimiter='\t')
			tsv_old = csv.reader(tsv_old, delimiter='\t')		
			row = next(tsv_old) # the purpose of this line is to skip the header in the csv file, so there is no need to iterate another 10K plus times through blast2go_table_walnut/walnut_interpro
			combined_row = row + ["Signature Description", "InterPro accession number", "InterPro description", "blast2go_process","blast2go_function","blast2go_component"]
			tsv_new.writerow(combined_row) # copy the header
			for row in tsv_old:
				id = str(row[0]).split("|")[0]
				#id = combined_anno_id_parser(row[1])
				#print (id)
				walnut_results_interpro = walnut_interpro_hashtable.get(id)
				#print (walnut_results_interpro)	
				if not walnut_results_interpro: #this is to create the blank spaces if there are no interpro results for corresponding IDs
					walnut_results_interpro = ["N/A","N/A","N/A"]
					#count_sequences_identification[1] += 1
				else:
					if len(walnut_results_interpro) > 3:
						check_cpf(walnut_results_interpro[3],"interpro")	
				walnut_results_blast2go = blast2go_hastable.get(id)
				
				if not walnut_results_blast2go: # if there is no match from blast2go use as filler
					walnut_results_blast2go = ["N/A","N/A","N/A"]
					count_sequences_identification[1] += 1
				else:
					count_sequences_identification[0] += 1	
					check_cpf(walnut_results_blast2go,"blast2go")
				combined_row = row + walnut_results_interpro + walnut_results_blast2go
				tsv_new.writerow(combined_row)	
	
	elif params == "interpro" or params == "interpro_custom_output":
		
		if len(args.interpro) == 1:
			#first make dicts for walnut interpro and blast2go for easier searching based off of IDs
			walnut_interpro_hashtable = make_walnut_interpro_dict(str(args.interpro[0]))
		else:
			walnut_interpro_hashtable = make_walnut_interpro_dict(str(args.interpro[0]))
			walnut_interpro_hashtable.update(make_walnut_interpro_dict(str(args.interpro[1])))

		#blast2go_hastable = make_blast2go_walnut_combined_dict(file_path_walnut_interpro)
		
		with open(args.input[0],'r') as tsv_old, \
		open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
			#print(os.path.realpath(__file__) + "//" + output)
			tsv_new = csv.writer(tsv_new, delimiter='\t')
			tsv_old = csv.reader(tsv_old, delimiter='\t')		
			row = next(tsv_old) # the purpose of this line is to skip the header in the csv file, so there is no need to iterate another 10K plus times through blast2go_table_walnut/walnut_interpro
			combined_row = row + ["Signature Description", "InterPro accession number", "InterPro description", "blast2go_process","blast2go_function","blast2go_component"]
			tsv_new.writerow(combined_row) # copy the header
			for row in tsv_old:
				#print (row)
				#id = combined_anno_id_parser(row[1])
				id = str(row[0]).split("|")[0]
				walnut_results_interpro = walnut_interpro_hashtable.get(id)
				#print (walnut_results_interpro)	
				if not walnut_results_interpro: #this is to create the blank spaces if there are no interpro results for corresponding IDs
					walnut_results_interpro = ["N/A","N/A","N/A"]
					count_sequences_identification[1] += 1
				else:
					count_sequences_identification[0] += 1
					if len(walnut_results_interpro) > 3:
						check_cpf(walnut_results_interpro[3],"interpro")	
				
				walnut_results_blast2go = None

				if not walnut_results_blast2go:
					walnut_results_blast2go = ["N/A","N/A","N/A"]
		
				combined_row = row + walnut_results_interpro + walnut_results_blast2go
				tsv_new.writerow(combined_row)		
		
	elif params == "blast2go" or params == "blast2go_custom_output":
		#first make dicts for walnut interpro and blast2go for easier searching based off of IDs
		blast2go_hastable = make_blast2go_walnut_combined_dict(str(args.blast2go))
		
		with open(args.input[0],'r') as tsv_old, \
		open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
			tsv_new = csv.writer(tsv_new, delimiter='\t')
			tsv_old = csv.reader(tsv_old, delimiter='\t')		
			row = next(tsv_old) # the purpose of this line is to skip the header in the csv file, so there is no need to iterate another 10K plus times through blast2go_table_walnut/walnut_interpro
			combined_row = row + ["Signature Description", "InterPro accession number", "InterPro description", "blast2go_process","blast2go_function","blast2go_component"]
			tsv_new.writerow(combined_row) # copy the header
			for row in tsv_old:
				#id = combined_anno_id_parser(row[0])
				id = row[0]
				id = str(row[0]).split("|")[0]
				
				walnut_results_interpro = ["N/A","N/A","N/A"]
				"""
				if not walnut_results_interpro:
					walnut_results_interpro = ["N/A","N/A","N/A"]
					count_sequences_identification[1] += 1
				else:
					count_sequences_identification[0] += 1
					if len(walnut_results_interpro) > 3:
						check_cpf(walnut_results_interpro[3],"interpro")	
				"""

				walnut_results_blast2go = blast2go_hastable.get(id)
				#print (walnut_results_blast2go)	
				if not walnut_results_blast2go: # if there is no match from blast2go use as filler
					walnut_results_blast2go = ["N/A","N/A","N/A"]
				else:
					check_cpf(walnut_results_blast2go,"blast2go")
	
				combined_row = row + walnut_results_interpro + walnut_results_blast2go
				tsv_new.writerow(combined_row)	


	#print (args.log)
	if args.log:	
		domain_ids = get_num_sequences_identification()
		at_least_1 = get_at_least_1() # C, P, F
		go_counts = get_go_counts()
		go_interpro_counts = get_go_interpro_counts()
		at_least_1_interpro = get_at_least_1_interpro()
 
		print ("appending to:\t" + args.log)
		if True:#os.path.exists(append_log + ".txt"):
			with open(os.path.dirname(os.path.realpath(__file__)) + "//" + args.log, 'a') as log:
				tsv_log = csv.writer(log, delimiter='\t')
				if not interpro_path is None:
					tsv_log.writerow(["Interpro File: "] + [str(interpro_path)]) # TODO put an if statement here to check for multiple interpro files and
				tsv_log.writerow(["Blast2GO File: "] + [str(blast2go_path)])
				tsv_log.writerow(["Number of sequences with Domain Identification: "] + [str(domain_ids[0])])
				tsv_log.writerow(["Number of sequences without Domain Identification: "] + [str(domain_ids[1])])
				tsv_log.writerow(["Blast2GO Gene Ontology Stats"])
				tsv_log.writerow(["Number of Components: "] + [str(go_counts[0])])
				tsv_log.writerow(["Number of Functions: "] + [str(go_counts[1])])
				tsv_log.writerow(["Number of Processes: "] + [str(go_counts[2])])
				tsv_log.writerow(["Number of Transcripts with at least 1 Component: "] + [str(at_least_1[0])])
				tsv_log.writerow(["Number of Transcripts with at least 1 Function: "] + [str(at_least_1[2])])
				tsv_log.writerow(["Number of Transcripts with at least 1 Process: "] + [str(at_least_1[1])])
				if not interpro_path is None:
					tsv_log.writerow(["Interpro Gene Ontology Stats (Totals)"])
					tsv_log.writerow(["Component: "] + [str(go_interpro_counts[0])])
					tsv_log.writerow(["Function: "] + [str(go_interpro_counts[1])])
					tsv_log.writerow(["Process: "] + [str(go_interpro_counts[2])])
					tsv_log.writerow(["Number of Transcripts with at least 1 Component: "] + [str(at_least_1_interpro[0])])
					tsv_log.writerow(["Number of Transcripts with at least 1 Function: "] + [str(at_least_1_interpro[2])])
					tsv_log.writerow(["Number of Transcripts with at least 1 Process: "] + [str(at_least_1_interpro[1])])



	print (str(time.clock() - start_time) + " seconds")
	print("complete -- annotation file now available")


if __name__ == '__main__':
	args = sys.argv	
	args = args[1:]
	print (args)
	main(args)

