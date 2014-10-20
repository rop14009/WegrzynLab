import datetime
import time
import os
import sys
import csv
import string

def parse_config_line(line):
	count = 0
	return_string = ""
	for char in line:
		if count >= 1 and char != "\n":
			return_string += char
		if char == ":":
			count += 1
	return return_string
		

'''
The following shows the index in list of settings
0)Path to query FASTA: something
1) Source Databases: 2
2) Sequence Search Application: d
3) Query Organism:sf
4) Path to formatted database1:asdf
5) Path to FASTA-version of database1:asdf
6) Path to search results from database1:asdf
7) Path to formatted database2 (organism-specific):asdf
8) Path to FASTA-version of database2:asdf
9) Path to search results from database2:asdf
10) Path to formatted database3 (full length): asdf
11) Path to FASTA-version of database3:asdf
12) Path to search results from database3:asdf
13) Full-length coverage requirement: asd
14) Minimum Evalue:asdf
15) Generate XML for Blast2GO:ASD
	'''
def parse_config_file(file_path):
	count = 0 # to keep track of which option is currently being parsed
	config_file = open(file_path, "r" )
	config_file_settings = []
	for line in config_file:
		config_file_settings.append(parse_config_line(line)) 
		count += 1
	print(config_file_settings)
	config_file.close()
	return config_file_settings
	
def get_gi_num_from_string(line):
	line = line[3:] # we know the first 3 characters are "gi|" so they can be removed
	return_string = ""
	for char in line:
		if char != "|":
			return_string += char
		else:
			return return_string
	return return_string
	
#helper method for multi fasta parse, just gets up until first space
def get_gi_string(line):
	temp_string = ""
	for char in line:
		if char == " ":
			#print (temp_string[1:])
			return temp_string[1:] # first char cut off since we know it to be ">"
		temp_string += char
	return temp_string[1:]		
	
# this method creates hashtable used for looking up fasta seq based off of gi
def multi_fasta_parse(file_name):
	fasta_db = dict()
	fasta_db_description = dict()
	fasta_db_species = dict()
	current_gi = ""
	current_protein_seq = ""
	current_desc = ""
	current_species = ""
	with open(file_name,'r') as file:
		for line in file:
			if line[:1] == ">": # line contains description of sequence
				if current_protein_seq != "":
					#print (current_gi)
					fasta_db[current_gi] = current_protein_seq
					fasta_db_description[current_gi] = current_desc
					fasta_db_species[current_gi] = current_species
					current_protein_seq = ""
					current_desc = ""
					current_species = ""
				
				current_gi =  get_gi_num_from_string(get_gi_string(line))
				
				current_desc = line[line.find(" "):].strip()
				current_species = line[line.find("[")+1:line.find("]")]
				
			else:
				current_protein_seq += line
	return [fasta_db, fasta_db_description, fasta_db_species]


'''
NCBI format
Field 1: query label
Field 2: target label
Field 3: percent identity
Field 4: alignment mismatches
Field 5: number of mismatches
Field 6: number of gap opens
Field 7: 1-based position of start in query
Field 8: 1-based position of end in query
Field 9: 1-based position of start in target
Field 10: 1-based position of end in target
Field 11: e value
Field 12: bit score
'''	
def ncbi_format_db_parse(file_name):
	ncbi_db = dict()
	global num_queries
	with open(file_name, "r") as file:
		file_tsv = csv.reader(file, delimiter='\t')
		for line in file_tsv:
			ncbi_db[ get_gi_num_from_string(line[1])] = line
			num_queries += 1
	return ncbi_db
	
	
'''
usearch_db format
find later
'''	
def usearch_format_db_parse(file_name):
	usearch_db = dict()
	global num_queries
	with open(file_name, "r") as file:
		file_tsv = csv.reader(file, delimiter='\t')
		for line in file_tsv:
			#print(get_gi_num_from_string(line[1]))
			num_queries += 1
			usearch_db[get_gi_num_from_string(line[1])] = line
	return usearch_db
	

#this function loads the list of 'uninformative' filter hits into memory for searching later on
def load_filter_list(file_name):
	return_list = []
	with open(file_name, "r") as file:
		for line in file:
			return_list.append(line[:-1]) #the last part of string is cut off since it is always the "/n" regex (due to formatting of the filter list)
	return return_list

	
def is_uninformative(fasta_db_element):
	for line in filter_list:
		if line in fasta_db_element:
			return True
	return False

def write_log(element, log_name):
	if number_db == 1:
		if not os.path.exists(log_name+".txt"): # if file doesnt exist create it
			print ("creating new logfile with name: " + log_name)
			file = open(log_name+".txt", "w")
			file.write(element+"\n")
			file.close()
		else:
			file = open(log_name+".txt", "a")
			file.write(element+"\n")
			file.close()
	else:
		counter -= 1
		if not os.path.exists(log_name+".txt"): # if file doesnt exist create it
			print ("creating log file with name: " + log_name)
			file = open(log_name+".txt", "w")
			if not counter > 0:
				file.write(element)
			else:
				file.write(element+"\n")
			file.close()
		else:
			file = open(log_name+".txt", "a")
			if not counter > 0:
				file.write(element)
			else:
				file.write(element+"\n")
			file.close()


	
def build_contaminants_db():
	contaminant_db = dict()
	bacteria_db = open("bacteria_db.txt", "r")
	fungi_db = open("fungi_db.txt","r")
	
	#bacteria first
	file_tsv = csv.reader(bacteria_db, delimiter='\t')
	for line in file_tsv:
		contaminant_db[line[0]] = line[0] #its faster to check if value exists in a hashtable than a regular list
		
	for line in fungi_db:
		contaminant_db[line] = line
	
	return contaminant_db
	
def write_contaminants_log(element,log_name):
	if number_db == 1:
		if not os.path.exists(log_name+".txt"): # if file doesnt exist create it
			file = open(log_name+".txt", "w")
			print ("creating contaminants log with name: " + log_name)
			file.write(element+"\n")
			file.close()
		else:
			file = open(log_name+".txt", "a")
			file.write(element+"\n")
			file.close()
	else:
		counter -= 1
		if not os.path.exists(log_name+".txt"): # if file doesnt exist create it
			file = open(log_name+".txt", "w")
			print ("creating contaminants log with name: " + log_name)
			if not counter > 0:
				file.write(element)
			else:
				file.write(element+"\n")
			file.close()
		else:
			file = open(log_name+".txt", "a")
			if not counter > 0:
				file.write(element)
			else:
				file.write(element+"\n")
			file.close()
	
	
def match_fasta(fasta, fasta_desc, fasta_species, db):
        global num_queries_informative_hit #1 or more informative hits
        global num_queries_no_hit #no hits
        global num_queries_uninformative

	for element in fasta:
		#check for contaminants
		if fasta_species[element] in contaminants:
			print("contaminants found: " + fasta_species[element])
			
			#append multiple contaminants per query
			if contaminants_found[element] is None and not db.get(element) is None:
				contaminants_found[element] = db.get(element)
			elif not contaminants_found[element] is None and not db.get(element) is None:
				contaminants_found[element] = contaminants_found[element] + db.get(element)
			else:
				print("this should not be reached")
			del db[element]
		print(element)
		print(db.get(element))			
		#make sure there is a match between fasta element and db
		if not db.get(element) is None:
			#returns true if "uninformative"
			#then if uninformative changes desc to uninformative
			if is_uninformative(fasta_desc[element]):
				fasta_desc[element] = "uninformative" + fasta_desc[element][fasta_desc[element].find("[")-1:]
				num_queries_uninformative += 1
			else:
				num_queries_informative_hit += 1

			get = db.get(element) + [fasta_desc[element][:fasta_desc[element].find("[")]] + [fasta_species[element]]
			if type(get) is list:
				tsv_new.writerow(get)
			else:
				tsv_new.writerow([get])
			#after match is found remove it from db
			del db[element]
		else:
			print ("no such element in db: " + element)
			num_queries_no_hit += 1
	

			
#entry point of script				
if __name__ == '__main__':
	start_time = time.clock()
	arguments_list = sys.argv
	settings = parse_config_file("configuration_file.txt") #sets up the settings
	output_log = "log"
	global filter_list
	global contaminants
	global number_db
	global counter
	global contaminants_found
	
	#The following global variables are used to record statistical data in order to generate a log
	
	#query length variables
	global longest_query_length
	global shortest_query_length
	global median_query_length
	
	longest_query_length = -2147483648
	shortest_query_length = 2147483647
	median_query_length = 0

	global total_query_sequences
	global avg_length_query_sequences
	
	total_query_sequences = 0
	avg_length_query_sequences = 0

	global num_queries # the number of queries
	global num_queries_informative_hit #1 or more informative hits
	global num_queries_no_hit #no hits
	global num_queries_uninformative
	
	num_queries = 0
	num_queries_informative_hit = 0
	num_queries_no_hit = 0
	num_queries_uninformative = 0
	print (num_queries)
	global db_type
	db_type = settings[2]
	filter_list = load_filter_list("filter_list.txt")
	contaminants = build_contaminants_db()
	contaminants_found = dict()
	print("contaminants db built")
	output = "default_output_annotation.tsv"
	number_db = int(settings[1]) # the number of databases being parsed
	counter = number_db
	db = dict()
	db2 = dict()
	db3 = dict()
	with open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
		tsv_new = csv.writer(tsv_new, delimiter='\t')
		#The first row
				
		row = ["Query","Subject_id","Identity(%)","Alignment_length","Mismatches","Number of gap opens","Query_start","Query_end" \
		,"Subject_start","Subject_end","E-value","Bit_score","Subject Description","Species"]
		
		for x in range(0, (number_db - 1)):
			row += ["Subject_description","Species","Subject_id","Alignment_length", \
			"Mismatches","Query_start","Query_end","Subject_start","Subject_end","E-value","Bit_score","Subject_description","Species"]
		row +=["has GFF3 data","GFF3 start","GFF3 stop","GFF3 ORF"]
		
		tsv_new.writerow(row)
		
		
		print(range (0, (number_db)))
		
		if number_db == 1: #database/fasta pair 1
			print("1 database/fasta pair")
			[fasta_db, fasta_db_description, fasta_db_species] = multi_fasta_parse(settings[6])
			
			db_time = time.clock()
			if settings[2] == "ncbi":
				db = ncbi_format_db_parse(settings[5])
			else:
				db = usearch_format_db_parse(settings[5])
			print (str(time.clock() - db_time) + " seconds")	
			match_fasta(fasta_db, fasta_db_description, fasta_db_species, db)
				
			#after parsing of all fasta elements add all missed hits to nohits file
			for key in db:
				write_log(db.get(key)[0],"nohits")
			
			for key in contaminants_found:
				print (contaminants_found.get(key))
				write_contaminants_log(contaminants_found.get(key),"contaminants")
			
		elif number_db == 2: #database/fasta pair 2
			print("2 database/fasta pair")
				
		elif number_db == 3: #database/fasta pair 3
			print("3 database/fasta pair")
				
		
	
	print (str(time.clock() - start_time) + " seconds")
	print("complete -- annotation file now available")


	if not os.path.exists(output_log+".txt"):
		print("creating log file: " + output_log)
		log_file = open(output_log + ".txt", "w")
		file.write("num_queries: " + num_queries + "\n")
		file.write("num_queries_informative_hit: " + num_queries_informative_hit + "\n")
		file.write("num_queries_no_hit: " + num_queries_no_hit + "\n")
                file.write("num_queries_uninformative: " + num_queries_uninformative + "\n")
		file.close()
		print ("log files complete -- exiting")
