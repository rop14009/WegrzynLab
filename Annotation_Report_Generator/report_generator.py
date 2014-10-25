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
[12:34:19 PM] Jill Wegrzyn: Yes for the e-value to make this decision on when to accept the "next best" hit where one would increment through 
contaminants and no hits in the top five list - I would set a default e-value to 1e-5 in the configuration file but the user could modify this if they would like in such a file.
[12:34:54 PM] Jill Wegrzyn: Within that e-value set, I would search through the list to find an informative, non-contaminant where possible.  
When not available - the top hit should be used and marked accordingly.
[12:35:47 PM] Sam Ginzburg: What should be done if there are no e values less than 1e-5?
[12:36:33 PM] Sam Ginzburg: and if there are two elements left with identical e values and one is contaminated and the other uninformative which one is chosen as the better hit?
[12:41:23 PM] Jill Wegrzyn: I would select the contaminant one as it is a better known in that case
[12:41:37 PM] Jill Wegrzyn: if no evalue are less than 1e-5, then we select the top hit
	
Hi Sam, a couple of other notes on the annotation approach - for assigning a contaminated hit, we want to first see if there is another hit that is a non-contam 
in the top five or so that has an evalue within acceptable limits and select this preferentially (the same is true for a noninformative hit).  
Some of this logic was built in for non-informative hits in the prior scripts and I wanted to see what it looks like presently	
'''	


def find_best_query_result(query1, query2):
	#note: query1 is the pre-existing value in db
	#these variables are needed to determine informative/contaminant hits
	global fasta_db
	global fasta_db_description
	global fasta_db_species

	global e_value # 1e-5 / 0.00001 
	
	#attempt 1) find a query that is not a contaminant and informative, with lowest e-value
	
	if not fasta_db_species[query1[1]] in contaminants and not fasta_db_species[query2[1]] in contaminants:
		if not is_uninformative(fasta_db_description[query1[1]] and not is_uninformative(fasta_db_description[query2[1]]:
			if parse_e_value(query1) < parse_e_value(query2):
				return query1
			else:
				return query2
	
	#attempt 2) find a contaminated query with lowest possible e-value and return that
	if not is_uninformative(fasta_db_description[query1[1]] and not is_uninformative(fasta_db_description[query2[1]]:
		if parse_e_value(query1) < parse_e_value(query2):
			return query1
		else:
			return query2
			
	#attempt 3) as a final resort, simply return the query with lowest possible e-value (this result will be both contaminated and uninformative)
	
	if parse_e_value(query1) < parse_e_value(query2):
		return query1
	else:
		return query2
	
	
	print ("error in find_best_query_result(query1, query2) this should not ever be reached")
	return ["",""] #this should not be reached, ever
	
def parse_e_value(e_val):
	number = e_val
	exp = ""
	e_val = "".join(e_val)

	base = number[:e_val.find("e")]
	exp = number[e_val.find("e")+2:]
	try:
		return (10 ** float(exp)) * float(base)
	except ValueError:
		print ("e_val ::: " + e_val)	
		print ("exp :::  " + exp)
		print ("base ::: " + base)	
		return float(e_val)	
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
	temp_query_group= dict()
	current_query = ""
	global num_queries
	
	first_row = 1	

	with open(file_name, "r") as file:
		file_tsv = csv.reader(file, delimiter='\t')
		
		for line in file_tsv:
			if first_row:
				temp_query_group[get_gi_num_from_string(line[1])] =  line
				first_row += 1
		
			if current_query != line[0]:
				#add the best query from the previous group of queries
				[best_query_gi, best_query] = find_best_query_result(temp_query_group)
				ncbi_db[best_query_gi] = best_query
				temp_query_group = dict() # start a new group of matching queries
				num_queries += 1
				current_query = line[0]
			else:
				temp_query_group[get_gi_num_from_string(line[1])] =  line
	return ncbi_db
	
	
'''
usearch_db format
find later
'''	
def usearch_format_db_parse(file_name):
	usearch_db = dict()
	temp_query_group = dict()
	current_query = ""
	first_row = 1
	
	global num_queries
	
	with open(file_name, "r") as file:
		file_tsv = csv.reader(file, delimiter='\t')
		for line in file_tsv:
			if not line[1] in usearch_db:
				usearch_db[line[1]] = line
			else:
				usearch_db[line[1]] = find_best_query_result(usearch_db[line[1]], line)
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
	
	
def match_fasta(db):
	global num_queries_informative_hit #1 or more informative hits
	global num_queries_no_hit #no hits
	global num_queries_uninformative

	global fasta_db
	global fasta_db_description
	global fasta_db_species
	
	
	for element in fasta_db:
		#check for contaminants
		if fasta_db_species[element] in contaminants:
			print("contaminants found: " + fasta_db_species[element])
			
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
			if is_uninformative(fasta_db_description[element]):
				fasta_db_description[element] = "uninformative" + fasta_db_description[element][fasta_db_description[element].find("[")-1:]
				num_queries_uninformative += 1
			else:
				num_queries_informative_hit += 1

			get = db.get(element) + [fasta_db_description[element][:fasta_db_description[element].find("[")]] + [fasta_db_species[element]]
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
	
	global fasta_db
	global fasta_db_description
	global fasta_db_species
	
	global e_value # this constant refers to the threshold to use for determining a "best match" in terms of queries
	e_value = 0.00001 
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
			print ("multi fasta file done")
			db_time = time.clock()
			if settings[2] == "ncbi":
				db = ncbi_format_db_parse(settings[5])
			else:
				db = usearch_format_db_parse(settings[5])
			print (str(time.clock() - db_time) + " seconds")	
			match_fasta(db)
				
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
