import datetime
import time
import os
import sys
import csv
import string
import subprocess
import re

import combine_annotations
import sync_repair


def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False



# this function tests the configuration options to ensure that all the file paths are valid, and the format is correct
def config_sanity_check(config_list):
	print ("Performing sanity check on config file before continuing")
	#check query fasta file
	if not os.path.exists(config_list[0]):
		print ("The specified query fasta file does not exist -- aborting execution")
		exit()

	# ensure number of databases is checked and that it is a number
	if config_list[1] != "" and config_list[1].isdigit():
		num_db = int(config_list[1])	

		if num_db > 3:
			print("The number of databases was set to a number higher than 3, since 3 is the max number allowed currently, this value will be reset to 3")
		
		if num_db == 1: # ensure that the data for the first db is entered correctly
			if (not config_list[4].isdigit()) or not os.path.exists(config_list[5]) or not os.path.exists(config_list[6]):
				print ("Either the score for database 1 was entered incorrectly, or the file paths for database 1 are invalid -- aborting execution")
				exit()
		if num_db == 2: # ensure that the data for the first and second db is entered correctly	
			if (not config_list[4].isdigit()) or not os.path.exists(config_list[5]) or not os.path.exists(config_list[6]) \
			       or (not config_list[7].isdigit()) or not os.path.exists(config_list[8]) or not os.path.exists(config_list[9]):
					print ("Either the score for database 1/2 was entered incorrectly or the file paths for databases 1/2 are invalid -- aborting execution")
					exit()
		if num_db == 3: # ensure that the data for all three databases has been entered correctly
			if (not config_list[4].isdigit()) or not os.path.exists(config_list[5]) or not os.path.exists(config_list[6]) \
				or (not config_list[7].isdigit()) or not os.path.exists(config_list[8]) or not os.path.exists(config_list[9]) \
				or (not config_list[10].isdigit()) or not os.path.exists(config_list[11]) or not os.path.exists(config_list[12]):
					print ("Either the score for database 1/2/3 was entered incorrectly or the file paths for databases 1/2/3 are invalid -- aborting execution")
					exit()
		if num_db <= 0:
			print ("The number of source databases must be at least 1 -- aborting execution")
			exit()
	else:
		print ("The number of databases was either not entered, or entered incorrectly -- aborting execution")
		exit()


	# check to ensure minimum e-value is entereted and minimum coverage requirement entered as well
	if not isfloat(config_list[13])	or not isfloat(parse_e_value(config_list[13])):
		print ("Either the e-value or the minimum coverage requirement has been entered incorrectly --  aborting execution")
		exit()
	
	print ("Sanity check complete --  now continuing with execution")

def parse_config_line(line):
	count = 0
	return_string = ""
	for char in line:
		if count >= 1 and char != "\n":
			return_string += char
		if char == ":":
			count += 1
	return return_string.strip()
		

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
	if not os.path.exists(file_path):
		print ("The configuration file is not in the current directory -- please check to see if configuration_file.txt is present - execution ")
		exit()	


	count = 0 # to keep track of which option is currently being parsed
	config_file = open(file_path, "r" )
	config_file_settings = []
	for line in config_file:
		config_file_settings.append(parse_config_line(line)) 
		count += 1
	print(config_file_settings)
	config_file.close()
	config_sanity_check(config_file_settings)
	return config_file_settings
	
def get_gi_num_from_string(line):
	#global is_tair
	#if is_tair:
		#print (line)
	#	return line

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


#gets the index of the nth occurance of an element in a list
def get_nth_index(n, element, l):
	#print (l)
	#print (n)
	#print ((l.find(element)+1))
	if n == 0:
		return l
	elif l.find(element) != -1:
		n -= 1
		return get_nth_index(n, element , l[l.find(element)+1:])
	else:
		print ("Error: there is no occurance of " + str(element) + " in list: " + str(l))
		return -1 


def chunk_line(line):
	temp_list = re.findall('(.*?)\]',line)
	return_list = list()

	for item in temp_list:
		new_item = item + ']'
		return_list.append(new_item)

	return return_list


	
# this method creates hashtable used for looking up fasta seq based off of gi
def multi_fasta_parse(file_name):
	global is_tair


	if not os.path.exists(file_name):
		print ("The specified file path to a multi-fasta file:\t"+file_name+" does not exist, please recheck file paths and names -- aborting execution")
		exit()

	gi_list = list()
	
	fasta_db = dict()
	fasta_db_description = dict()
	fasta_db_species = dict()
	
	gi_results_list = list()
	species_results_list = list()
	description_results_list = list()
	
	current_gi = ""
	current_protein_seq = ""
	current_desc = ""
	current_species = ""
	
	with open(file_name,'r') as file:
		for line in file:
			if ">" in line: # line contains description of sequence
				if current_protein_seq != "":
					for gi,species,desc in zip(gi_results_list, species_results_list, description_results_list):
						fasta_db[gi] = current_protein_seq
						fasta_db_description[gi] = desc
						fasta_db_species[gi] = species				
				
					gi_results_list = list()
					species_results_list = list()
					description_results_list = list()

					current_protein_seq = ""
				
				gi_list = chunk_line(line)
				#print (gi_list)
				for item in gi_list:
					gi_results_list.extend(re.findall('gi\|(\d*)\|',item))
					species_results_list.extend(re.findall('\[(.*?)\]',item))
					description_results_list.extend(re.findall('\|\s(.*?)\s\[',item))
				#print (zip(gi_results_list, species_results_list, description_results_list))
			else:
				current_protein_seq += line.strip().replace("\n","")

	#print (last_gi)
	#print (current_gi)

	#fasta_db[current_gi] = current_protein_seq
	#fasta_db_description[current_gi] = current_desc
	#fasta_db_species[current_gi] = current_species
	
	for gi,species,desc in zip(gi_results_list, species_results_list, description_results_list):
		fasta_db[gi] = current_protein_seq
		fasta_db_description[gi] = desc
		fasta_db_species[gi] = species

	#print ("test:  " + str(fasta_db.get("")))
	#print (fasta_db)
	print ("Finished building fasta db:\t" + str(file_name))
	print ("Fasta db contains: " + str(len(fasta_db)) + " elements")
	return [fasta_db, fasta_db_description, fasta_db_species]



def uninformative_debug_log(element):
	global date
	if not os.path.exists("uninformative_" + date + ".txt"):
		with open(os.path.dirname(os.path.realpath(__file__)) + "//" + "uninformative_" + date + ".txt", 'a') as tsv_log:
			tsv_log = csv.writer(tsv_log, delimiter='\t')
			for key,value in element.items():
				tsv_log.writerow([str(key)])





'''
Input: query1, query2 (each query is a list)

Output: The "better" query, based upon e-value, query coverage, contaminant, uninformative-ness and null queries (which should not happen but are still accounted for)
'''	

def find_best_query_result(query1, query2):


	# None detection

	if query1 is None and not query2 is None:
		return query2

	if not query1 is None and query2 is None:
		return query1

	if query1 is None and query2 is None:
		return None


        
	#debug var
	global debug_uninformative_list


	global contaminants
	#note: query1 is the pre-existing value in db
        #these variables are needed to determine informative/contaminant hits
        global fasta_db
        global fasta_db_description
        global fasta_db_species
	global fasta_no_gi
	global contaminants_found

	global min_coverage # 0.7 --> refers to 70% coverage
        global e_value # 1e-5 / 0.00001

        #print (query1)
	#print (query2)
	query1_gi = get_gi_num_from_string(query1[1])
        query2_gi = get_gi_num_from_string(query2[1])
	

	#print (fasta_no_gi)
	
	#print (query1[0].index("|"))


	#print (str(query1[0][0:query1[0].index("|")]))	
	#print (str(query2[0][0:query2[0].index("|")]))
	
	#print (query1)

	#print (query2)

	query1_coverage = float(query1[3]) / float(len(fasta_no_gi[query1[0]]))
	query2_coverage = float(query2[3]) / float(len(fasta_no_gi[query2[0]]))

	#print ("query1 coverage: " + str(query1_coverage))
	#print ("query2 coverage: " + str(query2_coverage))

	if is_uninformative(fasta_db_description[query1_gi]):
		#if fasta_db_description[query2_gi] != "uninformative":
		#	debug_uninformative_list[fasta_db_description[query1_gi]] = fasta_db_description[query1_gi]
		
		fasta_db_description[query1_gi] = "uninformative"					
	if is_uninformative(fasta_db_description[query2_gi]):
		
		#if fasta_db_description[query2_gi] != "uninformative":
		#	debug_uninformative_list[fasta_db_description[query2_gi]] = fasta_db_description[query2_gi]
		
		fasta_db_description[query2_gi] = "uninformative"


	# because of naming a double contradiction occurs here,
	# not uninformative = informative hit

	
	#print (fasta_db_description[query1_gi])
	#print (fasta_db_description[query2_gi])
	
	if not fasta_db_description[query1_gi] is "uninformative" and fasta_db_description[query2_gi] is "uninformative":
		return query1

	if fasta_db_description[query1_gi] is "uninformative" and not fasta_db_description[query2_gi] is "uninformative":
		return query2




	species = fasta_db_species[query1_gi]
	species = species.split(" ")
	if len(species) > 1:
		species = species[0] + " " + species[1]
	else:
		species = species[0]
	fasta_db_species[query1_gi] = species

	species2 = fasta_db_species[query2_gi]
	species2 = species2.split(" ")
	if len(species2) > 1:
		species2 = species2[0] + " " + species2[1]
	else:
		species2 = species2[0]
	fasta_db_species[query2_gi] = species2
	'''
	if species in contaminants:
		print ("contaminant:   " + str(species))
	if species2 in contaminants:
		print ("contaminant:   " + str(species2))
	'''


	if not species in contaminants and species2 in contaminants and query1_coverage < min_coverage:
		return query1
		#contaminants_found[query1_gi] = fasta_db_species[query1_gi]
		#print ("contam found")
	if not species2 in contaminants and species in contaminants and query2_coverage < min_coverage:
		return query2
		#contaminants_found[query2_gi] = fasta_db_species[query2_gi]
		#print ("contam found")

	# new algorithm --  if both queries are over the min cov requirement, chose the one with smaller e-value


	# if one is uninformative and the other a contaminant or vice-versa, pick the uninformative hit over the conta


	# if this point is reached, then compare the query coverages and e values to determine the better hit
	
	if query1_coverage >= min_coverage and query2_coverage >= min_coverage:
		if parse_e_value(query1[10]) >= parse_e_value(query2[10]):
	                del query2
			return query1
		else:
			del query1
			return query2
	elif query1_coverage >= min_coverage and not query2_coverage >= min_coverage:
		del query2
		return query1
	elif query2_coverage >= min_coverage and not query1_coverage >= min_coverage:	
		return query2
	else: # the scenario where both are below the min coverage
		del query2
		return query1


	print ("This statement has been reached, when it should never be reached")
	exit()

def parse_e_value(e_val):
	number = e_val
	exp = ""
	e_val = "".join(e_val)

	base = number[:e_val.find("e")]
	exp = number[e_val.find("e")+1:]
	#print (exp)
	
	if e_val.find("e") == -1:
		return float(e_val)

	try:
		return (10 ** float(exp)) * float(base)
	except ValueError:
		#print ("e_val ::: " + e_val)	
		#print ("exp :::  " + exp)
		#print ("base ::: " + base)	
		return float(e_val)	



def trim_query_name(query_name):
	if "|" in query_name:
		return query_name[:query_name.index("|")]
	else:
		return query_name

'''
WARNING THIS SECTION IS CURRENTLY NOT FUNCTIONING CORRECTLY - DO NOT USE
NCBI format
Field 1: query label
Field 2: target label
Field 3: percent identity
Field 4: alignment mismatches
Field 5: number of mismatches
Field 6: number of gap opens
Field 7: 1-based position of start in queryH promoter plus cDNA', 'gi|547235765|ref|WP_021971974.1|', '35.7', '42', '24', '1', '1250', '1134', '364', '405', '31', '33.5']
creating contaminants log with name: contaminants
Field 8: 1-based position of end in query
Field 9: 1-based position of start in target
Field 10: 1-based position of end in target
Field 11: e value
Field 12: bit score
'''	
#this method is currently depreicated (DO NOT USE UNTIL UPDATED)
def ncbi_format_db_parse(file_name):
	ncbi_db = dict()
	temp_query_group= dict()
	current_query = ""
		
	first_row = 1	
        
	global longest_query_length
	global shortest_query_length
	global median_query_length
	global avg_length_query_sequences
	global num_queries



	with open(file_name, "r") as file:
		file_tsv = csv.reader(file, delimiter='\t')
		for line in file_tsv:
			#num_queries += 1	
			[key, in_db] = is_query_in_db(line[0], ncbi_db)

			if in_db:
				ncbi_db[key] = find_best_query_result(ncbi_db[key], line)
			else:
				ncbi_db[get_gi_num_from_string(line[1])] = line
		
		'''
		length = int(line[4])
		avg_length_query_sequence = float((avg_length_query_sequences * (num_queries-1) + length) / num_queries)
		median_query_length.append(length)
		if length > longest_query_length:
			longest_query_length = length
		if length < shortest_query_length:
			shortest_query_length = length	
		'''
		return ncbi_db
	


# This is a helper method to diamond format db parse to ensure that the best mach for each query is made
def is_query_present(query):
	global query_gi_association_db_0
	global query_gi_association_db_1
	global query_gi_association_db_2
	
	global db_count

	return get_current_db().get(query)

def get_current_db():
	global db_count
	global query_gi_association_db_0
	global query_gi_association_db_1
	global query_gi_association_db_2
	if db_count == 0:
		return query_gi_association_db_0
	elif db_count == 1:
		return query_gi_association_db_1
	elif db_count == 2:
		return query_gi_association_db_2
	else:
		print ("error -- unrecognized db count")
		print (db_count)



'''
diamond_db format
find later
'''	
def diamond_format_db_parse(file_name):
	diamond_db = dict()
	temp_query_group = dict()
	current_query = ""
	first_row = 1
	line_count = 0
	global is_tair

	global contaminants	
	global contaminants_found
	global longest_query_length
        global shortest_query_length
        global median_query_length
	global avg_length_query_sequences
	global num_queries
	global num_queries_uninformative
	global fasta_db_species
	global fasta_no_gi
	global db_count

	global query_gi_association_db_0
	global query_gi_association_db_1
	global query_gi_association_db_2

	if not os.path.exists(file_name):
		print ("The diamond results file:\t"+file_name+" does not exist, please recheck name and filepath -- aborting execution")
		exit()


	with open(file_name, "r") as file:
		file_tsv = csv.reader(file, delimiter='\t')
		print ("Now parsing diamond DB:\t" + str(file_name))
		for line in file_tsv:
			line_count += 1
			if not fasta_db.get(str(get_gi_num_from_string(line[1]))) is None:
				#num_queries += 1
				#print (get_gi_num_from_string(line[1]))
				#[key, in_db] = is_query_in_db(line[0], diamond_db)		
				#print (fasta_db_description)

				#print (line)
				
				line[0] = str(line[0].split(" ")[0])
				
				if not fasta_no_gi.get(line[0]) is None and not diamond_db.get(line[0]) is None:
					diamond_db[line[0]] = find_best_query_result(diamond_db[line[0]], line)
					if diamond_db[line[0]][len(diamond_db[line[0]]) - 1] != file_name:
						diamond_db[line[0]].append(file_name)
				elif not fasta_no_gi.get(line[0]) is None:
					diamond_db[line[0]] = line
					diamond_db[line[0]].append(file_name)
			elif not fasta_no_gi.get(line[0]) is None:
				print ("A mismatch between the file: " + str(file_name) + " and its corresponding fasta db has occurred\n")
				print ("The ID: " + str(get_gi_num_from_string(line[1])) + " is present within the DB and not the fasta file, indicating that the files may potentially be out of sync")
				print ("Attempting to continue by downloading transcript data from NCBI...")
				gi_data = sync_repair.sync_repair(get_gi_num_from_string(line[1]))
				
				fasta_db[get_gi_num_from_string(line[1])] = gi_data[0]
				fasta_db_species[get_gi_num_from_string(line[1])] = gi_data[2]
				fasta_db_description[get_gi_num_from_string(line[1])] = gi_data[1]
				print ("Success - Managed to download missing transcript from NCBI and add them to FASTA Database, continuing normally")

				line[0] = str(line[0].split(" ")[0])

				if not fasta_no_gi.get(line[0]) is None and not diamond_db.get(line[0]) is None:
					diamond_db[line[0]] = find_best_query_result(diamond_db[line[0]], line)
					if diamond_db[line[0]][len(diamond_db[line[0]]) - 1] != file_name:
						diamond_db[line[0]].append(file_name)
				elif not fasta_no_gi.get(line[0]) is None:
					diamond_db[line[0]] = line
					diamond_db[line[0]].append(file_name)



		#print (len(get_current_db()))
		#print (len(diamond_db))
		#print ("db complete parsing")
		#print (get_current_db())
		#print (query_gi_association_db_0)
		db_count += 1
		#print ("diamond db #" + str(db_count) + " contains: " + str(len(diamond_db)) + " unique best hits and " + str(line_count) + " total elements")
		return diamond_db








def is_query_in_db(query, db):
	key = ""
	#print (query)
	for key,value in db.items():
		if query in value[0]:
			#print ("true")
			return key, True
	return key, False



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
	
	f = open(log_name,"a")
	f.write(str(element))


#this is a helper method for building the contaminants DB
def parse_contaminant_name(name):
	name = name.split(" ")
	if len(name) > 1:
		name = name[0] + " " + name[1]
	else:
		name = name[0]
	return name

def build_contaminants_db():
	global settings
	contaminant_db = dict()

	bacteria_db = open(settings[25], "r")
	fungi_db = open(settings[24],"r")
	insects_db = open(settings[23],"r")	
	
	#bacteria first
	#file_tsv = csv.reader(bacteria_db, delimiter='\t')
	
	if ("y" in settings[21] or "yes" in settings[21]) and os.path.exists(settings[25]): 
		# loads bacteria DB
		for line in bacteria_db:
			line = line[:-2] # remove the \r\n special characters at end of line
			line = line.strip() # remove excess spaces
			#print (line)
			line = parse_contaminant_name(line)	
			contaminant_db[line] = line #its faster to check if value exists in a hashtable than a regular list
	elif not os.path.exists(settings[25]) and ("y" in settings[21] or "yes" in settings[21]):
		print ("The bacteria contaminants database was selected to be included, however the file path is invalid -- exiting now")
		exit()	


	if ("y" in settings[20] or "yes" in settings[20]) and os.path.exists(settings[24]):	
		# loads fungi DB
		for line in fungi_db:
			line = line[:-2]
			#print(line)
			line = line.strip()
			line = parse_contaminant_name(line)
			contaminant_db[line] = line
	elif not os.path.exists(settings[24]) and ("y" in settings[20] or "yes" in settings[20]): 
		print ("The fungi contaminants database was selected to be included, however the file path is invalid -- exiting now")
		exit()

	if ("y" in settings[19] or "yes" in settings[19]) and os.path.exists(settings[23]):
		# loads insects DB
		for line in insects_db:	
			line = line[:-2]
			#print (line)
			line = line.strip()
			line = parse_contaminant_name(line)
			contaminant_db[line] = line
	elif ("y" in settings[19] or "yes" in settings[19]) and not os.path.exists(settings[23]):
		print ("The insects contaminants database was selected to be included, however the file path is invalid -- exiting now")
		exit()

	'''		
	print ("searching in contaminants for bacteria: ")
	print (contaminant_db.get("Escherichia coli"))
	print (contaminant_db.get("Agaporomorphus tambopatensis"))
	print (contaminant_db.get("Escherichia albertii"))
	print (contaminant_db.get("Denticollinae"))
	print (contaminant_db.get("Escherichia coli str. K-12 substr. MG1655"))	
	'''
	return contaminant_db
	
def write_contaminants_log(element,log_name):
	f = open(log_name, "a")
	f.write(str(element))

'''
Assumptions before matching
1) The best possible match has already been found for every term
2) IF a fasta file/wGI numbers is present then there then every element in the search results has a corresponding entry in the FASTA file (THE REVERSE IS NOT TRUE -- since only best matches remain in search results)

Steps
1) IF fasta file w/GI numbers is present, scan through
2) For each element in the FASTA file without GI numbers, write all matching search results to annotation.tsv

'''
	
	
def match_fasta(db):

	#print (len(db))
	global num_queries_informative_hit #1 or more informative hits
	global num_queries_no_hit #no hits
	global num_queries_uninformative
	global db_count
	global annotation_log_entries
	global temp_log_entries
	temp_log_entries = dict()	
	global fasta_db
	global fasta_db_description
	global fasta_db_species
	global fasta_no_gi
	global contaminants_found

	global query_gi_association_db_0
	global query_gi_association_db_1
	global query_gi_association_db_2

	global nohits_found

	global min_coverage
	nohits_found = dict()
	query_gi_final = dict()

	if db_count == 999:
		f_nohits = open(output_folder+"//nohits_"+ date + ".txt","a")
		f_contaminants = open(output_folder+"//contaminants_" + date + ".txt","a")		
	
	num_queries_no_hit = 0
	annotation_log_entries = dict()	
	print ("matching fasta")
	if fasta_no_gi != "no_file":
		#scan through every elemeny in fasta_no_gi and return the best match from DB
		for element in fasta_no_gi:
			# The query variable is the line from the DB containing the matching query
			
			query = db.get(element)

			# try to find a match with no isomer in DB
			if query is None:
				query = db.get(element.split("|")[0])
			# if its still None, then there is no hit possible	
			if not query is None:
				key = str(get_gi_num_from_string(query[1]))
				#print (query[0])
				#print (key)
				query_coverage = float(query[3]) / float(len(fasta_no_gi[query[0]]))

				if db_count == 999 and fasta_db_species.get(get_gi_num_from_string(query[1])) in contaminants:
					contaminants_found[element] = element
					f_contaminants.write(str(query + [fasta_db_description[key]] + [fasta_db_species[key]]) + "\n")
	
				if annotation_log_entries.get(key) is None:
					annotation_log_entries[element] = db[element] + [fasta_db_description[key]] + [fasta_db_species[key]]
					temp_log_entries[element] = db[element] + [fasta_db_description[key]] + [fasta_db_species[key]]
				else:
					#print ("hit from more than 1 DB detected")
					#print (db[key])
					#print ([fasta_db_description[key]] + [fasta_db_species[key]])
					annotation_log_entries[element] = find_best_query_result(annotation_log_entries[element], query)
					temp_log_entries[element] = find_best_query_result(temp_log_entries[element], query)
			else:
				num_queries_no_hit += 1	
				if db_count == 999:
					f_nohits.write(element+"\n")
					#print ("adding to nohits log")
					nohits_found[element] = element	
	print ("Finished matching db: " + str(get_db_name(db_count)))
	#db_count += 1			

def parse_fasta_no_gi(file_name):
	return_dict = dict()
	query = ""
	description = ""
	counter = 0
        with open(file_name,'r') as file:
                for line in file:
			if ">" in line:
				#counter += 1
				#print (query)
				#print (description)
				if description != "":
					query = str(query.split(" ")[0])
					#query = str(query.split("|")[0])
					#print (query)
					counter += 1
					#if not return_dict.get(query) is None:
					#	print ("repeat query name")
					#	print (return_dict[query])
					#	print (query)
					return_dict[query] = description
					description = ""
				query = str(line[1:-1])
				#print (query)
				#return_dict[query] = description
			else:
				if "\n" in line:
					description += str(line[:-1])
				else:			
					description += str(line)
	
	return_dict[query] = description
	print ("Finished parsing query fasta")
	print ("Query fasta contains: " + str(len(return_dict)) + " elements")
	#print (len(return_dict))
	#print (str(counter))
	return return_dict	


'''
Method of calculation of n50 statistic (n25, n75 could be added using a similar method)

step 1) Sum up all of the contig lengths to obtain sum of all contig lengths

step 2) set n50_lengths to 1/2 of sum of all contig lengths

step 3) while sum < n50_lengths, i++ (i is the index of n50 statistic)

step 4) return lengths[i-1]

'''

def get_n50_statistic(lengths):
	#print ("starting n50 calcs")
	n50_lengths = 0
	n50_index = 0
	sum = 0
	for element in lengths:
		n50_lengths += element	
	n50_lengths = 0.5 * n50_lengths
	
	while sum < n50_lengths:
		sum += lengths[n50_index]
		n50_index += 1
	return lengths[n50_index - 1]




def get_median(lengths):
	odd, halfway = len(lengths) % 2, len(lengths) / 2
	if odd:  # odd number
		return lengths[halfway]  # e.g. for length 3, returns lengths[1]
	else:
		return 0.5 * (lengths[halfway - 1] + lengths[halfway])




def write_xml(filename, results_db):
	global number_db
	global fasta_no_gi
	global fasta_db
	global is_tair
	global fasta_db_species
	global contaminants
	global min_coverage
	#print (filename)
	#print (number_db)
	#rite header first
	if not os.path.exists(filename+".xml"): # if file doesnt exist create it
		file = open(filename+".xml", "w")
		print ("creating xml output with name:\t" + filename + ".xml")
		#print (fasta_no_gi)
	
		for key in results_db:

			result = results_db[key]
			query_coverage = float(result[3]) / float(len(fasta_no_gi[result[0]])) 
			file.write("<?xml version=\"1.0\" ?>\n")
			file.write("<!DOCTYPE BlastOutput\n")
			file.write("\tPUBLIC \'-//NCBI//NCBI BlastOutput/EN\'\n")
			file.write("\t'NCBI_BlastOutput.dtd'>\n")
			file.write("<BlastOutput>\n")
			file.write("\t<BlastOutput_program>BLASTX</BlastOutput_program>\n")
			file.write("\t<BlastOutput_version>BLASTX 2.2.25+</BlastOutput_version>\n")
			file.write("\t<BlastOutput_reference>Altschul, et. al.</BlastOutput_reference>\n")
			file.write("\t<BlastOutput_db>" + str(get_db_name(number_db)) + "</BlastOutput_db>\n")
			file.write("\t<BlastOutput_query-ID>1</BlastOutput_query-ID>\n")
			if not is_tair:
				file.write("\t<BlastOutput_query-def>" + result[0] + "</BlastOutput_query-def>\n")
			else:
				file.write("\t<BlastOutput_query-def>" + result[1] + "</BlastOutput_query-def>\n")
			file.write("\t<BlastOutput_query-len>" + result[3] + "</BlastOutput_query-len>\n")
			file.write("\t<BlastOutput_param>\n")
			file.write("\t\t<Parameters>\n")
			file.write("\t\t\t<Parameters_expect>10</Parameters_expect>\n")
			file.write("\t\t\t<Parameters_include>0</Parameters_include>\n")
			file.write("\t\t\t<Parameters_sc-match>1</Parameters_sc-match>\n")
			file.write("\t\t\t<Parameters_sc-mismatch>-3</Parameters_sc-mismatch>\n")
			file.write("\t\t\t<Parameters_gap-open>5</Parameters_gap-open>\n")
			file.write("\t\t\t<Parameters_gap-extend>2</Parameters_gap-extend>\n")
			file.write("\t\t\t<Parameters_filter>D</Parameters_filter>\n")
			file.write("\t\t</Parameters>\n")
			file.write("\t</BlastOutput_param>\n")	
			file.write("\t<BlastOutput_iterations>\n")
			file.write("\t\t<Iteration>\n")		
			file.write("\t\t\t<Iteration_iter-num>1</Iteration_iter-num>\n")
			file.write("\t\t\t<Iteration_query-ID>1</Iteration_query-ID>\n")
			if not is_tair:
				file.write("\t\t\t<Iteration_query-def>" + result[0] + "</Iteration_query-def>\n")
			else:
				file.write("\t\t\t<Iteration_query-def>" + result[1] + "</Iteration_query-def>\n")
			file.write("\t\t\t<Iteration_query-len>" + result[3] + "</Iteration_query-len>\n")
			file.write("\t\t\t<Iteration_hits>\n")
			

			#depending on legnth of result loop through this later
		

			for count in range(0, number_db):
				if not fasta_db.get(get_gi_num_from_string(result[1])) is None:
					file.write("\t\t\t\t<Hit>\n")
					file.write("\t\t\t\t\t<Hit_num>" + str(count) + "</Hit_num>\n")
					file.write("\t\t\t\t\t<Hit_id>" + result[1] + "</Hit_id>\n")
					file.write("\t\t\t\t\t<Hit_def>" + result[13] + " [" + result[14] + "]" + "</Hit_def>\n")
					file.write("\t\t\t\t\t<Hit_accession>" + get_gi_num_from_string(result[1]) + "</Hit_accession>\n")
					file.write("\t\t\t\t\t<Hit_len>" + result[3] + "</Hit_len>\n")
					file.write("\t\t\t\t\t<Hit_hsps>\n")
					file.write("\t\t\t\t\t\t<Hsp>\n")
					file.write("\t\t\t\t\t\t\t<Hit_num>" + str(count) + "</Hit_num>\n")	
					file.write("\t\t\t\t\t\t\t<Hsp_bit-score>" + result[3] + "</Hsp_bit-score>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_score>0</Hsp_score>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_evalue>" + result[10] + "</Hsp_evalue>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_query-from>" + result[6] + "</Hsp_query-from>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_query-to>" + result[7] + "</Hsp_query-to>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_hit-from>0</Hsp_hit-from>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_hit-to>0</Hsp_hit-to>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_pattern-from>0</Hsp_pattern-from>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_pattern-to>0</Hsp_pattern-to>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_query-frame>0</Hsp_query-frame>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_hit-frame>0</Hsp_hit-frame>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_identity>0</Hsp_identity>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_positive>0</Hsp_positive>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_gaps>0</Hsp_gaps>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_align-len>" + result[3] + "</Hsp_align-len>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_density>0</Hsp_density>\n")
					
					#print (fasta_no_gi[str(result[0])])
					#print ( fasta_db[get_gi_num_from_string(result[1])].replace("\n", ""))
					file.write("\t\t\t\t\t\t\t<Hsp_qseq>" + fasta_no_gi[str(result[0])] + "</Hsp_qseq>\n")
					file.write("\t\t\t\t\t\t\t<Hsp_hseq>" +  fasta_db[get_gi_num_from_string(result[1])].replace("\n", "") + "</Hsp_hseq>\n")
					
					file.write("\t\t\t\t\t\t</Hsp>\n")
					file.write("\t\t\t\t\t</Hit_hsps>\n")
					file.write("\t\t\t\t</Hit>\n")
					


			file.write("\t\t\t</Iteration_hits>\n")
			file.write("\t\t</Iteration>\n")
			file.write("\t</BlastOutput_iterations>\n")
			file.write("</BlastOutput>\n")



	else: #the file already exists
		print ("The output for the xml file already exists -- check file directory to rename the file/remove it")


# input will always be annotation_log_entries
def calc_stats(results):
	#global contaminants_found
	global temp_log_entries
	global db_count
	global fasta_no_gi
	global fasta_db
	global fasta_db_species
	global longest_query_length
	global shortest_query_length
	global median_query_length
	global num_contaminants
	global avg_length_query_sequences
	global num_queries
	global num_queries_informative_hit
	global num_queries_no_hit
	global num_queries_uninformative
	global n50_statistic
	global top_ten_hits
	global top_ten_contaminants
	global contaminants
	global min_coverage
	#print ("db length:   " + str(len(results)))
	top_ten_hits = dict()
	top_ten_contaminants = dict()
	num_queries_uninformative = 0
	num_queries_informative_hit = 0
	num_queries = 0
	num_queries_no_contamaints = 0
	median_query_length = list()
	avg_length_query_sequences = 0
	num_contaminants = 0
	longest_query_length = 0
	shortest_query_length = sys.maxint

	
	#print ("calc stats called")
	#print ("len results: " + str(len(results)))
	for element in dict(results):
		temp = results[element]
		#print (temp)
		num_queries += 1
		gi = str(get_gi_num_from_string(temp[1]))
		query_length = len(fasta_no_gi[temp[0]])
		query_coverage = float(temp[3]) / float(query_length)
		#print (fasta_db_species[gi])
		#print (contaminants.get(fasta_db_species[gi]))
		
		if not contaminants.get(str(fasta_db_species[gi])) is None:
			num_contaminants += 1
			if not top_ten_contaminants.get(fasta_db_species[gi]) is None:
				top_ten_contaminants[fasta_db_species[gi]] += 1		
			else:
				top_ten_contaminants[fasta_db_species[gi]] = 1
			del results[element] # this deletes the transcript so it will not be included in the XML file or the annotation
		else:	
			num_queries_no_contamaints += 1
			if fasta_db_description[gi] == "uninformative":
				num_queries_uninformative += 1
			else:
				num_queries_informative_hit += 1
			
			if not top_ten_hits.get(str(temp[14])) is None:
				top_ten_hits[str(temp[14])] += 1
			else:
				top_ten_hits[str(temp[14])] = 1
			# only include non-contaminant hits in statistics calc
			median_query_length.append(query_length)
			avg_length_query_sequences = float((avg_length_query_sequences * (num_queries_no_contamaints-1) + query_length) / num_queries_no_contamaints)
		
			if query_length > longest_query_length:
				longest_query_length = query_length
			if query_length < shortest_query_length:
				shortest_query_length = query_length

	#num_queries_no_hit = len(fasta_no_gi) - len(results)
	#num_queries = len(results)
	median_query_length.sort() # sorts least to greatest
	median_query_length.reverse() # reverses the order to greatest to least (median remains identical)
	n50_statistic = get_n50_statistic(median_query_length)
	median_query_length = get_median(median_query_length)
	#num_queries = num_queries_uninformative + num_queries_informative_hit + num_queries_no_hit	
	#print (top_ten_hits)
	


def remove_smallest(top_list):
	#print (len(top_list))
	min_value = None
	k = None
	for key, value in top_list.iteritems():
		#print (key)
		#print (value)
		if min_value is None:
			min_value = value
			k = key
		elif value < min_value:
			min_value = value	
			k = key 
	
	#print (k)
	#print (min_value)	
	#print (top_list[k])
	if top_list:
		del top_list[k]
	return top_list



def get_top_ten(top_list):
	#ret_list = list()
	
	if len(top_list) < 10:
		return top_list
	else:
		while len(top_list) > 10:
			top_list = remove_smallest(top_list)
			#print (top_list)
	return top_list


# this is a helper method to allow the sorted() call to sort the tuples by second value
def getKey(item):
	return item[1]

def print_stats():
	global final_output_temp
	global fasta_no_gi
	global longest_query_length
	global shortest_query_length
	global median_query_length
	global num_contaminants
	global avg_length_query_sequences
	global num_queries
	global num_queries_informative_hit
	global num_queries_no_hit
	global num_queries_uninformative
	global db_count
	global n50_statistic	
	global top_ten_hits
	
	temp = list()

	temp.append(["DB: "] + [str(get_db_name(db_count))])
	temp.append(["Number of queries with an informative hit:"] + [str(num_queries_informative_hit)])
	temp.append(["Number of queries with an uninformative hit:"] + [str(num_queries_uninformative)])
	temp.append(["Number of contaminants:"] + [str(num_contaminants)])
	temp.append(["Number of queries with no hit:"] + [str(num_queries_no_hit)])
	temp.append(["The top 10 hits by species:"])
	if top_ten_hits:
		for key,value in sorted(get_top_ten(top_ten_hits).iteritems(), key=getKey, reverse=True):
			temp.append([str(key)+":"] + [str(value)])
	else:
		temp.append(["No hits from this DB (possible error)"])

	temp.append(["The top 10 contaminants by species:"])
	
	if top_ten_contaminants:
		for key,value in sorted(get_top_ten(top_ten_contaminants).iteritems(), key=getKey, reverse=True):
			temp.append([str(key)+": "]+ [str(value)])
	else:
		temp.append(["No contaminants present"])
	
	temp.append([])	
	
	final_output_temp.append(temp)
	db_count += 1

	



def get_db_name(db_number):
	if db_number == 0:
		return str(settings[6])
	elif db_number == 1:
		return str(settings[9])
	elif db_number == 2:
		return str(settings[12])
	elif db_number == 999:
		return "combined_db (sum of previous DBs)"
	else:
		return "error"

def print_summary_stats():
	global db_count
	global final_output_temp
	global annotation_log_entries
	global longest_query_length
	global shortest_query_length
	global median_query_length
	global num_contaminants
	global avg_length_query_sequences
	global num_queries
	global num_queries_informative_hit
	global num_queries_uninformative
	global top_ten_contaminants
	global num_queries_no_hit
	global final_output
	global temp_log_entries	
	global db
	global db2
	global db3

	db_combined = dict()
	

	db_combined.update(db)
	
	for element in db2:
		if db_combined.get(element) is None:
			db_combined[element] = db2[element]
		else:
			db_combined[element] = find_best_query_result(db_combined.get(element), db2[element])


	for element in db3:
		if db_combined.get(element) is None:
			db_combined[element] = db3[element]
		else:
			db_combined[element] = find_best_query_result(db_combined.get(element), db3[element])

	#db_combined.update(db2)
	#db_combined.update(db3)



	num_queries_no_hit = 0	
	#db_count = 999 
	match_fasta(db_combined)
	calc_stats(temp_log_entries)
	#print (len(fasta_no_gi) - num_queries_uninformative - num_queries_informative_hit)
	#num_queries_no_hit = len(fasta_no_gi) - len(temp_log_entries)
	temp = list()
	#print (num_queries_no_hit)
	
	temp.append(["Summary Statistics on the Transcriptome Assembly Input (Query)"])
	temp.append(["Total Number of Query Sequences: "] + [str(len(fasta_no_gi))])
	temp.append(["Median Query Length: "] + [str(median_query_length)])
	temp.append(["Average Query Length: "] + [str(round(avg_length_query_sequences,2))])
	temp.append(["Shortest Query Length: "] + [str(shortest_query_length)])
	temp.append(["Longest Query Length: "] + [str(longest_query_length)])
	temp.append(["N50 Statistic: "] + [str(n50_statistic)])
	temp.append(["Number of queries with an informative hit: "] + [str(num_queries_informative_hit)])
	temp.append(["Number of queries with an uninformative hit: "] + [str(num_queries_uninformative)])
	temp.append(["Number of queries with no hit: "] + [str(num_queries_no_hit)])
	temp.append(["Number of contaminants: "] + [str(num_contaminants)])
	temp.append(["The top 10 hits by species: "])
	if top_ten_hits:
		for key,value in sorted(get_top_ten(top_ten_hits).iteritems(), key=getKey, reverse=True):
			temp.append([str(key) + ":"] + [str(value)])
	else:
		temp.append(["No Hits from this DB (possible error)"])
	
	temp.append(["The top 10 contaminants: "])
	
	if top_ten_contaminants:
		for key,value in sorted(get_top_ten(top_ten_contaminants).iteritems(), key=getKey, reverse=True):
			temp.append([str(key)] + [str(value)])
	else:
		temp.append(["No contaminants present"])
	temp.append([])
	final_output_temp.append(temp)
		

def find_best_hit_between_databases(db1_weight, db2_weight, db1_hit, db2_hit):
	global min_coverage
	global e_value

	db1_gi = get_gi_num_from_string(db1_hit[1])
	db2_gi = get_gi_num_from_string(db2_hit[1])

	db1_coverage = float(db1_hit[3]) / float(len(fasta_no_gi[db1_hit[0]]))
	db2_coverage = float(db2_hit[3]) / float(len(fasta_no_gi[db2_hit[0]]))


	if db1_weight == db2_weight:
		return find_best_query_result(db1_hit, db2_hit)
	elif db1_weight < db2_weight:
		if fasta_db_description[db2_gi] != "uninformative" and db2_coverage >= min_coverage and parse_e_value(db2_hit[10]) >= e_value:
			return db2_hit
		else:
			return find_best_query_result(db1_hit, db2_hit)
	else: # db1_weight > db2_weight
		if fasta_db_description[db1_gi] != "uninformative" and db1_coverage >= min_coverage and parse_e_value(db1_hit[10]) >= e_value:
			return db1_hit
		else:
			return find_best_query_result(db1_hit, db2_hit)


def print_weighted_summary_stats():
	global db_count
	global final_output_temp
	global annotation_log_entries
	global longest_query_length
	global shortest_query_length
	global median_query_length
	global num_contaminants
	global avg_length_query_sequences
	global num_queries
	global num_queries_informative_hit
	global num_queries_uninformative
	global top_ten_contaminants
	global num_queries_no_hit
	global final_output
	global temp_log_entries
	global db
	global db2
	global db3
	global db_score
	global db2_score
	global db3_score

	db_combined = dict()
	
	db_combined.update(db)
	
	for element in db2:
		if db_combined.get(element) is None:
			db_combined[element] = db2[element]
		else:
			db_combined[element] = find_best_hit_between_databases(db_score, db2_score, db_combined.get(element), db2[element])


	for element in db3:
		if db_combined.get(element) is None:
			db_combined[element] = db3[element]
		else:
			if not db.get(element) is None:
				db_combined[element] = find_best_hit_between_databases(db_score, db3_score,db_combined.get(element), db3[element])
			else:
				db_combined[element] = find_best_hit_between_databases(db2_score, db3_score,db_combined.get(element), db3[element])
	

	num_queries_no_hit = 0
	#db_count = 999
	match_fasta(db_combined)
	calc_stats(temp_log_entries)	

	

	temp = list()
	temp.append(["Weighted Summary Statistics on the Transcriptome Assembly Input (Query)"])
	temp.append(["Total Number of Query Sequences: "] + [str(len(fasta_no_gi))])
	temp.append(["Median Query Length: "] + [str(median_query_length)])
	temp.append(["Average Query Length: "] + [str(round(avg_length_query_sequences,2))])
	temp.append(["Shortest Query Length: "] + [str(shortest_query_length)])
	temp.append(["Longest Query Length: "] + [str(longest_query_length)])
	temp.append(["N50 Statistic: "] + [str(n50_statistic)])
	temp.append(["Number of queries with an informative hit: "] + [str(num_queries_informative_hit)])
	temp.append(["Number of queries with an uninformative hit: "] + [str(num_queries_uninformative)])
	temp.append(["Number of queries with no hit: "] + [str(num_queries_no_hit)])
	temp.append(["Number of contaminants: "] + [str(num_contaminants)])
	temp.append(["The top 10 hits by species: "])
	if top_ten_hits:
		for key,value in sorted(get_top_ten(top_ten_hits).iteritems(), key=getKey, reverse=True):
			temp.append([str(key) + ":"] + [str(value)])
	else:
		temp.append(["No Hits from this DB (possible error)"])

	temp.append(["The top 10 contaminants: "])

	if top_ten_contaminants:
		for key,value in sorted(get_top_ten(top_ten_contaminants).iteritems(), key=getKey, reverse=True):
			temp.append([str(key)] + [str(value)])
	else:
		temp.append(["No contaminants present"])
	temp.append([])
	final_output_temp.append(temp)



#entry point of script				
if __name__ == '__main__':
	start_time = time.clock()
	global final_output_temp
	final_output_temp = list()
	global debug_uninformative_list
	debug_uninformative_list = dict()
	arguments_list = sys.argv
	global date	
	date = str(datetime.datetime.now()).replace(" ","_")
	print(date)
	
	global output_folder

	# ensure that the config file given as a parameter exists

	if len(arguments_list) < 2 or len(arguments_list) > 3:
		print ("No argument entered, this script must be ran with the configuration file as a parameter, of the form:")
		print ("python report_generator.py configuration_file.txt")
		print ("Execution of this script will now abort")
		exit()

	if os.path.exists(arguments_list[1]):
		settings = parse_config_file(str(arguments_list[1])) #sets up the settings
	else:
		print ("The configuration file specified does not exist -- aborting execution")
		exit()

	if settings[3] != "":	
		output_folder = str(settings[3]) + "_output_" + date
	else:
		output_folder = "output_" + date

	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	output_log = output_folder+"//log_" + date + ".txt"

	global final_output
	final_output = dict()	
	global top_ten_hits
	global top_ten_contaminants	
	global annotation_log_entries
	annotation_log_entries = dict()
	global temp_log_entries
	temp_log_entries = dict() # used for calculating statistics on DBs
	global filter_list
	global contaminants
	global number_db
	global counter
	global contaminants_found
	global nohits_found	
	global fasta_db
	global fasta_db_description
	global fasta_db_species
	
	nohits_found = dict()

	global fasta_no_gi

	global is_tair
	is_tair = False
	global db_count # this var keeps track of which db is currently being loaded into memory (in order to keep query/gi associations per db)
	db_count = 0

	# These three contain the query name gi number associations for the three DBs so matching can be done in O(1) speed, at the cost of additional ram usage
	global query_gi_association_db_0
	global query_gi_association_db_1
	global query_gi_association_db_2
	
	query_gi_association_db_0 = dict()
	query_gi_association_db_1 = dict()
	query_gi_association_db_2 = dict()

	global e_value # this constant refers to the threshold to use for determining a "best match" in terms of queries
	e_value = parse_e_value(settings[14]) 
	#print("min e val is :: " + str(e_value))
	#The following global variables are used to record statistical data in order to generate a log
	
	#query length variables
	global longest_query_length
	global shortest_query_length
	global median_query_length
	global num_contaminants
	global n50_statistic
	n50_statistic = 0
	num_contaminants = 0	
	longest_query_length = -2147483648
	shortest_query_length = 2147483647
	median_query_length = list()

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
	
	global db_type
	global min_coverage	
	global number_db	
	min_coverage = settings[13]


	db_type = settings[2]
	filter_list = load_filter_list("filter_list.txt")
	contaminants = build_contaminants_db()
	contaminants_found = dict()
	print("Contaminants DB built")
	output = output_folder+"//default_output_annotation_" + date +".tsv"
	number_db = int(settings[1]) # the number of databases being parsed
	counter = number_db

	global db
	global db2
	global db3


	db = dict()
	db2 = dict()
	db3 = dict()

	global db_score
	global db2_score
	global db3_score

	db_score = 0
	db2_score = 0
	db3_score = 0


	with open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
		tsv_new = csv.writer(tsv_new, delimiter='\t')
		#The first row
		if (db_type == "diamond") or (db_type == "blastx"):
			row = ["Query","Subject_id","Identity(%)","Alignment_length","Mismatches","Number of gap opens","Query_start","Query_end" \
			,"Subject_start","Subject_end","E-value","Bit_score","Origin Database","Subject Description","Species"]
			
			#for x in range(0, (number_db - 1)):
			#	row += ["Subject_description","Species","Subject_id","Alignment_length", \
			#	"Mismatches","Query_start","Query_end","Subject_start","Subject_end","E-value","Bit_score","Subject_description","Species"]
			#row +=["has GFF3 data","GFF3 start","GFF3 stop","GFF3 ORF"]
		else:
			print("ncbi - not supported in this version of the script - now aborting")
			exit()



		tsv_new.writerow(row)
		
		
		#print(range (0, (number_db)))
		
		if number_db == 1: #database/fasta pair 1
			print("1 database/fasta pair")
			db_score = int(settings[4])
			#print (settings[6] != "")
			if settings[0] != "":
				fasta_no_gi = parse_fasta_no_gi(settings[0])
			else:
				fasta_no_gi = "no_file"
			#print (fasta_no_gi["ptAAH promoter plus cDNA"])			
			[fasta_db, fasta_db_description, fasta_db_species] = multi_fasta_parse(settings[5])
			print ("multi fasta file done")
			
			if db_type == "ncbi":
				db = ncbi_format_db_parse(settings[6])
			else:
				db = diamond_format_db_parse(settings[6])
			print ("db complete -- now matching")
			print (str(time.clock() - start_time) + " ::  time to complete db")
			db_count = 0
			
			match_fasta(db)
			calc_stats(temp_log_entries)
			print_stats()
			final_output.update(temp_log_entries)

			uninformative_debug_log(debug_uninformative_list)

			if settings[15] == "yes" or settings[15] == "y":
				write_xml(output_folder+"//blastxml_" + "db1" + "_" + date, annotation_log_entries)
				#write_xml("blastxml_" + date, annotation_log_entries)

			print (str(time.clock() - start_time) + " :: time to match fasta")
			
		
		elif number_db == 2: #database/fasta pair 2
			print("2 database/fasta pair")
			db_score = int(settings[4])
			db2_score = int(settings[7])						
			if settings[0] != "":
				fasta_no_gi = parse_fasta_no_gi(settings[0])
			else:
				fasta_no_gi = "no_file"

			[fasta_db, fasta_db_description, fasta_db_species] = multi_fasta_parse(settings[5])
			[fasta_db2, fasta_db_description2, fasta_db_species2] = multi_fasta_parse(settings[8])

			


			fasta_db.update(fasta_db2)
			fasta_db_description.update(fasta_db_description2)
			fasta_db_species.update(fasta_db_species2)
			
			print ("multi fasta file done")

			
			if db_type == "ncbi":
				db = ncbi_format_db_parse(settings[6])
				db2 = ncbi_format_db_parse(settings[9])
			else:
				db = diamond_format_db_parse(settings[6])
				db2 = diamond_format_db_parse(settings[9])
		
			#print (len(db))
			#print (len(db2))
	
			print("db, db2 complete -- now matching")
			print(str(time.clock() - start_time) + " ::  time to complete db, db2")
			db_count = 0
			match_fasta(db)
			calc_stats(temp_log_entries)
			print_stats()

			final_output.update(temp_log_entries)
			#print (len(temp_log_entries))
			if settings[15] == "yes" or settings[15] == "y":
				write_xml(output_folder+"//blastxml_" + "db1_" + date, annotation_log_entries)

			match_fasta(db2)
			calc_stats(temp_log_entries)
			#print (len(temp_log_entries))
			print_stats()
			
			#final_output.update(temp_log_entries)
			#print (len(final_output))
			if settings[15] == "yes" or settings[15] == "y":
				write_xml(output_folder+"//blastxml_" + "db2_" + date, annotation_log_entries)
			print (str(time.clock() - start_time) + " :: time to match fasta")	


		elif number_db == 3: #database/fasta pair 3
			print("3 database/fasta pair")
			db_score = int(settings[4])
			db2_score = int(settings[7])
			db3_score = int(settings[10])
			if settings[0] != "":
				fasta_no_gi = parse_fasta_no_gi(settings[0])
			else:
				fasta_no_gi = "no_file"

			[fasta_db, fasta_db_description, fasta_db_species] = multi_fasta_parse(settings[5])
			[fasta_db2, fasta_db_description2, fasta_db_species2] = multi_fasta_parse(settings[8])
			[fasta_db3, fasta_db_description3, fasta_db_species3] = multi_fasta_parse(settings[11])

			fasta_db.update(fasta_db2)
			fasta_db.update(fasta_db3)
			fasta_db_description.update(fasta_db_description2)
			fasta_db_description.update(fasta_db_description3)
			fasta_db_species.update(fasta_db_species2)
			fasta_db_species.update(fasta_db_species3)

			print ("multi fasta file done")

			if db_type == "ncbi":
				db = ncbi_format_db_parse(settings[6])
				db2 = ncbi_format_db_parse(settings[9])
				db3 = ncbi_format_db_parse(settings[12])
			else:
				db = diamond_format_db_parse(settings[6])
				db2 = diamond_format_db_parse(settings[9])
				db3 = diamond_format_db_parse(settings[12])
			
			print("db, db2, db3 complete -- now matching")
			print(str(time.clock() - start_time) + " ::  time to complete db, db2, db3")
			db_count = 0
			match_fasta(db)
			calc_stats(temp_log_entries)
			print_stats()
			final_output.update(temp_log_entries)
			if settings[15] == "yes" or settings[15] == "y":
				write_xml(output_folder+"//blastxml_" + "db1_" + date, annotation_log_entries)
			match_fasta(db2)
			calc_stats(temp_log_entries)
			print_stats()
			final_output.update(temp_log_entries)
			if settings[15] == "yes" or settings[15] == "y":
				write_xml(output_folder+"//blastxml_" + "db2_" + date, annotation_log_entries)
			match_fasta(db3)
			calc_stats(temp_log_entries)
			print_stats()
			final_output.update(temp_log_entries)
			if settings[15] == "yes" or settings[15] == "y":
				write_xml(output_folder+"//blastxml_" + "db3_" + date, annotation_log_entries)
			print (str(time.clock() - start_time) + " :: time to match fasta")

	
		#print ("calculating statistics")

		#calc_stats(annotation_log_entries)

		
		#print ("writing annotation log entires")
		#for key in annotation_log_entries:
		#	tsv_new.writerow(final_output[key])
		
		#del annotation_log_entries
			
		#print ("writing no hits log")
		#after parsing of all fasta elements add all missed hits to nohits file
		'''
		f_nohits = open(output_folder+"//nohits_"+ date + ".txt","a")
		for item in nohits_found:
			f_nohits.write(item+"\n")
			#write_log([item],output_folder+"//nohits_"+ date + ".txt")
		f_nohits.close()	
		'''
	
		#print ("writing contaminants log")
		'''
		f_contaminants = open(output_folder+"//contaminants_" + date + ".txt","a")
		for key in contaminants_found:
			f_contaminants.write(key+"\n")
			#write_contaminants_log([contaminants_found.get(key)],output_folder+"//contaminants_" + date + ".txt")
		f_contaminants.close()
		#del contaminants_found	
		'''
		print (str(time.clock() - start_time) + " seconds")
		print("complete -- annotation file now available -- continuing on to create log file")
		#print ("printing db")
		#print (db)

		#num_queries_informative_hit = (num_queries - num_queries_uninformative)
		if db_score == db2_score == db3_score:	
			print_summary_stats()
		else:
			print_summary_stats()
			db_count = 999
			print_weighted_summary_stats()

		for key in temp_log_entries:
			tsv_new.writerow(temp_log_entries[key])	



	


	if settings[15] == "yes" or settings[15] == "y":
		write_xml(output_folder+"//blastxml_" + "combined_db" + "_" + date, temp_log_entries)

	#TODO when config file becomes a parameter => make it display the given filepath
	final_output_temp.append([["Path to configuration file: " + str(os.path.dirname(os.path.realpath(__file__))) + "/configuration_file.txt"]])
	#print (final_output_temp)
	if not os.path.exists(output_log + ".txt"):
		with open(os.path.dirname(os.path.realpath(__file__)) + "//" + output_log, 'a') as tsv_log:
			tsv_log = csv.writer(tsv_log, delimiter='\t')
			for count in range(0, len(final_output)):
				if final_output_temp:
					temp = final_output_temp.pop()
					for line in temp:
						#line = line.replace("\"","")
						tsv_log.writerow(line)
			
		
        
	go_counts = [0, 0, 0]
	go_interpro_counts = [0, 0, 0]
	at_least_1_interpro = [0, 0, 0]	
	domain_ids = [0, 0]
	at_least_1 = [0, 0, 0]
	


	if settings[16] != "" and os.path.exists(settings[16]):
                print ("option to include interpro output has been checked, filepath to interpro: " + settings[16])
		#print (settings[17])	
		if settings[17] == "" and os.path.exists(settings[17]):
                        combine_annotations.main(["--input"] + [output] + ["--interpro"] + [settings[16]] + ["--output"] + [output_folder+"//combined_annotation_"+date+".tsv"])
			go_interpro_counts = combine_annotations.get_go_interpro_counts()
			at_least_1_interpro = combine_annotations.get_at_least_1_interpro() # same order as other, (C,P,F)
		else:
                        combine_annotations.main(["--input"] + [output] + ["--blast2go"] + [settings[17]] + ["--interpro"] + [settings[16]] + ["--output"] + [output_folder+"//combined_annotation_"+date+".tsv"])
			domain_ids = combine_annotations.get_num_sequences_identification()
			at_least_1 = combine_annotations.get_at_least_1() # C, P, F
			go_counts = combine_annotations.get_go_counts()
			go_interpro_counts = combine_annotations.get_go_interpro_counts()
			at_least_1_interpro = combine_annotations.get_at_least_1_interpro() # same order as other, (C,P,F)
	elif settings[17] != "":
		combine_annotations.main(["--input"] + [output] + ["--blast2go"] + [settings[17]] + ["--output"] + [output_folder+"//combined_annotation_"+date+".tsv"]) 
		domain_ids = combine_annotations.get_num_sequences_identification()
		at_least_1 = combine_annotations.get_at_least_1() # C, P, F
		go_counts = combine_annotations.get_go_counts()
	else:
                print ("no interpro file has been specified, skipping combined annotation step")

	#potentially move this over combined_annotations.py eventually?
	if not os.path.exists(output_log + ".txt"):
		with open(os.path.dirname(os.path.realpath(__file__)) + "//" + output_log, 'a') as tsv_log:
			tsv_log = csv.writer(tsv_log, delimiter='\t')
			tsv_log.writerow(["Interpro File: "] + [str(settings[16])]) # TODO put an if statement here to check for multiple interpro files and
			if settings[17] != "" and not settings[17] is None:
				tsv_log.writerow(["Blast2GO File: "] + [str(settings[17])])
			tsv_log.writerow(["Number of sequences with Domain Identification: "] + [str(domain_ids[0])])
			tsv_log.writerow(["Number of sequences without Domain Identification: "] + [str(domain_ids[1])])
			tsv_log.writerow(["Blast2GO Gene Ontology Stats"])
			tsv_log.writerow(["Number of Components: "] + [str(go_counts[0])])
			tsv_log.writerow(["Number of Functions: "] + [str(go_counts[2])])
			tsv_log.writerow(["Number of Processes: "] + [str(go_counts[1])])
			tsv_log.writerow(["Number of Transcripts with at least 1 Component: "] + [str(at_least_1[0])])
			tsv_log.writerow(["Number of Transcripts with at least 1 Function: "] + [str(at_least_1[2])])
			tsv_log.writerow(["Number of Transcripts with at least 1 Process: "] + [str(at_least_1[1])])
			tsv_log.writerow(["Interpro Gene Ontology Stats (Totals)"])
			tsv_log.writerow(["Component: "] + [str(go_interpro_counts[0])])
			tsv_log.writerow(["Function: "] + [str(go_interpro_counts[2])])
			tsv_log.writerow(["Process: "] + [str(go_interpro_counts[1])])
			tsv_log.writerow(["Number of Transcripts with at least 1 Component: "] + [str(at_least_1_interpro[0])])
			tsv_log.writerow(["Number of Transcripts with at least 1 Function: "] + [str(at_least_1_interpro[2])])
			tsv_log.writerow(["Number of Transcripts with at least 1 Process: "] + [str(at_least_1_interpro[1])])


	print ("Completed everything -- Now exiting")
