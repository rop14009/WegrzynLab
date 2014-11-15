import datetime
import time
import os
import sys
import csv
import string
import multiprocessing

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
	last_gi = ""
	with open(file_name,'r') as file:
		for line in file:
			if line[:1] == ">": # line contains description of sequence
				if current_protein_seq != "":
					last_gi = current_gi
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

	#print (last_gi)
	#print (current_gi)

	fasta_db[current_gi] = current_protein_seq
	fasta_db_description[current_gi] = current_desc
	fasta_db_species[current_gi] = current_species

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
	global fasta_no_gi
	global contaminants_found

	global num_queries_informative_hit
	#global num_queries_no_hit #no hits
	global num_queries_uninformative

	global min_coverage # 0.7 --> refers to 70% coverage
        global e_value # 1e-5 / 0.00001

        query1_gi = get_gi_num_from_string(query1[1])
        query2_gi = get_gi_num_from_string(query2[1])
	

	# this is done to prevent the ptAAH cDNA etc... error where no match to fasta_no_gi can be made
	#query1[0] = query1[:-query1[0].find(" ")]
	#query2[0] = query2[:-query2[0].find(" ")]

	#int(line[7]) - int(line[6])
	#print (query1[0])
	#print (str(query1))
	#print (str(query2))
	#print (fasta_no_gi[str(query1[0])[2:-2]])

	#query1_search_length = int(query1[7]) - int(query1[6])
	#query2_search_length = int(query2[7]) - int(query2[6])
	
	#print (str(query1[0]))	
	#print (len(fasta_no_gi[str(query1[0])[2:-2]]))	
	
	#print (fasta_no_gi["ptAAH promoter plus cDNA"])
	

	#print ("query1 length: " + str(query1_search_length))
	#print ("query2 length: " + str(query2_search_length))


	#query1_coverage = float(query1_search_length) / float(len(fasta_no_gi[str(query1[0])]))
	#query2_coverage = float(query2_search_length) / float(len(fasta_no_gi[str(query2[0])]))

	query1_coverage = float(query1[3]) / float(len(fasta_no_gi[str(query1[0])]))
	query2_coverage = float(query2[3]) / float(len(fasta_no_gi[str(query2[0])]))

	print ("query1 coverage: " + str(query1_coverage))
	print ("query2 coverage: " + str(query2_coverage))

	


	'''
	print ("query1: " + str(query1))	
	print ("query2: " + str(query2))
	print ("fasta no gi lengths")
	print (len(fasta_no_gi[str(query1[0])[2:-2]]))
	print (len(fasta_no_gi[str(query2[0])[2:-2]]))
	print ("coverages")
	print (query1_coverage)
	print (query2_coverage)
	'''
	if is_uninformative(fasta_db_description[query1_gi]):
		fasta_db_description[query1_gi] = "uninformative"

	if is_uninformative(fasta_db_description[query2_gi]):
		fasta_db_description[query2_gi] = "uninformative"
	

	

	if not is_uninformative(fasta_db_description[query1_gi]) and is_uninformative(fasta_db_description[query2_gi]):
		return query1

	if is_uninformative(fasta_db_description[query1_gi]) and not is_uninformative(fasta_db_description[query2_gi]):
		return query2

	if is_uninformative(fasta_db_description[query1_gi]) and is_uninformative(fasta_db_description[query2_gi]):
		#print ("test --  both uninformative")
		if (query1_coverage >= min_coverage and query1_coverage > query2_coverage):
			return query1

		if (query2_coverage >= min_coverage and query2_coverage > query1_coverage):
			return query2
		return query1

	if not is_uninformative(fasta_db_description[query1_gi]) and not is_uninformative(fasta_db_description[query2_gi]):
		#print ("test --  both informative")
		if (query1_coverage >= min_coverage and query1_coverage > query2_coverage):
			return query1

		if (query2_coverage >= min_coverage and query2_coverage > query1_coverage):
			return query2
		return query1

	print ("this statement has been reached, when it should never be reached")
	

	'''
	if is_uf (query1_coverage > min_coverage and query1_coverage > query2_coverage)
                        return query1

                if (query2_coverage > min_coverage and query2_coverage > query1_coverage)
                        return query2
ninformative(fasta_db_description[query1_gi]):
		if is_uninformative(fasta_db_description[query2_gi]):
			return query1
		else:
			print ("bettery query than: " + str(query1) + " : " + str(fasta_db_description[query2_gi]))
			return query2
	elif not is_uninformative(fasta_db_description[query2_gi]):
		if (query1_coverage > min_coverage and query2_coverage > min_coverage):
			if query1_coverage > query2_coverage:
				return query1
			else:
				print ("bettery query than: " + str(query1) + " : " + str(fasta_db_description[query2_gi]))
				return query2
		else:
			return query1
	else:
		#print ("error in find_best_query_result(query1, query2) this should not ever be reached")
	        print ("last case scenario occured")
		return query2 #this should not be reached, ever
	'''	

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




'''
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
			num_queries += 1	
			[key, in_db] = is_query_in_db(line[0], ncbi_db)

			if in_db:
				ncbi_db[key] = find_best_query_result(ncbi_db[key], line)
			else:
				ncbi_db[get_gi_num_from_string(line[1])] = line
		

		length = int(line[8]) - int(line[7])
		avg_length_query_sequence = float((avg_length_query_sequences * (num_queries-1) + length) / num_queries)
		median_query_length.append(length)
		if length > longest_query_length:
			longest_query_length = length
		if length < shortest_query_length:
			shortest_query_length = length	

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

	global longest_query_length
        global shortest_query_length
        global median_query_length
	global avg_length_query_sequences
	global num_queries
	
	#median_query_length = 0

	#median_query_length = list()
	
	with open(file_name, "r") as file:
		file_tsv = csv.reader(file, delimiter='\t')
		for line in file_tsv:
			num_queries += 1
			#print (get_gi_num_from_string(line[1]))
			[key, in_db] = is_query_in_db(line[0], usearch_db)		

			if in_db:
				#print (usearch_db[key] )
				#print (line)
				usearch_db[str(key)] = find_best_query_result(usearch_db[key], line)
			else:
				usearch_db[str(get_gi_num_from_string(line[1]))] = line
		
			#length = int(line[7]) - int(line[6])
			length = float(line[3])
	
			#print (length)
			avg_length_query_sequences = float((avg_length_query_sequences * (num_queries-1) + length) / num_queries)
			#print (avg_length_query_sequences)
			#insert_in_order(median_query_length, length)
			median_query_length.append(length)	
			if length > longest_query_length:
				longest_query_length = length
			if length < shortest_query_length:
				shortest_query_length = length

		return usearch_db



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
			file.write("".join(element)+"\n")
			file.close()
		else:
			file = open(log_name+".txt", "a")
			file.write("".join(element)+"\n")
			file.close()
	else:
		counter -= 1
		if not os.path.exists(log_name+".txt"): # if file doesnt exist create it
			file = open(log_name+".txt", "w")
			print ("creating contaminants log with name: " + log_name)
			if not counter > 0:
				file.write(element)
			else:
				file.write("".join(element)+"\n")
			file.close()
		else:
			file = open(log_name+".txt", "a")
			if not counter > 0:
				file.write(element)
			else:
				file.write("".join(element)+"\n")
			file.close()

'''
Assumptions before matching
1) The best possible match has already been found for every term
2) IF a fasta file/wGI numbers is present then there then every element in the search results has a corresponding entry in the FASTA file (THE REVERSE IS NOT TRUE -- since only best matches remain in search results)

Steps
1) IF fasta file w/GI numbers is present, scan through
2) For each element in the FASTA file without GI numbers, write all matching search results to annotation.tsv


'''
	
	
def match_fasta(db):
	global num_queries_informative_hit #1 or more informative hits
	global num_queries_no_hit #no hits
	global num_queries_uninformative
	
	global annotation_log_entries
	
	global fasta_db
	global fasta_db_description
	global fasta_db_species
	global fasta_no_gi
	
	#print (fasta_no_gi)

	if fasta_no_gi != "no_file":
		#scan through every elemeny in fasta_no_gi and return the best match from DB
		for element in fasta_no_gi:
			#print (fasta_no_gi[element])
			[key, is_present] = is_query_in_db(element, db)
			#print (key)
			#print (db[key])
			#print (is_present)
			if is_present:
				if annotation_log_entries.get(key) is None:
					annotation_log_entries[key] = db[key] + [fasta_db_description[key]] + [fasta_db_species[key]]
				else:
					annotation_log_entries[key] = annotation_log_entries[key] + (db[key] + [fasta_db_description[key]] + [fasta_db_species[key]])
				del db[key]
			else:
				#print ("no hit += 1")
				#print (element)
				num_queries_no_hit += 1


def parse_fasta_no_gi(file_name):
	return_dict = dict()
	query = ""
	description = ""
        with open(file_name,'r') as file:
                for line in file:
			#print (line)
			if ">" in line:
				#print (query)
				#print (description)
				if description != "":
					#print (query)
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
	#print (return_dict)
	return return_dict	



'''
Method of calculation of n50 statistic (n25, n75 could be added using a similar method)

step 1) Sum up all of the contig lengths to obtain sum of all contig lengths

step 2) set n50_lengths to 1/2 of sum of all contig lengths

step 3) while sum < n50_lengths, i++ (i is the index of n50 statistic)

step 4) return lengths[i-1]

'''

def get_n50_statistic(lengths):
	print ("starting n50 calcs")
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



#entry point of script				
if __name__ == '__main__':
	start_time = time.clock()
	arguments_list = sys.argv
	
	date = str(datetime.datetime.now())
	print(date)
	settings = parse_config_file("configuration_file.txt") #sets up the settings
	output_log = "log_" + date + ".txt"
	
	global annotation_log_entries
	annotation_log_entries = dict()
	global filter_list
	global contaminants
	global number_db
	global counter
	global contaminants_found
	
	global fasta_db
	global fasta_db_description
	global fasta_db_species

	global fasta_no_gi
	
	global e_value # this constant refers to the threshold to use for determining a "best match" in terms of queries
	e_value = parse_e_value(settings[14]) 
	print("min e val is :: " + str(e_value))
	#The following global variables are used to record statistical data in order to generate a log
	
	#query length variables
	global longest_query_length
	global shortest_query_length
	global median_query_length
	
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
	min_coverage = settings[13]


	db_type = settings[2]
	filter_list = load_filter_list("filter_list.txt")
	contaminants = build_contaminants_db()
	contaminants_found = dict()
	print("contaminants db built")
	output = "default_output_annotation_" + date +".tsv"
	number_db = int(settings[1]) # the number of databases being parsed
	counter = number_db
	db = dict()
	db2 = dict()
	db3 = dict()
	with open(os.path.dirname(os.path.realpath(__file__)) + "//" + output, 'w') as tsv_new:
		tsv_new = csv.writer(tsv_new, delimiter='\t')
		#The first row
		if db_type == "usearch":
			row = ["Query","Subject_id","Identity(%)","Alignment_length","Mismatches","Number of gap opens","Query_start","Query_end" \
			,"Subject_start","Subject_end","E-value","Bit_score","Subject Description","Species"]
			
			for x in range(0, (number_db - 1)):
				row += ["Subject_description","Species","Subject_id","Alignment_length", \
				"Mismatches","Query_start","Query_end","Subject_start","Subject_end","E-value","Bit_score","Subject_description","Species"]
			row +=["has GFF3 data","GFF3 start","GFF3 stop","GFF3 ORF"]
		else:
			print("ncbi")




		tsv_new.writerow(row)
		
		
		#print(range (0, (number_db)))
		
		n50_statistic = None

		if number_db == 1: #database/fasta pair 1
			print("1 database/fasta pair")
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
				db = usearch_format_db_parse(settings[6])
			print ("db complete -- now matching")
			print (str(time.clock() - start_time) + " ::  time to complete db")
			match_fasta(db)
			print (str(time.clock() - start_time) + " :: time to match fasta")
			
		
		elif number_db == 2: #database/fasta pair 2
			print("2 database/fasta pair")
						
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
				db = usearch_format_db_parse(settings[6])
				db2 = usearch_format_db_parse(settings[9])
			
			print("db, db2 complete -- now matching")
			print(str(time.clock() - start_time) + " ::  time to complete db, db2")
			match_fasta(db)
			match_fasta(db2)
			print (str(time.clock() - start_time) + " :: time to match fasta")	


		elif number_db == 3: #database/fasta pair 3
			print("3 database/fasta pair")

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
				db = usearch_format_db_parse(settings[6])
				db2 = usearch_format_db_parse(settings[9])
				db3 = usearch_format_db_parse(settings[12])

			print("db, db2, db3 complete -- now matching")
			print(str(time.clock() - start_time) + " ::  time to complete db, db2, db3")
			match_fasta(db)
			match_fasta(db2)
			match_fasta(db3)
			print (str(time.clock() - start_time) + " :: time to match fasta")


	
		print ("calculating n50 statistic")
		median_query_length.sort() # sorts least to greatest
		median_query_length.reverse() # reverses the order to greatest to least (median remains identical)
		n50_statistic = get_n50_statistic(median_query_length)
		print (str(time.clock() - start_time) + " :: time to calculate n50 statistic")
		
		median_query_length = get_median(median_query_length)

		print ("writing annotation log entires")
		for key in annotation_log_entries:
			tsv_new.writerow(annotation_log_entries[key])
		
			
		print ("writing no hits log")
		#after parsing of all fasta elements add all missed hits to nohits file
		for key in db:
			write_log(db.get(key)[0],"nohits_"+ date)


		print ("writing contaminants log")
		for key in contaminants_found:
			#print (contaminants_found.get(key))
			write_contaminants_log("\t".join(contaminants_found.get(key)),"contaminants_" + date)
		
		
	
	print (str(time.clock() - start_time) + " seconds")
	print("complete -- annotation file now available")


	if not os.path.exists(output_log+".txt"):
		with open(os.path.dirname(os.path.realpath(__file__)) + "//" + output_log, 'w') as tsv_log:
			tsv_log = csv.writer(tsv_log, delimiter='\t')
			tsv_log.writerow(["num_queries: "] + [str(num_queries)])
			tsv_log.writerow(["median query length: "] + [str(median_query_length)])
			tsv_log.writerow(["average query length: "] + [str(round(avg_length_query_sequences,2))])
			tsv_log.writerow(["shortest query length: "] + [str(shortest_query_length)])
			tsv_log.writerow(["longest query length: "] + [str(longest_query_length)])
			tsv_log.writerow(["n50 statistic: "] + [str(n50_statistic)])
			tsv_log.writerow(["num_queries_informative_hit: "] + [str(num_queries_informative_hit)])
			tsv_log.writerow(["num_queries_no_hit: "] + [str(num_queries_no_hit)])
			tsv_log.writerow(["num_queries_uninformative: "] + [str(num_queries_uninformative)])
			print ("log files complete -- exititing")
			exit()
