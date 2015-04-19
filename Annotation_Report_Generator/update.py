'''
This script updates the contmainant databases from their source, the NCBI database
Author: Sam Ginzburg
'''

import urllib2
import os
import re
import sys
import multiprocessing
from multiprocessing import Pipe

def download(url_info, file_name, download_status):
	base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy"
	final_url = base_url + "&WebEnv=" + url_info[0] + "&query_key=" + str(url_info[1]) + "&retstart=0" + "&retmax=10000" + "&rettype=taxon&retmode=text"
	#print (url_info[2])	
	f = open(os.getcwd() + file_name, 'wb')
	for x in range(1, (url_info[2]/10000 + 2)):
		try:
			response = urllib2.urlopen(final_url)
		except:
			print ("Unable to get response -- error")
		f.write(response.read())
		final_url = base_url + "&WebEnv=" + url_info[0] + "&query_key=" + str(url_info[1]) + "&retstart=" + str((x*10000)) + "&retmax=" + str(10000+(x*10000)) + "&rettype=taxon&retmode=text"
		download_status.send(round(100.00*float(((x)*10000)) / (float(url_info[2])+20000), 2))
	#print (file_name+" is complete")
	download_status.send(100.00)


def get_ncbi_search_results(query):
	base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
	query_url = "esearch.fcgi?db=taxonomy&term=" + query + "&usehistory=y"
	total_url = base_url + query_url

	web_env = ""
	key = 0
	count = 0 # number of terms in results

	try:
		response = urllib2.urlopen(total_url, timeout=5)
		#print (response.info())
		search_result = response.read()
		#print (search_result)
		
		web_env_re = re.compile('<WebEnv>(\S+)<\/WebEnv>')
		web_env = web_env_re.search(search_result).group()[8:-9]
	
		count_re = re.compile('<Count>(\S+)<\/Count>')		
		count = int(count_re.search(search_result).group()[7:-8])

		key_re = re.compile('<QueryKey>(\S+)<\/QueryKey>')
		key = int(key_re.search(search_result).group()[10:-11])
	except:
		print ("Failure to reach NCBI database for search query: \t" + query)

	return [web_env, key, count]

if __name__ == '__main__':

	# Check file paths, and make them if they don't exist
	if not os.path.exists("contaminant_databases"):
		os.makedirs("contaminant_databases")



	
	bacteria_parent, bacteria_child = Pipe()
	fungi_parent, fungi_child = Pipe()
	insects_parent, insects_child = Pipe()
	print ("Downloading all the databases from NCBI now...")

	bacteria_info = get_ncbi_search_results("txid2[Subtree]")
	fungi_info = get_ncbi_search_results("txid4751[Subtree]")
	insects_info = get_ncbi_search_results("txid6960[Subtree]")

	base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy"

	url_insects = base_url + "&WebEnv=" + insects_info[0] + "&query_key=" + str(insects_info[1]) + "&retstart=0" + "&retmax=" + str(insects_info[2]) + "&rettype=taxon&retmode=text"
	url_fungi = base_url + "&WebEnv=" + fungi_info[0] + "&query_key=" + str(fungi_info[1]) + "&retstart=0" + "&retmax=" + str(fungi_info[2]) + "&rettype=taxon&retmode=text"
	url_bacteria = base_url + "&WebEnv=" + bacteria_info[0] + "&query_key=" + str(bacteria_info[1]) + "&retstart=0" + "&retmax=" + str(bacteria_info[2]) + "&rettype=taxon&retmode=text"
	# check if contaminant_databases folder exists, if not create it

	#if not os.path.exists(os.getcwd()+"contaminant_databases"):
	#	os.makedirs(os.getcwd()+"contaminant_databases")
	'''
	download(bacteria_info, '/contaminant_databases/bacteria_db.txt')
	download(insects_info, '/contaminant_databases/insects_db.txt')
	download(fungi_info, '/contaminant_databases/fungi_db.txt')
	'''

	p_bacteria = multiprocessing.Process(target=download, args=(bacteria_info,'/contaminant_databases/bacteria_db.txt',bacteria_child,))
	p_fungi = multiprocessing.Process(target=download, args=(fungi_info,'/contaminant_databases/fungi_db.txt',fungi_child,))
	p_insects = multiprocessing.Process(target=download, args=(insects_info,'/contaminant_databases/insects_db.txt',insects_child,))
	
	p_bacteria.start()
	p_fungi.start()
	p_insects.start()

	
	bacteria_status = bacteria_parent.recv()
	fungi_status = fungi_parent.recv()
	insects_status = insects_parent.recv()

	while bacteria_status < 100 or fungi_status < 100 or insects_status < 100:
		if bacteria_status < 100:
			bacteria_status = bacteria_parent.recv()
		if fungi_status < 100:
			fungi_status = fungi_parent.recv()
		if insects_status < 100:
			insects_status = insects_parent.recv()
		'''
		sys.stdout.write("Bacteria:\t"+str(bacteria_status)+"")
		sys.stdout.write("Fungi:\t"+str(fungi_status)+"")
		sys.stdout.write("Insects:\t"+str(insects_status)+"")
		'''
		sys.stdout.write("Percent Completion:\t" + str(round((bacteria_status + fungi_status + insects_status)/3,2)) + "%\r")
		sys.stdout.flush()
	p_bacteria.join()
	p_fungi.join()
	p_insects.join()

	print ("Downloads Complete")

