import urllib2
import re

# download fasta format attempt
def download_nucleotides(url_info, num_retry):
        base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
        final_url = base_url + "&WebEnv=" + url_info[0] + "&query_key=" + str(url_info[1]) + "&rettype=fasta&retmode=text"
	nucleotides = "null"
	
	if num_retry == 0:
		print ("All retries for download_nucleoties failed")
		exit()

	try:
		response = urllib2.urlopen(final_url, timeout=15)
		nucleotides = response.read()
		#print (nucleotides)

		nucleotides = nucleotides.split("\n")
		#print (nucleotides[0])
		header = str(nucleotides[0])
		del nucleotides[0]
		nucleotides = ''.join(nucleotides)
	except:
		print ("unable to get response -- net error")
		return download_nucleotides(url_info,num_retry-1)
	#print ([str(nucleotides),description_text,species_text])	
	return str(nucleotides)


def download(gi, url_info, num_retry):
        description_text = "null"
        species_text = "null"
        nucleotides = "null"
	
	if num_retry == 0:
		print ("All retries for download failed")
		exit()

        try:
		response2 = urllib2.urlopen("http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=protein&dopt=genpept&sort=&val="+gi+"&from=begin&to=end&extrafeat=984&maxplex=3", timeout=15)
		r2 = response2.read()
		#print ("http://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=protein&dopt=genpept&sort=&val="+gi+"&from=begin&to=end&extrafeat=984&maxplex=3")

		species = re.compile('ORGANISM  (.*?)\n')
		species_text = species.search(r2).group()[10:-1]
		#print (species_text)
		description = re.compile("DEFINITION  (.*?)ACCESSION", re.DOTALL)
		description_text = description.search(r2).group()[12:-10]
		description_text = description_text.replace("            ", " ")
		description_text = description_text.replace("\n","")	

	except Exception as e:
                print ("unable to get response -- net error")
		print (str(gi) + "\t:"  + str(e))
		print ("desc:\t"+str(description_text))
		return download(gi,url_info, num_retry-1)
        #print ([str(nucleotides),description_text,species_text])
        return [description_text,species_text]



def query_retry_loop(query):
	print ("Now starting retry attempts...")
	result = None
	for attempts in range(1,15): # retry 15 times
		if result is None:
			result = get_ncbi_protein_backup(query, attempts)		
	if result is None:
		print ("All retry attempts for ncbi query data failed")
	else:
		print ("Retry attempts succeeded! Now continuing as usual ...")
	return result


def get_ncbi_protein_backup(query, attempt_num):
	base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        query_url = "esearch.fcgi?db=protein&term=" + query + "&usehistory=y"
        total_url = base_url + query_url
        #print (total_url)
        web_env = ""
        key = 0
        count = 0 # number of terms in results
        try:
                response = urllib2.urlopen(total_url, timeout=15)
                #print (response.info())
                search_result = response.read()
                #print (search_result)
		web_env_re = re.compile('<WebEnv>(\S+)<\/WebEnv>')
                web_env = web_env_re.search(search_result).group()[8:-9]
                #print ("web_env\t"+web_env)
                key_re = re.compile('<QueryKey>(\S+)<\/QueryKey>')
                key = int(key_re.search(search_result).group()[10:-11])
                #print ("key:\t"+str(key))
		return [web_env, key]
        except:
		print ("retry #" + str(attempt_num) + " failed")
		return None


def get_ncbi_protein_query_data(query):
        base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        query_url = "esearch.fcgi?db=protein&term=" + query + "&usehistory=y"
        total_url = base_url + query_url
	#print (total_url)
        web_env = ""
        key = 0
        count = 0 # number of terms in results
	try:
		response = urllib2.urlopen(total_url, timeout=15)
		#print (response.info())
		search_result = response.read()
		#print (search_result)
	except:
		print ("error -- unable to obtain get ncbi protein query data")
		return query_retry_loop(query)

	try:
		web_env_re = re.compile('<WebEnv>(\S+)<\/WebEnv>')
		web_env = web_env_re.search(search_result).group()[8:-9]
		#print ("web_env\t"+web_env)
		key_re = re.compile('<QueryKey>(\S+)<\/QueryKey>')
		key = int(key_re.search(search_result).group()[10:-11])
		#print ("key:\t"+str(key))
	except:
		print ("error in regex")
		print (search_result)
		return query_retry_loop(query)		

        return [web_env, key]

# return peptide sequence, description, species
def sync_repair(gi_id):
	#print (gi_id)
	data = get_ncbi_protein_query_data(gi_id)
	#print (data)

	nucleotides = download_nucleotides(data, 15)
	[desc, species] = download(gi_id, data, 15)

	return [nucleotides, desc, species]

"""
print (sync_repair("228412"))
print (sync_repair("12321310"))
print (sync_repair("15892919"))
print (sync_repair("116492625"))
print (sync_repair("334187585"))
print (sync_repair("74745129"))
print (sync_repair("364506258"))
print (sync_repair("675969702"))
print (sync_repair("491668487"))
"""

