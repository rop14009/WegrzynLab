import urllib2
import re

def download(url_info):
        base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
        final_url = base_url + "&WebEnv=" + url_info[0] + "&query_key=" + str(url_info[1]) + "&rettype=fasta&retmode=text"
	
	description_text = "null"
	species_text = "null"
	nucleotides = "null"


	try:
		response = urllib2.urlopen(final_url)
		nucleotides = response.read()
		nucleotides = nucleotides.split("\n")
		del nucleotides[0]
		nucleotides = ''.join(nucleotides)
		

		response_genpept = urllib2.urlopen(base_url +  "&WebEnv=" + url_info[0] + "&query_key=" + str(url_info[1]) + "&rettype=genpept&retmode=text")
		response_gen_pep = response_genpept.read()
		
		#print (response_gen_pep)
		description = re.compile('title\s"(.*?)\"', re.DOTALL)
		species = re.compile('taxname\s"(.*?)\"', re.DOTALL)
		description_text = description.search(response_gen_pep).group()#[7:-1]
		species_text = species.search(response_gen_pep).group()[9:-1]
		
		
		#print (response_genpept.read())
	except:
		print ("unable to get response -- net error")
	
	return [str(nucleotides),description_text,species_text]

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
		print ("error")

	web_env_re = re.compile('<WebEnv>(\S+)<\/WebEnv>')
	web_env = web_env_re.search(search_result).group()[8:-9]
	#print ("web_env\t"+web_env)
	key_re = re.compile('<QueryKey>(\S+)<\/QueryKey>')
	key = int(key_re.search(search_result).group()[10:-11])
	#print ("key:\t"+str(key))

        return [web_env, key]

# return peptide sequence, description, species
def sync_repair(gi_id):
	#print (gi_id)
	data = get_ncbi_protein_query_data(gi_id)
	#print (data)
	return download(data)

#sync_repair("228412")
#print (sync_repair("255546175"))

