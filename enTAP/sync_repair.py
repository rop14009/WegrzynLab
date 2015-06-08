import urllib2
import re

# download fasta format attempt
def download_nucleotides(url_info):
        base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
        final_url = base_url + "&WebEnv=" + url_info[0] + "&query_key=" + str(url_info[1]) + "&rettype=fasta&retmode=text"

	nucleotides = "null"


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
	#print ([str(nucleotides),description_text,species_text])	
	return str(nucleotides)


def download(gi):
        description_text = "null"
        species_text = "null"
        nucleotides = "null"


        try:
                response = urllib2.urlopen("http://www.ncbi.nlm.nih.gov/protein/"+gi+"?report=genpept", timeout=15)
                r = response.read()
		header = re.compile('<div class="rprtheader">\s\s\s\<h1\>(.*?)\<\/', re.DOTALL)
		header_text = header.search(r).group()[31:-2]
		#print (header_text)
	
		species = re.compile('\s\[(.*?)\]')
		species_text = species.search(header_text).group()[2:-1]
	
		description = re.compile('(.*?)\s\[')
		description_text = description.search(header_text).group()[:-2]


	except:
                print ("unable to get response -- net error")
        #print ([str(nucleotides),description_text,species_text])
        return [description_text,species_text]




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

	nucleotides = download_nucleotides(data)
	[desc, species] = download(gi_id)

	return [nucleotides, desc, species]

print (sync_repair("228412"))
print (sync_repair("12321310"))
print (sync_repair("15892919"))
print (sync_repair("116492625"))
print (sync_repair("356517118"))
print (sync_repair("334187585"))
