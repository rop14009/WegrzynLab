import os
import time
import csv


def write_file(element):
	if not os.path.exists("fungi_db.txt"): # if file doesnt exist create it
		file = open("fungi_db.txt", "w")
		file.write(element+"\n")
		file.close()
	else:
		file = open("fungi_db.txt", "a")
		file.write(element+"\n")
		file.close()



if __name__ == '__main__':
	print()
	counter = 0
	with open("fungi.csv", "r") as file:
		file_tsv = csv.reader(file, delimiter='\t')
		try:
			print("testestst".encode("utf-8"))
			for line in file_tsv:
				newstring = "".join(line).encode("utf-8")
				newstring = newstring.decode("ascii","ignore")
	
				if newstring.find("<i>") >= 0 and newstring.find("</i>") >= 0:
					newstring = newstring[newstring.find("<i>")+3:] # cut off first half
					newstring = newstring[:newstring.find("</i>")] # cut off second half
					write_file(newstring)
					counter += 1
				
		except UnicodeDecodeError:
			print("Unicode decode error")
	
	print(counter)
	
	
	