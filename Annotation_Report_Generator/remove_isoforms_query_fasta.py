import csv
import os
import sys


if __name__ == '__main__':
	arguments_list = sys.argv
	
	print (arguments_list)
	file_name = arguments_list[1]
	write_file = open(file_name[:-6] + "_new.fasta","w")

	with open(file_name, "r") as file:
		for line in file:
			if ">" in line:
				temp = line
				temp = temp[1:]
				temp = temp.split(" ")[0]
				temp = temp.split("|")[0]
				write_file.write(">" + temp + "\n")
			else:
				write_file.write(line)






