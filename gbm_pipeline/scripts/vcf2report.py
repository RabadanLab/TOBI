#!/usr/bin/env python

import sys
import re

# this function gets a list, and makes it lowercase	
def list_lower(l):
	y = [0] * len(l)
	for i, item in enumerate(l):
		if item[0:2] == "db":
			y[i] = item
		else:
			y[i] = item.lower()
	return y



# dict for entries in INFO field (all keys in the vcf get an entry)
d_info = {}

# list for entries in FORMAT field
format = []
# list for entries in header
header=[]
# number of samples in vcf file
number_samples = 1 

# boolean isfirst? if yes, we print the header once and turn off
isfirst = 1
# boolean verbose header
# verboseheader = 1  

# read file
contents = sys.stdin.read().split("\n")

# the idea: get all the possible keys in d_info with one run thro the file
# then run thro the file again and use d_info
# going thro the file twice is perhaps not the most efficient way, 
# but the final report should be small
# for line in contents: 
# 
# 	# skip header lines with double # 
# 	if ( line[0:2] == "##" ):
# 		pass;
# 	# header
# 	elif ( line[0:1] == "#" ):
# 		header = line.split("\t")
# 		# get the number of samples
# 		number_samples = len(header) - 9
# 	# non-empty line
# 	elif ( line ):
# 		# for elts in INFO field, make a dict
# 		for x in line.split("\t")[7].split(";"): 
# 			# pattern must be: blob=blog or we wont consider
# 			if (re.search(r'(\w+)=(\w+)', x)):
# 				d_info[x.split("=")[0]] = "."
# 		format = line.split("\t")[8].split(":")

# loop thro the file a second time
for line in contents: 

	# double # header 
	if ( line[0:2] == "##" ):
		# if argument, dont print double # header
		if (len(sys.argv) == 1): print(line)

		# grab INFO subfield ID
		match = re.search(r'##INFO=<ID=(\w+)(.*)', line)
		# define d_info
		if match: 
			d_info[match.group(1)] = "."

	# single # header
	elif ( line[0:1] == "#" ):
		header = line.split("\t")
		# get the number of samples
#		number_samples = len(header) - 9
	# non-empty line
	elif ( line ):
		# if isfirst bool, print the header
		if (isfirst):
			# get list of elts in format field
			format = line.split("\t")[8].split(":")

			# print first part of header
			print("\t".join(list_lower(header[0:7])) + "\t"),
			
			# print indel
			print("indel" + "\t"),

			# print INFO part of header -  sorted keys
			print("\t".join(list_lower(sorted(d_info))) + "\t"),

			# print FORMAT part of header
			for i in range(1,1+number_samples): 
				myjoiner = "_" + str(i) + "\t"
				print(myjoiner.join(list_lower(format)) + "_" + str(i) + "\t"),

			print("\n"),

			# turn off flag
			isfirst = 0

		# reset dict
		for x in d_info: d_info[x]="."

		# split line in tab
		linelist = line.split("\t")

		# print first part of line
		print("\t".join(linelist[0:7]) + "\t"),
		
		# add indel
		indel = "."
		
		# loop thro INFO keys - if key found, add val
		for x in linelist[7].split(";"):
			# pattern must be: blob=blog or we wont consider
			if (x == "INDEL"):
				indel = "1"
			elif (re.search(r'(\w+)=(\w+)', x)):
				d_info[x.split("=")[0]] = x.split("=")[1]
		
		# print INFO part of line - sorted values
		print(indel + "\t"),
		print("\t".join([d_info[x] for x in sorted(d_info)])),
		print("\t"),

		# print FORMAT part of header
		for i in range(9,9+number_samples): 
			print("\t".join(linelist[i].split(":")) + "\t"),

		print("\n"),
	