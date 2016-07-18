#!/usr/bin/env python
import sys
import re
import argparse
import os

# this function gets a list, and makes it lowercase	
def list_lower(l):
	y = [0] * len(l)
	for i, item in enumerate(l):
		if item[0:2] == "db":
			y[i] = item
		else:
			y[i] = item.lower()
	return y
def get_arg():
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'--debug',
		default = False,
		action = 'store_true',
		help = "Debug/verbose flag. Default: False"                
	)
	parser.add_argument(
		'--inputdir',
		help = "input"
	)
	parser.add_argument(
		'--output',
		help = "output"
	)
	parser.add_argument(
        '--case_name',
        help = "case name"
        )
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	args = parser.parse_args()
	# add key-value pairs to the args dict
	vars(args)['cwd'] = os.getcwd()
	
	return args

def to_tsv(input_str,case_name):
	isfirst = 1
	myheader = []
	myfields = []
	output = ""
	input_str = input_str.split("\n")
	
	for line in input_str:
		if line:
			line = line.replace("\n", "")
			if (isfirst):
				myheader = line.split("\t")
				isfirst = 0;
				# print CASE header
				output += "case\t"
				for field in myheader:
					if (field == "is"):
						output += "is_1\tis_2\t"
					elif (field == "dp4"):
						output += "dp4_1\tdp4_2\tdp4_3\tdp4_4\t"
					elif (field == "pv4"):
						output += "pv4_1\tpv4_2\tpv4_3\tpv4_4\t"
					elif (field == "eff"):
						output += "effect\teffect_impact\tfunctional_class\tcodon_change\tamino_acid_change\tamino_acid_length\tgene_name\ttranscript_biotype\tgene_coding\ttranscript_id\texon_rank\tgenotype_number\twarnings\t"
					else:
						output += field + "\t"
				output += "\n"
			else:
				# print case name
				output += case_name + "\t"
				myfields = line.split("\t") 
				for i,field in enumerate(myfields):
					if (myheader[i] == "is"):
						if (field == "."):
							output += ".\t.\t"
						else:
							output += field.replace(",", "\t") + "\t"
					elif (myheader[i] == "dp4"):
						if (field == "."):
							output += ".\t.\t.\t.\t"
						else:
							output += field.replace(",", "\t") + "\t"
					elif (myheader[i] == "pv4"):
						if (field == "."):
							output += ".\t.\t.\t.\t"
						else:
							output += field.replace(",", "\t") + "\t"
					elif (myheader[i] == "eff"):
						field = field.replace("||","|.|")
						field = field.replace("||","|.|")
						field = field.replace("||","|.|")
						if (field == "."):
							output += ".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t"
						elif ("WARNING" in field):
							field = field.replace(")", "")
							field = field.replace("(", "|")
							output += field.replace("|", "\t") + "\t"
						else:
							field = field.replace(")", "")
							field = field.replace("(", "|")
							output += field.replace("|", "\t") + "\t.\t"
					else:
						output += field + "\t"
				output += "\n"
	return output

def to_report(vcf_input):
# dict for entries in INFO field (all keys in the vcf get an entry)
	#output string
	output = ""
	d_info = {}
	
	# list for entries in FORMAT field
	#format_field = []
	# list for entries in header
	header=[]
	# number of samples in vcf file
	#number_samples = 1 
	
	# boolean isfirst? if yes, we print the header once and turn off
	isfirst = 1
	
	# read file
	contents = str(vcf_input).split("\n")
	
	# loop thro the file a second time
	for line in contents: 
		# double # header 
		if ( line[0:2] == "##" ):
			# if argument, dont print double # header
			if (len(sys.argv) == 1): output +=line
	
			# grab INFO subfield ID
			match = re.search(r'##INFO=<ID=(\w+)(.*)', line)
			# define d_info
			if match: 
				d_info[match.group(1)] = "."
	
		# single # header
		elif ( line[0:1] == "#" ):
			header = line.split("\t")
			# get the number of samples
			#number_samples = len(header) - 9
		# non-empty line
		elif ( line ):
			# if isfirst bool, print the header
			if (isfirst):
				# get list of elts in format field
				#format_field = line.split("\t")[8].split(":")
	
				# print first part of header
				output += "\t".join(list_lower(header[0:7])) + "\t"
				
				# print indel
				output +="indel" + "\t"
	
				# print INFO part of header -  sorted keys
				output +="\t".join(list_lower(sorted(d_info))) + "\t"
	
				# print FORMAT part of header
				#for i in range(1,1+number_samples): 
				#	myjoiner = "_" + str(i) + "\t"
				#	output +=myjoiner.join(list_lower(format_field)) + "_" + str(i) + "\t"
	
				output +="\n"
	
				# turn off flag
				isfirst = 0
	
			# reset dict
			for x in d_info: d_info[x]="."
	
			# split line in tab
			linelist = line.split("\t")
	
			# print first part of line
			output += "\t".join(linelist[0:7]) + "\t"
			
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
			output +=indel + "\t"
			output +="\t".join([d_info[x] for x in sorted(d_info)])
			output +="\t"
	
			# print FORMAT part of header
			#for i in range(9,9+number_samples): 
			#	output +="\t".join(linelist[i].split(":")) + "\t"
	
			output +="\n"
	return output

def main():
	args = get_arg()
	case_name = args.case_name
	print("[Converting " + case_name +" to .tsv]")
	inputfile = open(args.inputdir+"/"+case_name+".vcf","r")
	case_file = inputfile.read()
	inputfile.close()
	#find GERP++ and replace with GERP
	case_file = case_file.replace('GERP++','GERP')
	#vcf2report
	case_file = to_report(case_file)
	#parse_tsv case_name
	case_file = to_tsv(case_file,case_name)
	#get rid of # and ' characters
	case_file = case_file.replace("#","")
	case_file = case_file.replace("'","")
	#write to new file
	outputfile = open(args.output+"/filter/"+case_name+"_not_filt.tsv","w")
	outputfile.write(case_file)
	outputfile.close()

if __name__ == "__main__":
	main()
