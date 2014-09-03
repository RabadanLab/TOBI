#!/usr/bin/env python

# stdin:	tsv annotated file
# stdout:	formatted columns

import sys
import re

isfirst = 1
myheader = []
myfields = []

args = sys.argv[1:]
if args:
	case_name = args[0]

for line in sys.stdin:
	line = line.replace("\n", "")
	if (isfirst):
		myheader = line.split("\t")
		isfirst = 0;
		# print CASE header
		if args:
			print("case\t"),
		for field in myheader:
			if (field == "is"):
				print("is_1\tis_2\t"),
			elif (field == "dp4"):
				print("dp4_1\tdp4_2\tdp4_3\tdp4_4\t"),
			elif (field == "pv4"):
				print("pv4_1\tpv4_2\tpv4_3\tpv4_4\t"),
			elif (field == "eff"):
				print("effect\teffect_impact\tfunctional_class\tcodon_change\tamino_acid_change\tamino_acid_length\tgene_name\ttranscript_biotype\tgene_coding\ttranscript_id\texon_rank\tgenotype_number\twarnings\t"),
			else:
				print(field + "\t"),
		print("\n"),
	else:
		# print case name
		if args:
			print(case_name + "\t"),
		myfields = line.split("\t") 
		for i,field in enumerate(myfields):
			if (myheader[i] == "is"):
				if (field == "."):
					print(".\t.\t"),
				else:
					print(field.replace(",", "\t") + "\t"),
			elif (myheader[i] == "dp4"):
				if (field == "."):
					print(".\t.\t.\t.\t"),
				else:
					print(field.replace(",", "\t") + "\t"),
			elif (myheader[i] == "pv4"):
				if (field == "."):
					print(".\t.\t.\t.\t"),
				else:
					print(field.replace(",", "\t") + "\t"),
			elif (myheader[i] == "eff"):
				field = field.replace("||","|.|")
				field = field.replace("||","|.|")
				field = field.replace("||","|.|")
				if (field == "."):
					print(".\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t"),
				elif ("WARNING" in field):
					field = field.replace(")", "")
					field = field.replace("(", "|")
					print(field.replace("|", "\t") + "\t"),
				else:
					field = field.replace(")", "")
					field = field.replace("(", "|")
					print(field.replace("|", "\t") + "\t.\t"),
			else:
				print(field + "\t"),
		print("\n"),


