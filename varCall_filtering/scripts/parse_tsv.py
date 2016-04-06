# stdin:	tsv annotated file
# stdout:	formatted columns


def convert(input_str,case_name):
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

