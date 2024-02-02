####################################################################################### 
#
# Remove LTR candidates without a protein match (filter from LTRDdigest gff3 output) and 
# returns the filteret gff and a fasta files with the complete LTR
#
# Usage: LTRdigest_parse.py -f <complete.fasta> -g <ltrdigest.gff3> -o <out-index>
#
# Where:
# complete.fasta = the fasta file produced by LTRDigest containing all LTR sequences 
# 					complete in the same order as the sorted input gff
# ltrdigest.gff3 = output from LTRDigest in GFF 3 format
# output = Name index for the output fasta and gff files
#
# 
#######################################################################################

#!/usr/bin/python

import getopt
import sys
import os
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

ltr_seq = tuple()
genes_seq = dict()
duplicated = list()

fasta = ""
gff = ""

outgff = ""
tempout = ""
outname_temp = ""
fastaout = ""

seq_ID_old = ""
seq_id = ""

new_seq = False
first_line = True

count_ltr = 0

# Check for the arguments and print useful help messages

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:g:o:",["fasta=","gff=","output="])
except getopt.GetoptError:
    print '\n', '####     Invalid use     ####', '\n'
    print 'Usage: LTRdigest_parse.py -f <complete.fasta> -g <ltrdigest.gff3> -o <out-index>'
    print 'For help use LTRdigest_parse.py -h'
    sys.exit(99)

for opt, arg in opts:
    if opt == '-h':
        print ('\nRemove LTR candidates without a protein match (filter fro LTRDdigest gff3 output) and returns the filteret gff and a fasta files with the complete LTR\n')
        print ('Usage: LTRdigest_parse.py -f <complete.fasta> -g <ltrdigest.gff3> -o <out-index>\n')
        print ('Where: complete.fasta = the fasta file produced by LTRDigest containing all LTR sequences complete in the same order as the sorted input gff')
        print ('ltrdigest.gff3 = output from LTRDigest in GFF 3 format')
        print ('output = Name index for the output fasta and gff files\n')
        sys.exit()
    elif opt in ("-f", "--fasta"):
        fasta = open(arg)
    elif opt in ("-g", "--gff"):
        gff_file = arg
	gff = open(arg)
    elif opt in ("-o", "--output"):
        outname_temp = arg + ".tmp"
        temp_file = open(outname_temp,"w")
        outname_fasta = arg + ".fasta"
        fastaout = open(outname_fasta,"w")
        outname_gff = arg + ".gff3"
        outgff = open(outname_gff,"w")
    else:
        assert False, "unhandled option"
        sys.exit(99)



## Open the fasta file and save the sequences in numeric order because it corresponds to the ID in Parent=repeat_regionID
for ltr_record in SeqIO.parse(fasta, "fasta"):
    ltr_header = str(ltr_record.description)
    ltr_seq = ltr_seq + (ltr_header,)
    ## to print the output later
    bases = (ltr_record.seq)
    genes_seq[ltr_header]= bases


tempout = open(outname_temp, "w") # start writing in a temp output with all entries with a protein match

## Open GFF file, filter interest sequences and print outputs
for line in gff:
	
	## Deal with comments
	if line.startswith("###"): ## seq dividers
		new_seq = True
		continue
	elif line.startswith("##") and new_seq is False : ## comment, only the first line will be kept because seq tag will be updated in GFF
		if first_line is True:
			first_line = False
			outgff.write(line)
		continue
	
	## Deal with LTR annotation
	else: 
		line = line.rstrip()
		itens = line.split("\t")
		feature = str(itens[2:3]) # It will be 'protein_match' if it was annotated with LTRDigest, sequences with that is what we want to keep
		seq_ID = str(itens[0]) # sequence ID tag, internal from LTRdigest
        
        if "protein_match" in feature: # check if had a protein match
            tempout.write(line + "\n")

tempout.close()

######### Writing outputs

### In the fasta file
temp_file = open(outname_temp)
for entry in temp_file:
	entry = entry.rstrip()
	entry = entry.split("\t")
				
	# Get the Parent=repeat_regionID number to find the sequence in the fasta file
	entry_ID = str(entry[len(entry)-1]) # last entry
	entry_ID = entry_ID.split(";")[0]
        entry_number = entry_ID.split("retrotransposon")
	entry_number = int(entry_number[1])
        seq_id = ltr_seq[entry_number-1] # number -1 to index properly (starts at 0)
        if seq_id not in  duplicated:
            duplicated.append(seq_id)
	    count_ltr = count_ltr+1
    
            ## write the fasta output
	    better_seq_header = seq_id.split("_")
	    better_seq_id = str(better_seq_header[0] + "_" + better_seq_header[1])
	    pos = str("repeat_region" + str(entry_number) + "[" + better_seq_header[2] + "," + better_seq_header[3] + "]")
	    fasta_format_string = SeqRecord(genes_seq[seq_id], id=better_seq_id, description=pos)
	    fastaout.write(fasta_format_string.format("fasta"))

fastaout.close()


## In the gff file
temp_file = open(outname_temp)
entry_list = list()
for entry in temp_file:
	entry = entry.rstrip()
	entry = entry.split("\t")
    ## Get the Parent=repeat_regionID number to find the sequence in the fasta file
	entry_ID = str(entry[len(entry)-1]) # last entry
	entry_ID = entry_ID.split(";")[0]
    	entry_number = entry_ID.split("retrotransposon")
	entry_number = entry_number[1]
   	if entry_number not in entry_list:
		entry_list.append(entry_number)

 		## write new entry separator in gff
    		outgff.write("###\n")
    		gff = open(gff_file)
    
    	### Check all annotations for the LTR selected
    		for line in gff:
        		line = line.split("\t")
			description = line[-1]
			description = description.replace("=", " ")
			description = description.replace(";", " ")

			repeat_ID = str("repeat_region" + entry_number)
			ltr_ID = str("LTR_retrotransposon" + entry_number)
		
			if re.search(r"\b"+repeat_ID+r"\b", description) or re.search(r"\b"+ltr_ID+r"\b", description):
            			seq_id = ltr_seq[int(entry_number)-1] # number -1 to index properly (starts at 0)
                        	better_seq_header = seq_id.split("_")
                        	better_seq_id = str(better_seq_header[0] + "_" + better_seq_header[1])
                        	better_entry = line
                        	better_entry[0] = better_seq_id  ## update seq id
                        	better_entry = "\t".join(better_entry)
                        	outgff.write(better_entry)

outgff.close()

## remove temp file
os.remove(outname_temp)


## Useful info
useful_info = "\nThere are in total " + str(count_ltr) + " LTRretrotransposon annotated.\n"
print (useful_info)
