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
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

ltr_seq = tuple()
genes_seq = dict()

fasta = ""
gff = ""

outgff = ""
tempout = ""
outname_temp = ""
fastaout = ""

seq_ID_old = ""
seq_id = ""

new_seq = False
include_seq = False
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
		
######## Writing outputs
		### Print sequence entry in outputs they had been filtered
		if new_seq is True and include_seq is True: ## if new_seq is True is because all annotations from last sequence entry was analysed, so it should be save if filtered
			include_seq = False #reset check variable
			count_ltr = count_ltr+1 
			n = 1
			tempout.close() #close temp out
			outgff.write("###\n")
			#save all notation for this entry
			temp_file = open(outname_temp)
			for entry in temp_file:
				entry = entry.rstrip()
				entry = entry.split("\t")
				
				## print the fasta output for the entire seq, which is the first entry annotation/ position
				if n == 1:
					n = 2
					# Get the Parent=repeat_regionID number to find the sequence in the fasta file
					entry_ID = str(entry[len(entry)-1]) # last entry
					entry_number = entry_ID.split("region")
					entry_number = int(entry_number[1])
					seq_id = ltr_seq[entry_number-1] # number -1 to index properly (starts at 0)
					## write the fasta output
					better_seq_header = seq_id.split("_")
					better_seq_id = better_seq_header[0]
					pos = str(entry_ID + " " + "[" + better_seq_header[1] + "," + better_seq_header[2] + "]")
					fasta_format_string = SeqRecord(genes_seq[seq_id], id=better_seq_id, description=pos)
					fastaout.write(fasta_format_string.format("fasta"))
				
				## print the gff output, but better because the original has a especial sequence tag and it's cofunsing.
				# this one will have the original sequence ID (from fast input in LTRHarvest)
				seq_tag = seq_id.split("_")
				better_entry = entry
				better_entry[0] = seq_tag[0] # update seq name
				better_entry = "\t".join(better_entry)
				outgff.write(better_entry + "\n")							
		
		## Just don't save last entry because it was not annotated		
		if new_seq is True and include_seq is False: 
			try:
				tempout.close() #close temp out
			except:
				pass
			
########
		
		### Filter LTR with annotation
		## Save a new sequence entry in the temp output
		if (seq_ID != seq_ID_old) and new_seq is True: # If it is a new sequence
			new_seq = False
			seq_ID_old = seq_ID
			tempout = open(outname_temp, "w") # start writing in a temp output with all entries for this seqID
			tempout.write(line + "\n")
			if "protein_match" in feature: # check if had a protein match
				include_seq = True
			continue
		# Continue saving the same sequence entry in the temp output
		elif (seq_ID == seq_ID_old) and new_seq is False:
			tempout.write(line + "\n")
			if "protein_match" in feature: # check if had a protein match
				include_seq = True


## close output files
outgff.close()
fastaout.close()	

## remove temp file
os.remove(outname_temp)			


## Useful info
useful_info = "\nThere are in total" + count_ltr + "LTRretrotransposon annotated.\n")
print ("useful_info)

	
