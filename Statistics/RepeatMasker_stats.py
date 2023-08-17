#!/usr/local/bin/python


"""
Summary stats from all TE annotations based on RepeatMasker '.out' output file 
        
Usage = RepeatMasker_stats.py -G <genome_size(bp)> -i <RepeatMasker.out> -o <output_prefix>

Where: 
RepeatMasker.out = '.out' file output from RepeatMasker containing all repeats annotated in the genome
output_prefix = prefix name for the output file

Options: 
-h for usage help

"""

import sys, getopt, re 

# Check for the arguments, open the inputs and print useful help messages

try:
    opts, args = getopt.getopt(sys.argv[1:],"hG:i:o:",["output"])
except getopt.GetoptError:
    print ('\n', '####     Invalid use     ####', '\n')
    print ('Usage = RepeatMasker_stats.py <Options> -i <RepeatMasker.out> -o <output_prefix>')
    print ('For help use RepeatMasker_stats.py -h')
    sys.exit(99)
    
for opt, arg in opts:
    if opt == '-h':
        print ('\n', 'Summary stats from all TE annotations based on RepeatMasker \'.out\' output file.', '\n')
        print ('Usage = RepeatMasker_stats.py -G <genome_size(bp)> -i <RepeatMasker.out> -o <output_prefix>')
        print ('Where:\nRepeatMasker.out = \'.out\' file output from RepeatMasker containing all repeats annotated in the genome')
        print ('output_prefix = prefix name for the output file.')
        sys.exit()
    elif len(opts) == 3:
        if opt in ("-i"):
            TE_file = open(arg)
            # save the last line
            # lines = open(arg).readlines()
            # last = lines[-1]
        if opt in ("-G"):
            genome_size = int(arg)
        if opt in ("-o", "--output"):
            n_out1 = arg + "_size.txt"
            out1 = open(n_out1, "w")
            n_out2 = arg + "_stats.txt"
            out2 = open(n_out2, "w")
    elif len(opts) < 3:
        print ('\n', '###    Argument missing   ###', '\n', '\n' 'Use -h option for help\n')
        sys.exit(1)
    else:
        assert False, "unhandled option"
        
        
## Create variables
TE_size = {}
group_size = {}
n_lines = 0
query = ''
beg = int()
end = int()
TE_cat = ''

header = True


## Output headers
out1.write("class/family\t" + "total_size\(bp)\n")

## Get the total size of all matches
for line in TE_file:    
    # skip header
    if header is True:
        if n_lines < 2:
            n_lines = n_lines + 1
            continue
        else:
            header = False
            continue    
    
    # save the info needed for stats
    else:        
        # format line
        line = re.sub(r"^\s+", "", line, flags=re.UNICODE)
        line = " ".join(re.split("\s+", line, flags=re.UNICODE))
        line = line.split(" ")
        query = str(line[4])
        beg = int(line[5])
        end = int(line[6])
        TE_cat = str(line[10])
                
        # TE size
        total_l = end - beg
        
        # Save in out1
        out1.write(TE_cat + "\t" + str(total_l) + "\t" + "\n")
    
        ## Add totals to dictionary
        if TE_cat in TE_size:
            TE_size[TE_cat] = TE_size[TE_cat] + total_l
            
        if TE_cat not in TE_size:
            TE_size[TE_cat] = total_l
                


## Get the group totals
n_interspersed = ("Low_complexity","Simple_repeat","Satellite","rRNA","snRNA")
#rolling_circles = "RC/Helitron"
#unclassified = "Unknown"
#class1 (SINE, LTR, LINE)
#class2 (DNA)

group_size = {'non interspersed':0,
             'rolling circles':0,
             'unclassified': 0,
             'class1':0,
             'class2':0}

for key, value in TE_size.items():
    if key in n_interspersed:
        group_size['non interspersed'] = group_size['non interspersed'] + value
    elif key == "RC/Helitron":
        group_size['rolling circles'] = group_size['rolling circles'] + value
    elif key == "Unknown":
        group_size['unclassified'] = group_size['unclassified'] + value
    elif key.startswith("DNA"):
        group_size['class2'] = group_size['class2'] + value
    elif key.startswith("SINE") or key.startswith("LTR") or key.startswith("LINE"):
        group_size['class1'] = group_size['class1'] + value
        
        
## Print stats output
out2.write("\nGenome Size = " + str(genome_size) +  "\n")
out2.write("TE (bp) = " + str(sum(group_size.values())) +  "\n")
out2.write("TE % = " + str(round((sum(group_size.values())*100)/genome_size,2)) +  "\n\n")
out2.write("====================== Class totals ============================\n")

# Print per group totals
out2.write("class\t" + "total_size (bp)\t" + "genome (%)\n")
for key, value in sorted(group_size.items()):
    p_genome = (value*100) / genome_size
    out2.write(key + "\t" + str(value) + "\t" + str(round(p_genome,2)) + "%" + "\n")

    
out2.write("\n====================== Class/Family totals ============================\n")
# Print list with all
out2.write("class\t" + "total_size (bp)\t" + "genome (%)\n")
for key, value in sorted(TE_size.items()):
    p_genome = (value*100) /genome_size
    out2.write(key + "\t" + str(value) + "\t" + str(round(p_genome,2)) + "%" + "\n")
    

out1.close()
out2.close()
