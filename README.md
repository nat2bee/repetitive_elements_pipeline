# Pipeline to identify repetitive elements on non-model species genomes.

This pipeline is adapted from `Avrilomics` posts (http://avrilomics.blogspot.com/2015/09/creating-merged-repeatmodeler-and.html) and from `Matthew Berriman, Avril Coghlan, Isheng Jason Tsai et al.` (Creation of a comprehensive repeat library for a newly sequenced parasitic worm genome, 15 May 2018, PROTOCOL [Version 1] available at Protocol Exchange [+https://doi.org/10.1038/protex.2018.054+]). The pipeline reasoning is described in more detail there. 

Herein is provided a detailed and update version (Dez 2023) of this protocol using a bee species (*Tetrapedia diversipes*) as example.

## Step by Step

### 1- Create multiple repeat libraries for your genome

When we look for repeats in non-model species using only one approach may be problematic because available databases are restrictive. There are multiple programs to identify repeatitive elements in the genome, some of them are more general and some of them are tailored to identify certain repeat types. In this pipeline I have used the RepeatModeler<sup>1</sup> (v open-1.0.11), TransposonPSI<sup>2</sup> and LTRharvest<sup>3</sup> (included in GenomeTools v 1.5.8). So first we need to build one repeat library for each one of these software.

**PS UPDATE** *Recently (January 2024) I have updates the program versions and used the same pipeline, comments were included to incorporate this update. Tested program versions were RepeatModeler Version 2.0.5 (Dependencies: TRF 4.09, RECON , RepeatScout 1.0.6, RepeatMasker 4.1.6), TransposonPSI Version 1.0.0-3 (conda install), LTRharvest (included in GenomeTools v1.6.5 with ncbi-blast-2.15.0)*

**RepeatModeler** combine different approaches to identify the repeats, like RECON - good at identifying less conserved repetitive elements, RepeatScout - good for finding highly conserved repetitive elements, and Tandem Repeats Finder. 

To build the repeat library with RepeatModeler:
```
- First build the genome database
$ BuildDatabase -name tdiv -engine ncbi TDIV.fa

- Then run the program
$ RepeatModeler -engine ncbi -pa 11 -database tdiv >& run2.out

>> consensi.fa.classified is the final output
```
**PS UPDATE** *`-pa` was replaced by `-threads` in the updated version of RepeatModeler*

**TransposonPSI** searchs the genome for LTR retroelements (e.g. gypsy and copia), non-LTR retroelements (e.g. LINE retrotransposon ORFs), DNA transposons (e.g. cryptons and other families), and helitrons. It used PSI-BLAST to identify proteins enconded by transposable elements, thus it is useful to find degenerate elements.

To build the repeat library with TransposonPSI:
```
$ transposonPSI.pl TDIV.fa nuc

- Outputted files will be
>>
genome.fasta.TPSI.allHits : all HSPs reported by PSITBLASTN in btab format.
genome.fasta.TPSI.allHits.chains : collinear HSPs are chained together into larger chains (more complete element regions).
genome.fasta.TPSI.allHits.chains.gff3 : the chains in gff3 format
genome.fasta.TPSI.allHits.chains.bestPerLocus : a DP scan is performed, extracting the best scoring chain per genomic locus.
genome.fasta.TPSI.allHits.chains.bestPerLocus.gff3 : the best chains in gff3 format.
```

**LTRharvest** is optimized to identify LTR retrotransposons based on their characteristic structure. LTRharvest perform multiple filterings: it computes boundary positions, filters for LTR length and distance, and TSD length and motifs. Filters can be turned on/off and the pipeline can be tailored for species-specific need (if known).

To build the repeat library with LTRharvest:
```
- First generate the sufix array
$ gt suffixerator -db TDIV.fa -indexname TDIV.fsa -tis -suf -lcp -des -ssp -sds -dna

- Then run the program
$ nohup gt ltrharvest -index TDIV.fsa -v -out pred-TDIV.fsa -outinner pred-inner-TDIV.fsa -gff3 pred-TDIV.gff &

- Understand the output
`Each comment line starts with the comment symbol #. Each non-comment line denotes a LTR retrotransposon prediction with starting and ending positions of the whole LTR retrotransposon, the left LTR instance and the right LTR instance, respectively. Furthermore, for each of these elements, the corresponding element length is reported as well as a percentage similarity of the two LTRs. The last integer of each line denotes the number of the input sequence, the LTR retrotransposon prediction occurs in. The input sequence numbers are counted from 0.`
```

### 2- Improve and format libraries

**TransposonPSI**</br>
First get the repeats in fasta format based on the gff output from TransposonPSI. I have used samtools and bedtools for that.
```
- Generate the reference genome index
$ samtools faidx TDIV.fa

- Run bedtools
$ bedtools getfasta -fo TDIV.fa.TPSI.allHits.chains.bestPerLocus.fasta -s -fi TDIV.fa -bed TDIV.fa.TPSI.allHits.chains.bestPerLocus.gff3
```

Next classify the repeats with RepeatClassifier (included in RepeatModeler). 
```
- First build database
$ BuildDatabase -name TDIV.fa.TPSI.allHits.chains.bestPerLocus -engine ncbi TDIV.fa.TPSI.allHits.chains.bestPerLocus.fasta

- Then run the program
$ RepeatClassifier -engine ncbi -stockholm ../tdiv-families.stk -consensi TDIV.fa.TPSI.allHits.chains.bestPerLocus.fasta
```


Check if sequences of <50 are present in the library. Here I did it using SeqKit<sup>4</sup> (v0.11.0):
```
$ SeqKit_v0.11.0 stats TDIV.fa.TPSI.allHits.chains.bestPerLocus.fasta.classified.filtered
```
If they are, remove them from your library:
```
$ SeqKit_v0.11.0 seq -m 50 TDIV.fa.TPSI.allHits.chains.bestPerLocus.fasta.classified > TDIV.fa.TPSI.allHits.chains.bestPerLocus.fasta.classified.filtered
```

**LTRharvest**<br/>
Use LTRdigest (also included in GenomeTools v 1.5.8) to search the LTRharvest library for functional annotation in the LTRs. It search for similarities with protein HMMs at Pfam databases to make the annotation. <br/>

**PS UPDATE** *The Pfam HMM database did work with the updated version of LTRdigest because it was an older HMM format. I tried to convert the formats with no success, so I didn't use this database. The new version of the Gypsy database worked with no problem.* <br/>

I have used two HMM databses in this step:<br/>
1- From Pfam, PF03732 = Retrotrans_gag
```
$ wget http://pfam.xfam.org/family/PF03732/hmm
$ mv hmm HMM1.hmm 
```
2- From the REPET database 'ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm'
```
$ mv ProfilesBankForREPET_Pfam27.0_GypsyDB.hmm HMM2.hmm
```

Running LTRdigest with these databases
```
$ gt gff3 -sort pred-TDIV.gff > sorted_pred-TDIV.gff
$ gt ltrdigest -hmms HMM*.hmm -outfileprefix TDIV-ltrs sorted_pred-TDIV.gff TDIV.fsa > TDIV-ltrs_ltrdigest.gff3

- Outputed files will be
>>
TDIV-ltrs_ppt.fas
TDIV-ltrs_3ltr.fas
TDIV-ltrs_5ltr.fas
TDIV-ltrs_complete.fas 

- Understanding the output
`
• protein domains,
• polypurine tracts (PPT) and 
• primer binding sites (PBS)

LTRdigest computes the boundaries and attributes of the features that fit the user-supplied model and outputs them in GFF3 format [2] (in addition to the existing LTR retrotransposon annotation), as well as the corresponding sequences in multiple FASTA format. In addition, a tab-separated summary file is created that can conveniently and quickly browsed for results.

No screen output (except possible error messages) is produced since the GFF3 output on stdout is redirected to a file. Additionally, the files mygenome-ltrs conditions.csv, mygenome- -ltrs 3ltr.fas, mygenome-ltrs 5ltr.fas, mygenome-ltrs ppt.fas, mygenome- -ltrs pbs.fas, mygenome-ltrs tabout.csv and one FASTA file for each of the HMM models will be created and updated during the computation. As the files are buffered, it may take a while before first output to these files can be observed.`
```

Then filter the output from LTRDigest to keep only LTRs that matched the database, i.e. LTRs that contain functional elements.<br/>
For this I have written my on scripts 'LTRdigest_parse.py'<br/>
*Usage: LTRdigest_parse.py -f <complete.fasta> -g <ltrdigest.gff3> -o <out-index>*<br/>

**PS UPDATE** *For the new version of LTRdigest use the updated version of the script (LTRdigest_parse_new.py, also in this repository). The usage is the same but it now runs in python3 instead of python2.*

```
$ python ../LTRdigest_parse.py -f TDIV-ltrs_complete.fas -g TDIV-ltrs_ltrdigest.gff3 -o TDIV-ltrs_ltrdigest_filtered
```

***PS**: There is a trick with LTRDigest; since we use the off sorted it changed the index of the sequence tags. So for example, seq2574 from LTRharvest is not the same as seq2574 in LTRdigest output. Because of that the input must the fasta generated by LTRdigest ("TDIV-ltrs_complete.fas"), it has the correct repeat_region order number (Parent=repeat_region2994, the number represent the number of the parent fasta read in numeric order).*<br/>

**PS UPDATE** *The option `-engine ncbi` is no longer necessary in the new version of RepeatClassifier, and thus should be omitted.* <br/>

Finally, use RepeatClassifier to classify the repeats in the right format
```
- First build database
$ BuildDatabase -name TDIV-ltrs_ltrdigest_filtered -engine ncbi TDIV-ltrs_ltrdigest_filtered.fasta

- Then run the program
$ RepeatClassifier -engine ncbi -stockholm tdiv-families.stk -consensi TDIV-ltrs_ltrdigest_filtered.fasta
```

### 3- Merge all repeat libraries

Format library fasta files to have one read per line (unbreake line)
```
- RepeatModeler library
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < consensi.fa.classified > RepeatModeler_library_classified.fasta

- TransposonPSI library
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < TDIV.fa.TPSI.allHits.chains.bestPerLocus.fasta.classified > TPSI_library_classified.fasta

- LTRharvest
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < TDIV-ltrs_ltrdigest_filtered.fasta.classified > ltrharvest_library_classified.fasta
```

Merge the three libraries
```
$ cat RepeatModeler_library_classified.fasta  TPSI_library_classified.fasta ltrharvest_library_classified.fasta > merged.fa
```

Remove comments from the fasta file (necessary to run the merging program)
```
$ sed -e 's/>* .*$//' merged.fasta > merged.fa
```

Use USEARCH<sup>5</sup> (v 11.0.667) to create a non-redundant library.
```
- sort the sequences by length
$ usearch -sortbylength merged.fa -fastaout merged.sorted.fa --log usearch.log

- cluster sequences in fasta that are >=80% identical
$ usearch -cluster_fast merged.sorted.fa --id 0.8 --centroids my_centroids.fa --uc result.uc -consout final.nr.consensus.fa -msaout aligned.fasta --log usearch2.log

- Output of interest
>> final.nr.consensus.fa 
```

Again, use RepeatClassifier to classify the combined repeat library in the right format
```
- Build the database
$ BuildDatabase -name Tdiversipes_TEs-combined -engine ncbi final.nr.consensus.fa 

- Classify
$ RepeatClassifier -engine ncbi -stockholm tdiv-families.stk -consensi final.nr.consensus.fa
```


### 4- Merge your customized repeat library to other libraries of possible interest already available at RepeatMasker

**PS UPDATE** *For the new version of RepeatMasker this step was ignored because queryRepeatDatabase.pl is not available and famdb.py failed to compile the hymenoptera partition.*

RepeatMasker contains repeat libraries from other species (mostly model species). You can search this library for certain taxonomic groups to use their repeats if you think that is reaseonable in your case.<br/>
For (*Tetrapedia diversipes*) I have searched for (*Apis*) (honeybee) repeats with:
```
$ RepeatMasker/util/queryRepeatDatabase.pl -species apis -stat
```

Then I've combined these repeats with my costume database. The output was my final repeat library.
```
$ /RepeatMasker/util/queryRepeatDatabase.pl -species apis > apis_repeats_filter.fa
$ cat final.nr.consensus.fa.classified apis_repeats_filter.fa > Tdiversipes_TEs.fasta
```

### 5- Use RepeatMasker with your final library
Finally get the repetitive elements on your genome based on a search with your costumized library using RepeatMasker<sup>6</sup> (v 4.1.0). 
```
$ RepeatMasker -pa 3 -s -lib Tdiversipes_TEs.fasta -xsmall -gff -html -gccalc TDIV_genome_05_2019.fasta
```

## Repository content

### Others
- **LTRdigest_parse.py**: Script to remove LTR candidates without a protein match (filter from LTRDdigest gff3 output).


### Statistics
- **TEs_around_genes.R**: Script to identify all repetitive elements occuring upstream/downstream a gene reagion (5kb window).

- **TEsfamily_enriched.R**: Script to test whether certain repetitive elements categories/families are enriched in a certain dataset in comparison to the entire genome.


### Figures
- **TEsInGenome.R**: Script to create the figure illutrating the proportion of TEs class/families in the genome.

- **Plot-TEs_around_genes.R**: Script to create the plot showing the position of TEs around genes. 





## License
```
This work is distributed under the GPLv3 license. Reuse of code derived from this repository is permitted under two conditions:

Proper attribution (i.e., citation of the associated publication; see CITATION.cff and above).
Publication of reused scripts on an open-access platform, such as Github.
```


## References
<sup>1</sup> Robert Hubley, Arian Smit - Institute for Systems Biology. RepeatModeler; Jullien M. Flynn - Cornell University. LTR Pipeline Extensions <http://www.repeatmasker.org/RepeatModeler/> <br/>
<sup>2</sup> Brian J. Haas, TransposonPSI, 2007-2011 <http://transposonpsi.sourceforge.net><br/>
<sup>3</sup> Ellinghaus, D. et al. LTRharvest, an efficient and flexible software for de novo detection of LTR retrotransposons. BMC bioinformatics 9, 18 (2008)<br/>
<sup>4</sup> Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11(10): e0163962. https://doi.org/10.1371/journal.pone.0163962<br/>
<sup>5</sup> Edgar,RC (2010) Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19), 2460-2461.
doi: 10.1093/bioinformatics/btq461<br/>
<sup>6</sup> Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015 <http://www.repeatmasker.org>.

