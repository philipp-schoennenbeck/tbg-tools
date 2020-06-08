# tbg-tools v0.2

##If you have version 0.1 download v0.2! There is a huge bug still in v0.1!
The purpose of this program is to quickly look up positions on scaffolds and get some further informations about these positions.
It is not in its finished version and there are surely still some bugs to be found. 

python version = 3.6+
### Usage:

python3 main_program.py --help

-> shows the different options of the program

Note: if you do not have a lot of RAM available, use the -w option

More options for each program can be found with python3 main_program.py create -h, python3 main_program.py search -h ...

#### create: this has to be used first to create a tbg file
e.g. : python3 main_program.py create -g example.gff -f example.fa -o example.tbg -t 4 -w -v -hro example.tsv

tbg files are binary files and can not be read by other programs but you can also
create a human readable tsv file which will be much larger.
The program itself uses the binary tbg files



#### search:

e.g. :  python3 main_program.py search -n example.tbg -b snps.tsv -o positions_found.tsv -t 4 -w

Search for positions within the tbg file. When using a bed/tsv file for positions use this pattern:

scaffold1   123551

scaffold1   662234

scaffold2   100123  100160

If a third column exists the region between the second and the third columns will be searched for.
The result file is a small portion of the human readable tsv file which only contains the found positions.
If some positions are not found within a gene, they can be written to a rest file (option -r)


#### convert:
e.g. : python3 main_program.py convert -n example.tbg -o example.tsv

Converts tbg files to human readable tsv files

#### pop_gen:
e.g. : pop_gen -n example.tbg -o example_results_vcf.tsv -vcf example_vcf_file.vcf -sf example_stats_vcf.tsv -t 4 -v

Uses vcf- or sync-files to search for specific positions. tbg-tools can only lookup SNPs not indels.


## human readable file columns:
ScaffID pos base_in_reference   coding_base triplett_position   coded_AA    geneID  4ds >A  >T  >C  >G 

scaffID
The Id of the scaffold

position
The coordinate of the current nucleotide

base in reference
The base in the reference genome at this position

coding base
The base found in the gene at this position. Is equal to base in reference if the gene is on + strand.
Otherwise it is the complementary base.

triplett position
Position in the current amino acid. Possible positions are 1, 2 and 3

coded amino acid
The amino acid which is coded at this position

gene ID
Id of the current gene

4ds
True if a change of the nucleotide would not lead to a change of the amino acid. False otherwise

A, T, C, G
Last four columns show what would happen if a different amino acid (in order A, T, C, G) would be in this spot.
