# tbg 

The purpose of this program is to quickly look up positions on scaffolds and get some further informations about these positions.
It is not in its finished version and there are surely still some bugs to be found. 

### Usage:

python3 main_program.py --help

-> shows the different options of the program

Note: if you do not have a lot of RAM available, use the -w option


####create: this has to be used first to create a tbg file
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

still working on that


## human readable file columns:
ScaffID pos base_in_reference   coding_base triplett_position   coded_AA    geneID  4ds >A  >T  >C  >G 

\>A, \>T , \>C, \>G : Which amino acid will be coded if the nucleotide changes to this 

4ds : True if  \>A, \>T , \>C, \>G show all the same amino acid
