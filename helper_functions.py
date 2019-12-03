import re
import struct
import binascii
bases = ["A", "C", "G", "U"]
aa_codes = {
"AAA":["Lys","K"],
"AAC":["Asn","N"],
"AAG":["Lys","K"],
"AAU":["Asn","N"],
"ACA":["Thr","T"],
"ACC":["Thr","T"],
"ACG":["Thr","T"],
"ACU":["Thr","T"],
"AGA":["Arg","R"],
"AGC":["Ser","S"],
"AGG":["Arg","R"],
"AGU":["Ser","S"],
"AUA":["Ile","I"],
"AUC":["Ile","I"],
"AUG":["Met","M"],
"AUU":["Ile","I"],
"CAA":["Gln","Q"],
"CAC":["His","H"],
"CAG":["Gln","Q"],
"CAU":["His","H"],
"CCA":["Pro","P"],
"CCC":["Pro","P"],
"CCG":["Pro","P"],
"CCU":["Pro","P"],
"CGA":["Arg","R"],
"CGC":["Arg","R"],
"CGG":["Arg","R"],
"CGU":["Arg","R"],
"CUA":["Leu","L"],
"CUC":["Leu","L"],
"CUG":["Leu","L"],
"CUU":["Leu","L"],
"GAA":["Glu","E"],
"GAC":["Asp","D"],
"GAG":["Glu","E"],
"GAU":["Asp","D"],
"GCA":["Ala","A"],
"GCC":["Ala","A"],
"GCG":["Ala","A"],
"GCU":["ALa","A"],
"GGA":["Gly","G"],
"GGC":["Gly","G"],
"GGG":["Gly","G"],
"GGU":["Gly","G"],
"GUA":["Val","V"],
"GUC":["Val","V"],
"GUG":["Val","V"],
"GUU":["Val","V"],
"UAA":["Stop","stop"],
"UAC":["Tyr","Y"],
"UAG":["Stop","stop"],
"UAU":["Tyr","Y"],
"UCA":["Ser","S"],
"UCC":["Ser","S"],
"UCG":["Ser","S"],
"UCU":["Ser","S"],
"UGA":["Stop","stop"],
"UGC":["Cys","C"],
"UGG":["Trp","W"],
"UGU":["Cys","C"],
"UUA":["Leu","L"],
"UUC":["Phe","F"],
"UUG":["Leu","L"],
"UUU":["Phe","F"]
}


def load_aa_codes(code="default"):
    aa_code = {}
    with open("aa_codes.txt", "r") as f:
        found_code = False
        for line in f:
            line = line.strip()
            if found_code:
                if line[0] == "q":
                    f.close()
                    return aa_code
                line = line.split(";")
                aa_code[line[0]] = [line[1], line[2]]
            elif line == f"@{code}":
                found_code = True
    f.close()
    return aa_code


def get_code(code, aa_codes, three_letter_code = False):
    """Returns the amino acid code out of a string with 3 characters.
    By default the one letter code is returned. If three_letter_code is True, the three letter code will be returned"""
    code = code.upper()
    try:
        if three_letter_code:
            return aa_codes[code][0]
        return aa_codes[code][1]
    except Exception:
        raise Exception(f"\"{code}\" has no corresponding amino acid")


def get_all_aa(code, position, three_letter_code=False):
    """Returns all possible amino acids with changes at the position"""
    try:
        triplet = [code[0], code[1], code[2]]
        amino_acids = []
        for base in bases:
            triplet[position] = base
            amino_acids.append(get_code("".join(triplet), three_letter_code))
        return amino_acids
    except Exception:
        raise Exception(f"Cannot find all possible outcomes when changing one base. (Base: {code},"
                        f" position: {position})")


def get_other_strand(strand):
    """"Returns the corresponding bases of the other DNA strand"""
    strand = strand.upper()
    strand = strand.replace("A", "u").replace("T", "a").replace("C", "g").replace("G", "c")[::-1].upper()
    return strand


def load_gff(path, types=None, verbose=False):
    """"Loads a gff file and searches for all types of the given type.
     Returns a list of all the entries with that type."""
    with open(path, "r") as f:
        if verbose:
            print("Starting to load gff!")
        # try:
        #     gff_data = np.asarray(f)
        #     if type == None:
        #         return gff_data
        #     else:
        #         selected_data = gff_data[gff_data[:,2] == type]
        #         return selected_data
        # except:

        if types is not None:
            gff_data = {i: [] for i in types}
        else:
            gff_data = {"all": []}
        for line in f:
            if line[0] == "#":
                next
            line_split = line.split("\t")
            if types is None:
                if len(line_split) >= 3:
                    gff_data["all"].append(line)
            else:
                if len(line_split) >= 3:
                    for i in types:
                        if line_split[2] == i:
                            gff_data[i].append(line)
        if verbose: print("Loading gff completed!")
        return gff_data


def load_fasta(path, seperator=None, verbose=False):
    """Loads a fasta file. Returns a dictionary with the scaffold names as the key and the sequence as the value."""
    if verbose:
        print("Starting to load fasta!")
    try:
        with open(path, "r") as f:
            scaffolds = {}
            current_scaffold = ""
            for line in f:

                line = line.strip()
                if line[0] == ">":
                    line = line[1:]
                    if seperator is not None:
                        line = line.split(seperator)
                        line = line[0]
                    scaffolds[line] = []
                    current_scaffold = line
                    if verbose: print(f"scaffold {line} loading")
                else:
                    scaffolds[current_scaffold].append(line)
        f.close()
        for scaffold in scaffolds.keys():
            scaffolds[scaffold] = "".join(scaffolds[scaffold])
        if verbose: print("Loading fasta completed")
        return scaffolds
    except:
        raise Exception(f"Error with opening the fasta file \"{path}\"")


def get_sequence_from_fasta(scaffold, positions, forward=True):
    """Gets a sequence of a gene from a fasta file. positions can be a list with a start and end point or a list of
    lists with multiple start and endpoints"""
    sequence = ""
    try:
        if type(positions[0]) != list:
            positions = [positions]
    except:
        raise Exception("positions is of type " + str(type(positions)))
    for position in positions:
        sequence = f"{sequence}{scaffold[int(position[0])-1:int(position[1])]}"
    if not forward:
        sequence = get_other_strand(sequence)
    return sequence


def coding_strand_to_rna(strand):
    """returns the coding strand to the rna strand (T --> U)"""
    strand = strand.upper().replace("T","U")
    return strand


def rna_to_protein(rna, aa_codes):
    """translates the rna sequence to the amino acid sequence"""
    rna = rna.upper()
    protein = []
    for position in range(0,len(rna),3):
        triplet = rna[position:position+3]
        if len(triplet) == 3:
            protein.append(aa_codes[triplet][1])
    return "".join(protein[:-1])


def sort_first(element):
    return element[0]


def get_genes_and_cds(gff_path, verbose=False):
    """loads a gff file and finds all the genes, CDS and mRNA and returns the genes with their corresponding CDS"""
    if verbose:print("Starting to structure gff data!")
    gff_data = load_gff(gff_path, ["gene","CDS", "mRNA"],verbose=verbose)
    gene_and_cds = {}
    mrna = {}
    for gene in gff_data["gene"]:
        gene = gene.strip()
        gene = gene.split("\t")
        gene_split = gene[8].split(";")
        for i in gene_split:
            if "ID=" == i[0:3]:
                gene_id = i[3:]
                if gene[6] == "+":
                    gene_and_cds[gene_id] = [[], True, gene[0]]
                else:
                    gene_and_cds[gene_id] = [[], False, gene[0]]
    for rna in gff_data["mRNA"]:
        rna = rna.strip()
        rna  = rna.split("\t")
        rna_split = rna[8].split(";")
        for i in rna_split:
            if i[0:7] == "Parent=":
                rna_parent = i[7:]
            if i[0:3] == "ID=":
                rna_id = i[3:]
        mrna[rna_id] = rna_parent
    for cds in gff_data["CDS"]:
        cds = cds.strip()
        cds = cds.split("\t")
        cds_split = cds[8].split(";")
        for i in cds_split:
            if i[0:7] == "Parent=":
                if i[7:] in gene_and_cds:
                    gene_and_cds[i[7:]][0].append([cds[3], cds[4]])
                elif i[7:] in mrna:
                    gene_and_cds[mrna[i[7:]]][0].append([cds[3], cds[4]])
    delete_list = []
    for key in gene_and_cds.keys():
        if gene_and_cds[key][0] == []:
            delete_list.append(key)
    for key in delete_list:
        del gene_and_cds[key]
    for key in gene_and_cds.keys():
        gene_and_cds[key][0].sort(key=sort_first)
    return gene_and_cds


def get_position_in_scaffold(exon_list, position):
    """calculates the position of the nucleotide on the scaffold"""
    exon_length = [int(i[1]) - int(i[0])+1 for i in exon_list]
    sum_of_exons = 0
    for exon in range(len(exon_length)):
        if position < sum_of_exons + exon_length[exon]:
            return position + int(exon_list[exon][0]) - sum_of_exons
        sum_of_exons += exon_length[exon]
    return -1


def get_current_triplett(sequence, position, forward=True):
    """returns the triplett of the current position and the position within the triplett"""
    if not forward:
        sequence = sequence[::-1]
        position = len(sequence)-position-1
    triplett_position = position % 3
    start_of_triplett = int(position/3)*3
    triplett = sequence[start_of_triplett:start_of_triplett+3]
    return triplett, triplett_position


def format_string(length_of_string=0):
    return "IIccccI?cccc"

    # fmt = ["c" * length_of_string[0]]
    # fmt.append("IccIc")
    # fmt.append("c" * length_of_string[1])
    # fmt.append("?cccc")

    return "".join(fmt)

def byte_size_of_fmt(fmt):
    sizes = {"I" : 4, "c" : 1, "?" : 1}
    size = 0
    for i in fmt:
        size += sizes[i]
    return size

def write_binary(data, path):
    """binary file:
    int i(how many scaffolds are there),
    i ints (how long are the strings of the scaffolds),
    characters for the scaffolds
    same thing for genes
    after that: format string: IIccccI?cccc
    --> scaffold number, position, base on ref, base in gene, position in triplet, amino acid, gene number,
        4ds, to A, to C, to G, to T"""
    
    
    genes = {}
    genes_list = []
    scaffolds = {}
    scaffolds_list = []
    scaffolds_counter = 0
    genes_counter = 0
    for line in data:
        data_split = line.strip().split("\t")
        if data_split[0] not in scaffolds:
            scaffolds[data_split[0]] = scaffolds_counter
            scaffolds_counter += 1
            scaffolds_list.append(data_split[0])
        if data_split[6] not in genes:
            genes[data_split[6]] = genes_counter
            genes_counter += 1
            genes_list.append(data_split[6])

    with open(path, "wb") as f:
        f.write(struct.pack("I", len(genes)))
        f.write(struct.pack("I"*len(genes), *[len(i) for i in genes_list]))
        for gene in genes_list:
            f.write(struct.pack("c" * len(gene), *[bytes(i, "utf-8") for i in gene]))

        f.write(struct.pack("I", len(scaffolds)))
        f.write(struct.pack("I"*len(scaffolds), *[len(i) for i in scaffolds_list]))
        for scaffold in scaffolds_list:
            f.write(struct.pack("c" * len(scaffold), *[bytes(i, "utf-8") for i in scaffold]))
        for line in data:
            data_split = line.strip().split("\t")
            for i in [5, 8, 9, 10, 11]:
                if data_split[i] == "stop":
                    data_split[i] = "-"
            data_split[0] = scaffolds[data_split[0]]
            data_split[1] = int(data_split[1])
            # data_split[4] = int(data_split[4])
            data_split[6] = genes[data_split[6]]
            if data_split[7] == "True":
                data_split[7] = True
            elif data_split[7] == "False":
                data_split[7] = False
            single_char_list = []
            for i in data_split:
                if type(i) == str:
                    a = [bytes(b, "utf-8") for b in i]
                    single_char_list.extend(a)

                else:
                    single_char_list.append(i)
            fmt = format_string(0)
            f.write(struct.pack(fmt, *single_char_list))

    f.close()


def read_file(path):
    data = {}
    with open(path, "r") as f:
        for line in f:
            line_split = line.strip().split("\t")
            line_split[1] = int(line_split[1])
            if line_split[7] == "False":
                line_split[7] = False
            elif line_split[7] == "True":
                line_split[7] = True
            if line_split[0] in data:
                data[line_split[0]][line_split[1]] = line_split[2:]
            else:
                data[line_split[0]] = {}
                data[line_split[0]][line_split[1]] = line_split[2:]
    f.close()
    return data
            #todo: read file, figure out what to search
            
def read_binary_file(path):
    data = {}
    with open(path, "rb") as f:
        genes = {}
        number_of_genes = struct.unpack("I", f.read(4))[0]
        genes_str_length = struct.unpack("I" * number_of_genes, f.read(4 * number_of_genes))
        for i in range(number_of_genes):
            scaffold_name = struct.unpack("c" * genes_str_length[i], f.read(genes_str_length[i]))
            scaffold_name = "".join([str(letter, "utf-8") for letter in scaffold_name])
            genes[i] = scaffold_name
        scaffolds = {}
        number_of_scaffolds = struct.unpack("I", f.read(4))[0]
        scaffolds_str_length = struct.unpack("I" * number_of_scaffolds, f.read(4 * number_of_scaffolds))
        for i in range(number_of_scaffolds):
            scaffold_name = struct.unpack("c" * scaffolds_str_length[i], f.read(scaffolds_str_length[i]))
            scaffold_name = "".join([str(letter, "utf-8") for letter in scaffold_name])
            scaffolds[i] = scaffold_name
        while True:
            scaffold_and_position = f.read(8)
            rest = f.read(byte_size_of_fmt(format_string()) - 8)
            if not scaffold_and_position or not rest:
                # EOF
                break
            scaffold_and_position =struct.unpack("II", scaffold_and_position)
            if scaffold_and_position[0] in data:
                data[scaffold_and_position[0]][scaffold_and_position[1]] = rest
            else:
                data[scaffold_and_position[0]] = {}
                data[scaffold_and_position[0]][scaffold_and_position[1]] = rest

            # print(line_translated)
        # while True:
        #     bytes_line = f.read(byte_size_of_fmt(format_string()))
        #     if not bytes_line:
        #         # EOF
        #         break
        #     line = struct.unpack(format_string(), bytes_line)
        #     line_translated = [scaffolds[line[0]], line[1], str(line[2], "utf-8"), str(line[3], "utf-8"),
        #                        str(line[4], "utf-8"), str(line[5], "utf-8"), genes[line[6]], line[7],
        #                        str(line[8], "utf-8"), str(line[9], "utf-8"), str(line[10], "utf-8"),
        #                        str(line[11], "utf-8")]
        #     for i in [5, 8, 9, 10, 11]:
        #         if line_translated[i] == "-":
        #             line_translated[i] = "stop"
        #
        #     if line_translated[0] in data:
        #         data[line_translated[0]][line_translated[1]] = line_translated[2:]
        #     else:
        #         data[line_translated[0]] = {}
        #         data[line_translated[0]][line_translated[1]] = line_translated[2:]
        #
        #     # print(line_translated)
    return data


if __name__ == "__main__":
    pass