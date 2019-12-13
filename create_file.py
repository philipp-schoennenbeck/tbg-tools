from helper_functions import *
from datetime import datetime
from time import sleep

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
        line = line.strip()
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
        f.write(struct.pack("I" * len(genes), *[len(i) for i in genes_list]))
        for gene in genes_list:
            f.write(struct.pack("c" * len(gene), *[bytes(i, "utf-8") for i in gene]))

        f.write(struct.pack("I", len(scaffolds)))
        f.write(struct.pack("I" * len(scaffolds), *[len(i) for i in scaffolds_list]))
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
            fmt = format_string()
            try:
                f.write(struct.pack(fmt, *single_char_list))
            except:
                pass

    f.close()


def create_the_file(gff_file, fasta_file, outfile_hr="default.tsv", outfile_bin="defailt.bin", verbose=False,
                    create_binary=False, write_tsv=True, aa_code="default"):
    starting_time = datetime.now()
    timestamps = [0]
    outlist = []
    aa_codes = load_aa_codes(aa_code)
    time = datetime.now()
    gff_data = get_genes_and_cds(gff_file, verbose=verbose)
    timestamps.append((datetime.now() - time).total_seconds())
    time = datetime.now()
    fasta_data = load_fasta(fasta_file, seperator=" ", verbose=verbose)
    timestamps.append((datetime.now() - time).total_seconds())

    gene_sequences = {}
    time = datetime.now()
    for gene in gff_data.keys():
        sequence = get_sequence_from_fasta(fasta_data[gff_data[gene][2]], gff_data[gene][0])
        gene_sequences[gene] = coding_strand_to_rna(sequence)
    timestamps.append((datetime.now() - time).total_seconds())
    time_generate = datetime.now()

    for i in range(10):
        timestamps.append(0)
    #
    with open("protein_test.fa", "w") as f:
        for gene in gene_sequences.keys():
            if gene == "maker-scaffold1079_size269654-augustus-gene-0.59":
                a = gene_sequences[gene]
                b = gff_data[gene]
                c = get_other_strand(gene_sequences[gene])
                d = gene_sequences[gene][::-1]
            f.write(gene + "\n")
            if gff_data[gene][1]:
                f.write(rna_to_protein(gene_sequences[gene], aa_codes) + "\n")
            else:
                f.write(rna_to_protein(get_other_strand(gene_sequences[gene]), aa_codes) + "\n")

    f.close()

    # with open("radix_whole.gff") as f:
    #     a = f.read()
    #     a = a.split("\n")
    #     while True:
    #         sleep(20)


    if verbose:
        print("Starting to calculate the nucleotides!")
        number_of_genes = len(gff_data)
        counter = 0.0
        printer = 0
    for gene in gene_sequences.keys():
        if verbose:
            counter += 1.0
            if counter / number_of_genes >= printer:
                print(f"{round(printer * 100)}% done")
                printer += 0.05
        forward = gff_data[gene][1]
        for nucleotide in range(len(gene_sequences[gene])):
            scaff_id = gff_data[gene][2]
            time = datetime.now()
            position = get_position_in_scaffold(gff_data[gene][0], nucleotide)
            timestamps[4] += (datetime.now() - time).total_seconds()
            time = datetime.now()
            triplett, position_in_triplett = get_current_triplett(gene_sequences[gene], nucleotide, forward)
            timestamps[5] += (datetime.now() - time).total_seconds()
            time = datetime.now()
            base = gene_sequences[gene][nucleotide].replace("U", "T")
            if forward:
                base_in_gene = base
            else:
                base_in_gene = str(get_other_strand(base)).replace("U", "T")
            timestamps[6] += (datetime.now() - time).total_seconds()
            time = datetime.now()
            amino_acid = get_code(triplett, aa_codes)
            if amino_acid is None:
                amino_acid = "+"
                # continue
                # raise Exception(position, scaff_id, gene)
            gene_id = gene
            timestamps[7] += (datetime.now() - time).total_seconds()
            time = datetime.now()
            amino_acids = get_all_aa(triplett, position_in_triplett, aa_codes)
            timestamps[8] += (datetime.now() - time).total_seconds()

            four_ds = all([i == amino_acids[0] for i in amino_acids])
            A, C, G, T = amino_acids
            outlist.append("\t".join((str(i) for i in [scaff_id, position, base, base_in_gene,
                                                       position_in_triplett +1 , amino_acid, gene_id,
                                                       four_ds, A, C, G, T, "\n"])))
    if write_tsv:
        with open(outfile_hr, "w") as outf:
            time = datetime.now()
            for i in outlist:
                outf.write(i)
            timestamps[9] += (datetime.now() - time).total_seconds()
        outf.close()
    if create_binary:
        time = datetime.now()
        write_binary(outlist, outfile_bin)
        timestamps[10] = (datetime.now() - time).total_seconds()

    timestamps[3] = (datetime.now() - time_generate).total_seconds()
    timestamps[0] = (datetime.now() - starting_time).total_seconds()
    print(timestamps)




def write_human_readable(path_bin, path_hr=None):
    if path_hr is None:
        path_hr = path_bin[:-3] + "tsv"
    with open(path_hr, "w") as outf:
        with open(path_bin, "rb") as f:
            genes = {}
            number_of_genes = struct.unpack("I", f.read(4))[0]
            genes_str_length = struct.unpack("I" * number_of_genes, f.read(4 * number_of_genes))
            for i in range(number_of_genes):
                gene_name = struct.unpack("c" * genes_str_length[i], f.read(genes_str_length[i]))
                gene_name = "".join([str(letter, "utf-8") for letter in gene_name])
                genes[i] = gene_name
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
                scaffold_and_position = struct.unpack("II", scaffold_and_position)
                rest = [str(i) for i in  decode_line(rest, genes)]
                outf.write(scaffolds[scaffold_and_position[0]] + "\t" + str(scaffold_and_position[1]) + "\t" +
                               "\t".join(rest) + "\n")




if __name__ == "__main__":
    pass
