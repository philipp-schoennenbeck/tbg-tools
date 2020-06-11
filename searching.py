from helper_functions import *
from copy import deepcopy
import create_file


def check_snps(nucleotide_file, snp_file=None, snps=None, binary=False, outfile=["snps.tsv"], rest_file=None, threads=1,
               low_ram=False, verbose=False):
    if snp_file is not None:
        if verbose:
            print("Load in snp_file(s)!")
        snps_dict = {}
        snps_dict_per_file = {file : [] for file in snp_file}
        for snp_file_single in snp_file:

            with open(snp_file_single, "r") as f:

                file = f.read().split("\n")
                for i in range(len(file)):

                    dummy = file[i].split("\t")
                    if len(file[i].strip("\t\n ")) == 0:
                        continue
                    if len(dummy) == 2:
                        snps_dict[(dummy[0], int(dummy[1]))] = True
                        snps_dict_per_file[snp_file_single].append((dummy[0], int(dummy[1])))
                    elif len(dummy) == 3:
                        for j in range(int(dummy[1]), int(dummy[2])):
                            snps_dict[(dummy[0], j)] = True
                            snps_dict_per_file[snp_file_single].append((dummy[0], j))
        snps = []
        for snp in snps_dict.keys():
            snps.append(snp)
    if rest_file is None:
        rest_file = [None for _ in snp_file]

    if low_ram:
        if verbose:
            print("Searching in tbg file for SNPs...")
        rest,found = check_snps_low_ram(nucleotide_file, snps, verbose=verbose)
    else:
        rest,found = check_snps_normal(nucleotide_file, snps, threads, verbose=verbose)
    if snp_file is not None:
        if len(outfile) != 1:
            for outfile_single, snp_file_single, rest_file_single in zip(outfile, snp_file, rest_file):
                with open(outfile_single, "w") as outf:
                    if rest_file_single is not None:
                        rest_file_out = open(rest_file_single, "w")
                    for snp in snps_dict_per_file[snp_file_single]:
                        if snp in found:

                            for line in found[snp]:
                                line = "\t".join([str(i) for i in snp]) + "\t" + "\t".join([str(i) for i in line])
                                outf.write(line)
                        else:
                            if rest_file_single is not None:
                                rest_file_out.write("\t".join([str(i) for i in snp]) + "\n")
        else:
            with open(outfile[0], "w") as outf:
                for snp_file_single, rest_file_single in zip(snp_file, rest_file):
                    if rest_file_single is not None:
                        rest_file_out = open(rest_file_single, "w")
                    for snp in snps_dict_per_file[snp_file_single]:

                        if snp in found:
                            for line in found[snp]:
                                line = "\t".join([str(i) for i in snp]) + "\t" + "\t".join([str(i) for i in line])
                                outf.write(line)
                        else:
                            if rest_file_single is not None:
                                rest_file_out.write("\t".join([str(i) for i in snp]) + "\n")
                    if rest_file_single is not None:
                        rest_file_out.close()

    else:
        with open(outfile[0], "w") as outf:
            for keys in found.keys():
                for line in found[keys]:
                    line = "\t".join([str(i) for i in keys]) + "\t" + "\t".join([str(i) for i in line])
                    outf.write(line)
        if rest_file[0] is not None:
            r = open(rest_file[0], "w")
            for snp in rest:
                r.write(snp)


def check_snps_normal(nucleotide_file, snps=None, threads=1, verbose=False):

    snps_found = {}
    if verbose:
        print("Starting to load the tbg file!")

    try:
        data, genes, scaffolds = read_binary_file(nucleotide_file, threads=threads)
    except Exception:
        raise Exception("Error with loading the tbg file!")


    # check for every SNP if it is in the tbg file
    rest = []
    if verbose:
        print("Looking for SNPs in the tbg file!")
    for snp in snps:

        if snp[0] in scaffolds:
            scaffold = scaffolds[snp[0]]
        else:
            scaffold = None

        if scaffold in data:
            if snp[1] in data[scaffold]:
                snps_found[(snp[0], snp[1])] = []
                for line in data[scaffold][snp[1]]:
                    snp_information = decode_line(line, genes)
                    snp_information.append("\n")
                    snps_found[(snp[0], snp[1])].append(snp_information)
            else:
                rest.append(f"{snp[0]}\t{snp[1]}\n")
        else:
            rest.append(f"{snp[0]}\t{snp[1]}\n")
    return (rest, snps_found)

def check_snps_low_ram(nucleotide_file, snps=None, verbose=False):

    snps_dic = {i:None for i in snps}


    with open(nucleotide_file, "rb") as f:
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
            rest.append("\n")
            if (scaffolds[scaffold_and_position[0]], scaffold_and_position[1]) in snps_dic:
                if snps_dic[ (scaffolds[scaffold_and_position[0]], scaffold_and_position[1])] is None:

                    snps_dic[scaffolds[scaffold_and_position[0]], scaffold_and_position[1]] = [rest]
                else:
                    snps_dic[scaffolds[scaffold_and_position[0]], scaffold_and_position[1]].append(rest)

    rest=[]
    return_snps = {}
    for i in snps_dic.keys():
        if snps_dic[i] is None:
            rest.append(f"{i[0]}\t{i[1]}\n")
        else:
            return_snps[i] = snps_dic[i]
    return (rest,return_snps)





def check_gene(tbg_file, genes_to_find, outfile, verbose, rest):
    if verbose:
        print("Searching for genes!")
    outfile = outfile[0]
    with open(tbg_file, "rb") as f:
        genes = {}
        genes_found = {}
        number_of_genes = struct.unpack("I", f.read(4))[0]
        genes_str_length = struct.unpack("I" * number_of_genes, f.read(4 * number_of_genes))
        for i in range(number_of_genes):
            gene_name = struct.unpack("c" * genes_str_length[i], f.read(genes_str_length[i]))
            gene_name = "".join([str(letter, "utf-8") for letter in gene_name])
            genes[i] = gene_name
            genes_found[gene_name] = True
        scaffolds = {}
        number_of_scaffolds = struct.unpack("I", f.read(4))[0]
        scaffolds_str_length = struct.unpack("I" * number_of_scaffolds, f.read(4 * number_of_scaffolds))
        for i in range(number_of_scaffolds):
            scaffold_name = struct.unpack("c" * scaffolds_str_length[i], f.read(scaffolds_str_length[i]))
            scaffold_name = "".join([str(letter, "utf-8") for letter in scaffold_name])
            scaffolds[i] = scaffold_name

        rest_file = None
        if rest:
            rest_file = open(rest, "w")
        for gene in genes_to_find:
            if gene not in genes_found:
                if rest:
                    rest_file.write(gene)
                if verbose:
                    print(f"{gene} was not found!")
        if rest:
            rest_file.close()
        if verbose:
            print("Writing the result file...")
        fmt_size_before_gene = byte_size_of_fmt("IIcccc")
        with open(outfile, "w") as outfile:
            while True:
                line = f.read(byte_size_of_fmt(format_string()))
                gene = line[fmt_size_before_gene:fmt_size_before_gene + 4]

                if not line:
                    # EOF
                    break
                gene = struct.unpack("I", gene)
                if genes[gene[0]] in genes_to_find:
                    line = struct.unpack(format_string(), line)
                    line_translated = [scaffolds[int(line[0])], str(int(line[1])), str(line[2], "utf-8"),
                                       str(line[3], "utf-8"),
                                       str(line[4], "utf-8"), str(line[5], "utf-8"), genes[line[6]], str(line[7]),
                                       str(line[8], "utf-8"), str(line[9], "utf-8"), str(line[10], "utf-8"),
                                       str(line[11], "utf-8"), "\n"]
                    for i in [5, 8, 9, 10, 11]:
                        if line_translated[i] == "-":
                            line_translated[i] = "stop"
                    outfile.write("\t".join(line_translated))


def check_scaffold(tbg_file, scaffolds_to_find, outfile, verbose, rest):
    if verbose:
        print("Searching for scaffolds!")
    outfile = outfile[0]
    with open(tbg_file, "rb") as f:
        genes = {}
        scaffolds_found = {}
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
            scaffolds_found[scaffold_name] = True

        rest_file = None
        if rest:
            rest_file = open(rest, "w")
        for scaffold in scaffolds_to_find:
            if scaffold not in scaffolds_found:
                if rest:
                    rest_file.write(scaffold)
                if verbose:
                    print(f"{scaffold} was not found!")
        if rest:
            rest_file.close()
        if verbose:
            print("Writing the result file...")
        with open(outfile, "w") as outfile:
            while True:
                line = f.read(byte_size_of_fmt(format_string()))
                scaffold = line[0:4]

                if not line:
                    # EOF
                    break
                scaffold = struct.unpack("I", scaffold)
                if scaffolds[scaffold[0]] in scaffolds_to_find:
                    line = struct.unpack(format_string(), line)
                    line_translated = [scaffolds[int(line[0])], str(int(line[1])), str(line[2], "utf-8"),
                                       str(line[3], "utf-8"),
                                       str(line[4], "utf-8"), str(line[5], "utf-8"), genes[line[6]], str(line[7]),
                                       str(line[8], "utf-8"), str(line[9], "utf-8"), str(line[10], "utf-8"),
                                       str(line[11], "utf-8"), "\n"]
                    for i in [5, 8, 9, 10, 11]:
                        if line_translated[i] == "-":
                            line_translated[i] = "stop"
                    outfile.write("\t".join(line_translated))


if __name__ == "__main__":
    pass
