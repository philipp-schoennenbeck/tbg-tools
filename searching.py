from helper_functions import *


def load_file(path, binary=False):
    data = {}
    if not binary:
        data = read_file(path)
    if binary:
        data, genes, scaffolds = read_binary_file(path)
    return data


def check_snps(nucleotide_file, snp_file=None, snps=None, binary=False, outfile="snps.tsv", rest_file=None, threads=1):
    # load in all the snps
    if snp_file is not None:
        snps = []
        with open(snp_file, "r") as f:

            file = f.read().split("\n")
            for i in range(len(file)):

                dummy = file[i].split("\t")
                if len(file[i].strip("\t\n ")) == 0 or len(dummy) != 2:
                    continue
                snps.append([dummy[0], int(dummy[1])])

    # load in the tbg file
    scaffolds = {}
    genes = {}
    if binary:
        try:
            if threads > 1:
                data, genes, scaffolds = read_binary_file(nucleotide_file, threads=threads)
            else:
                data, genes, scaffolds = read_binary_file_no_threads(nucleotide_file)
        except Exception:
            raise Exception("Error with loading the tbg file!")
    else:
        try:
            data = read_file(nucleotide_file)
        except Exception:
            raise Exception("Error with loading in the hr file!")

    # check for every SNP if it is in the tbg file
    real_snps = []
    rest = []
    for snp in snps:
        if binary:
            if snp[0] in scaffolds:
                scaffold = scaffolds[snp[0]]
            else:
                scaffold = None
        else:
            scaffold = snp[0]
        if scaffold in data:
            if snp[1] in data[scaffold]:
                for line in data[scaffold][snp[1]]:
                    snp_dummy = snp
                    if binary:
                        snp_information = decode_line(line, genes)
                        snp_dummy.extend(snp_information)
                        snp_dummy.append("\n")
                        snp_dummy = "\t".join([str(i) for i in snp_dummy])

                    else:

                        snp_information = line
                        snp_dummy.extend(snp_information)
                        snp_dummy.append("\n")
                        snp_dummy = "\t".join([str(i) for i in snp_dummy])

                    real_snps.append(snp_dummy)
            else:
                if rest_file is not None:
                    rest.append(f"{snp[0]}\t{snp[1]}\n")
        else:
            if rest_file is not None:
                rest.append(f"{snp[0]}\t{snp[1]}\n")

    # output
    if rest_file is not None:
        r = open(rest_file, "w")
        for snp in rest:
            r.write(snp)
    f = open(outfile, "w")
    for snp in real_snps:
        f.write(snp)


def check_gene(tbg_file, genes_to_find, outfile, verbose, rest):
    if verbose:
        print("Searching for genes!")
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
