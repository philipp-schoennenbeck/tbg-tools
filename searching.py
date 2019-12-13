from helper_functions import *
from datetime import datetime

def load_file(path, binary=False):
    if not binary:
        data = read_file(path)
    if binary:
        data, genes, scaffolds = read_binary_file(path)
    return data


def check_snps(nucleotide_file, snp_file=None, snps=None, binary=False, outfile="snps.tsv", rest_file=None):

    if snp_file is not None:
        snps = []
        with open(snp_file, "r") as f:
            snps = f.read().split("\n")
            for i in range(len(snps)):
                dummy = snps[i].split("\t")
                snps[i] = [dummy[0], int(dummy[1])]

    if binary:
        data, genes, scaffolds = read_binary_file(nucleotide_file)
    else:
        data = read_file(nucleotide_file)
    if rest_file is not None:
        r = open(rest_file, "w")
    with open(outfile, "w") as f:
        for snp in snps:
            if binary:
                if snp[0] in scaffolds:
                    scaffold = scaffolds[snp[0]]
            else:
                scaffold = snp[0]
            if scaffold in data:
                if snp[1] in data[scaffold]:
                    if binary:
                        snp_information = decode_line(data[scaffold][snp[1]], genes)
                        snp.extend(snp_information)
                        snp.append("\n")
                        snp = "\t".join([str(i) for i in snp])

                    else:
                        snp_information = data[scaffold][snp[1]]
                        snp.extend(snp_information)
                        snp.append("\n")
                        snp = "\t".join([str(i) for i in snp])

                    f.writelines(snp)
                else:
                    if rest_file is not None:
                        r.write(f"{scaffold}\t{snps[0]}\n")
            else:
                if rest_file is not None:
                    r.write(f"{scaffold}\t{snps[0]}\n")


if __name__ == "__main__":
    pass