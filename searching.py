from helper_functions import *

def load_file(path, binary=False):
    if not binary:
        data = read_file(path)
    if binary:
        data, genes, scaffolds = read_binary_file(path)
    return data


def check_snps(nucleotide_file, snp_file=None, snps=[], binary=False, outfile="snps.tsv", rest_file=None, threads=1):
    #load in all the snps
    if snp_file is not None:
        with open(snp_file, "r") as f:

            file = f.read().split("\n")
            for i in range(len(file)):

                dummy = file[i].split("\t")
                if len(file[i].strip("\t\n ")) == 0 or len(dummy) != 2:
                    continue
                snps.append([dummy[0], int(dummy[1])])

    #load in the tbg file
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

    #check for every SNP if it is in the tbg file
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

    #output
    if rest_file is not None:
        r = open(rest_file, "w")
        for snp in rest:
            r.write(snp)
    f = open(outfile, "w")
    for snp in real_snps:
        f.write(snp)



if __name__ == "__main__":
    pass