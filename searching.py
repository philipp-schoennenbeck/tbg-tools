from helper_functions import *
import datetime

def load_file(path, binary=False):
    if not binary:
        data = read_file(path)
    if binary:
        data = read_binary_file(path)
    return data

if __name__ == "__main__":
    gff_file = "E_coli.gff"
    fasta_file = "E_coli.fna"
    outfile = "E_coli"
    # gff_file = "radix.gff"
    # fasta_file = "radix.fa"
    # outfile = "radix"

    time = datetime.now()
    radix_data = load_file(f"{outfile}.tsv", binary=False)
    print((datetime.now()-time).total_seconds())
    # print(radix_data)
    time = datetime.now()
    radix_data = load_file(f"binary_{outfile}.bin", binary=True)
    print((datetime.now() - time).total_seconds())
    # print(radix_data)