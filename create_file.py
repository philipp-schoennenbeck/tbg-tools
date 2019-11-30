from helper_functions import *
from datetime import datetime

if __name__ == "__main__":
    def create_the_file():
        starting_time = datetime.now()
        timestamps = [0]
        outlist = []

        gff_file = "E_coli.gff"
        fasta_file = "E_coli.fna"
        outfile = "E_coli"
        # gff_file = "radix.gff"
        # fasta_file = "radix.fa"
        # outfile = "radix.tsv"
        verbose = True
        time = datetime.now()
        gff_data = get_genes_and_cds(gff_file,verbose=verbose)
        timestamps.append((datetime.now() - time).total_seconds())
        time = datetime.now()
        fasta_data = load_fasta(fasta_file,seperator=" ",verbose=verbose)
        timestamps.append((datetime.now()-time).total_seconds())

        gene_sequences = {}
        time = datetime.now()
        for gene in gff_data.keys():
            sequence = get_sequence_from_fasta(fasta_data[gff_data[gene][2]], gff_data[gene][0])
            gene_sequences[gene] = coding_strand_to_rna(sequence)
        timestamps.append((datetime.now() - time).total_seconds())
        time_generate = datetime.now()
        for i in range(10):
            timestamps.append(0)

        with open(outfile + ".tsv","w") as outf:
            if verbose:
                print("Starting to create the file!")
                number_of_genes = len(gff_data)
                counter = 0.0
                printer = 0
            for gene in gene_sequences.keys():
                if verbose:
                    counter += 1.0
                    if counter/number_of_genes >= printer:
                        print(f"{round(printer*100)}% done")
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
                        base_in_gene = str(get_other_strand(base)).replace("U","T")
                    timestamps[6] += (datetime.now() - time).total_seconds()
                    time = datetime.now()
                    amino_acid = get_code(triplett)
                    gene_id = gene
                    timestamps[7] += (datetime.now() - time).total_seconds()
                    time = datetime.now()
                    amino_acids = get_all_aa(triplett, position_in_triplett)
                    timestamps[8] += (datetime.now() - time).total_seconds()

                    four_ds = all([i == amino_acids[0] for i in amino_acids])
                    A, C, G, T = amino_acids
                    outlist.append("\t".join((str(i) for i in [scaff_id, position, base, base_in_gene,
                                                               position_in_triplett, amino_acid, gene_id,
                                                               four_ds, A, C, G, T, "\n" ])))
            outf.write("".join(outlist))
            timestamps[9] += (datetime.now() - time).total_seconds()
        time = datetime.now()
        write_binary(outlist, "binary_" + outfile + ".bin")
        timestamps[10] = (datetime.now() - time).total_seconds()
        outf.close()
        timestamps[3] = (datetime.now() - time_generate).total_seconds()
        timestamps[0] = (datetime.now() - starting_time).total_seconds()
        print(timestamps)

    def load_file(path, binary=False):
        if not binary:
            data = read_file(path)
        if binary:
            pass


    create_the_file()