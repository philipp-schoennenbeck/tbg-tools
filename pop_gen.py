from helper_functions import *

def analyze_sync_file(nucleotide_file, sync_file, output, rest, threads):
    sync_data = {}
    order = "ATCGND"
    data, genes, scaffolds = read_binary_file(nucleotide_file, threads=threads)
    with open(sync_file, "r") as sync:
        for line in sync:
            line = line.strip()
            columns = line.split("\t")
            sync_data[(columns[0], int(columns[1]), columns[2])] = columns[3:]

    outfile = open(output, "w")
    for keys in sync_data.keys():
        if scaffolds[keys[0]] in data:
            if keys[1] in data[scaffolds[keys[0]]]:
                for line in data[scaffolds[keys[0]]][keys[1]]:
                    line = decode_line(line, genes)
                    outfile.write("\t".join([str(keys[0]), str(keys[1]), str(keys[2])]))
                    for population in sync_data[keys]:
                        amino_acids = []
                        current_pop = population.split(":")
                        current_pop = [int(x) for x in current_pop]
                        for base in range(4):
                            if current_pop[base] == 0:
                                amino_acids.append("0")
                            else:
                                amino_acids.append(str(current_pop[base]) + line[6 + base])
                        if all([x == line[6] for x in line[7:]]):
                            amino_acids.append(str(current_pop[4]) + line[6])
                        else:
                            amino_acids.append(str(current_pop[4]) + "-")

                        amino_acids = ":".join(amino_acids)
                        outfile.write("\t" + amino_acids)
                    outfile.write("\n")

    outfile.close()

if __name__ == "__main__":
    analyze_sync_file("E_coli.tbg", "E_coli.sync", "sync_results.tsv",0,2)
