from helper_functions import *

def analyze_sync_file(nucleotide_file, sync_file, output, rest, threads,low_ram=False, stat_file=None):
    sync_data = {}

    data, genes, scaffolds = read_binary_file(nucleotide_file, threads=threads)
    populations = 0
    with open(sync_file, "r") as sync:
        for line in sync:
            line = line.strip()
            columns = line.split("\t")
            if len(columns) <= 3:
                continue
            populations = len(columns[3:])
            sync_data[(columns[0], int(columns[1]), columns[2])] = columns[3:]
    if stat_file is not None:
        four_ds = [0 for i in range(populations)]
        syn_changes = [0 for i in range(populations)]
        non_syn_changes = [0 for i in range(populations)]
        ATCG_dict = {"A": 6, "T": 7, "C": 8, "G": 9}
    outfile = open(output, "w")
    for keys in sync_data.keys():
        if scaffolds[keys[0]] in data:
            if keys[1] in data[scaffolds[keys[0]]]:
                for line in data[scaffolds[keys[0]]][keys[1]]:
                    line = decode_line(line, genes)
                    outfile.write("\t".join([str(keys[0]), str(keys[1]), str(keys[2])]))
                    pop_counter = 0
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
                        if stat_file is not None:
                            if all([line[7+i] == line[6] for i in range(3)]):
                                four_ds[pop_counter] += 1
                            for i in range(4):
                                if line[ATCG_dict[keys[2]]] !=i and line[ATCG_dict[keys[2]]] == line[6+i]:
                                    syn_changes[pop_counter] += current_pop[i]
                                else:
                                    non_syn_changes[pop_counter] += current_pop[i]
                        pop_counter += 1
                    outfile.write("\n")
    if stat_file is not None:
        with open(stat_file, "w") as f:
            f.write("4ds_positions\t" + "\t".join([str(i) for i in four_ds]) + "\n")
            f.write("synonymous_changes\t" + "\t".join([str(i) for i in syn_changes]) + "\n")
            f.write("non_synonymous_changes\t" + "\t".join([str(i) for i in non_syn_changes]) + "\n")
    outfile.close()

if __name__ == "__main__":
    analyze_sync_file("E_coli.tbg", "E_coli.sync", "sync_results.tsv",0,2)
