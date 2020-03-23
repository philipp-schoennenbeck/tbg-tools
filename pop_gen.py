from helper_functions import *
import searching
import copy

def analyze_sync_file(nucleotide_file, sync_file, output, restfile=None, threads=1,low_ram=False, stat_file=None, verbose=False):
    sync_data = {}
    snps = []

    populations = 0
    if verbose:
        print("Loading sync file!")
    with open(sync_file, "r") as sync:
        for line in sync:

            line = line.strip()
            columns = line.split("\t")

            if len(columns) <= 3:
                continue

            populations = len(columns[3:])
            snps.append((columns[0], int(columns[1])))
            sync_data[(columns[0], int(columns[1]), columns[2])] = columns[3:]

    if stat_file is not None:
        four_ds = [0 for i in range(populations)]
        syn_changes = [0 for i in range(populations)]
        non_syn_changes = [0 for i in range(populations)]
        ATCG_dict = {"A": 6, "T": 7, "C": 8, "G": 9}
    if verbose:
        print("Loading tbg file!")
    if low_ram:
        rest, found = searching.check_snps_low_ram(nucleotide_file, snps, verbose=False)
    else:
        rest, found = searching.check_snps_normal(nucleotide_file, snps, threads, verbose=False)
    if verbose:
        print("Write new sync file!")
    outfile = open(output, "w")
    if restfile is not None:
        rest_out = open(restfile,"w")
    for keys in sync_data.keys():

        if (keys[0], keys[1]) in found:
            for line in found[(keys[0],keys[1])]:
                scaf_pos = [keys[0], keys[1]]
                scaf_pos.extend(line)
                line = scaf_pos
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
                            amino_acids.append(str(current_pop[base]) + line[8 + base])
                    if all([x == line[8] for x in line[8:]]):
                        amino_acids.append(str(current_pop[4]) + line[8])
                    else:
                        amino_acids.append(str(current_pop[4]) + "-")
                    amino_acids.append(str(current_pop[5]))
                    amino_acids = ":".join(amino_acids)
                    outfile.write("\t" + amino_acids)
                    if stat_file is not None:
                        if all([line[9+i] == line[8] for i in range(3)]):
                            four_ds[pop_counter] += 1
                        for i in range(4):
                            if line[ATCG_dict[keys[2]]] !=i and line[ATCG_dict[keys[2]]] == line[6+i]:
                                syn_changes[pop_counter] += current_pop[i]
                            else:
                                non_syn_changes[pop_counter] += current_pop[i]
                        non_syn_changes[pop_counter] += current_pop[4]
                    pop_counter += 1
                outfile.write("\n")
        else:
            if restfile is not None:

                out = [str(i) for i in keys]
                out.extend(sync_data[keys])
                rest_out.write("\t".join(out))
                rest_out.write("\n")

    if stat_file is not None:
        with open(stat_file, "w") as f:
            f.write("4ds_positions\t" + "\t".join([str(i) for i in four_ds]) + "\n")
            f.write("synonymous_changes\t" + "\t".join([str(i) for i in syn_changes]) + "\n")
            f.write("non_synonymous_changes\t" + "\t".join([str(i) for i in non_syn_changes]) + "\n")
    outfile.close()

def analyze_vcf_file(nucleotide_file, vcf_file, output, restfile=None, threads=1,low_ram=False, stat_file=None, verbose=False):
    vcf_data = {}
    snps = []
    ATCG_dict = {"A": 6, "T": 7, "C": 8, "G": 9}

    if verbose:
        print("Loading vcf file!")
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line[0] == "#":
                continue
            line = line.strip()
            columns = line.split("\t")

            snps.append((columns[0], int(columns[1])))
            vcf_data[(columns[0], int(columns[1]), columns[3], columns[4],columns[2])] = columns[5:]

    if stat_file is not None:
        snps_not_in_gene = 0
        syn_changes = 0
        non_syn_changes = 0
        snps_in_gene = 0
        not_snps_variants = 0
    if verbose:
        print("Loading tbg file!")
    if low_ram:
        rest, found= searching.check_snps_low_ram(nucleotide_file,snps,verbose=False)
    else:
        rest, found = searching.check_snps_normal(nucleotide_file, snps, threads, verbose=False)
    if verbose:
        print("Write new vcf file!")
    outfile = open(output, "w")
    if restfile is not None:
        rest_output = open(restfile, "w")
    for keys in vcf_data.keys():

        if (keys[0], keys[1]) in found:

            for line in found[(keys[0],keys[1])]:

                if len(keys[2]) == 1 and keys[2] in "ATGC":
                    alt_aa = []
                    ref_aa = line[ATCG_dict[keys[2]]]

                    for alt in keys[3].split(","):

                        if len(alt) != 1 or alt not in "ATGC":
                            alt_aa.append("-")
                            not_snps_variants += 1
                            continue
                        snps_in_gene += 1
                        alt_aa.append(line[ATCG_dict[alt]])
                        if stat_file is not None:
                            # print(line[ATCG_dict[alt]], ref_aa)
                            if line[ATCG_dict[alt]] == ref_aa:
                                syn_changes += 1
                            else:
                                non_syn_changes += 1
                    alt_aa = ",".join(alt_aa)
                    outfile.write("\t".join([str(keys[0]),str(keys[1])]) + f"\t{line[4]}\t{ref_aa}\t{alt_aa}\t" + "\t".join(vcf_data[keys]) + "\n")
                else:
                    not_snps_variants += 1

        else:
            snps_not_in_gene += 1
            if restfile is not None:
                out = [str(i) for i in [keys[0], keys[1], keys[4], keys[2], keys[3]]]
                out.extend(vcf_data[keys])
                rest_output.write("\t".join(out))
                rest_output.write("\n")
    outfile.close()

    if stat_file is not None:
        with open(stat_file, "w") as f:
            f.write("SNPs found in coding sequences\t" +str(snps_in_gene)+ "\n")
            f.write("SNPs not found in coding sequences\t" + str(snps_not_in_gene) + "\n")
            f.write("variants which are not SNPs\t" + str(not_snps_variants) + "\n")
            f.write("synonymous_changes\t" + str(syn_changes) + "\n")
            f.write("non_synonymous_changes\t" + str(non_syn_changes) + "\n")



if __name__ == "__main__":
    pass
    # analyze_vcf_file("crip/crip.tbg", "crip/example_vcf_file.vcf", "crip/crip_vcf_results.tsv",restfile="crip/crip_vcf_rest.tsv",threads=4,low_ram=True,stat_file="stats_test_vcf.tsv")
    analyze_sync_file("crip/crip.tbg", "test.sync", "crip/crip_sync_results.tsv",
                     restfile="crip/crip_sync_rest.tsv", threads=4, low_ram=True, stat_file="stats_test_sync.tsv")
