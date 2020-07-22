from helper_functions import *
import searching
import copy
import os

def analyze_sync_file(nucleotide_file, sync_file, output, restfile=None, threads=1,low_ram=False, stat_file=None, verbose=False):
    sync_data = {i:{} for i in sync_file}
    snps_dict = {}

    populations = 0
    if verbose:
        print("Loading sync file!")


    for sync_f in sync_file:
        with open(sync_f, "r") as sync:
            for line in sync:
                if line[0] == "#":
                    continue
                line = line.strip()
                columns = line.split("\t")

                if len(columns) <= 3 or columns[2] == "N":
                    continue

                populations = len(columns[3:])
                snps_dict[(columns[0], int(columns[1]))] = True
                sync_data[sync_f][(columns[0], int(columns[1]), columns[2])] = columns[3:]
    snps = []
    for snp in snps_dict.keys():
        snps.append(snp)

    if stat_file is not None:
        if len(stat_file) > 1:
            four_ds = {j:[0 for i in range(populations)] for j in sync_file}
            syn_changes = {j:[0 for i in range(populations)] for j in sync_file}
            non_syn_changes = {j:[0 for i in range(populations)] for j in sync_file}
            ATCG_dict = {"A": 8, "T": 9, "C": 10, "G": 11}
        else:

            four_ds = [0 for i in range(populations)]
            syn_changes = [0 for i in range(populations)]
            non_syn_changes = [0 for i in range(populations)]
            ATCG_dict = {"A": 8, "T": 9, "C": 10, "G": 11}

    if verbose:
        print("Loading tbg file!")

    if low_ram:
        rest, found = searching.check_snps_low_ram(nucleotide_file, snps, verbose=False)
    else:
        rest, found = searching.check_snps_normal(nucleotide_file, snps, threads, verbose=False)

    if verbose:
        print("Write new sync file!")
    for file in output:
        if os.path.isfile(file):
            os.remove(file)

    if restfile is not None:
        for file in restfile:
            if os.path.isfile(file):
                os.remove(file)
    for sync_f, counter in zip(sync_file, range(len(sync_file))):
        if len(output) > 1:
            outfile = open(output[counter], "w")
        else:
            outfile = open(output[0], "a")
        if restfile is not None:
            if len(restfile) >1:
                rest_output = open(restfile[counter], "w")
            else:
                rest_output = open(restfile[0], "a")
        for keys in sync_data[sync_f].keys():
            if (keys[0], keys[1]) in found:
                for line in found[(keys[0],keys[1])]:
                    scaf_pos = [keys[0], keys[1]]
                    scaf_pos.extend(line)
                    line = scaf_pos
                    outfile.write("\t".join([str(keys[0]), str(keys[1]), line[5]]))

                    pop_counter = 0

                    for population in sync_data[sync_f][keys]:
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
                                if len(stat_file) > 1:
                                    four_ds[sync_f][pop_counter] += 1
                                else:
                                    four_ds[pop_counter] += 1
                            for i,j in zip(range(4), "ATCG"):
                                if keys[2] != j and line[ATCG_dict[keys[2]]] == line[8+i]:

                                    if len(stat_file) > 1:
                                        syn_changes[sync_f][pop_counter] += current_pop[i]
                                    else:
                                        syn_changes[pop_counter] += current_pop[i]
                                elif keys[2] != j:

                                    if len(stat_file) > 1:
                                        non_syn_changes[sync_f][pop_counter] += current_pop[i]
                                    else:
                                        non_syn_changes[pop_counter] += current_pop[i]
                            if len(stat_file) > 1:
                                non_syn_changes[sync_f][pop_counter] += current_pop[4]
                            else:
                                non_syn_changes[pop_counter] += current_pop[4]
                        pop_counter += 1
                    outfile.write("\n")
            else:
                if restfile is not None:

                    out = [str(i) for i in keys]
                    out.extend(sync_data[sync_f][keys])
                    rest_output.write("\t".join(out))
                    rest_output.write("\n")

    if stat_file is not None:
        if len(stat_file) == 1:
            with open(stat_file[0], "w") as f:
                f.write("4ds_positions\t" + "\t".join([str(i) for i in four_ds]) + "\n")
                f.write("synonymous_changes\t" + "\t".join([str(i) for i in syn_changes]) + "\n")
                f.write("non_synonymous_changes\t" + "\t".join([str(i) for i in non_syn_changes]) + "\n")
        else:
            for sync_f, counter in zip(sync_file, range(len(sync_file))):
                with open(stat_file[counter], "w") as f:
                    f.write("4ds_positions\t" + "\t".join([str(i) for i in four_ds[sync_f]]) + "\n")
                    f.write("synonymous_changes\t" + "\t".join([str(i) for i in syn_changes[sync_f]]) + "\n")
                    f.write("non_synonymous_changes\t" + "\t".join([str(i) for i in non_syn_changes[sync_f]]) + "\n")
    outfile.close()

def analyze_vcf_file(nucleotide_file, vcf_file, output, restfile=None, threads=1,low_ram=False, stat_file=None, verbose=False, individual_stat_file=False):
    vcf_data = {i:{} for i in vcf_file}
    snps_dict = {}
    ATCG_dict = {"A": 6, "T": 7, "C": 8, "G": 9}

    if individual_stat_file:
        ind_names = {}
    if verbose:
        print("Loading vcf file!")

    for vcf_f in vcf_file:
        with open(vcf_f, "r") as vcf:
            for line in vcf:
                if line[:2] == "##":
                    continue
                if line[0] == "#" and individual_stat_file is not None:
                    line = line.strip()
                    line = line.split()
                    ind_names[vcf_f] = line[9:]
                    snps_not_in_gene_ind = {i: [0 for _ in ind_names[vcf_f]] for i in vcf_file}
                    syn_changes_ind = {i: [0 for _ in ind_names[vcf_f]] for i in vcf_file}
                    non_syn_changes_ind = {i: [0 for _ in ind_names[vcf_f]] for i in vcf_file}
                    snps_in_gene_ind = {i: [0 for _ in ind_names[vcf_f]] for i in vcf_file}
                    not_snps_variants_ind = {i: [0 for _ in ind_names[vcf_f]] for i in vcf_file}
                    no_information_about_variants = {i: [0 for _ in ind_names[vcf_f]] for i in vcf_file}
                    positions_with_ref_nucl = {i: [0 for _ in ind_names[vcf_f]] for i in vcf_file}
                    continue
                line = line.strip()
                columns = line.split("\t")
                if len(columns) <= 3:
                    continue
                snps_dict[(columns[0], int(columns[1]))] = True
                vcf_data[vcf_f][(columns[0], int(columns[1]), columns[3], columns[4],columns[2])] = columns[5:]
    snps = []
    for snp in snps_dict.keys():
        snps.append(snp)
    if stat_file is not None:
        if len(stat_file) > 1:
            snps_not_in_gene = {i:0 for i in vcf_file}
            syn_changes = {i:0 for i in vcf_file}
            non_syn_changes = {i:0 for i in vcf_file}
            snps_in_gene = {i:0 for i in vcf_file}
            not_snps_variants = {i:0 for i in vcf_file}
        else:

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

    for file in output:
        if os.path.isfile(file):
            os.remove(file)

    if restfile is not None:
        for file in restfile:
            if os.path.isfile(file):
                os.remove(file)

    for vcf_f, counter in zip(vcf_file, range(len(vcf_file))):
        if len(output) > 1:
            outfile = open(output[counter], "w")
        else:
            outfile = open(output[0], "a")
        if restfile is not None:
            if len(restfile) > 1:
                rest_output = open(restfile[counter], "w")
            else:

                rest_output = open(restfile[0], "a")
        for keys in vcf_data[vcf_f].keys():
            if (keys[0], keys[1]) in found:
                for line in found[(keys[0],keys[1])]:
                    if len(keys[2]) == 1 and keys[2] in "ATGC":
                        alt_aa = []
                        ref_aa = line[ATCG_dict[keys[2]]]

                        for alt in keys[3].split(","):

                            if len(alt) != 1 or alt not in "ATGC":
                                alt_aa.append("-")
                                if stat_file is not None:
                                    if len(stat_file) > 1:
                                        not_snps_variants[vcf_f] += 1
                                    else:

                                        not_snps_variants += 1
                                continue
                            if stat_file is not None:
                                if len(stat_file) > 1:
                                    snps_in_gene[vcf_f] += 1
                                else:
                                    snps_in_gene += 1
                            alt_aa.append(line[ATCG_dict[alt]])
                            if stat_file is not None:
                                if line[ATCG_dict[alt]] == ref_aa:
                                    if len(stat_file) > 1:
                                        syn_changes[vcf_f] += 1
                                    else:
                                        syn_changes += 1
                                else:
                                    if len(stat_file) > 1:
                                        non_syn_changes[vcf_f] += 1
                                    else:
                                        non_syn_changes += 1
                        alt_aa_string = ",".join(alt_aa)

                        outfile.write("\t".join([str(keys[0]),str(keys[1])]) + "\t" + line[4] + f"\t{line[4]}\t{ref_aa}\t{alt_aa_string}\t" + "\t".join(vcf_data[vcf_f][keys]) + "\n")
                    else:
                        if stat_file is not None:
                            if len(stat_file) > 1:
                                not_snps_variants[vcf_f] += 1
                            else:
                                not_snps_variants += 1
                if individual_stat_file:
                    line = vcf_data[vcf_f][keys]
                    gt_position = None
                    formats = line[3].split(":")
                    for pos, format in enumerate(formats):
                        if format == "GT":
                            gt_position = pos
                            break
                    if gt_position is None:
                        continue
                    for ind_counter, ind in enumerate(line[4:]):

                        information = ind.split(":")
                        information[gt_position] = information[gt_position].replace("|", "/")
                        allels = information[gt_position].split("/")
                        allel_counter = 0
                        for allel in allels:
                            if allel == "." or allel == "*":
                                no_information_about_variants[vcf_f][ind_counter] += 1
                                continue
                            if allel_counter == 0:
                                snps_in_gene_ind[vcf_f][ind_counter] += 1
                            allel_counter += 1

                            if allel == "0":
                                positions_with_ref_nucl[vcf_f][ind_counter] += 1
                                continue
                            else:
                                allel = int(allel)
                                if ref_aa == alt_aa[allel-1]:
                                    syn_changes_ind[vcf_f][ind_counter] += 1
                                elif alt_aa[allel-1] == "-":
                                    not_snps_variants_ind[vcf_f][ind_counter] += 1
                                else:
                                    non_syn_changes_ind[vcf_f][ind_counter] += 1

            else:
                if individual_stat_file:
                    for ind_counter, ind in enumerate(ind_names[vcf_f]):
                        snps_not_in_gene_ind[vcf_f][ind_counter] += 1

                if stat_file is not None:
                    if len(stat_file) > 1:
                        snps_not_in_gene[vcf_f] += 1
                    else:
                        snps_not_in_gene += 1
                if restfile is not None:
                    out = [str(i) for i in [keys[0], keys[1], keys[4], keys[2], keys[3]]]
                    out.extend(vcf_data[vcf_f][keys])
                    rest_output.write("\t".join(out))
                    rest_output.write("\n")
        outfile.close()

    if stat_file is not None:
        if len(stat_file) == 1:

            with open(stat_file[0], "w") as f:
                f.write("SNPs found in coding sequences\t" +str(snps_in_gene)+ "\n")
                f.write("SNPs not found in coding sequences\t" + str(snps_not_in_gene) + "\n")
                f.write("variants which are not SNPs\t" + str(not_snps_variants) + "\n")
                f.write("synonymous_changes\t" + str(syn_changes) + "\n")
                f.write("non_synonymous_changes\t" + str(non_syn_changes) + "\n")
        else:
            for vcf_f, counter in zip(vcf_file, range(len(vcf_file))):
                with open(stat_file[counter], "w") as f:
                    f.write("SNPs found in coding sequences\t" + str(snps_in_gene[vcf_f]) + "\n")
                    f.write("SNPs not found in coding sequences\t" + str(snps_not_in_gene[vcf_f]) + "\n")
                    f.write("variants which are not SNPs\t" + str(not_snps_variants[vcf_f]) + "\n")
                    f.write("synonymous_changes\t" + str(syn_changes[vcf_f]) + "\n")
                    f.write("non_synonymous_changes\t" + str(non_syn_changes[vcf_f]) + "\n")

    if individual_stat_file:
        for vcf_f in vcf_file:

            if "\\" in vcf_f:
                file = vcf_f.split("\\")
            else:
                file = vcf_f.split("/")
            with open(file[-1] + "_tbg_ind_stats.tsv", "w") as f:
                indeces = ["Variants found in coding sequences with information",
                           "Variants not found in coding sequences",
                            "variants which are not SNPs",
                            "synonymous_changes",
                            "non_synonymous_changes",
                            "no_information_about_variant",
                            "positions_with_ref_nucleotide",
                            "total dN/dS"
                            ]
                stats = [snps_in_gene_ind[vcf_f],
                         snps_not_in_gene_ind[vcf_f],
                         not_snps_variants_ind[vcf_f],
                         syn_changes_ind[vcf_f],
                         non_syn_changes_ind[vcf_f],
                         no_information_about_variants[vcf_f],
                         positions_with_ref_nucl[vcf_f],
                         ]
                f.write("stats\t")
                f.write("\t".join(ind_names[vcf_f]) + "\n")

                for i in range(len(stats)):
                    f.write(indeces[i] + "\t" + "\t".join([str(stat) for stat in  stats[i]]) + "\n")
                f.write(indeces[-1] + "\t" + "\t".join(str(a/b) for a,b in zip(non_syn_changes_ind[vcf_f], syn_changes_ind[vcf_f])))



if __name__ == "__main__":
    pass
