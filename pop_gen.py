from helper_functions import *
import searching
import copy

def analyze_sync_file(nucleotide_file, sync_file, output, restfile=None, threads=1,low_ram=False, stat_file=None, verbose=False):
    sync_data = {i:{} for i in sync_file}
    snps = []

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

                if len(columns) <= 3:
                    continue

                populations = len(columns[3:])
                snps.append((columns[0], int(columns[1])))
                sync_data[sync_f][(columns[0], int(columns[1]), columns[2])] = columns[3:]

    if stat_file is not None:
        if len(stat_file) > 1:
            four_ds = {j:[0 for i in range(populations)] for j in sync_file}
            syn_changes = {j:[0 for i in range(populations)] for j in sync_file}
            non_syn_changes = {j:[0 for i in range(populations)] for j in sync_file}
            ATCG_dict = {"A": 6, "T": 7, "C": 8, "G": 9}
        else:

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
                    outfile.write("\t".join([str(keys[0]), str(keys[1]), str(keys[2])]))
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
                            for i in range(4):
                                if line[ATCG_dict[keys[2]]] !=i and line[ATCG_dict[keys[2]]] == line[6+i]:
                                    if len(stat_file) > 1:
                                        syn_changes[sync_f][pop_counter] += current_pop[i]
                                    else:
                                        syn_changes[pop_counter] += current_pop[i]
                                else:
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

def analyze_vcf_file(nucleotide_file, vcf_file, output, restfile=None, threads=1,low_ram=False, stat_file=None, verbose=False):
    vcf_data = {i:{} for i in vcf_file}
    snps = []
    ATCG_dict = {"A": 6, "T": 7, "C": 8, "G": 9}

    if verbose:
        print("Loading vcf file!")
    for vcf_f in vcf_file:
        with open(vcf_f, "r") as vcf:
            for line in vcf:
                if line[0] == "#":
                    continue
                line = line.strip()
                columns = line.split("\t")
                snps.append((columns[0], int(columns[1])))
                vcf_data[vcf_f][(columns[0], int(columns[1]), columns[3], columns[4],columns[2])] = columns[5:]

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
                                if len(stat_file) > 1:
                                    not_snps_variants[vcf_f] += 1
                                else:

                                    not_snps_variants += 1
                                continue
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
                        alt_aa = ",".join(alt_aa)
                        outfile.write("\t".join([str(keys[0]),str(keys[1])]) + f"\t{line[4]}\t{ref_aa}\t{alt_aa}\t" + "\t".join(vcf_data[vcf_f][keys]) + "\n")
                    else:
                        if len(stat_file) > 1:
                            not_snps_variants[vcf_f] += 1
                        else:
                            not_snps_variants += 1

            else:
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



if __name__ == "__main__":
    pass
    # analyze_vcf_file("crip/crip.tbg", ["crip/example_vcf_file.vcf","crip/example_vcf_file.vcf"], ["crip/crip_vcf_results.tsv", "crip/crip_vcf_results2.tsv"],restfile=["crip/crip_vcf_rest.tsv","crip/crip_vcf_rest2.tsv"],threads=4,low_ram=True,stat_file=["stats_test_vcf.tsv","stats_test_vcf2.tsv"])
    analyze_sync_file("crip/crip.tbg", ["test.sync","test.sync"], ["crip/crip_sync_results.tsv","crip/crip_sync_results2.tsv"],
                     restfile=["crip/crip_sync_rest.tsv","crip/crip_sync_rest2.tsv"], threads=4, low_ram=True, stat_file=["stats_test_sync.tsv","stats_test_sync2.tsv"])

# pop_gen -n crip/crip.tbg -o crip/crip_results_vcf.tsv -vcf crip/example_vcf_file.vcf crip/example_vcf_file2.vcf -sf crip/crip_stats_vcf.tsv -r crip/crip_vcf_rest.vcf crip/crip_vcf_rest2.vcf -w -t 4 -v
# pop_gen -n crip/crip.tbg -o crip/crip_results_sync.tsv crip/crip_results_sync2.tsv -s test.sync test2.sync -sf crip/crip_stats_sync.tsv crip/crip_stats_sync2.tsv -r crip/crip_sync_rest.sync crip/crip_sync_rest2.sync -w -t 4 -v