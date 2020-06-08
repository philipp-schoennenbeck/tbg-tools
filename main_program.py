import create_file
import searching
import pop_gen
import argparse
import sys
import os.path
from datetime import datetime


if __name__ == "__main__":
    # time = datetime.now()
    description = "This is a program to work with tbg files.\n" \
                  "tbg files are files to fast and easily find out if specific positions on scaffolds are within a gene " \
                  "or not and if they are within a gene, you get additional information about this position\n" \
                  "create\tcreate the tbg file with a gff and a fasta file\n" \
                  "search\tsearches in the tbg file for specific SNPs\n" \
                  "convert\tconvert the tbg file to a human readable file"
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-#", "--version", help="print program version", action="store_true")

    if len(sys.argv) >= 2:
        if sys.argv[1] == "create":
            parser.add_argument("-f", "--fasta", help="path to fasta file", required=True)
            parser.add_argument("-g", "--gff", help="path to gff file", required=True)
            parser.add_argument("-o", "--outfile", help="path to the created file, default is the gff name + .tbg",
                                default=None)
            parser.add_argument("-hr", "--human_readable",const="outfile_hr.tsv", default=None, help="creates a human readable file aswell,"
                                                                "default is outfile_hr.tsv",nargs="?", type=str,
                                action="store")
            parser.add_argument("-aa", "--amino_acid_codes", help="selects the amino acid code from aa_codes.txt,"
                                                                  " default is default", default="default")
            parser.add_argument("-p", "--protein_file", help="path to an extra protein fasta file", default=None)
            parser.add_argument("-n", "--gene_sequence_file", help="path to an extra gene sequence file (nucleotide "
                                                                   "sequence)", default=None)

            parser.add_argument("-v", "--verbose", help="increases verbosity", action="store_true", default=False)
            parser.add_argument("-t", "--threads", help="number of threads to be used", default=1, type=int)
            parser.add_argument("-w", "--low_ram", help="option for systems with low RAM, will create some intermediate"
                                                        "files, if used with --threads", action="store_true", default=False)
            parser.add_argument("-i", "--ignore", help="discard genes with length which is not dividable by 3, default"
                                                       " is adding \"A's\" until it is dividable by 3 (poly-A tail)", default=False,
                                action="store_true")
            parser.add_argument( "--filter", help="list of source programs (column two of the gff) to use seperated by "
                                                  "a whitespace", nargs="+", default=None)
            args = parser.parse_args(sys.argv[2:])


            if args.outfile:
                outfile = args.outfile
            else:
                outfile = args.gff[:-3] + "tbg"

            if not os.path.isfile(args.gff):
                raise FileNotFoundError(f"gff file was not found (\"{args.gff}\")")
            if not os.path.isfile(args.fasta):
                raise FileNotFoundError(f"fasta file was not found (\"{args.fasta}\")")
            create_file.create_the_file(args.gff, args.fasta, outfile_bin=outfile, outfile_hr=args.human_readable,
                                        verbose=args.verbose, create_binary=True,
                                        aa_code=args.amino_acid_codes, threads=args.threads, protein=args.protein_file,
                                        low_ram=args.low_ram, write_gene=args.gene_sequence_file,filter=args.filter, ignore=args.ignore)
        elif sys.argv[1] == "search":
            parser.add_argument("-n", "--tbg_file", help="path to the tbg file created with \"create\"", required=True)
            snp_group = parser.add_mutually_exclusive_group(required=True)
            snp_group.add_argument("-b", "--bed", help="path to tab seperated bed file with the SNPs, \"scaffold  position\", "
                                                       "it is also possible to search entire regions with \"scaffold  start   end\","
                                                       "multiple files can be seperated by space", nargs="+")
            snp_group.add_argument("-s", "--snps", help="list of SNPs separated by space e.g. \"scaffold1,position1 "
                                                     "scaffold2,position2 \" ", nargs="+")
            parser.add_argument("-o", "--outfile", help="path to the resulted outfile,"
                                                        " default is \"default.tsv\","
                                                        "multiple files can be seperated by space", nargs="+")
            parser.add_argument("-r", "--rest", help="path to the file with all none relevant SNPs, default is no file,"
                                                     "multiple files can be seperated by space, has to be the same number as the bed files or none", nargs="+")
            parser.add_argument("-v", "--verbose", help="increases verbosity", action="store_true")
            parser.add_argument("-t", "--threads", help="number of threads to be used", default=1, type=int)
            snp_group.add_argument("-g", "--genes", help="list of genes for a human readable file with only these genes,"
                                                      "seperated by space e.g. \"gene1 gene2 gene3\"", nargs="+")
            snp_group.add_argument("-c", "--scaffolds", help="list of scaffolds for a human readable file with only these scaffolds,"
                                                      "seperated by space e.g. \"scaffold1 scaffold2 scaffold3\"", nargs="+")
            parser.add_argument("-w", "--low_ram", help="option for systems with low RAM, will create some intermediate"
                                                        "files, if used with --threads", action="store_true", default=False)
            args = parser.parse_args(sys.argv[2:])

            if not os.path.isfile(args.tbg_file):
                raise FileNotFoundError(f"tbg file was not found (\"{args.tbg_file}\")")
            if args.bed:
                for bedfile in args.bed:
                    if not os.path.isfile(bedfile):
                        raise FileNotFoundError(f"bed file was not found (\"{bedfile}\")")
                number_of_files = len(args.bed)
                if args.rest is not None:
                    rest_files = args.rest
                    if len(rest_files) != number_of_files:
                        raise ValueError("Number of restfiles must match bed files!")
                else:
                    rest_files = None

                if args.outfile is not None:
                    out_files = args.outfile
                    if len(out_files) != 1 and len(out_files) != number_of_files:
                        raise ValueError("Number of outfiles must match bed files or be 1!")
                else:
                    out_files = "default.tsv"
                searching.check_snps(args.tbg_file, snp_file=args.bed, binary=True, outfile=out_files,
                                     rest_file=args.rest, threads=args.threads, low_ram=args.low_ram, verbose=args.verbose)
            elif args.snps:
                snps = [i.split(",") for i in args.snps]
                snps = [(i[0], int(i[1])) for i in snps]
                searching.check_snps(args.tbg_file, snps=snps, binary=True, outfile=args.outfile,
                                     rest_file=args.rest, threads=args.threads, low_ram=args.low_ram, verbose=args.verbose)
            elif args.genes:
                searching.check_gene(args.tbg_file, args.genes, args.outfile, args.verbose, args.rest)
            elif args.scaffolds:
                searching.check_scaffold(args.tbg_file, args.scaffolds, args.outfile, args.verbose, args.rest)
        elif sys.argv[1] == "convert":
            parser.description = "Converts the tbg file to a human readable tsv file." \
                                 " These files can get very big"
            parser.add_argument("-n", "--tbg_file", help="path to the tbg file.", required=True)
            parser.add_argument("-o", "--outfile", help="path the to human readable file, default ist the tbg"
                                                       " file with .tsv")
            parser.add_argument("-v", "--verbose", help="increases verbosity", action="store_true", default=False)
            args = parser.parse_args(sys.argv[2:])
            if not os.path.isfile(args.tbg_file):
                raise FileNotFoundError(f"tbg file was not found (\"{args.tbg_file}\")")
            create_file.write_human_readable(args.tbg_file, path_hr=args.outfile)
        elif sys.argv[1] == "pop_gen":
            parser.description = "Some tools to work with population genetic data files"
            parser.add_argument("-n", "--tbg_file", help="path to the tbg file.", required=True)
            parser.add_argument("-o", "--outfile", help="path the to human readable file, default ist the tbg"
                                                        " file with .tsv, multiple files can be seperated by space"
                                                        "(number of files has to match number of files in -s or -vcf or just a single one for an combined file)",nargs="+")
            parser.add_argument("-v", "--verbose", help="increases verbosity", action="store_true", default=False)
            parser.add_argument("-s", "--sync_file", help="Path to a sync file, multiple files can be seperated by space", default=None, nargs="+")
            parser.add_argument("-vcf", "--vcf_file", help="Path to a vcf file, multiple files can be seperated by space", default=None,nargs="+")

            parser.add_argument("-sf", "--stats_file", help="Path to an optioal stats file to write, multiple files can be seperated by space (number of files has to match number of files in -s or -vcf or just a single one for an combined file)", default=None, nargs="+")
            parser.add_argument("-t", "--threads",help="number of threads to be used", default=1, type=int)
            parser.add_argument("-r", "--rest", help="Path to a potential file with not found positions,  multiple files can be seperated by space (number of files has to match number of files in -s or -vcf or just a single one for an combined file)", default=None,nargs="+")
            parser.add_argument("-w", "--low_ram", help="option for systems with low RAM", action="store_true",
                                default=False)
            args = parser.parse_args(sys.argv[2:])
            if not os.path.isfile(args.tbg_file):
                raise FileNotFoundError(f"tbg file was not found (\"{args.tbg_file}\")")

            number_of_files = 0
            if args.sync_file is not None:
                number_of_files = len(args.sync_file)
            elif args.vcf_file is not None:
                number_of_files = len(args.vcf_file)

            if args.stats_file is not None:
                stats_files = args.stats_file
                if len(stats_files) != 1 and len(stats_files) != number_of_files:
                    raise ValueError("Number of statsfiles must match vcf/sync files or be 1!")
            else:
                stats_files = None
            if args.rest is not None:
                rest_files = args.rest
                if len(rest_files) != 1 and len(rest_files) != number_of_files:
                    raise ValueError("Number of restfiles must match vcf/sync files or be 1!")
            else:
                rest_files = None
            if args.outfile is not None:
                out_files = args.outfile
                if len(out_files) != 1 and len(out_files) != number_of_files:
                    raise ValueError("Number of outfiles must match vcf/sync files or be 1!")
            else:
                out_files = [args.tbg_file[:-4] + "_results.tsv"]

            if args.sync_file is not None:
                sync_files = args.sync_file
                for file in sync_files:
                    if not os.path.isfile(file):
                        raise FileNotFoundError(f"sync file was not found (\"{file}\")")
                pop_gen.analyze_sync_file(args.tbg_file, sync_files, out_files, rest_files, threads=args.threads,
                                          stat_file=stats_files, low_ram=args.low_ram, verbose=args.verbose)
            if args.vcf_file is not None:
                vcf_files = args.vcf_file

                for file in vcf_files:
                    if not os.path.isfile(file):
                        raise FileNotFoundError(f"vcf file was not found (\"{file}\")")
                pop_gen.analyze_vcf_file(args.tbg_file, vcf_files, out_files, rest_files, threads=args.threads,
                                          stat_file=stats_files, low_ram=args.low_ram, verbose=args.verbose)
        else:
            args = parser.parse_args()
            if args.version:
                print("TBG v0.2")
            else:
                parser.print_help()
    else:
        parser.print_help()
    # print((datetime.now() - time).total_seconds())
# c
# search -n E_coli.tbg -o E_coli_specific_gene_tbg.tsv -v -t 4 -g gene-b3268