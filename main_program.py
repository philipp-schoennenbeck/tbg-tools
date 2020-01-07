import create_file
import searching
import argparse
import sys
import os.path


if __name__ == "__main__":
    description = "This is a program to work with tbg files.\n" \
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
            parser.add_argument("-hr", "--human_readable", help="creates a human readable file aswell",
                                action="store_true")
            parser.add_argument("-aa", "--amino_acid_codes", help="selects the amino acid code from aa_codes.txt,"
                                                                  " default is default", default="default")
            parser.add_argument("-p", "--protein_file", help="path to an extra protein fasta file", default=None)
            parser.add_argument("-hro", "--human_readable_outfile", help="path to the human readable file,"
                                                                         "default is \"outfile\"_hr.tsv")
            parser.add_argument("-v", "--verbose", help="increases verbosity", action="store_true")
            parser.add_argument("-t", "--threads", help="number of threads to be used", default=1, type=int)
            args = parser.parse_args(sys.argv[2:])


            if args.outfile:
                outfile = args.outfile
            else:
                outfile = args.gff[:-3] + "tbg"
            if args.human_readable:
                hr = True
                if args.human_readable_outfile:
                    hro = args.human_readable_outfile
                else:
                    hro = outfile[:-4] + "_hr.tsv"
            else:
                hr = False
                hro = ""
            if not os.path.isfile(args.gff):
                raise FileNotFoundError(f"gff file was not found (\"{args.gff}\")")
            if not os.path.isfile(args.fasta):
                raise FileNotFoundError(f"fasta file was not found (\"{args.fasta}\")")
            create_file.create_the_file(args.gff, args.fasta, outfile_bin=outfile, outfile_hr=hro,
                                        verbose=args.verbose, create_binary=True, write_tsv=hr,
                                        aa_code=args.amino_acid_codes, threads=args.threads, protein=args.protein_file)
        elif sys.argv[1] == "search":
            parser.add_argument("-f", "--tbg_file", help="path to the tbg file created with \"create\"", required=True)
            snp_group = parser.add_mutually_exclusive_group(required=True)
            snp_group.add_argument("-b", "--bed", help="path to tab seperated bed file with the SNPs, \"scaffold  position\"")
            snp_group.add_argument("-s", "--snps", help="list of SNPs separated by space e.g. \"scaffold1,position1"
                                                     "scaffold2,position2 \" ", nargs="+")
            parser.add_argument("-o", "--outfile", help="path to the outfile for relevant SNPs,"
                                                        " default is \"SNPs.tsv\"", default="SNPs.tsv")
            parser.add_argument("-r", "--rest", help="path to the file with all none relevant SNPs, default is no file")
            parser.add_argument("-v", "--verbose", help="increases verbosity", action="store_true")
            parser.add_argument("-t", "--threads", help="number of threads to be used", default=1, type=int)
            args = parser.parse_args(sys.argv[2:])

            if not os.path.isfile(args.tbg_file):
                raise FileNotFoundError(f"tbg file was not found (\"{args.tbg_file}\")")
            if args.bed:
                if not os.path.isfile(args.bed):
                    raise FileNotFoundError(f"bed file was not found (\"{args.bed}\")")
                searching.check_snps(args.tbg_file, snp_file=args.bed, binary=True, outfile=args.outfile,
                                     rest_file=args.rest, threads=args.threads)
            elif args.snps:
                snps = [i.split(",") for i in args.snps]
                searching.check_snps(args.tbg_file, snps=snps, binary=True, outfile=args.outfile,
                                     rest_file=args.rest, threads=args.threads)
        elif sys.argv[1] == "convert":
            parser.description = "Converts the tbg file to a human readable tsv file." \
                                 " These files can get very big"
            parser.add_argument("-n", "--tbg_file", help="path to the tbg file.", required=True)
            parser.add_argument("-o", "--outfile", help="path the to human readable file, defailt ist the tbg"
                                                       " file with .tsv")
            parser.add_argument("-v", "--verbose", help="increases verbosity", action="store_true")
            args = parser.parse_args(sys.argv[2:])
            if not os.path.isfile(args.tbg_file):
                raise FileNotFoundError(f"tbg file was not found (\"{args.tbg_file}\")")
            create_file.write_human_readable(args.tbg_file, path_hr=args.outfile)
        else:
            args = parser.parse_args()
            if args.version:
                print("TBG v0.1")
            else:
                parser.print_help()
    else:
        parser.print_help()
