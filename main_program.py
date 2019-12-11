import create_file
import helper_functions
import searching
import argparse
import sys

if __name__ == "__main__":
    print(sys.argv)
    description = "This is a program to work with ??? files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-#", "--version", help="print program version", action="store_true")
    parser.add_argument("-v", "--verbose", help="increases verbosity", action="store_true")

    if len(sys.argv) >= 2:
        if sys.argv[1] == "create":
            parser.add_argument("-f", "--fasta", help="path to fasta file", default=False)
            parser.add_argument("-g", "--gff", help="path to gff file", default=False)
            parser.add_argument("-o", "--outfile", help="path to the created file, default is the gff name",
                                default=None)
            parser.add_argument("-hr", "--human_readable", help="reates a human readable file aswell",
                                action="store_true")
            parser.add_argument("-aa", "--amino_acid_codes", help="selects the amino acid code from aa_codes.txt,"
                                                                  " default is default", default="default")

            parser.add_argument("-hro", "--human_readable_outfile", help="path to the human readable file")
            args = parser.parse_args(sys.argv[2:])

            if args.fasta and args.gff:
                if args.outfile:
                    outfile = args.outfile
                else:
                    outfile = args.gff[:-3] + "bin"
                if args.human_readable:
                    hr = True
                    if args.human_readable_outfile:
                        hro = args.human_readable_outfile
                    else:
                        hro = outfile[:-4] + "_hr.tsv"
                else:
                    hr = False
                    hro = ""

                create_file.create_the_file(args.gff, args.fasta, outfile_bin=outfile, outfile_hr=hro,
                                            verbose=args.verbose, create_binary=True, write_tsv=hr,
                                            aa_code=args.amino_acid_codes )
        if sys.argv[1] == "search":
            parser.add_argument("-n", "--nucleotide_file", help="path to the binary nucleotide file create with"
                                                                " \"create\"", required=True)
            parser.add_argument("-b", "--bed", help="path to bed file")
            parser.add_argument("-s", "--snps", help="list of SNPs separated by space e.g. \"scaffold1,position1"
                                                     "scaffold2,position2 \" ", nargs="+")
            parser.add_argument("-o", "--outfile", help="path to the outfile for relevant SNPs,"
                                                        " default is \"SNPs.tsv\"", default="SNPs.tsv")
            parser.add_argument("-r", "--rest", help="path to the file with all none relevant SNPs, default is no file")
            args = parser.parse_args(sys.argv[2:])
            if args.bed:
                searching.check_snps(args.nucleotide_file, snp_file=args.bed, binary=True, outfile=args.outfile,
                                     rest_file=args.rest)
            elif args.snps:
                snps = [i.split(",") for i in args.snps]
                searching.check_snps(args.nucleotide_file, snps=snps, binary=True, outfile=args.outfile,
                                     rest_file=args.rest)
        if sys.argv[1] == "convert":
            pass
