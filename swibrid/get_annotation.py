"""\
prepare shortened annotation file for insert annotation. input is a GTF file (e.g., from GENCODE), output is a bed file
"""


def setup_argparse(parser):
    parser.add_argument("-i", "--input", help="""input gtf""")
    parser.add_argument("-o", "--out", help="""output bed""")


def run(args):
    with open(args.out, "w") as outf:
        for line in open(args.input):
            if line.startswith("#"):
                continue
            ls = line.strip().split("\t")
            info = dict((l.split()[0], l.split()[1].strip('" ')) for l in ls[8].split(";")[:-1])
            if ls[2] == "gene":
                gene = info["gene_name"]
                outf.write("{0}\t{1}\t{2}\t{3}\n".format(ls[0], int(ls[3]) - 1, ls[4], gene))
            elif ls[2] == "exon":
                gene = info["gene_name"]
                exon_number = info["exon_number"]
                outf.write(
                    "{0}\t{1}\t{2}\t{3}.exon{4}\n".format(
                        ls[0], int(ls[3]) - 1, ls[4], gene, exon_number
                    )
                )
