import os
import sys
import pandas as pd
import glob
import re
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    "-i",
    "--indir",
    dest="indir",
    help="""input directory with results.tsv and files""",
)
parser.add_argument("-o", "--outf", dest="outf", help="""output file (excel)""")
parser.add_argument("--study", dest="study", help="""name of study""")
parser.add_argument("--genome", dest="genome", help="""genome assembly""")
parser.add_argument(
    "--no_links",
    dest="no_links",
    action="store_true",
    default=False,
    help="""do not create links in excel""",
)

args = parser.parse_args()

dfs = [f for f in os.listdir(".") if f.endswith(".tsv")]

if not args.no_links:
    hub_url = "https://bimsbstatic.mdc-berlin.net/hubs/collab_bobermay_kdelaro/{0}/hub.txt".format(
        args.study
    )
    session_url = "https://bimsbstatic.mdc-berlin.net/hubs/collab_bobermay_kdelaro/{0}_session.txt".format(
        args.genome
    )
    ucsc_link = (
        "http://genome.mdc-berlin.net/cgi-bin/hgTracks?"
        "db={0}&"
        "hubClear={1}&"
        "hgS_loadUrlName={2}&"
        "hsS_doLoadUrl=submit"
    ).format(args.genome, hub_url, session_url)


def widen_coords(coords, w):
    chrom, start, end = re.split("[:-]", coords)
    return "{0}:{1}-{2}".format(chrom, int(start) - w, int(end) + w)


with pd.ExcelWriter(args.outf) as writer:
    for rf in glob.glob(os.path.join(args.indir, "*.tsv")):
        name = rf.split("/")[-1].split("_results.tsv")[0]
        try:
            tmp = pd.read_csv(rf, sep="\t", header=0)
            if not args.no_links:
                tmp["switch_coords"] = tmp["switch_coords"].apply(
                    lambda x: '=HYPERLINK("{0}","{1}")'.format(
                        ucsc_link + "&position=" + widen_coords(x, 100), x
                    )
                )
                tmp["insert_coords"] = tmp["insert_coords"].apply(
                    lambda x: '=HYPERLINK("{0}","{1}")'.format(
                        ucsc_link + "&position=" + widen_coords(x, 100), x
                    )
                )
            tmp.to_excel(writer, sheet_name=name, index=False)
            sys.stderr.write("adding inserts for {0}\n".format(name))
        except pd.errors.EmptyDataError:
            sys.stderr.write("no inserts for {0}\n".format(name))
            pass
