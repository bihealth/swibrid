"""\
collect stats from different samples: this will collect all summary files for the listed samples and produce one big table.
optionally, this will also collect all inserts and produce a multi-sheet excel file.
clustering analysis results can also be collected into an excel file.
"""


def setup_argparse(parser):
    parser.add_argument("--samples", dest="samples", required=True, help="""required: comma-separated list of samples""")
    parser.add_argument(
        "--sample_stats",
        dest="sample_stats",
        required=True,
        help="""required: output file with sample_stats""",
    )
    parser.add_argument(
        "--summary_path",
        dest="summary_path",
        default="pipeline/{sample}/{sample}_summary.csv",
        help="""path pattern for summary files [%(default)s]""",
    )
    parser.add_argument(
        "--inserts_path",
        dest="inserts_path",
        default="pipeline/{sample}/{sample}_inserts.tsv",
        help="""path pattern for insert tables [%(default)s]""",
    )
    parser.add_argument(
        "--cluster_analysis_path",
        dest="cluster_analysis_path",
        default="pipeline/{sample}/{sample}_cluster_analysis.csv",
        help="""path pattern for cluster_analysis files [%(default)s]""",
    )   
    parser.add_argument(
        "--cluster_stats",
        dest="cluster_stats",
        help="""output file with cluster stats""",
    )
    parser.add_argument("--inserts", dest="inserts", help="""output file with inserts""")
    parser.add_argument(
        "--clean_column_names",
        dest="clean_column_names",
        action="store_true",
        help="""re-name columns to make features a bit more explicit""",
    )

column_name_mapping = {
    "nclusters_initial":		"clusters_initial",
    "nclusters_final":			"clusters_raw",
    "nclusters_eff":			"clusters_eff_raw",
    "nreads_final":		        "reads",
    "mean_length":		        "read_length_mean",
    "std_length":		        "read_length_std",
    "mean_GC":		            "GC_content_mean",
    "std_GC":		            "GC_content_std",
    "mean_cluster_size":		"cluster_size_mean_raw",
    "std_cluster_size":			"cluster_size_std_raw",
    "cluster_gini":			    "cluster_gini_raw",
    "cluster_entropy":			"cluster_entropy_raw",
    "cluster_inverse_simpson":	"cluster_inverse_simpson_raw",
    "top_clone_occupancy":		"occupancy_top_cluster_raw",
    "big_clones_occupancy":		"occupancy_big_clusters_raw",
    "nclusters_eff_downsampled":			"clusters_eff",
    "nclusters_final_downsampled":		    "clusters",
    "mean_cluster_size_downsampled":		"cluster_size_mean",
    "std_cluster_size_downsampled":			"cluster_size_std",
    "cluster_gini_downsampled":			    "cluster_gini",
    "cluster_entropy_downsampled":			"cluster_entropy",
    "cluster_inverse_simpson_downsampled":	"cluster_inverse_simpson",
    "top_clone_occupancy_downsampled":		"occupancy_top_cluster",
    "big_clones_occupancy_downsampled":		"occupancy_big_clusters",
    "median_inversion_size":	"inversion_size",
    "median_duplication_size":	"duplication_size",
    "mean_length_SA":		    "read_length_SA_mean",
    "std_length_SA":		    "read_length_SA_std",
    "mean_length_SG":		    "read_length_SG_mean",
    "std_length_SG":		    "read_length_SG_std",
    "mean_length_SM":		    "read_length_SM_mean",
    "std_length_SM":		    "read_length_SM_std",
    "frac_clusters_SA1":		"pct_clusters_SA1",
    "frac_clusters_SA2":		"pct_clusters_SA2",
    "frac_clusters_SG1":		"pct_clusters_SG1",
    "frac_clusters_SG2":		"pct_clusters_SG2",
    "frac_clusters_SG3":		"pct_clusters_SG3",
    "frac_clusters_SM":		    "pct_clusters_SM",		
    "frac_clusters_SG4":		"pct_clusters_SG4",
    "frac_clusters_Sa":		    "pct_clusters_Sa",
    "frac_clusters_Sg1":		"pct_clusters_Sg1",
    "frac_clusters_Sg2b":		"pct_clusters_Sg2b",
    "frac_clusters_Sg2c":		"pct_clusters_Sg2c",
    "frac_clusters_Sg3":		"pct_clusters_Sg3",
    "frac_clusters_Sm":		    "pct_clusters_Sm",		
    "alpha_ratio_reads":		"alpha_ratio_reads",
    "alpha_ratio_clusters":		"alpha_ratio_clusters",
    "insert_frequency":		    "pct_templated_inserts",
    "mean_insert_frequency":    "mean_cluster_insert_frequency",
    "n_untemplated_switch":		"untemplated_inserts",
    "n_untemplated_switch_SA":	"untemplated_inserts_SA",
    "n_untemplated_switch_SG":	"untemplated_inserts_SG",
    "n_untemplated_switch_SM":	"untemplated_inserts_SM",
    "n_homology_switch":		"homology",
    "n_homology_switch_SA":		"homology_SA",
    "n_homology_switch_SG":		"homology_SG",
    "n_homology_switch_SM":		"homology_SM",
    "frac_blunt":	        	"pct_blunt",
    "frac_blunt_SA":		    "pct_blunt_SA",
    "frac_blunt_SG":	     	"pct_blunt_SG",
    "frac_blunt_SM":    		"pct_blunt_SM",
    "frac_breaks_inversions":	"pct_inversions",
    "frac_breaks_duplications":	"pct_duplications",
    "mean_intra_break_size":    "intraswitch_size_mean",
    "std_intra_break_size":	    "intraswitch_size_std",
    "frac_breaks_single":		"pct_direct_switch",
    "frac_breaks_sequential":	"pct_sequential_switch",
    "frac_breaks_intra":		"pct_intraswitch_deletion",
    "frac_breaks_inversions_intra":		"pct_intraswitch_inversion",
    "frac_breaks_duplications_intra":	"pct_intraswitch_duplication",
    "spread_SA1":	            "break_dispersion_SA1",
    "spread_SA2":           	"break_dispersion_SA2",
    "spread_SG1":            	"break_dispersion_SG1",
    "spread_SG2":             	"break_dispersion_SG2",
    "spread_SG3":            	"break_dispersion_SG3",
    "spread_SG4":             	"break_dispersion_SG4",
    "spread_SG2B":	         	"break_dispersion_SG2B",
    "spread_SG2C":	         	"break_dispersion_SG2C",
    "spread_SM":            	"break_dispersion_SM",
    "spread_SA":	         	"break_dispersion_SA",
    "spread_SG":	         	"break_dispersion_SG",
    "spread_SE":                "break_dispersion_SE",
    "frac_breaks_SM_SM":		"pct_breaks_SM_SM",
    "frac_breaks_SM_SE":		"pct_breaks_SM_SE",
    "frac_breaks_SE_SE":		"pct_breaks_SE_SE",
    "frac_breaks_SM_SG":		"pct_breaks_SM_SG",
    "frac_breaks_SG_SG":		"pct_breaks_SG_SG",
    "frac_breaks_SG_SE":		"pct_breaks_SG_SE",
    "frac_breaks_SM_SA":		"pct_breaks_SM_SA",
    "frac_breaks_SG_SA":		"pct_breaks_SG_SA",
    "frac_breaks_SE_SA":		"pct_breaks_SE_SA",
    "frac_breaks_SA_SA":		"pct_breaks_SA_SA",
    "homology_fw":				"homology_score_fw",
    "homology_rv":				"homology_score_rv",
    "homology_fw_SM_SG":		"homology_score_fw_SM_SG",
    "homology_rv_SM_SG":		"homology_score_rv_SM_SG",
    "homology_fw_SM_SA":		"homology_score_fw_SM_SA",
    "homology_rv_SM_SA":		"homology_score_rv_SM_SA",
    "homology_fw_SM_SE":		"homology_score_fw_SM_SE",
    "homology_rv_SM_SE":		"homology_score_rv_SM_SE",
    "c_opt":                    "clustering_cutoff",
    "size_GC_bias":             "PCR_bias_GC",
    "size_length_bias":         "PCR_bias_length",
}

def run(args):
    import pandas as pd
    from logzero import logger
    import glob
    from openpyxl import Workbook

    assert args.samples is not None, "no list of samples given!"
    
    samples = args.samples.split(",")

    dfs = dict(
        (
            sample,
            pd.read_csv(gl, header=None, index_col=0).squeeze().dropna(),
        )
        for sample in samples
        for gl in glob.glob(args.summary_path.format(sample=sample))
    )
    dfs = dict((k, v[v.index.notnull()]) for k, v in dfs.items())
    df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1).T

    if args.clean_column_names:
        logger.info("cleaning up column names")
        df.rename(columns=column_name_mapping, inplace=True)
        for c in df.columns:
            if c.startswith('pct_'):
                df[c] = 100*df[c]
                
    logger.info("saving sample stats to " + args.sample_stats)
    df.to_csv(args.sample_stats)

    if args.inserts is not None:
        wb = Workbook()
        with pd.ExcelWriter(args.inserts, engine="openpyxl") as writer:
            writer.workbook = wb
            for sample in samples:
                for gl in glob.glob(args.inserts_path.format(sample=sample)):
                    try:
                        tmp = pd.read_csv(
                            gl,
                            sep="\t",
                            header=0,
                        )
                    except pd.errors.EmptyDataError:
                        pass
                    tmp.to_excel(writer, sheet_name=sample, index=False)
                    logger.info("adding inserts for {0}".format(sample))

    if args.cluster_stats is not None:
        wb = Workbook()
        with pd.ExcelWriter(args.cluster_stats, engine="openpyxl") as writer:
            writer.workbook = wb
            for sample in samples:
                for gl in glob.glob(args.cluster_analysis_path.format(sample=sample)):
                    tmp = pd.read_csv(
                        gl,
                        header=0,
                        index_col=0,
                    ).sort_values("size", ascending=False)
                    tmp.to_excel(writer, sheet_name=sample, index=True)
                    logger.info("adding clustering stats for {0}".format(sample))
