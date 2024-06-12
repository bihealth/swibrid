import numpy as np
from logzero import logger

ncodes = dict(zip(list("ACGTacgt"), [1, 2, 3, 4, -1, -2, -3, -4]))

IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "K": "[GT]",
    "S": "[GC]",
    "Y": "[CT]",
    "M": "[AC]",
    "W": "[AT]",
    "B": "[CGT]",
    "H": "[ACT]",
    "D": "[AGT]",
    "V": "[ACG]",
    "N": "[ACGT]",
}

RC = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

isotype_colors = {
    "SM": "#000000",
    "SG3": "#FF96BC",
    "SG1": "#FF629B",
    "SA1": "#13C203",
    "SG2": "#C81E5C",
    "SG4": "#76002B",
    "SE": "#E3D913",
    "SA2": "#006709",
    "Sm": "#000000",
    "Sg3": "#FF96BC",
    "Sg1": "#FF629B",
    "Sa": "#13C203",
    "Sg2b": "#C81E5C",
    "Sg2c": "#76002B",
    "Se": "#E3D913",
}


def parse_switch_coords(switch_coords):
    switch_chrom = switch_coords.split(":")[0]
    switch_start, switch_end = map(int, switch_coords.split(":")[1].split("-"))
    switch_orientation = switch_coords.split(":")[2]
    return switch_chrom, switch_start, switch_end, switch_orientation


def read_switch_anno(switch_annotation):
    logger.info("using switch annotation from " + switch_annotation)
    switch_anno = []
    for line in open(switch_annotation):
        ls = line.strip("\n").split("\t")
        chrom = ls[0]
        start = int(ls[1])
        end = int(ls[2])
        switch_anno.append((chrom, start, end) + tuple(ls[3:]))
    switch_anno = sorted(switch_anno, key=lambda x: (x[0], x[1]))
    return switch_anno


def get_switch_coverage(switch_anno, switch_chrom, switch_start, switch_end):
    switch_coverage = np.zeros(switch_end - switch_start, dtype=int)
    anno_recs = []
    for rec in intersect_intervals(
        [(switch_chrom, switch_start, switch_end)], switch_anno, loj=True
    ):
        anno_recs.append(rec)
        start = int(rec[3][1])
        end = int(rec[3][2])
        switch_coverage[(start - switch_start) : (end - switch_start)] += 1

    cov_int = list(nonzero_intervals(switch_coverage, switch_start))
    Ltot = sum(e - s for s, e in cov_int)
    eff_start = cov_int[0][0]
    eff_end = eff_start + Ltot
    return cov_int, Ltot, eff_start, eff_end, anno_recs


def merge_intervals(intervals, add_count=False):
    """interval merging function from here: http://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
    but extended to include chromosome as first entry of tuple, and to optionally count number of overlaps
    """

    sorted_by_lower_bound = sorted(intervals, key=lambda tup: (tup[0], tup[1]))
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher + ((1,) if add_count else ()))
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[1] <= higher[1]
            if higher[0] == lower[0] and higher[1] <= lower[2]:
                upper_bound = max(lower[2], higher[2])
                merged[-1] = (
                    lower[0],
                    lower[1],
                    upper_bound,
                ) + (
                    (merged[-1][3] + 1,) if add_count else ()
                )  # replace by merged interval

            else:
                merged.append(higher + ((1,) if add_count else ()))

    return merged


def intersect_intervals(intervals1, intervals2, loj=False):
    """borrowed from https://codereview.stackexchange.com/questions/178427/given-2-disjoint-sets-of-intervals-find-the-intersections"""

    def getIntersection(interval_1, interval_2):
        # if on different chromosomes
        chrom = interval_1[0]
        if chrom != interval_2[0]:
            return None
        # if on the same chromosome
        start = max(interval_1[1], interval_2[1])
        end = min(interval_1[2], interval_2[2])
        if start < end:
            ret = (chrom, start, end)
            if len(interval_2) > 3:
                ret += (interval_2,)
            return ret
        return None

    iter1 = iter(intervals1)
    iter2 = iter(intervals2)

    try:
        interval1 = next(iter1)
        interval2 = next(iter2)
    except StopIteration:
        return

    while True:
        intersection = getIntersection(interval1, interval2)
        if intersection:
            yield intersection
            try:
                if not loj and intersection[2] == interval1[2]:
                    interval1 = next(iter1)
                else:
                    interval2 = next(iter2)
            except StopIteration:
                return

        try:
            while interval1[0] > interval2[0] or interval1[1] >= interval2[2]:
                interval2 = next(iter2)
            while interval2[0] > interval1[0] or interval2[1] >= interval1[2]:
                interval1 = next(iter1)
        except StopIteration:
            return


def interval_length(intervals):
    return sum((max(i[2], i[1]) - min(i[2], i[1])) for i in intervals)


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing
    (taken and modified for nan values from here: http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python)
    """

    p = np.asfarray(p)
    ok = np.isfinite(p)
    by_descend = p[ok].argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p[ok])) / np.arange(len(p[ok]), 0, -1)
    q = np.zeros_like(p)
    q[ok] = np.minimum(1, np.minimum.accumulate(steps * p[ok][by_descend]))[by_orig]
    q[~ok] = np.nan
    return q


def parse_LAST_pars(pars):
    with open(pars) as inf:
        line = inf.readline()
        while line:
            if line.startswith("# mean delete size"):
                mean_del_size = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# mean insert size"):
                mean_ins_size = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# delOpenProb"):
                p_open_del = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# insOpenProb"):
                p_open_ins = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# delExtendProb"):
                p_extend_del = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# insExtendProb"):
                p_extend_ins = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# delExistCost"):
                cost_open_del = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# insExistCost"):
                cost_open_ins = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# delExtendCost"):
                cost_extend_del = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# insExtendCost"):
                cost_extend_ins = float(line.split(":")[1].strip())
                line = inf.readline()
            if line.startswith("# substitution percent identity"):
                p_sub = float(line.split(":")[1].strip()) / 100.0
                line = inf.readline()

            if line.startswith("# count matrix"):
                ll = inf.readline()
                p_c = []
                for k in range(4):
                    ll = inf.readline()
                    p_c.append(list(map(float, ll.split()[2:])))
                p_c = np.array(p_c)

            if line.startswith("# probability matrix"):
                ll = inf.readline()
                p_mat = []
                for k in range(4):
                    ll = inf.readline()
                    p_mat.append(list(map(float, ll.split()[2:])))
                p_mat = np.array(p_mat)

            if line.startswith("# score matrix"):
                ll = inf.readline()
                if not ll.startswith("#"):
                    break
                cost_mat = []
                for k in range(4):
                    ll = inf.readline()
                    cost_mat.append(list(map(int, ll.split()[2:])))
                cost_mat = np.array(cost_mat)

            line = inf.readline()

    return {
        "mean_del_size": mean_del_size,
        "mean_ins_size": mean_ins_size,
        "p_open_del": p_open_del,
        "p_open_ins": p_open_ins,
        "p_extend_del": p_extend_del,
        "p_extend_ins": p_extend_ins,
        "p_sub": p_sub,
        "p_c": p_c,
        "p_mat": p_mat,
        "cost_open_del": cost_open_del,
        "cost_open_ins": cost_open_ins,
        "cost_extend_del": cost_extend_del,
        "cost_extend_ins": cost_extend_ins,
        "cost_mat": cost_mat,
    }


def decode_match(match_string):
    import re

    tmp = re.split("[@:-]", match_string)
    return (tmp[0], int(tmp[1]), int(tmp[2]))


def decode_coords(coord_string):
    tmp = coord_string.split(":")
    return (
        tmp[0],
        int(tmp[1].split("-")[0]),
        int(tmp[1].split("-")[1]),
        "." if len(tmp) <= 2 else tmp[2],
    )


def decode_insert(insert_string):
    import re

    m = re.match(
        r"insert_"
        r"(?P<switch_chroml>.+):(?P<switch_left>\d+)_"
        r"(?P<gapl>[-\d]+)_"
        r"(?P<istart>\d+)-(?P<iend>\d+)_"
        r"(?P<orientation>[\+\-])_"
        r"(?P<insert_chrom>.+):(?P<insert_start>\d+)-(?P<insert_end>\d+)_"
        r"(?P<gapr>[-\d]+)_"
        r"(?P<switch_chromr>.+):(?P<switch_right>\d+)",
        insert_string,
    )
    return m


def nonzero_intervals(vec, start):
    """
    Find islands of non-zeros in the vector vec
    """
    if len(vec) == 0:
        return []
    elif not isinstance(vec, np.ndarray):
        vec = np.array(vec)

    (edges,) = np.nonzero(np.diff((vec == 0) * 1))
    edge_vec = [edges + 1]
    if vec[0] != 0:
        edge_vec.insert(0, [0])
    if vec[-1] != 0:
        edge_vec.append([len(vec)])
    edges = np.concatenate(edge_vec)
    return zip(start + edges[::2], start + edges[1::2])


def shift_coord(coord, cov_int, ignore=False):
    """shift a genomic coordinate into a relative coordinate over the interval set cov_int"""
    shift = cov_int[0][0]
    for k in range(1, len(cov_int)):
        if coord >= cov_int[k][0]:
            shift += cov_int[k][0] - cov_int[k - 1][1]
        elif coord < cov_int[k][0] and coord > cov_int[k - 1][1]:
            return cov_int[k - 1][1] - shift if ignore else np.nan
        elif coord < cov_int[k][1]:
            break
    return coord - shift


def get_leader_height(Z, C):
    """get leader heights of nodes in linkage for a given clustering"""
    import scipy.cluster.hierarchy
    import pandas as pd

    L, M = scipy.cluster.hierarchy.leaders(Z, C.astype(np.int32))
    h = pd.concat(
        [
            pd.Series(Z[:, 2], index=Z[:, 0].astype(int)),
            pd.Series(Z[:, 2], index=Z[:, 1].astype(int)),
        ],
        axis=0,
    )
    if np.isin(L, h.index).all():
        return pd.Series(h.loc[L].values, index=M)
    else:
        return pd.Series([h.max()], index=M)


def filter_clustering(Z, C, p=0.95, min_size=0):
    """filter a clustering such that only the largest clusters (and those with the smallest leader height)
    containing at least a fraction p of reads are retained"""
    clusts, cinv, csize = np.unique(C, return_inverse=True, return_counts=True)
    heights = get_leader_height(Z, C)
    o = np.lexsort([heights.max() - heights.loc[clusts].values, csize])[::-1]
    o_rev = np.zeros_like(o)
    o_rev[o] = np.arange(len(o))
    remove = ((csize < min_size) | (np.cumsum(csize[o])[o_rev] >= p * len(C)))[cinv]
    C_filtered = np.copy(C)
    C_filtered[remove] = -1
    return C_filtered


def f2(p, x):
    """helper function for fitting a double exponential trend"""
    return p[0] * np.exp(-p[1] * x) + p[2] * np.exp(-p[3] * x)


def res2(p, x, y):
    """helper function for the residual (of the log values) of f2 and y"""
    return np.log(f2(p, x)) - np.log(y)


def get_gap_positions(msa, max_gap=25):
    """find gap positions in MSA (sparse matrix)
    returns a tuple of vectors of length # gaps,  containing
    the read (=row) index, the left and right (=col) positions and the size
    """
    # get indices of nonzero elements
    ii, jj = msa.nonzero()
    # find all "normal" breakpoint positions (between empty and nonzero)
    bps = np.where((np.diff(jj) > 1) & (np.diff(ii) == 0))[0]
    gap_left = jj[bps] + 1
    gap_right = jj[bps] + np.diff(jj)[bps]
    gap_read = ii[bps]
    gap_size = gap_right - gap_left

    assert np.all(gap_size >= 0), "negative gap sizes!"

    return gap_read, gap_left, gap_right, gap_size


def vrange(starts, stops):
    """Create concatenated ranges of integers for multiple start/stop

    Parameters:
        starts (1-D array_like): starts for each range
        stops (1-D array_like): stops for each range (same shape as starts)

    Returns:
        numpy.ndarray: concatenated ranges

    For example:

        >>> starts = [1, 3, 4, 6]
        >>> stops  = [1, 5, 7, 6]
        >>> vrange(starts, stops)
        array([3, 4, 4, 5, 6])

    """
    import numpy as np

    stops = np.asarray(stops)
    ll = stops - starts  # Lengths of each range.
    return np.repeat(stops - ll.cumsum(), ll) + np.arange(ll.sum())


def remove_gaps(msa, gaps=None, max_gap=75, return_sparse=True):
    """removes gaps smaller than max_gap from MSA"""
    import numpy as np
    import scipy.sparse

    if gaps is None:
        read_idx, pos_left, pos_right, gap_size = get_gap_positions(msa)
    else:
        read_idx = gaps["read_idx"]
        gap_left = gaps["gap_left"]
        gap_right = gaps["gap_right"]
        gap_size = gaps["gap_size"]

    remove = gap_size <= max_gap
    gaps_to_remove = gap_size[remove]
    read_idx_to_remove = np.repeat(read_idx[remove], gaps_to_remove)
    pos_idx_to_remove = vrange(gap_left[remove], gap_right[remove])

    if return_sparse:
        msa_cleaned = scipy.sparse.csr_matrix(
            (msa.data // 10, (msa.row, msa.col)), shape=msa.shape, dtype=np.int8
        ).tolil()
        vals_replace = np.repeat(
            np.max(
                np.array(
                    [
                        msa_cleaned[(read_idx, gap_left - 1)].todense().A1,
                        msa_cleaned[(read_idx, gap_right)].todense().A1,
                    ]
                ),
                0,
            )[remove],
            gaps_to_remove,
        )
        msa_cleaned[(read_idx_to_remove, pos_idx_to_remove)] = vals_replace
        return msa_cleaned.tocsr()
    else:
        msa_cleaned = np.zeros(msa.shape, dtype=np.int8)
        msa_cleaned[(msa.row, msa.col)] = msa.data // 10
        vals_replace = np.repeat(
            np.max(
                np.array(
                    [msa_cleaned[(read_idx, gap_left - 1)], msa_cleaned[(read_idx, gap_right)]]
                ),
                0,
            )[remove],
            gaps_to_remove,
        )

        msa_cleaned[(read_idx_to_remove, pos_idx_to_remove)] = vals_replace
        return msa_cleaned


def parse_range(numString: str):
    """parses a range given like so: 1,3,5-8 into a set of integers"""
    import itertools

    def expand_range(s: str):
        z = s.split("-")
        ll = len(z)
        if ll == 1:
            return [int(z[0])]
        elif ll == 2:
            return list(range(int(z[0]), int(z[1]) + 1))
        else:
            raise IndexError("too many values in range!")

    return sorted(set(itertools.chain(*map(expand_range, numString.split(",")))))


def construct_mut_matrix(mutations, nreads, npos, max_bias=0.25):
    """returns a sparse matrix of dimensions nreads x npos
    containing alternative nucleotides encoded as integers
    (uses only positions with strand bias < max_bias
    """
    import scipy.sparse
    import pandas as pd

    take = (mutations["strand_bias"] - 0.5).abs() < max_bias
    if sum(take) > 0:
        i = []
        j = []
        d = []
        for _, row in mutations[take].iterrows():
            for k, x in enumerate("ACGT"):
                c = "alt_" + x
                if not pd.isnull(row[c]):
                    ridx = np.array(list(map(int, row[c].split(","))))
                    use = ridx < nreads
                    i.append(ridx[use])
                    j.append(np.repeat(row["rel_pos"], sum(use)))
                    d.append(np.repeat(k + 1, sum(use)))

        i = np.concatenate(i)
        j = np.concatenate(j)
        d = np.concatenate(d)

        mut = scipy.sparse.csc_matrix((d, (i, j)), shape=(nreads, npos), dtype=np.int8)
    else:
        mut = scipy.sparse.csc_matrix(([], ([], [])), shape=(nreads, npos), dtype=np.int8)

    return mut


def nmf_consistency(U_max, U_all):
    """assesses consistency over list of NMF-induced clusterings"""
    from sklearn.metrics.cluster import contingency_matrix

    gmax = np.argmax(U_max, axis=1)
    c = []
    for u in U_all:
        g = np.argmax(u, axis=1)
        mm = np.argmax(contingency_matrix(gmax, g), axis=1)
        try:
            c.append(mm[g])
        except:
            pass
    return np.mean(np.array(c), 0)


def get_switch_iis(anno_recs, cov_int, binsize):
    """makes vector of switch region annotations over region bins of size binsize"""
    switch_iis = []
    for rec in anno_recs:
        start = shift_coord(int(rec[3][1]), cov_int)
        end = shift_coord(int(rec[3][2]), cov_int)
        switch_iis.append([rec[3][3].upper()] * (end // binsize - start // binsize))
    return np.concatenate(switch_iis)


def calculate_gini(x, w=None):
    """calculate Gini coefficient of array x; see https://stackoverflow.com/questions/48999542/more-efficient-weighted-gini-coefficient-in-python"""
    import numpy as np

    x = np.asarray(x)
    if w is not None:
        w = np.asarray(w)
        sorted_indices = np.argsort(x)
        sorted_x = x[sorted_indices]
        sorted_w = w[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = np.cumsum(sorted_w, dtype=float)
        cumxw = np.cumsum(sorted_x * sorted_w, dtype=float)
        return np.sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) / (cumxw[-1] * cumw[-1])
    else:
        sorted_x = np.sort(x)
        n = len(x)
        cumx = np.cumsum(sorted_x, dtype=float)
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    taken from https://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    """
    if weights.sum() <= 0:
        return (np.nan, np.nan)
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values - average) ** 2, weights=weights)
    return (average, np.sqrt(variance))
