import concurrent.futures
import os
import subprocess
import threading
import tempfile

from Bio import SeqIO, AlignIO, Seq, Align
from scipy import stats
import pandas as pd

from . import file_util
from . import tree

""" Various alignment generation and manipulation tools. Alignments are generally performed with mafft via calls to the
command line. Biopython objects are used to represent sequences and alignments.
 Note: many functions currently assume alignment files are in Fasta format."""


def custom_mafft_aln(
        in_file,
        out_file,
        structural=True,
        quiet=False
):
    """ Customisable alignment using command line MAFFT. Exposure of different parameters is added here as required. """

    os.system(f"mafft --reorder{' --quiet' if quiet else ''}{' --dash --originalseqonly' if structural else ''}"
              f" {in_file} > {out_file}")


def read_aln(aln, format='fasta'):
    """ Check if aln is already a Bio.Align.MSA object. If not, assume aln is a file and read in. """
    if isinstance(aln, Align.MultipleSeqAlignment):
        return aln
    else:
        return AlignIO.read(aln, format)

def ava_aligned_distance(aln, metric='id'):
    """ Compute pairwise distances based on a given sequence alignment and returns of a matrix of said distances.
    Default distance metric is p-distance.

    TODO: this is currently pairwise identity not distance. """

    aln = read_aln(aln)
    ids = [seq.id for seq in aln]
    dists = pd.DataFrame(index=ids, columns=ids)

    if metric == 'id':
        for i in range(len(ids)):
            dists[aln[i].id][aln[i].id] = 1
            for j in range(i + 1, len(ids)):
                match = 0  # n of non-gap cols with character match
                non_gap = 0  # n of cols with both seqs non-gap
                for col in range(aln.get_alignment_length()):
                    char1 = aln[i][col]
                    char2 = aln[j][col]
                    if char1 != '-' and char2 != '-':
                        non_gap += 1
                        if char1 == char2:
                            match += 1
                dists[aln[i].id][aln[j].id] = dists[aln[j].id][aln[i].id] = match / non_gap


    # Implement further distance metrics here

    else:
        raise RuntimeError(f"Metric '{metric}' is not implemented.")

    return dists


def ava_pairwise_distance(
        seq_file,
        aligner='mafft',
        structural=False,
        metric='id'
):
    """ Perform all-verse-all pairwise alignment of a set of sequences independent of any computed alignment. Returns a
    matrix of distances by a given metric.

    NOTE: the values actually represent similarities, not distances, currently
    """

    if structural and aligner != 'mafft':
        raise RuntimeError(f"Structural alignment is only available using mafft.")

    ids = [seq.id for seq in SeqIO.parse(seq_file, 'fasta')]
    dists = pd.DataFrame(index=ids, columns=ids)

    for i in range(len(ids)):
        dists[ids[i]][ids[i]] = 1
        for j in range(i+1, len(ids)):
            file_util.extract_fasta(seq_file, 'temp.fa', [ids[i], ids[j]], ungap=True)
            if structural:
                struct_aln('temp.fa', 'temp.aln', quiet=True)
            else:
                default_aln('temp.fa', 'temp.aln', quiet=True)

            aln = AlignIO.read('temp.aln', 'fasta')
            pw_dist_matrix = ava_aligned_distance('temp.aln', metric)
            dists[aln[i].id][aln[j].id] = dists[aln[j].id][aln[i].id] = pw_dist_matrix.iloc[0][1]

            os.remove('temp.fa')
            os.remove('temp.aln')

    return dists


# def
#
#
# def split_cols():
#     return
#
# def filter_seqs():
#     """ Filter sequences  """


def remove_gap_only_cols(in_file, out_file):
    """ Remove all columns where all sequences contain a gap. """

    og_aln = read_aln(in_file)

    # Find cols with all gap chars
    gap_cols = []
    for col in range(og_aln.get_alignment_length()):
        characters = [og_aln[i,col] for i in range(len(og_aln))]
        non_gap = False
        for char in characters:
            if char != '-':
                non_gap = True
                break
        if not non_gap:
            gap_cols.append(col)

    new_aln = trim_aln(og_aln, gap_cols)
    AlignIO.write(new_aln, out_file, 'fasta')


def filter_cols(aln, threshold=0.01):
    """ Simple column filtering based purely on occupancy """
    aln = read_aln(aln)
    to_trim = []
    for i in range(aln.get_alignment_length()):
        col = aln[:, i]
        if 1 - col.count('-') / aln.get_alignment_length() < threshold:
            to_trim.append(col)

    return trim_aln(aln, to_trim)


def run_trimal(in_aln, out_aln, gt=0.1):
    """ Run TrimAl """

    subprocess.run(["trimal", "-in", in_aln, "-out", out_aln, "-gt", str(gt)])


def trim_aln(aln, to_trim, out_aln=None):
    """ Trim indices in cols from a multiple sequence alignment. """

    if not isinstance(aln, Align.MultipleSeqAlignment):  # Assume alignment file name provided
        aln = AlignIO.read(aln, "fasta")
    
    # Check that all col indices are in range of sequence length
    if not isinstance(to_trim, list): to_trim = [to_trim]
    to_trim = list(set(to_trim))  # Remove duplicate column indices
    for col in to_trim:
        if not 0 <= col < aln.get_alignment_length():
            raise RuntimeError(f"Invalid column index: {col}")

    # Work out which columns we want
    to_retain = [col for col in range(aln.get_alignment_length()) if col not in to_trim]

    trimmed = []
    for seq in aln:
        trimmed_seq = ''.join([seq.seq[i] for i in to_retain])
        trimmed_record = seq[:]
        trimmed_record.seq = Seq.Seq(trimmed_seq)
        trimmed.append(trimmed_record)
        
    if out_aln:
        SeqIO.write(trimmed, out_aln, "fasta")
        # Remove gap only cols
        os.system(f"trimal -in {out_aln} -out {out_aln} -noallgaps")
    else:
        # TODO: If returning Bio.MSA, currently maintains gap only cols
        return Align.MultipleSeqAlignment(trimmed)


def default_aln(
        in_file,
        out_file,
        quiet=False,
        nice=None,
        max_threads=8
):
    """ Conduct sequence alignment using MAFFT with default parameters. """

    args = []

    if nice:
        args.extend(["nice", "-n", str(nice)])

    args.extend(["mafft", "--reorder", "--anysymbol", "--thread", str(max_threads)])

    if quiet:
        args.extend(["--quiet"])

    args.append(in_file)

    with open(out_file, 'w') as out_file:
        subprocess.run(args, stdout=out_file)


def struct_aln(in_file, out_file, quiet=False):
    """Conduct a structure-guided alignment via MAFFT-DASH."""
    # aln_name = (in_file.split('.fa')[0] if '.fa' in in_file else in_file) if not out_file else out_file.split('.')[0]
    os.system(f"mafft --dash --reorder --originalseqonly {'--quiet ' if quiet else ''}{in_file} > {out_file}")

def linsi_aln(in_file, out_file, quiet=False, nice=0):
    """ Conduct a sequence alignment using MAFFT with the --localpair option. Other parameters are currently default
     only. """
    os.system(f"nice -n {nice} mafft --reorder --localpair {'--quiet ' if quiet else ''}{in_file} > {out_file}")


def sub_aln(
        in_file,
        out_file,
        target_seqs,
        realign=False,
        structural=False,
        quiet=True,
        missing_error=True,
        missing_warning=True
):
    """ Extract and optionally realign a subset of sequences from an existing alignment.

    #TODO: Verify this works as expected when realign set to False """

    new_seqs = []
    for seq in SeqIO.parse(in_file, "fasta"):
        if seq.name in target_seqs:
            new_seqs.append(seq)

    # If required, raise error/print warning if any target sequences are not found
    if len(new_seqs) != len(target_seqs):
        if missing_error:
            raise RuntimeError(f"Sequences {', '.join(list(set(target_seqs)-set([seq.name for seq in new_seqs])))} not found in original "
                               f"alignment.")
        elif missing_warning:
            print(f"Sequences {', '.join(list(set(target_seqs)-set([seq.name for seq in new_seqs])))} not found in original "
                               f"alignment.")
            

    if realign:
        SeqIO.write(new_seqs, f"temp_{out_file}.aln", "fasta")
        custom_mafft_aln(f"temp_{out_file}.aln", out_file, structural=structural, quiet=quiet)
        os.remove(f"temp_{out_file}.aln")

    else:

        # Gap only cols may be present
        retain = []
        for i in range(len(new_seqs[0])):
            content = False
            for seq in new_seqs:
                if seq[i] != '-':
                    content = True
                    break
            if content:   # Retain if at least one sequence with content here
                retain.append(i)
        # Trim all sequences to remove gap only cols
        for seq in new_seqs:
            seq.seq = Seq.Seq(''.join([seq[i] for i in retain]))

        SeqIO.write(new_seqs, out_file, "fasta")


def clade_sub_aln(
        in_file,
        out_file,
        full_tree,
        target_clade_node=None,
        target_clade_bounds=None,
        out_leaf=None,
        reroot=False,
        exclude=None,
        include=None,
        realign=False,
        structural=False,
        quiet=True,
        missing_error=True
):
    """ Construct a sub-alignment from a clade in a tree.  """

    if not exclude:
        exclude = []
    if not include:
        include = []

    full_tree = tree.load_tree(full_tree)  # If a treefile is provided, load tree object
    if reroot:
    # Root tree so clade of interest is monophyletic
        full_tree = tree.root_tree(full_tree, mode="outgroup", og_bounds=target_clade_bounds, in_leaf=out_leaf,
                                   og_node=target_clade_node)

    # Get required leaves
    clade_leaves = tree.get_subtree_leaves(full_tree, anc=target_clade_node, ext_nodes=target_clade_bounds)
    if include:  # Explicitly provided 'include' leaves override 'exclude' leaves
        profile_leaves = [leaf for leaf in clade_leaves if leaf in include]
    else:
        profile_leaves = [leaf for leaf in clade_leaves if leaf not in exclude]

    # Make sub-alignment
    sub_aln(in_file, out_file, profile_leaves, realign=realign, structural=structural, quiet=quiet,
            missing_error=missing_error)


def merge_cluster_alns(
        ref_aln,
        cluster_alns,
        flag_missing=True,
        out_file=None,
        no_return=True
):
    """ Given a reference alignment containing cluster representatives, merge alignments for each cluster. """

    if no_return and not out_file:
        raise RuntimeError("An output alignment file handle must be provided if no_return is True.")

    # Read in alns as BioAlign objects if not already
    if isinstance(ref_aln, str):
        ref_aln = AlignIO.read(ref_aln, "fasta")
    if isinstance(cluster_alns[0], str):
        cluster_alns = [AlignIO.read(aln, "fasta") for aln in cluster_alns]

    print(f"CHECK FIVE: {len(ref_aln)}")
    print(f"CHECK SIX: {len(cluster_alns)}")

    # Should be exactly one reference sequence in each cluster alignment, else assume that it has been filtered
    # and ignore that cluster
    ref_names = [seq.name for seq in ref_aln]  # Order all lists below by order of ref_names
    cluster_reps = [None for ref in ref_names]
    cluster_rep_idxs = [None for ref in ref_names]  # Index for rep in each cluster alignment
    cluster_alns_filt = [None for ref in ref_names]
    rep_filt_idxs = [None for ref in ref_names]  # Index in ref aln of reps which map to clusters
    for aln in cluster_alns:
        cnt = 0  # Check how many reference aln seqs are in this cluster
        for i in range(len(aln)):
            this_name = aln[i].name
            if this_name in ref_names:
                this_ref = this_name
                cluster_idx = i  # Index of rep in cluster aln
                ref_idx = ref_names.index(this_name)
                cnt += 1
        if cnt > 1:
            raise RuntimeError(f"Cluster alignment with first sequence {aln[0].name} contains multiple reference sequences.")
        elif cnt == 1:  # Add data to lists with order matching reference aln
            cluster_reps[ref_idx] = this_ref
            cluster_rep_idxs[ref_idx] = cluster_idx
            cluster_alns_filt[ref_idx] = aln
            rep_filt_idxs[ref_idx] = ref_idx
        elif flag_missing:
            print(f"Note: cluster alignment with first sequence {aln[0].name} does not contain a reference sequence.")

    # Remove cluster reps not present in any cluster from consideration - then they are all consistently ordered
    cluster_reps = [rep for rep in cluster_reps if rep != None]
    cluster_rep_idxs = [idx for idx in cluster_rep_idxs if idx != None]
    cluster_alns_filt = [aln for aln in cluster_alns_filt if aln != None]
    rep_filt_idxs = [idx for idx in rep_filt_idxs if idx != None]

    merged_aln = []  # Store SeqRecords for merged alignment

    for cluster in cluster_alns_filt:  # Create SeqRecord objects with all relevant data
        for seq in cluster:
            merged_aln.append(SeqIO.SeqRecord(Seq.Seq(''), name=seq.name, id=seq.id, description=seq.description))

    col_tracker = [0 for rep in cluster_reps]  # Track which index we're up to in each cluster alignment

    # Determine the last content position for each rep
    rep_ends = []
    for i in range(len(cluster_reps)):
        for j in range(len(ref_aln[0])-1, 0, -1):
            if ref_aln[rep_filt_idxs[i]][j] != '-':
                rep_ends.append(j)
                break

    for i in range(len(ref_aln[0])):  # For each column in reference alignment
        to_append = []  # Character directly mapping to this ref column to append to each seq in merged aln
        for j in range(len(cluster_reps)):  # For rep (and by extension cluster) in ref aln which maps to a cluster
            if ref_aln[rep_filt_idxs[j]][i] != '-':  # There is content in this column
                # Insert content for columns directly preceding this col in cluster aln that don't exist in ref aln
                while True:  # Move col by col through this cluster aln until we find new content in rep seq
                    if cluster_alns_filt[j][cluster_rep_idxs[j]][col_tracker[j]] == '-':  # Col not present in ref aln
                        # Introduce a new column in merged alignment
                        l = 0  # Merged aln seqs index
                        for k in range(len(cluster_reps)):
                            if k == j:  # Extract content for new col from this cluster aln
                                for m in range(len(cluster_alns_filt[k])):
                                    merged_aln[l].seq += cluster_alns_filt[j][m][col_tracker[j]]
                                    l += 1
                            else:  # Add gap everywhere else
                                for m in range(len(cluster_alns_filt[k])):
                                    merged_aln[l].seq += '-'
                                    l += 1
                        col_tracker[j] += 1
                    else:  # This is col that maps to col i of ref aln
                        to_append.extend([cluster_alns_filt[j][n][col_tracker[j]]
                                          for n in range(len(cluster_alns_filt[j]))])
                        col_tracker[j] += 1
                        break  # Found col that maps to ref column i so stop here

            else:  # No content to map for any sequences in this cluster for reference col i
                to_append.extend(['-' for seq in cluster_alns_filt[j]])

        # Append mapped column to merged aln
        for l in range(len(merged_aln)):
            merged_aln[l].seq += to_append[l]

        # For any reps which have no content left, append C-terminal insertions of their cluster aln to merged aln
        for j in range(len(cluster_reps)):
            if rep_ends[j] == i:  # This was last residue of this reference seq
                for k in range(col_tracker[j], len(cluster_alns_filt[j][0])):  # Append these cols to merged aln
                    l = 0  # Index of merged aln
                    for m in range(len(cluster_reps)):
                        if m == j:  # Append content from remainder of cluster alignment
                            for n in range(len(cluster_alns_filt[m])):
                                merged_aln[l].seq += cluster_alns_filt[m][n][k]
                                l += 1
                        else:  # Append a gap everywhere else
                            for n in range(len(cluster_alns_filt[m])):
                                merged_aln[l].seq += '-'
                                l += 1

        # TODO: For testing - REMOVE
        if i%100 == 0:
            print(f'Finished column {i}')

    # Write and/or return merged aln
    if out_file:
        SeqIO.write(merged_aln, out_file, "fasta")
    if not no_return:
        return merged_aln


def indel_split(in_file, out_file=None, format='fasta'):
    """ Process indel-rich alignments to split likely-misaligned positions. This creates a more gappy alignment but
    potentially improves overall quality and simplifies the resolution of indels during ancestral reconstruction.
    TODO: First implementing smart_trim below as a simpler approach to test improvement """

    # Only fasta currently implemented
    if format == 'fasta':
        # TODO: Implement fasta format
        pass

    else:
        raise RuntimeError(f'indel_split is not currently implemented for {format} format.')

def smart_trim(
        in_file,
        out_file=None,
        format='fasta',
        base_threshold=0.05,
        p_threshold=0.05,
        max_threads=10
):
    """ Trim columns informed by column occupancy and pairwise sequence similarity of sequences with characters in
     low-occupancy positions. Trimming in this manner is likely to remove erroneously aligned columns without a large
     loss of phylogenetic signal, and reduce complexity for indel resolution in ASR. Note that if ASR is being performed
     using the resulting alignment, an assessment should be made as to whether trimmed columns are likely to be present
     in any targetted node. """

    full_aln = AlignIO.read(in_file, format)

    pw_coverage_full = {}

    # Multi-thread pairwise coverage calculations
    pw_coverage_executor = concurrent.futures.ThreadPoolExecutor(max_workers=max_threads)
    pw_coverage_lock = threading.Lock()

    # Method for pairwise coverage calculation of sequence j with all sequences from j+1 to len(aln)
    def pw_coverage(j):
        # Calculate pairwise alignment coverages (as aligned proportion of shorter sequence)
        for k in range(j + 1, len(full_aln)):
            seqs = {full_aln[j].name: full_aln[j], full_aln[k].name: full_aln[k]}
            # Determine shortest ungapped seq
            shorter = min(seqs.keys(), key=lambda seq: len(str(seqs[seq].seq).replace('-', '')))
            longer = [key for key in seqs.keys() if key != shorter][0]
            # Find proportion of shorter seq positions aligned w/ longer
            aligned_cnt = 0  # Positions with content in short seq with content also in long seq
            for pos in range(full_aln.get_alignment_length()):
                if seqs[shorter].seq[pos] != '-' and seqs[longer].seq[pos] != '-':
                    aligned_cnt += 1
            aligned_proportion = aligned_cnt / len(str(seqs[shorter].seq).replace('-', ''))
            with pw_coverage_lock:
                pw_coverage_full[tuple(seqs.keys())] = aligned_proportion
        print(f"PW coverages done for seq {j}")

    future_store = [pw_coverage_executor.submit(pw_coverage, j) for j in range(len(full_aln) - 1)]
    concurrent.futures.wait(future_store)

    # TODO: Now in nested function for multi-threading - remove this if it works
    # # Calculate pairwise alignment coverage (as aligned proportion of shorter sequence) for all seqs in full_aln
    # for j in range(len(full_aln) - 1):
    #     for k in range(j + 1, len(full_aln)):
    #         seqs = {full_aln[j].name: full_aln[j], full_aln[k].name: full_aln[k]}
    #         # Determine shortest ungapped seq
    #         shorter = min(seqs.keys(), key=lambda seq: len(str(seqs[seq].seq).replace('-', '')))
    #         longer = [key for key in seqs.keys() if key != shorter][0]
    #         # Find proportion of shorter seq positions aligned w/ longer
    #         aligned_cnt = 0  # Positions with content in short seq with content also in long seq
    #         for pos in range(full_aln.get_alignment_length()):
    #             if seqs[shorter].seq[pos] != '-' and seqs[longer].seq[pos] != '-':
    #                 aligned_cnt += 1
    #         aligned_proportion = aligned_cnt / len(str(seqs[shorter].seq).replace('-', ''))
    #         pw_coverage_full[tuple(seqs.keys())] = aligned_proportion

    to_trim = []  # Indices of columns to trim
    under_co_cnt = []

    # Multi-thread column assessment for trimming
    col_assess_executor = concurrent.futures.ThreadPoolExecutor(max_workers=max_threads)
    to_trim_lock = threading.Lock()  # For appending cols to trim
    under_co_lock = threading.Lock()  # For keeping track of # of cols assessed for trimming

    # Method for assessing columns for potential trimming
    def assess_col(i):

        # Check if column occupancy is below base_threshold and assessable for trimming
        if 1 - (full_aln[:, i].count('-') / len(full_aln)) < base_threshold:
            with under_co_lock:
                under_co_cnt.append(i)

            # Collect seqs with content in column i
            sub_aln = [full_aln[j] for j in range(len(full_aln)) if full_aln[j][i] != '-']

            # Sort coverage comparisons by whether both seqs have content in column of interest
            pw_coverage_pos = {}  # Both seqs have content in column of interest
            pw_coverage_neg = {}  # One or both do not
            sub_pairs = []  # List of immutable sets representing pairs of seqs in sub_aln
            for j in range(len(sub_aln) - 1):
                for k in range(j, len(sub_aln)):
                    # TODO: Check if below line is correct (was throwing error, have made a change but not sure its right)
                    sub_pairs.append(tuple([sub_aln[j].name, sub_aln[k].name]))
            for seqs, coverage in pw_coverage_full.items():
                if seqs in sub_pairs:
                    pw_coverage_pos[seqs] = coverage
                else:
                    pw_coverage_neg[seqs] = coverage

            # Maintain column only if positive pairwise coverages are significantly higher than negatives
            t, p = stats.ttest_ind(list(pw_coverage_pos.values()), list(pw_coverage_neg.values()), equal_var=False,
                                   alternative='greater')
            if p >= p_threshold:  # Trim if NOT significantly higher coverage in aligned seqs
                with to_trim_lock:
                    to_trim.append(i)

            print(f'Assessed col {i}')

        else:
            print(f'Col {i} above CO threshold')

    future_store = [col_assess_executor.submit(assess_col, i) for i in range(full_aln.get_alignment_length())]
    concurrent.futures.wait(future_store)

    # TODO: Now in nested function for multi-threading - remove if all works
    # for i in range(full_aln.get_alignment_length()):   # For each column

        # # Check if column occupancy is below base_threshold and assessable for trimming
        # if 1 - (full_aln[:, i].count('-') / len(full_aln)) < base_threshold:
        #     under_co_cnt += 1
        #
        #
        #     # Collect seqs with content in column i
        #     sub_aln = [full_aln[j] for j in range(len(full_aln)) if full_aln[j][i] != '-']
        #
        #     # Sort coverage comparisons by whether both seqs have content in column of interest
        #     pw_coverage_pos = {}  # Both seqs have content in column of interest
        #     pw_coverage_neg = {}  # One or both do not
        #     sub_pairs = []  # List of immutable sets representing pairs of seqs in sub_aln
        #     for j in range(len(sub_aln)-1):
        #         for k in range(j, len(sub_aln)):
        #             # TODO: Check if below line is correct (was throwing error, have made a change but not sure its right)
        #             sub_pairs.append(tuple([sub_aln[j].name, sub_aln[k].name]))
        #     for seqs, coverage in pw_coverage_full.items():
        #         if seqs in sub_pairs:
        #             pw_coverage_pos[seqs] = coverage
        #         else:
        #             pw_coverage_neg[seqs] = coverage
        #
        #     # Maintain column only if positive pairwise coverages are significantly higher than negatives
        #     t, p = stats.ttest_ind(list(pw_coverage_pos.values()), list(pw_coverage_neg.values()), equal_var=False,
        #                            alternative='greater')
        #     if p >= p_threshold:  # Trim if NOT significantly higher coverage in aligned seqs
        #         to_trim.append(i)
        #
        #         trimmed_cnt += 1
        #
        #     print(f'Assessed col {i}')
        #
        # else:
        #     print(f'Col {i} above CO threshold')

    # Trim appropriate columns
    trimmed = trim_aln(full_aln, to_trim)

    print(f'Trimmed {len(to_trim)} of {len(under_co_cnt)} under occupancy threshold.')

    # Trimmed cols may be out of order due to multi-threading
    to_trim.sort()

    print(to_trim)

    # Write trimmed aln to file if specified, else return MSA object
    if out_file:
        AlignIO.write(trimmed, out_file, format)
        return to_trim   # Return trimmed cols as list
    else:
        return trimmed


def phylo_smart_trim(
        in_file,
        out_file=None,
        log=None,
        format='fasta',
        occ_threshold=0.05,
        max_threads=10
):
    return


def map_cols(aln, ref_name, positions):
    """ Identify alignment columns corresponding to specific (ungapped) site indices of a reference sequence. """

    ref_seq = SeqIO.index(aln, 'fasta')[ref_name]  # TODO: This is saving whole file in memory - change to just parsing

    pos_map = {}   # Un-gapped to gapped position for ref_seq
    ungapped_idx = 0   # Ungapped index of next non-gap character

    for i in range(len(ref_seq)):
        if not ref_seq[i] == '-':
            pos_map[ungapped_idx] = i
            ungapped_idx += 1

    gapped_pos = [pos_map[pos] for pos in positions]

    return gapped_pos


def map_pos_to_aln(aln, ungapped_positions):
    """ Map ungapped sequence positions to alignment columns. Ungapped_positions is a dictionary mapping
    {seq_id : [list_of_pos_nums]. Ouptuts a dictionary in same format with mapped columns. """

    mappings = {}

    for seq in SeqIO.parse(aln, "fasta"):

        if seq.name in ungapped_positions:

            pos_map = {}  # Un-gapped to gapped position for ref_seq
            ungapped_idx = 0  # Ungapped index of next non-gap character

            for i in range(len(seq)):
                if not seq[i] == '-':
                    pos_map[ungapped_idx] = i
                    ungapped_idx += 1

            # Note positions arguments start from 1 --> need to -1 to get idx
            mappings[seq.name] = [pos_map[pos-1] for pos in ungapped_positions[seq.name]]

    return mappings


def extract_mapped_cols(
        aln,
        ref_name,
        positions,
        out_file=None,
        no_return=False
):
    """ Return a dictionary mapping sequence ID to characters in columns corresponding to a reference sequence and
    associated un-gapped positions. Optionally, write sub-sequences containing cols of interest to alignment file. """

    mapped_cols = map_cols(aln, ref_name, positions)

    aln_seqs = list(SeqIO.parse(aln, 'fasta'))
    sub_aln_seqs = []
    for seq in aln_seqs:
        subseq = ""
        for col in mapped_cols:
            subseq += seq[col]
        seq.seq = Seq.Seq(subseq)

    if out_file:
        SeqIO.write(aln_seqs, out_file, 'fasta')

    if not no_return:
        subseqs = []
        for seq in aln_seqs:
            subseqs.append((seq.name, str(seq.seq)))
        return subseqs


def tag_aln_by_identity(
        alignment_file,
        reference_file,
        output_file,
        identity_threshold=0.9,
        coverage_threshold=0.9,
        reference_name_map=None
):
    """
    Tag sequences in an alignment by their best match to a reference set using MAFFT --localpair.

    Parameters:
    - alignment_file: path to FASTA alignment file
    - reference_file: path to reference FASTA file
    - output_file: where to write tagged output FASTA
    - identity_threshold: minimum identity required for tagging (0-1)
    - coverage_threshold: minimum coverage of reference required (0-1)
    - reference_name_map: optional dict mapping reference indices (0-based) to custom tag strings
    """

    def run_mafft_pairwise(seq1, seq2):
        """
        Run MAFFT --localpair between two SeqRecords and return aligned pair.
        """
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as input_fasta, \
             tempfile.NamedTemporaryFile(mode="r", delete=False) as output_fasta:
            SeqIO.write([seq1, seq2], input_fasta, "fasta")
            input_fasta.flush()

            cmd = ["mafft", "--quiet", "--localpair", input_fasta.name]
            with open(output_fasta.name, "w") as out_f:
                subprocess.run(cmd, stdout=out_f, stderr=subprocess.DEVNULL, check=True)

            aligned = list(SeqIO.parse(output_fasta.name, "fasta"))

        os.remove(input_fasta.name)
        os.remove(output_fasta.name)

        return aligned[0], aligned[1]

    def calc_identity_and_coverage(aln1, aln2):
        """
        Calculate sequence identity and reference coverage from aligned sequences.
        - Identity = matches / aligned positions (excluding gaps in either)
        - Coverage = aligned positions / ungapped reference length
        """
        s1 = str(aln1.seq)
        s2 = str(aln2.seq)

        aligned_positions = 0
        matches = 0
        ref_len = len(s2.replace("-", ""))

        for a, b in zip(s1, s2):
            if a == "-" or b == "-":
                continue
            aligned_positions += 1
            if a == b:
                matches += 1

        identity = matches / aligned_positions if aligned_positions > 0 else 0
        coverage = aligned_positions / ref_len if ref_len > 0 else 0
        return identity, coverage

    # Read alignment and reference sequences
    alignment = list(SeqIO.parse(alignment_file, "fasta"))
    references = list(SeqIO.parse(reference_file, "fasta"))

    tagged_records = []

    for query in alignment:
        best_tag = None
        best_identity = 0

        for i, ref in enumerate(references):
            try:
                aligned_query, aligned_ref = run_mafft_pairwise(query, ref)
                identity, coverage = calc_identity_and_coverage(aligned_query, aligned_ref)

                if identity >= identity_threshold and coverage >= coverage_threshold:
                    if identity > best_identity:
                        # Use user-defined tag if available; fallback to reference ID
                        best_tag = reference_name_map[i] if reference_name_map and i in reference_name_map else ref.id
                        best_identity = identity

            except subprocess.CalledProcessError as e:
                print(f"MAFFT failed for pair {query.id} vs {ref.id}: {e}")

        if best_tag:
            query.id += f"_{best_tag}"

        query.description = ""  # Remove FASTA description for clean output
        tagged_records.append(query)

    # Write updated alignment to file
    SeqIO.write(tagged_records, output_file, "fasta")


# Pre-defined metal sites of some MBL superfamily proteins for mapping alignment columns
# Ordering of columns is: 3 x alpha site, 3 x beta site, bridging residue (note col for bridging residue is included
# regardless of whether Aspartate appears there. Mapped as Uniprot ID : list of positions
METAL_SITES = { # NOTE: These are indices (i.e. aln_position - 1)

    # TODO: Incorporate new information about evolution of metal site usage between sub-families

    # E coli ZipD - TODO: Need to check - does a different distal get used in the beta site between classes?
    # "P0A8V0" : [64, 66, 141, 68, 69, 270, 212],

    # AIM-1
    "B5DCA0" : [103, 105, 181, 107, 108, 249, 206],  # Have verified that these are now correct
    "UPI00017F685B" : [103, 105, 181, 107, 108, 249, 206],

    # Asp-bridging sequence from B3 exp iter4 (bug in priority for cluster rep selection)
    "UPI002FBC470E" : [70, 72, 155, 74, 75, 216, 176],

    # Yeast Trz1 - NOTE: most distal Histidine is not position-homologous to that of B3 MBLs
    # 3 x alpha, 3 x beta, bridging, phosphate co-ordinating Histidine (8 positions total)
    "P36159" : [539, 541, 669, 543, 544, 758, 698, 736]

}


