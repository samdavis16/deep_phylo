import concurrent.futures
import threading

from scipy import stats
from Bio import Align, AlignIO


def trim_aln(aln : Align.MultipleSeqAlignment, cols):
    """ Trim indices in cols from a multiple sequence alignment. """

    # Check that all col indices are in range of sequence length
    if not isinstance(cols, list): cols = [cols]
    cols = list(set(cols))  # Remove duplicate column indices
    for col in cols:
        if not 0 <= col < aln.get_alignment_length():
            raise RuntimeError(f"Invalid column index: {col}")

    # Trim from highest to lowest index
    cols.sort()
    for col in cols[::-1]:
        if col == aln.get_alignment_length() - 1:  # Trimming last col of current aln
            aln = aln[:, :-1]
        elif col == 0:    # Trimming first column
            aln = aln[:, 1:]
        else:   # Trimming non-end column and need to combine two sub-alns
            aln = aln[:, :col] + aln[:, col+1:]

    return aln


def smart_trim(in_file, out_file=None, format='fasta', base_threshold=0.05, p_threshold=0.05, max_threads=10):
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
