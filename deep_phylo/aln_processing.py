import os

from Bio import SeqIO

""" Functions for analysis and automated QC of MSAs. """


def seq_lens(seqs, descending=True):
    """ Return a list of doubles (seq_name : len). Input may be either a fasta file or list of Bio.SeqRecords. """

    if isinstance(seqs, str):  # Read in seqs from file
        seqs = list(SeqIO.parse(seqs, 'fasta'))

    lens = [(seq.name, len(seq)) for seq in seqs]
    reverse = True if descending else False
    lens.sort(reverse=reverse, key=lambda x:x[1])

    return lens


def len_filter(seqs, out_file=None, no_return=False, min=None, max=None):
    """ Filter sequences by length. Input can be fasta or a list of Bio.SeqRecords. Output can be to fasta, as a list,
    or both."""

    if no_return and not out_file:
        raise RuntimeError("If no return, an output file name must be specified.")

    if not (min or max):
        raise RuntimeError("At least one of minimum or maximum length must be specified.")

    if isinstance(seqs, str):  # Read in seqs from file
        seqs = list(SeqIO.parse(seqs, 'fasta'))

    lens = seq_lens(seqs)  # List of doubles: (name, len)

    if min:
        lens = [length for length in lens if length[1] >= min]
    if max:
        lens = [length for length in lens if length[1] <= max]

    seqs = [seq for seq in seqs if seq.name in [length[0] for length in lens]]

    if out_file:
        SeqIO.write(seqs, out_file, 'fasta')

    if not no_return:
        return seqs


def trim_retention(aln_file, gt, rm_temp=True):
    """ Return a dictionary of seq : seq% retained when aln is trimmed at given column occupancy threshold.
    Set rm_temp=False to retain trimmed aln files. """

    trimmed_file = f"{aln_file.split('.')[0]}_t{int(gt*100)}.aln"
    os.system(f"trimal -in {aln_file} -out {trimmed_file} -gt {gt}")

    ret_dict = {}  # Seq ID mapped to content retention rate
    untrimmed_seqs = SeqIO.index(aln_file, 'fasta')
    trimmed_seqs = SeqIO.index(trimmed_file, 'fasta')

    for name in trimmed_seqs.keys():
        # Find number of non-gap characters in trimmed aln for this seq
        trimmed_seq = trimmed_seqs[name].seq
        non_gaps = len(trimmed_seq) - trimmed_seq.count('-')

        # Find ungapped len of original sequence
        untrimmed_seq = untrimmed_seqs[name].seq
        ungapped_len = len(untrimmed_seq) - untrimmed_seq.count('-')

        try:
            # Return as percentage
            ret_dict[name] = non_gaps / ungapped_len * 100
        except ZeroDivisionError:  # Edge case - input seq contains no content - needs to be filtered anyway so set tr=0
            ret_dict[name] = 0

    if rm_temp:
        os.system(f"rm {trimmed_file}")

    return ret_dict


def trim_gap_pc(aln_file, gt, rm_temp=True):
    """ Return a dictionary of seq : gap proportion when aln is trimmed at given threshold
    Set rm_temp=False to retain trimmed aln files."""

    trimmed_file = f"{aln_file.split('.')[0]}_t{int(gt*100)}.aln"
    os.system(f"trimal -in {aln_file} -out {trimmed_file} -gt {gt}")

    gap_dict = {}
    trimmed_seqs = SeqIO.index(trimmed_file, 'fasta')

    for name in trimmed_seqs.keys():
        # Find percentage of gap characters in trimmed aln for this seq
        trimmed_seq = trimmed_seqs[name].seq
        gap_dict[name] = trimmed_seq.count('-') / len(trimmed_seq) * 100

    if rm_temp:
        os.system(f"rm {trimmed_file}")

    return gap_dict


def mean_col_agreement(aln, min_occupancy=0.5):
    """ Return a dictionary of sequence IDs mapped to the mean of frequencies of the sequence's character in each column
     for which it has one of the 20 amino acids, if said column's overall occupancy is above a given threshold.
     Ambiguous characters are ignored. """

    aln = list(SeqIO.parse(aln, "fasta"))
    col_freqs = {seq.name : [] for seq in aln}

    for i in range(len(aln[0])):
        aa_cnts = {aa : 0 for aa in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']}
        for seq in aln:
            if seq[i] in aa_cnts.keys():
                aa_cnts[seq[i]] += 1
        total = sum(aa_cnts.values())
        if total == 0 or total/len(aln) < min_occupancy:
            # Column only has ambiguous characters and/or is below occupancy threshold
            continue
        aa_freqs = {aa : cnt/total for aa, cnt in aa_cnts.items()}
        for seq in aln:
            if seq[i] in aa_cnts.keys():
                col_freqs[seq.name].append(aa_freqs[seq[i]])

    mean_agreement = {seq : sum(freqs) / len(freqs) for seq, freqs in col_freqs.items()}
    return mean_agreement
