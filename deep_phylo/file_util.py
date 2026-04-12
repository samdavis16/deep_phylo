import gzip

from Bio import SeqIO, Seq


def records_from_fasta(
        in_file,
        targets=None,
        exclude=None,
        flag_missing='warning'
):
    """ Parse sequences from a fasta file as SeqRecords. If targets are not specified, load all sequences from the file
    except a specified list. target sequence IDs missing from the input file can be handled by 'warning' (default),
     'error', or 'ignore'. """

    if targets and exclude:
        raise RuntimeError("At most one one of target or exclude sequences may be specified.")

    all_records = SeqIO.index(in_file, 'fasta')

    if targets:
        if isinstance(targets, str):
            targets = [targets]
        retain_records = []
        for target in targets:
            if target in all_records.keys():
                retain_records.append(all_records[target])
            else:
                if flag_missing == "error":
                    raise RuntimeError(f"Target sequence: {target} not found in {in_file}.")
                elif flag_missing == "ignore":
                    pass
                else:
                    raise RuntimeWarning(f"Target sequence: {target} not found in {in_file}.")
        return retain_records

    elif exclude:
        return [record for record in all_records if record.name not in exclude]

    else:
        return all_records


def raw_extract_fasta(
        seq_db,
        out_fa,
        target_seqs,
        boundaries=None,
        id_format='uniref50'
):
    """ Extract a list of target_seqs from a Fasta database and write to a new Fasta file. This parses the db file
    line-by-line and is useful when the db is too large to be read in as a collection of sequence objects (e.g. via
     BioPython). The full sequence can be extracted, or else boundaries (e.g. of specific domains) can be provided if
     desired

     TODO: could allow for several fasta outputs based on different sequence sets to save reading file multiple times"""

    if boundaries:
        if len(target_seqs) != len(boundaries):
            raise RuntimeError("If boundaries are provided, exactly one set must be provided for each target sequence.")
        else:
            boundary_map = {seq:bounds for seq, bounds in zip(target_seqs, boundaries)}

    with open(seq_db, 'r') as in_file:
        with open(out_fa, 'w') as out_file:
            in_seq = False
            seq = ''
            while True:
                line = in_file.readline()
                if len(line) == 0:  # EOF
                    if in_seq:
                        # Extract sub-sequence if boundaries provided
                        subseq = seq[boundary_map[target_id][0]-1:boundary_map[target_id][1]] if boundaries else seq
                        out_file.write(subseq+'\n')
                    break
                line = line.strip()
                if line.startswith('>'):
                    line = line.split('>')[1]
                    # If in_seq, write completed sequence (or sub-sequence) to file
                    if in_seq:
                        subseq = seq[boundary_map[target_id][0] - 1:boundary_map[target_id][1]] if boundaries else seq
                        out_file.write(subseq + '\n')
                    # New ID based on format indicated TODO: only Uniref50 currently implemented
                    if id_format == 'uniref50':
                        seq_id = line.split()[0].split('UniRef50_')[1]
                    else:
                        seq_id = line.split()[0]
                    if seq_id in target_seqs:
                        out_file.write('>' + seq_id + '\n')
                        seq = ''
                        in_seq = True  # Extracting this seq
                        target_id = seq_id
                    else:
                        in_seq = False
                else:
                    if in_seq:
                        seq += line


def extract_subset_fasta(
        full_fa,
        out_fa,
        target_seqs,
        boundaries=None
):
    """ Extract subset from Fasta file. Uses BioPython.SeqIO.parse with optional decompression on the fly. """

    if boundaries:
        if len(target_seqs) != len(boundaries):
            raise RuntimeError("If boundaries are provided, exactly one set must be provided for each target sequence.")
        else:
            boundary_map = {seq:bounds for seq, bounds in zip(target_seqs, boundaries)}

    with open(out_fa, 'w') as out_file:
        # Check if file is gzipped
        if ".gz" in full_fa:
            handle = gzip.open(full_fa, "rt")
        else:
            handle = open(full_fa)

        for seq in SeqIO.parse(handle, "fasta"):
            if seq.name in target_seqs:
                subseq = seq[boundary_map[seq.name][0] - 1:boundary_map[seq.name][1]] if boundaries else seq
                out_file.write(f">{seq.name}\n{str(seq.seq)}\n")

        handle.close()


def extract_fasta(
        seq_fa,
        out_fa,
        target_seqs,
        coords=None,
        ungap=False,
        idx_file=None
):
    """ Using BioPython SeqIO to parse sequence DB to a dictionary of records then extracting seqs in target_seqs
    """  """
    TODO: currently is directly writing sequences to output fasta as Uniref50 annotations are breaking stuff
     TODO: could allow for several fasta outputs based on different sequence sets to save reading file multiple times"""

    if coords and len(coords) != len(target_seqs):
        raise RuntimeError("If specified, coordinates must be provided for each target sequence.")

    if idx_file:
        # Create or load fasta idx
        fa_idx = SeqIO.index_db(
            idx_file,
            seq_fa,
            "fasta"
        )

        with open(out_fa, 'w') as out_f:
            # Extract (optionally with coordinates) via index
            SeqIO.write(
                (fa_idx[target_seqs[i]][coords[i][0]-1 : coords[i][1]]
                 if coords else fa_idx[target_seqs[i]]
                 for i in range(len(target_seqs))
                 if target_seqs[i] in fa_idx),
                out_f,
                "fasta"
            )

        fa_idx.close()

    else:

        seq_subset = []
        for seq in SeqIO.parse(seq_fa, "fasta"):
            if seq.name in target_seqs:
                if ungap:
                    seq.seq = Seq.Seq(str(seq.seq).replace("-",""))
                if coords:  # Extract only a subsequence according to coordinates if specified
                    this_coords = coords[target_seqs.index(seq.name)]
                    seq_subset.append(seq[this_coords[0]-1 : this_coords[1]])
                else:
                    seq_subset.append(seq)

        SeqIO.write(seq_subset, out_fa, "fasta")



def ungap_fasta(in_file, out_file):
    """ Output a Fasta file with ungapped versions of input sequences. """
    seq_names = list(SeqIO.index(in_file, 'fasta').keys())  # List of all seq names
    extract_fasta(in_file, out_file, seq_names, ungap=True)


def merge_fastas(in_files, out_file):
    """ Merge multiple fasta files, removing any duplicate sequences by name. """

    records = []  # SeqRecords
    seq_names = []  # Check for duplicate seq IDs between input files

    for file in in_files:
        for seq_record in SeqIO.parse(file, 'fasta'):
            if seq_record.name not in seq_names:
                seq_names.append(seq_record.name)
                records.append(seq_record)

    SeqIO.write(records, out_file, 'fasta')


def merge_seqs(
        in_fastas=None,
        in_records=None,
        exclude=None,
        out_file=None,
        no_return=False
):
    """ Merge sequences from any number of fasta and/or Bio.SeqRecord inputs, removing duplicates by name between (not
    within) input sets/files. SeqRecord inputs can be a list of multiple separate SeqRecord lists, or a single list.
    Sequences for exclusion can be specified by name. Output can either be a new fasta file, or a list of SeqRecords,
    or both."""

    if no_return and not out_file:
        raise RuntimeError("If no return is specified, an output file name must be specified.")

    if not in_fastas:
        in_fastas = []
    if not in_records:
        in_records = []
    if not exclude:
        exclude = []

    merged_records = []
    exclude_names = set(exclude)  # Check that these seqs are not added / re-added. Update as new sets added

    # Merge SeqRecord inputs
    if in_records:

        # Check for multiple collections of SeqRecord lists
        if isinstance(in_records[0], list):  # Multiple collections
            for collection in in_records:
                for record in collection:
                    if record.name not in exclude_names:
                        merged_records.append(record)
                exclude_names.update([record.name for record in collection])

        else:  # Single collection
            merged_records = [record for record in in_records if record.name not in exclude_names]
            exclude_names.update([seq.name for seq in in_records])

    # Merge fasta inputs
    if not isinstance(in_fastas, list):
        in_fastas = [in_fastas]

    for fasta in in_fastas:
        seqs = list(SeqIO.parse(fasta, 'fasta'))
        for record in seqs:
            if record.name not in exclude_names:
                merged_records.append(record)
        exclude_names.update([seq.name for seq in seqs])

    if out_file:
        SeqIO.write(merged_records, out_file, 'fasta')

    if not no_return:
        return merged_records


def tab_del_file_to_dict(file):
    """ Parse a generic key:val tab-delimited file to a dictionary. Lines are assumed to contain a key:value pair only
     and all tabs after the first are ignored. """

    with open(file) as in_file:
        out_dict = {}
        for line in in_file:

            data = line.strip().split('\t', maxsplit=1)
            if len(data) > 1:  # If val is not null
                out_dict[data[0]] = data[1]
            else:
                out_dict[data[0]] = None

    return out_dict


########## TREE FILE PARSERS ##########

def specific_node_labels(
        in_nwk,
        out_nwk,
        retain_nodes,
        relabel=None,
        ignore_missing=False
):
    """ Write new Newick file with only specific internal nodes labeled. Internal nodes must already be annotated in the
     original file and correspond to labelled provided in retain_nodes. These can be relabeled in the output file
     by providing a list of new labels in an order corresponding to the original label list.
     Note: minimal checks currently performed for integrity of input tree file. """

    # Check for duplicate labels
    node_set = set(retain_nodes)
    if len(retain_nodes) != len(node_set):
        raise RuntimeError("Duplicate internal nodes provided.")

    # Keep track of found node labels to check for any missing ones
    labels_found = []

    # Read in Newick string assuming single tree in file
    in_file = open(in_nwk, 'r')
    in_str = in_file.read()
    in_file.close()

    with open(out_nwk, 'w') as out_file:
        i = 0
        found_label = False
        while i < len(in_str):
            if in_str[i] == ')':  # Found internal node label
                out_file.write(')')
                found_label = True
                label = ''
                i += 1
                while found_label:
                    if in_str[i] == ':' or in_str[i] == ';':
                        found_label = False
                    else:
                        label += in_str[i]
                        i += 1
                i += 1
                try:  # Check if this is a label to retain
                    idx = retain_nodes.index(label)
                    labels_found.append(label)
                    if relabel:
                        out_file.write(f"{relabel[idx]}{':' if in_str[i-1] == ':' else ';'}")
                    else:
                        out_file.write(f"{label}{':' if in_str[i-1] == ':' else ';'}")
                except ValueError:  # Don't write this label to output Newick
                    out_file.write(':' if in_str[i-1] == ':' else ';')
            else:      # Not in node label
                out_file.write(in_str[i])
                i += 1

    # Check for provided node labels not encountered
    missed = []
    for node in retain_nodes:
        if node not in labels_found:
            missed.append(node)
    if missed:
        missed_str = ', '.join(missed)
        if ignore_missing:
            print(f"Nodes: [{missed_str}] not found in tree file and were ignored.")
        else:
            raise RuntimeError(f"Nodes: [{missed_str}] not found in tree file.")

