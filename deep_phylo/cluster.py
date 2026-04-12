"""
Functions for sequence clustering, mostly utilising and extending mmseqs2. File naming conventions used by
mmseqs2 are used here.

This module implements a hierarchical clustering framework used predominantly in the clade representation/expansion
workflow. It's functionality has extended far beyond its initial intended purpose and needs significant re-structuring.

Internal database identifiers are referred to as db_id(s), while the external (original)
accessions are referred to as acc(s).
"""

import glob
import math
import multiprocessing
import os
import random
import shutil
import subprocess

import numpy as np
from Bio import SeqIO, AlignIO
from sklearn import metrics

from . import aln
from . import file_util


def create_seq_db(fasta_names, db_name, quiet=True):
    """ Create an mmseqs database from one or several (optional gzipped) fasta files. """

    if not isinstance(fasta_names, list):
        fasta_names = [fasta_names]

    # Database name should conform with convention for this module
    db_name = db_name.split('DB')[0] + 'DB'

    # Remove any existing DB files associated with the seq_db or clust_db names
    for file in glob.glob(f"{db_name}*"):
        os.remove(file)

    args = ["mmseqs", "createdb"]

    for seq_file in fasta_names:
        args.append(seq_file)

    args.append(db_name)

    subprocess.run(args,
                   check=True,
                   stdout=subprocess.DEVNULL if quiet else None)


def get_internal_ids(db_name):
    """ Extract the internal IDs from a sub-database or clustering database (via the corresponding index file). Note
    that this will only extract the representative sequences for each cluster for a cluster DB. """

    db_name = db_name.split('DB')[0] + 'DB'

    int_ids = []
    with open(db_name + '.index') as index_file:
        for line in index_file:
            int_ids.append(line.strip().split('\t')[0])

    return int_ids


def get_external_ids(db_name):
    """ Extract the external IDs from a given database. Note: for sub-databases for which .lookup files are sym-linked,
    this function will extract all external IDs from the full database. This function is not suitable for cluster
     databases. """

    db_name_split = db_name.split('DB')[0]
    ext_ids = []

    with open(db_name_split + ".lookup") as in_file:
        for line in in_file:
            ext_ids.append(line.strip().split('\t')[1])

    return ext_ids


def get_cluster_members(db_name):
    """ Return a list of all members of all clusters from a cluster DB. """

    db_name = db_name.split('DB')[0] + 'DB'

    members = set()
    for clust in cluster_map(db_name).values():
        members.update(clust)

    return list(members)


def db_id_map(db_name, int2ext=True):
    """ For a given mmseqs sequence database, return a dictionary mapping identifier types (internal to external by
    default - make int2ext False for the reverse. """

    db_name = db_name.split('DB')[0] + 'DB'

    id_map = {}

    with open(db_name + '.lookup') as file:
        while True:
            line = file.readline()
            if not line.strip():
                break
            split_line = line.split()
            if int2ext:
                id_map[split_line[0]] = split_line[1]
            else:
                id_map[split_line[1]] = split_line[0]

    return id_map


def db_idx_dict(db_name):
    """ Read in index file for sequence database. Internal IDs are stored as strings (not ints) in all workflows.
    Index positions are stored as ints. """

    db_name = db_name.split('DB')[0] + 'DB'

    db_idx = {}
    with open(db_name + '.index') as idx_file:
        for line in idx_file:
            data = line.split()
            db_idx[data[0]] = int(data[1])  # str : int

    return db_idx


def db_seq_lens(seq_db, ext_ids=False):
    """ Return a dictionary mapping {ID : sequence_len} for a sequence database. Note: output for cluster databases will
    be meaningless. """

    filename = seq_db.split('DB')[0] + 'DB.index'

    if ext_ids:   # Need to map to external IDs
        id_map = db_id_map(seq_db)

    len_dict = {}
    with open(filename) as idx_file:
        for line in idx_file:
            data = line.strip().split('\t')
            # Internal ID is in 1st column - map to external if required
            len_dict[id_map[data[0]] if ext_ids else data[0]] = int(data[2]) - 2   # Index length includes \0 and \n

    return len_dict


def cluster_map(cluster_db):
    """ Return a dictionary mapping (internal) cluster representative identifiers to all members (including the
     representative). """

    cluster_db = cluster_db.split('DB')[0] + 'DB'

    map = {}

    with open(cluster_db) as db_file:
        rep = db_file.readline().strip()  # First cluster rep
        clust = [rep]   # Begin first cluster
        for line in db_file:
            if line.startswith("\0"):    # End of cluster
                map[rep] = clust    # Record previous cluster
                if line != "\0":  # NOT the last cluster
                    # Initiate new cluster
                    rep = line.strip().strip('\0')
                    clust = [rep]
            else:
                clust.append(line.strip())

    return map


def rep_map_reconstruction(
        og_cluster_db,
        master_seq_db=None,
        rep_fasta=None,
        sub_seq_db=None,
        sub_cluster_db=None
):
    """ Generates a map of {original:sampled} cluster representatives for a (possibly re-sampled) cascade clustering or
    re-sampled extracted representatives. The cascaded cluster DB
    or a sub-DB containing the resampled representatives can be provided as sub_db. If the provided databases do not
    represent a perfect cascade clustering (i.e. the smaller clustering/sub-database does not contain exactly one member
    from each cluster in the original clustering, None is returned.

    TODO: Had to make changes debugging this so now only works when sub_db is a cascade clustering of og_cluster_db
    TODO: Also make work for sequence sub-DB """

    og_cluster_db = og_cluster_db.split('DB')[0] + 'DB'
    if rep_fasta:
        if not master_seq_db:
            raise RuntimeError("The master sequence database must be provided if external representative IDs are "
                               "provided.")
    elif sub_seq_db:
        sub_seq_db = sub_seq_db.split('DB')[0] + 'DB'
    elif not sub_cluster_db:
        raise RuntimeError("Sub-cluster must be indicated either in the form of a fasta file or database of cluster "
                           "representatives, or a sub-clustering database.")
    else:
        sub_cluster_db = sub_cluster_db.split('DB')[0] + 'DB'

    # Map og cluster members to their rep
    og_cluster_map = cluster_map(og_cluster_db)
    og_member2rep = {}
    for rep, members in og_cluster_map.items():
        for member in members:
            og_member2rep[member] = rep

    # Get internal IDs of sub-clustering (possibly re-sampled) reps

    if rep_fasta:
        id_map = db_id_map(master_seq_db, int2ext=False)
        sub_reps = [id_map[seq.name] for seq in SeqIO.parse(rep_fasta, "fasta")]

    elif sub_seq_db:
        sub_reps = get_internal_ids(sub_seq_db)

    else:
        sub_reps = get_cluster_members(sub_cluster_db)  # cascade cluster members

    try:
        rep_map = {og_member2rep[rep] : rep for rep in sub_reps}
    except KeyError:   # At least one ID in sub-DB not present in original clustering
        print(f"At least one ID in subDB not present in original clustering {og_cluster_db}.")
        return None

    if set(rep_map.keys()) != set(og_cluster_map.keys()):
        print(f"The provided sub-cluster representatives is not comprised of exactly one representative from each "
              f"cluster in {og_cluster_db}.")
        return None

    return rep_map


def write_rep_map(filename, rep_map):
    """ Writes a representative map of {original : re-sampled} representatives to a tab-delimited file. """

    filename = filename.split('.repmap')[0] + '.repmap'

    with open(filename, 'w') as file:
        for og, new in rep_map.items():
            file.write(f"{og}\t{new}\n")


def read_rep_map(filename):
    """ Reads a representative map of {original : re-sampled} representatives from a tab-delimited file. """

    filename = filename.split('.repmap')[0] + '.repmap'
    rep_map = {}

    with open(filename) as file:
        for line in file:
            ids = line.strip().split('\t')
            rep_map[ids[0]] = ids[1]

    return rep_map


def create_sub_db(
        db_name,
        sub_db_name,
        sub_db_ids=None,
        sub_accs=None,
        report_missing=True,
        quiet=True
):
    """ Use mmseqs createsubdb to create a sub-database on either internal or external identifiers. """

    db_name = db_name.split('DB')[0] + 'DB'
    sub_db_name = sub_db_name.split('DB')[0] + 'DB'

    if sub_db_ids:
        sub_db_ids = set(sub_db_ids)  # For lookup efficiency
        if sub_accs:
            raise RuntimeError("Sub-database must be defined on either internal or external identifiers, not both.")
    elif not sub_accs:
        raise RuntimeError("A subset of (internal OR external) identifiers must be provided to define a sub-database.")
    else:
        sub_accs = set(sub_accs)

    db_lookup = db_name + ".lookup"

    # Create temporary file listing db_ids to retain
    with open(db_lookup) as lookup:
        with open('tmp_list', 'w') as tmp_list:

            found = []   # Keep track of IDs found in original DB lookup file
            lookup_line = lookup.readline()
            while lookup_line:

                if sub_db_ids:
                    this_id = lookup_line.strip().split('\t')[0]
                    if this_id in sub_db_ids:
                        found.append(this_id)
                        tmp_list.write(f"{this_id}\n")
                else:
                    this_id = lookup_line.strip().split('\t')[1]
                    if this_id in sub_accs:
                        found.append(this_id)
                        internal_id = lookup_line.strip().split('\t')[0]
                        tmp_list.write(f"{internal_id}\n")

                lookup_line = lookup.readline()

    # Create subDB
    args = ["mmseqs", "createsubdb", "tmp_list", db_name, sub_db_name]
    subprocess.run(args,
                   check=True,
                   stdout=subprocess.DEVNULL if quiet else None)

    # Report IDs not found in original DB if required
    if report_missing and len(found) != len(sub_db_ids if sub_db_ids else sub_accs):
        missing = [acc for acc in (sub_db_ids if sub_db_ids else sub_accs) if acc not in found]
        if missing:
            print(f"IDs not found in original database lookup: {missing}")

    os.remove('tmp_list')


def create_intersection_db(
        master_db,
        target_db,
        int_db,
        id_type='external'
):
    """ Creates an mmseqs database containing entries in the intersection of two existing databases. A sub-db of the
    first provided db (master_db) is created. By default, external identifiers are used to determine matching sequences.
    Use acc_type='internal' to strictly check for matching internal IDs (only suitable if target_db is a sub-database of
    master_db). If one or both databases are very large, efficiency will be increased when target_db is the smaller. """

    master_db = master_db.split('DB')[0] + 'DB'
    target_db = target_db.split('DB')[0] + 'DB'
    int_db = int_db.split('DB')[0] + 'DB'

    with open(target_db+'.index') as target_index:
        target_ids_int = set([line.strip().split('\t')[0] for line in target_index])  # Internal target ids

    with open(target_db+'.lookup') as target_lookup:  # Note: this could still take some time if target_db is a sub-db of a larger db
        target_id_map = {line.strip().split('\t')[0] : line.strip().split('\t')[1] for line in target_lookup
                         if line.strip().split('\t')[0] in target_ids_int}


    # Search for matching IDs in master DB
    sub_ids = []   # Internal or external IDs in master DB which intersect with target
    with open(master_db+'.lookup') as master_lookup:

        if id_type == 'external':
            target_ids = set(target_id_map.values())  # Check for matching of external IDs
            for line in master_lookup:
                if line.strip().split('\t')[1] in target_ids:
                    sub_ids.append(line.strip().split('\t')[0])   # Need the internal ID of master for subDB


        elif id_type == 'internal':
            target_ids = set(target_id_map.keys())  # Check for matching of internal IDs
            for line in master_lookup:
                if line.strip().split('\t')[0] in target_ids:
                    sub_ids.append(line.strip().split('\t')[0])

        else:
            raise RuntimeError("ID type for DB intersection must be either 'internal' or 'external'.")

    # Create intersection DB as sub-DB from master
    create_sub_db(master_db, int_db, sub_db_ids=sub_ids)


def create_profile_db(
        profile_name,
        msa_file=None,
        msa_db=None,
        msa_db_name=None,
        match_ratio=0.5
):
    """ From a FASTA-formatted MSA, create a mmseqs profile database. match_ratio defines the column occupancy
     required for inclusion in profile (using match-mode=1). Use match_ratio=None to use first sequence as centre
     sequence for assigning columns (i.e. use mmseqs default of match-mode=0). """

    if msa_file:

        if msa_db:
            raise RuntimeError("Either an MSA file or database (not both) should be provided for profile construction.")

        else:
            # Need to convert FASTA aln to Stockholm aln
            msa_seqs = SeqIO.parse(msa_file, "fasta")
            stockholm_name = msa_file.split('.')[0] + '.stockholm'  # temp stockholm-formatted alignment
            SeqIO.write(msa_seqs, stockholm_name, 'stockholm')

            # Create alignment database
            if not msa_db_name:
                msa_db_name = profile_name.split('DB')[0] + '_msaDB'
            os.system(f"mmseqs convertmsa {stockholm_name} {msa_db_name}")

            # Remove temp stockholm aln but retain new MSA DB
            os.remove(stockholm_name)

    else:
        if msa_db:
            msa_db_name = msa_db
        else:
            raise RuntimeError("Either an MSA file or database should be provided for profile construction.")

    # Build PSSM using match-mode 1 (column inclusion is on basis of occupancy %)
    profile_db_name = profile_name.split('DB')[0] + 'DB'  # Keeps naming conventions consistent
    match_ratio_option = f" --match-mode 1 --match_ratio {match_ratio}" if match_ratio else ""
    os.system(f"mmseqs msa2profile {msa_db_name} {profile_db_name} --match-{match_ratio_option}")


def getDBLen(seq_db):
    """ Read the database index file to determine the number of sequences in a DB. Assumes index file is
    {DB_name}.index.
    TODO: Could reduce memory burden here by not reading whole index file at once. """
    db_name = seq_db.split('DB')[0] + 'DB'
    with open(db_name+'.index') as index:
        return len(index.readlines())


def removeDB(db_name):
    """ Remove an mmseqs database and associated files. To be used (with caution!) for cleaning up unnecessary files.
     Assumes that naming conventions imposed by mmseqs are followed. """
    db_name = db_name.split('DB')[0]+'DB'
    for file in glob.glob(f"{db_name}*"):
        os.remove(file)


def profile_score(
        profile_db,
        target_file=None,
        target_db=None,
        target_db_name=None,
        prefilter_name=None,
        align_name=None,
        max_seqs=None,
        prefilter_threshold=0,
        sens=9,
        e_threshold=None,
        return_scores=False,
        delete_prefilter=True
):
    """ Runs a profile search against a target sequence database. The prefilter and alignment steps are run separately
     allowing for maximum parameter flexibility. NOTE: Currently, prefilter is run with highly sensitive parameters
     and no e-value threshold is applied to align by default. If a sequence file in FASTA format is provided, a target
     database will be created.
     TODO: turns out you can run the search workflow and still have full flexibility with individual modules's params
     TODO: need to integrate option(s) for alignment coverage """

    profile_db = profile_db.split('DB')[0] + 'DB'

    if target_file:

        if target_db:
            raise RuntimeError("Either a sequence file in FASTA format or a sequence database (not both) muse be supplied.")

        else:
            # Create target seq DB
            if not target_db_name:
                target_db_name = target_file.split('.')[0] + 'DB'
            else:
                target_db_name = target_db_name.split('DB')[0] + 'DB'
            os.system(f"mmseqs createdb {target_file} {target_db_name}")

    else:
         if not target_db:
             raise RuntimeError("Either a sequence file in FASTA format or a sequence database muse be supplied.")
         else:
             target_db_name=target_db

    # Run mmseqs profile -> seqs search

    # Prefilter params
    max_seqs_option = f"--max-seqs {max_seqs if max_seqs else getDBLen(target_db_name)}"  # Allow all to pass by default
    pref_threshold_option = f"--min-ungapped-score {prefilter_threshold}"
    sens_option = f"-s {sens}"
    prefilter_name = prefilter_name.split('DB')[0]+'DB' if prefilter_name else \
        f"{profile_db.split('DB')[0]}_{target_db_name.split('DB')[0]}_prefDB"

    # Run prefilter
    os.system(f"mmseqs prefilter {profile_db} {target_db_name} {prefilter_name} {max_seqs_option} "
              f"{pref_threshold_option} {sens_option}")

    # Align params
    e_threshold_option = f"-e {e_threshold if e_threshold else 'inf'}"
    align_name = align_name.split('DB')[0]+'DB' if align_name else \
        f"{profile_db.split('DB')[0]}_{target_db_name.split('DB')[0]}_alignDB"

    # Run align
    os.system(f"mmseqs align {profile_db} {target_db_name} {prefilter_name} {align_name} {e_threshold_option}")

    # Create TSV file of search results
    tsv_name = align_name.split('DB')[0]+'.tsv'
    os.system(f"mmseqs createtsv {profile_db} {target_db_name} {align_name} {tsv_name}")

    if delete_prefilter:
        os.system(f"rm {prefilter_name}*")

    if return_scores:   # Read tsv to get dictionary of {target_seq : aln_score}
        with open(tsv_name) as tsv:
            return {line.strip().split('\t')[1] : line.strip().split('\t')[2] for line in tsv}


def cluster_custom(
        seq_db,
        clust_db,
        min_seq_id,
        linclust=False,
        max_hits_per_query=None,
        tmp_dir=None,
        min_coverage=0.8,
        coverage_type=1,
        clust_mode=2,
        quiet=True,
        max_threads=20,
        nice=None
):
    """ Use mmseqs cluster or linclust on a pre-prepared mmseqs sequence database. Note: only some parameters are controllable
     through this function. See docs for justifications of fixed params - this workflow is customised for
     analysing relationships in protein superfamilies by clustering at different levels/hierarchies. This function uses
     the 'cluster' workflow by default. Specify linclust=True for linear time (less rigorous) clustering. """

    # Ensure database names follow convention
    seq_db = seq_db.split('DB')[0] + 'DB'
    clust_db = clust_db.split('DB')[0] + 'DB'

    for file in glob.glob(f"{clust_db}*"):
        os.remove(file)

    # Make tmp directory if required
    if not tmp_dir and not os.path.exists('tmp'):
        os.mkdir('tmp')

    if linclust:

        args = []

        # Optional nice
        if nice:
            args.extend(["nice", "-n", str(nice)])

        args.extend([
            "mmseqs", "linclust",
            seq_db,
            clust_db,
            tmp_dir if tmp_dir else "tmp",
            "--min-seq-id", str(min_seq_id),
            "-c", str(min_coverage),
            "--cov-mode", str(coverage_type),
            "--cluster-mode", str(clust_mode),
            "--threads", str(max_threads),
        ])

        if quiet:
            args.extend(["-v", '0'])

        # TODO: For debugging - REMOVE try/except around the subprocess.run
        try:
            subprocess.run(args, check=True)
        except subprocess.CalledProcessError as e:
            print(e.returncode)
            print(e.output)
            print(e.stdout)
            print(e.stderr)
            raise

    else:

        args = []

        # Optional nice
        if nice:
            args.extend(["nice", "-n", str(nice)])

        args.extend([
            "mmseqs", "cluster",
            seq_db,
            clust_db,
            tmp_dir if tmp_dir else "tmp",
            "--max-seqs", str(max_hits_per_query if max_hits_per_query else getDBLen(seq_db)),
            "--min-seq-id", str(min_seq_id),
            "-c", str(min_coverage),
            "--cov-mode", str(coverage_type),
            "--cluster-mode", str(clust_mode),
            "--cluster-reassign",
            "--threads", str(max_threads),
        ])

        if quiet:
            args.extend(["-v", '0'])

        subprocess.run(args, check=True)

    # TODO: ##################################### OLD CODE, REMOVE ONCE VERIFIED ABOVE CHANGES WORK ##################

    ### Prefilter params ###

    # # Allow all comparisions for which to pass prefilter by default - no constraint on cluster sizes
    # max_seqs_option = f"--max-seqs {max_hits_per_query if max_hits_per_query else getDBLen(seq_db)}"
    #
    # # TODO: perhaps these files are all put in tmp by default? Check and remove naming options if so (same with align)
    # prefilter_name = prefilter_name.split('DB')[0] + 'DB' if prefilter_name else \
    #     f"{seq_db.split('DB')[0]}_clust_{min_seq_id}_prefDB"
    #
    # ### Align params ###
    #
    # min_seq_id_option = f"--min-seq-id {min_seq_id}"
    # min_coverage_option = f"-c {min_coverage}"  # Note: bi-directional (mmseqs default) is fixed coverage type
    #
    # # Default for --cov-mode and --cluster-mode is Target coverage, greedy (by length) respectively
    # # Allows fragments to be clustered appropriately without risk of fragments becoming cluster reps
    # coverage_type_option = f"--cov-mode {coverage_type}"  # "Target coverage" is default
    # clust_mode_option = f"--cluster-mode {clust_mode}"    # In conjunction with target coverage, suitable for fragments
    #
    # # TODO: See above
    # align_name = align_name.split('DB')[0] + 'DB' if align_name else \
    #     f"{seq_db.split('DB')[0]}_clust_{min_seq_id}_alignDB"
    #
    # # Cluster reassignment is enabled (flag in command below)
    #
    # # Max threads to use
    # threads_option = f"--threads {max_threads}"
    #
    # # Optionally use nice if running on shared compute resources
    # nice_option = f"nice -n {nice} " if nice else ""
    #
    # if linclust:  # Linear time clustering
    #
    #
    #
    # else:  # More rigorous clustering
    #     # Run mmseqs cluster with above parameters + some fixed options
    #
    #     os.system(f"{nice_option}mmseqs cluster {seq_db} {clust_db} {tmp_dir if tmp_dir else 'tmp'} {max_seqs_option} "
    #               f"{min_seq_id_option} {min_coverage_option} {coverage_type_option} {clust_mode_option} "
    #               f"--cluster-reassign {threads_option}")

    #TODO: ################################# End of old code #########################################################

    # TODO: Troubleshooting cluster DB merging for MMseqs - This works on test dataset, need to verify it generalises
    # Check if merged cluster file was successfully produced
    if clust_db not in os.listdir():  # If not, manually merge segments
        seg_db_files = [file for file in os.listdir() if clust_db+"." in file and file.split(".")[-1][0] in "0123456789"]
        seg_db_files.sort(key=lambda x:int(x.split('.')[-1]))
        with open(clust_db, 'w') as full_clust_file:
            for seg_file in seg_db_files:
                with open(seg_file) as in_file:
                    for line in in_file:
                        full_clust_file.write(line)
                os.remove(seg_file)

    # TODO: See above re: temp files
    # if remove_temp_db:
    #     pass # TODO: Add removal of temporary databases (prefilter + alignment)


def rr_seqs(
        in_file,
        rr_file,
        min_seq_id,
        min_coverage,
        linclust=False,
        max_threads=8,
        nice=None,
        quiet=True
):
    """ Run redundancy reduction by clustering input (fasta) sequences and extracting cluster representatives. Writes
    a fasta file of redundancy reduced sequence. Note: most customisability of cluster_custom not currently exposed.

    # TODO: This is only being used to create temporary databases / clusterings used within other workflows
    """

    # Create seq DB
    create_seq_db(in_file, "temp_rr_fullDB")

    # Cluster
    cluster_custom("temp_rr_fullDB", "temp_rr_clustDB", min_seq_id, linclust=linclust, min_coverage=min_coverage,
                   max_threads=max_threads, nice=nice, quiet=quiet)

    # Extract cluster reps
    extract_cluster_reps("temp_rr_fullDB", "temp_rr_clustDB", fasta_name=rr_file)

    # Remove temp files
    for file in glob.glob(f"temp_rr_fullDB*"):
        os.remove(file)
    for file in glob.glob("temp_rr_clustDB*"):
        os.remove(file)


def extract_cluster_reps(
        seq_db,
        cluster_db,
        sub_db_name=None,
        fasta_name=None,
        priority_resample=None,
        quiet=True
):
    """ Extraction of cluster representatives from an mmseqs cluster database as a sub-database and/or Fasta file. """

    if not (sub_db_name or fasta_name):
        raise RuntimeError("A sub-database or fasta file name must be provided.")

    # Ensure database naming conventions
    seq_db = seq_db.split('DB')[0] + 'DB'
    cluster_db = cluster_db.split('DB')[0] + 'DB'

    # If sub_db is not to be retained, add temp tag so all generated files can be removed
    sub_db_arg = sub_db_name if sub_db_name else f"temp_{cluster_db}"

    # Create sub-DB with cluster reps
    if priority_resample:  # If cluster reps should be re-sampled

        # If external IDs provided, need to convert
        try:
            int(priority_resample[0][0])  # If no error, already have internal IDs
        except ValueError:  # Convert ext to int IDs
            id_map = db_id_map(seq_db, int2ext=False)
            priority_resample = [[id_map[ext_id] for ext_id in level if ext_id in id_map.keys()]
                                 for level in priority_resample]

            rep_map = priority_resample_cluster_reps(cluster_db, 1, priority_resample)[0]

            create_sub_db(seq_db, sub_db_arg, sub_db_ids=list(rep_map.values()))

    else:
        args = ["mmseqs", "createsubdb", cluster_db, seq_db, sub_db_arg]
        subprocess.run(args,
                       check=True,
                       stdout=subprocess.DEVNULL if quiet else None)

    if fasta_name:
        args = ["mmseqs", "convert2fasta", sub_db_arg, fasta_name]
        subprocess.run(args,
                       check=True,
                       stdout=subprocess.DEVNULL if quiet else None)

    if not sub_db_name:  # Remove sub-database files if not to be retained
        for file in glob.glob(f"temp_{cluster_db}*"):
            os.remove(file)


def single_cluster_subdb(
        seq_db,
        cluster_db,
        cluster_rep,
        id_type='int',
        sub_db_name=None,
        id_map=None,
        cluster_idx=None,
        quiet=True
):
    """ Extract all members of a given cluster from a cluster database to a sub-database. Either the internal or external
    cluster representative can be provided to identify the cluster of interest. If external ID type is provided, a
    mapping of external to internal IDs can be provided or computed here. Cluster reps index dictionary can be provided
    or created from .index file. """

    # Ensure database naming conventions
    seq_db = seq_db.split('DB')[0] + 'DB'
    cluster_db = cluster_db.split('DB')[0] + 'DB'
    sub_db_name = sub_db_name.split('DB')[0] + 'DB'

    # Map external to internal ID if required
    if id_type != 'int':

        if not id_map: # Create id_map - should already be read in for repetitive tasks!
            id_map = db_id_map(seq_db, int2ext=False)

        # Get internal ID
        try:
            cluster_rep = id_map[cluster_rep]
        except KeyError:
            raise RuntimeError(f"External ID {cluster_rep} not found in database {seq_db}.")

    # If cluster index not available, read it in
    if not cluster_idx:
        cluster_idx = {}
        with open(cluster_db + '.index') as idx_file:
            for line in idx_file:
                data = line.split()
                cluster_idx[data[0]] = int(data[1])  # str : int

    # Get file offset of cluster
    try:
        rep_idx = cluster_idx[cluster_rep]
    except KeyError:
        raise RuntimeError(f"Internal ID {cluster_rep} is not a cluster representative of {cluster_db}.")

    # Find cluster of interest
    with open(cluster_db) as cluster_file:
        cluster_file.seek(rep_idx)
        with open('tmp_list', 'w') as list_file:  # Write temp file containing cluster seqs
            while True:
                line = cluster_file.readline()
                if line.startswith('\0'):  # End of current cluster
                    break
                list_file.write(line)

    # Create sub-DB
    args = ["mmseqs", "createsubdb", "tmp_list", seq_db, sub_db_name]
    subprocess.run(args,
                   check=True,
                   stdout=subprocess.DEVNULL if quiet else None)
    os.remove('tmp_list')


def single_cluster_fasta(
        seq_db,
        cluster_db,
        cluster_rep,
        fasta_name,
        id_type='int',
        sub_db_name=None,
        id_map=None,
        cluster_idx=None,
        quiet=True
):
    """ Extract a cluster as a fasta file via a sub-database using single_cluster_db. The sub-database can be retained
    or deleted. """

    if not sub_db_name:  # Sub-database to be deleted at end
        sub_db_name = 'temp_DB'

    # Create subDB
    single_cluster_subdb(seq_db, cluster_db, cluster_rep, id_type=id_type, sub_db_name=sub_db_name, id_map=id_map,
                         cluster_idx=cluster_idx)

    # Create fasta file
    args = ["mmseqs", "convert2fasta", sub_db_name, fasta_name]
    subprocess.run(args,
                   check=True,
                   stdout=subprocess.DEVNULL if quiet else None)

    if sub_db_name == 'temp_DB':
        for file in glob.glob(f"temp_DB*"):
            os.remove(file)


def resample_cluster_reps(
        og_cluster_db,
        n_resamplings,
        cluster_map=None
):
    """ Sample representatives from each cluster, mapping all back to the original cluster rep. Returns a list of dicts
     mapping (internal) original cluster rep identifiers to the resampled rep identifier. """

    og_cluster_db = og_cluster_db.split('DB')[0] + 'DB'

    if not cluster_map:  # Make map of original reps to cluster members if not already available
        cluster_map = cluster_map(og_cluster_db)

    resamples_list = []

    for i in range(n_resamplings):
        resample = {}

        for og_rep, clust in cluster_map.items():
            # Randomly sample a new representative for each cluster
            resample[og_rep] = clust[random.randint(0, len(clust)-1)]

        resamples_list.append(resample)

    return resamples_list


def len_resample_cluster_reps(
        og_cluster_db,
        n_resamplings,
        seq_lens_dict,
        upper_len=None,
        lower_len=None,
        cluster_map=None
):
    """ For clusters with representatives outside some defined length range, attempt to re-sample a new rep from the
    same cluster within the bounds. If none is available, the original representative is retained.
    TODO: Exclusion of non-conforming clusters option """

    og_cluster_db = og_cluster_db.split('DB')[0] + 'DB'

    if not (upper_len or lower_len):
        raise RuntimeError("At least one of upper length bound or lower length bound must be provided for length-based "
                           "cluster representative re-sampling.")

    # If only one of upper or lower length bound is provided, make the other +inf or 0 (respectively)
    upper_len = upper_len if upper_len else math.inf
    lower_len = lower_len if lower_len else 0

    if not cluster_map:  # Make map of original reps to cluster members if not already available
        cluster_map = cluster_map(og_cluster_db)

    resamples_list = []

    for i in range(n_resamplings):
        resample = {}

        for og_rep, clust in cluster_map.items():
            # Check if og_rep is within preferred length range
            if lower_len < seq_lens_dict[og_rep] < upper_len:  # No need to re-sample
                resample[og_rep] = og_rep
            else:    # Re-sample from cluster reps within length range
                allowable = [seq for seq in clust if lower_len < seq_lens_dict[seq] < upper_len]
                if allowable:  # At least one member in desired range - randomly sample from these
                    resample[og_rep] = allowable[random.randint(0, len(allowable)-1)]
                else:   # None in desired range - stick with original
                    resample[og_rep] = og_rep

        resamples_list.append(resample)

    return resamples_list


def priority_resample_cluster_reps(
        og_cluster_db,
        n_resamplings,
        priority,
        cluster_map=None
):
    """ Re-sample cluster representatives if a member with higher priority is present in a given cluster. Priority
    is provided as a list of lists with the order of the outer list determining priority. Returns a list of
    cluster re-sampling maps ( {old_rep : new_rep} ).
    """

    og_cluster_db = og_cluster_db.split('DB')[0] + 'DB'

    if not cluster_map:  # Make map of original reps to cluster members if not already available
        cluster_map = cluster.cluster_map(og_cluster_db)  # TODO: Should rename this variable to avoid shadowing

    resamples_list = []

    for i in range(n_resamplings):

        resample = {}

        for og_rep, clust in cluster_map.items():
            found = False
            for level in priority:  # Find top priority level (if any) for which at least one cluster rep is present
                for member in clust:
                    if member in level:  # Highest priority present in cluster
                        if og_rep in level:  # If original rep is also in this priority level
                            resample[og_rep] = og_rep
                        else:  # Resample this member
                            resample[og_rep] = member
                        found = True
                        break
                if found:
                    break
            if not found:
                resample[og_rep] = og_rep

        resamples_list.append(resample)

    return resamples_list


def reconcile_rep_resample(og_clustering, rep_fasta=None, sub_clustering=None):
    """ Given a set of representatives resulting from a resampling, map the resampled  """


def assert_cluster_reps():
    return


def cascade_cluster_single(
        seq_db,
        og_cluster_db,
        new_cluster_db,
        min_seq_id,
        min_coverage=0.8,
        linclust=False,
        rep_map=None,
        sub_db_name=None,
        max_hits_per_query=None,
        tmp_dir=None,
        max_threads=20,
        nice=None,
        quiet=True
):
    """ Perform a single step cascade clustering of representatives from an existing cluster database. The original
     representatives can be used (default) or a mapping to alternative representatives can be provided. """

    seq_db = seq_db.split('DB')[0] + 'DB'
    og_cluster_db = og_cluster_db.split('DB')[0] + 'DB'
    new_cluster_db = new_cluster_db.split('DB')[0] + 'DB'

    if not sub_db_name:
        sub_db_name = "temp_subDB"
    else:
        sub_db_name = sub_db_name.split('DB')[0] + 'DB'

    # Get sub-db for clustering
    if not rep_map:  # Use original representatives
        args = ["mmseqs", "createsubdb", og_cluster_db, seq_db, sub_db_name]
        subprocess.run(args,
                       check=True,
                       stdout=subprocess.DEVNULL if quiet else None)
    else:
        create_sub_db(seq_db, sub_db_name, sub_db_ids=list(rep_map.values()), report_missing=True)

    # Cluster sub-db
    cluster_custom(sub_db_name, new_cluster_db, min_seq_id, min_coverage=min_coverage, linclust=linclust,
                max_hits_per_query=max_hits_per_query, tmp_dir=tmp_dir,
                   max_threads=max_threads, nice=nice, quiet=quiet)

    if sub_db_name == "temp_subDB":
        for file in glob.glob(f"temp_subDB*"):
            os.remove(file)


def cascade_cluster_multisample(
        seq_db,
        og_cluster_db,
        new_cluster_dbs,
        n_samples,
        min_seq_id,
        min_coverage=0.8,
        linclust=False,
        rep_maps=None,
        include_og_reps=False,
        priority=None,
        upper_len=None,
        lower_len=None,
        seq_lens_dict=None,
        og_cluster_map=None,
        sub_db_names=None,
        prefilter_names=None,
        align_names=None,
        max_hits_per_query=None,
        tmp_dir=None,
        remove_temp_db=True,
        max_threads=20,
        nice=None
):
    """ From an existing clustering of a sequence database, use n resamplings of cluster representatives and
     cascade cluster at a lower minimum sequence identity for each resampling. Specific mappings of cluster
     representatives can be asserted, or can be generated randomly (keep rep_maps=None), in which case the original
     representatives can optionally be enforced as the first 'sample'. """

    seq_db = seq_db.split('DB')[0] + 'DB'
    og_cluster_db = og_cluster_db.split('DB')[0] + 'DB'
    new_cluster_dbs = [new_cluster_db.split('DB')[0] + 'DB' for new_cluster_db in new_cluster_dbs]

    if not sub_db_names:
        sub_db_names = [None for i in range(n_samples)]
    if not prefilter_names:
        prefilter_names = [None for i in range(n_samples)]
    if not align_names:
        align_names = [None for i in range(n_samples)]

    if not rep_maps:  # Need to generate resampled cluster reps

        if include_og_reps:  # Use original reps as first sample

            if priority:  # Re-sampling based on priority sets required

                first_sample = priority_resample_cluster_reps(og_cluster_db, 1, priority, cluster_map=og_cluster_map)[0]

            elif upper_len or lower_len:  # Length-based cluster rep re-sampling required

                if not seq_lens_dict:   # Compute sequence length dictionary if not available
                    seq_lens_dict = db_seq_lens(seq_db)

                # Get single length-based resampling of cluster reps
                first_sample = len_resample_cluster_reps(og_cluster_db, 1, seq_lens_dict, upper_len, lower_len,
                                                         cluster_map=og_cluster_map)[0]

            else:  # No re-sampling required

                if og_cluster_map:
                    og_reps = list(og_cluster_map.keys())
                else:  # Need to extract from database
                    og_reps = get_internal_ids(og_cluster_db)

                first_sample = {rep : rep for rep in og_reps}

            # First sample is og_reps (possibly length-resampled), remainder are completely random samples
            rep_maps = [first_sample] + resample_cluster_reps(og_cluster_db, n_samples-1, cluster_map=og_cluster_map)

        else:  # Randomly sample all
            rep_maps = resample_cluster_reps(og_cluster_db, n_samples, cluster_map=og_cluster_map)

    # Ensure relevant parameters supplied for all resamplings
    if not (n_samples == len(rep_maps) == len(new_cluster_dbs) == len(sub_db_names) == len(prefilter_names) ==
            len(align_names)):
        raise RuntimeError("Parameter list length mismatch.")

    # Perform cascade clusterings for each resampling of representatives
    for i in range(n_samples):
        cascade_cluster_single(
            seq_db,
            og_cluster_db,
            new_cluster_dbs[i],
            min_seq_id,
            min_coverage=min_coverage,
            linclust=linclust,
            rep_map=rep_maps[i],
            sub_db_name=sub_db_names[i],
            max_hits_per_query=max_hits_per_query,
            tmp_dir=tmp_dir,
            max_threads=max_threads,
            nice=nice
        )


def cluster_comparison_single(
        clustering_1,
        clustering_2,
        rep_maps=None,
        metric="NMI"
):
    """ Compare two clusterings of the same original database by an available metric. If the clusterings are cascade
    clustering and the original cluster representatives have not been used, mappings of {original_rep : resampled_rep}
    must be provided. Clusterings can be specified as either the database name or a pre-computed cluster map
    {rep : [members]}.

    TODO: Add option to provide original database and compute the rep_maps in case of comparing cascade clusterings from same og database"""

    # Note possibly confusing terminology: 'representative' could be reps of the provided (new) clusterings or in the
    # case of cascade clustering, differently sampled members of clusters which now represent ALL members of new clusters

    if type(clustering_1) is type(clustering_2) is str:
        # Computer cluster maps
        clustering_1 = cluster_map(clustering_1.split('DB')[0] + 'DB')
        clustering_2 = cluster_map(clustering_2.split('DB')[0] + 'DB')
    elif not (type(clustering_1) is type(clustering_2) is dict):
        raise RuntimeError("Clusterings must both be provided as either mmseqs databases or cluster maps.")

    if not rep_maps:  # Assume that both clusters contain the same sequence IDs (using same reps if cascade clustering)
        int_ids_1 = clustering_1.keys()
        int_ids_2 = clustering_2.keys()
        if set(int_ids_1) != set(int_ids_2):
            raise RuntimeError("Internal ID sets of clusterings are not consistent - supply a mapping if different "
                               "representative sequences were used to cascade cluster.")
        rep_maps = [{int_id: int_id for int_id in int_ids_1} for i in range(2)]  # Format matches user-provided rep_maps

    else:  # Need to reverse the old:new rep maps to be new:old rep maps
        rep_maps = [{new: old for old, new in rep_map.items()} for rep_map in rep_maps]

    # Replace IDs in clusterings with mapped og IDs

    mapped_clustering_1 = {}
    for rep, members in clustering_1.items():
        mapped_rep = rep_maps[0][rep]
        mapped_members = []
        for member in members:
            mapped_members.append(rep_maps[0][member])
        mapped_clustering_1[mapped_rep] = mapped_members

    mapped_clustering_2 = {}
    for rep, members in clustering_2.items():
        mapped_rep = rep_maps[1][rep]
        mapped_members = []
        for member in members:
            mapped_members.append(rep_maps[1][member])
        mapped_clustering_2[mapped_rep] = mapped_members

    # Check that the set of clustered members is identical between clusterings
    members_1 = set()
    for clust in mapped_clustering_1.values():
        members_1.update(clust)
    members_2 = set()
    for clust in mapped_clustering_2.values():
        members_2.update(clust)
    if members_1 != members_2:
        raise RuntimeError("Total sets of cluster members not identical between provided clusterings.")

    # Define ordering for input to cluster comparison metrics
    members_list = list(members_1)
    members_list.sort()

    # Get map of member : cluster_rep for both clusterings
    member2rep_1 = {}
    for rep, members in mapped_clustering_1.items():
        for member in members:
            member2rep_1[member] = rep

    member2rep_2 = {}
    for rep, members in mapped_clustering_2.items():
        for member in members:
            member2rep_2[member] = rep

    # Get cluster labels (mapped back to og reps) for all members for each clustering
    labels_1 = [member2rep_1[member] for member in members_list]
    labels_2 = [member2rep_2[member] for member in members_list]

    # Compare clusterings with specified metric
    if metric == "MI":
        return metrics.mutual_info_score(labels_1, labels_2)
    elif metric == "NMI":
        return metrics.normalized_mutual_info_score(labels_1, labels_2)
    elif metric == "AMI":
        return metrics.adjusted_mutual_info_score(labels_1, labels_2)
    elif metric == "ARI":
        return metrics.adjusted_rand_score(labels_1, labels_2)
    else:
        raise RuntimeError(f"Clustering comparison metric {metric} is not currently implemented.")


def cluster_comparison_multi(
        clusterings,
        rep_maps=None,
        metric="NMI",
        clustering_names=None
):
    """ Make single (pairwise) comparisons for all combinations of provided clusterings. Determine the mean metric
    across all pairwise comparisons, and the clustering with the highest mean for comparisons in which it is involved.
    Returns (pairwise_metric_dict, mean_pw_metric, individual_clustering_means, best_clustering). """

    if not isinstance(clusterings, list):
        raise RuntimeError("Multiple clusterings must be provided as a list of database names or cluster maps.")
    else:
        if rep_maps:
            if len(rep_maps) != len(clusterings):
                raise RuntimeError("If rep_maps are provided, they must be given as a list of equal length to clusterings.")
        if clustering_names:
            if len(clustering_names) != len(clusterings):
                raise RuntimeError("If clustering names are provided, they must be given as a list of equal length to clusterings.")
        else:   # Just name them numerically
            clustering_names = [i for i in range(len(clusterings))]

    score_map = {}   # {tuple(clustering1,clustering2) : score}

    # Get score for all pairwise comparisons
    for i in range(len(clusterings)-1):
        for j in range(i+1, len(clusterings)):
            score_map[tuple([clustering_names[k] for k in [i,j]])] = \
                cluster_comparison_single(clusterings[i], clusterings[j],
                                          rep_maps=[rep_maps[i], rep_maps[j]] if rep_maps else None, metric=metric)

    # Compute total metric mean
    total_mean = np.mean(list(score_map.values()))

    # Compute the mean of each clustering's pairwise comparisons to each other clustering
    individual_means = []
    for name in clustering_names:
        individual_means.append((name, np.mean([score for pair, score in score_map.items() if name in pair])))

    # Extract best clustering
    best = max(individual_means, key=lambda x: x[1])[0]

    return score_map, total_mean, individual_means, best


class HierarchicalClustering:
    """ A multi-level Hierarchical clustering of a sequence database. Clusterings at lower sequence identities should
    be cascade clusters from a clustering at a higher sequence identity. While the default setting for most automatic
    workflows is multi-sampling of cluster representatives for cascade clusterings, followed by selection of the
    clustering with the highest mean pairwise agreement for further clustering, tree-like structures of relationships
    between clusterings is not enforced. Currently supports a single clustering at the highest sequence identity. """


    def __init__(
            self,
            master_seq_db,
            top_clustering=None,
            parents=None,
            sample_scores=None,
            level_names=None,
            clustering_scores=None,
            metrics=None,
            parent_default=False,
            priority=None,
            upper_lens=None,
            lower_lens=None,
            hc_name=None
    ):

        self.master_seq_db = master_seq_db.split('DB')[0] + 'DB'

        if top_clustering:
            self.top_clustering = top_clustering.split('DB')[0] + 'DB'
        else:
            self.top_clustering = None

        self.hc_name = hc_name  # Tag to be added to DB names to identify this analysis

        # List of metrics for clustering comparison
        self.metrics = metrics if metrics else []

        # Mean of scores for samples at corresponding level of hierarchy
        self.clustering_scores = clustering_scores if clustering_scores else []
        self.parent_default = parent_default    # True when the first sample of a level is always used as next parent

        # Attributes associated with levels of the hierarchy - lists should all be the same length and indexed equivalently
        self.level_names = level_names if level_names else []

        # Mapping of {sample:[scores]} for each level of the hierarchy
        self.sample_scores = sample_scores if sample_scores else []
        # [{sample1.1:[scores], sample1.2:[scores]....}, {sample2.1:[scores], sample2.2:[scores]...}.......]
        self.parents = parents

        self.priority = priority

        # If parent_default is True, caveats can be placed on length of representatives (where alternatives are available)

        if not upper_lens:
            self.upper_lens = [None for level in self.level_names]
        else:
            self.upper_lens = upper_lens

        if not lower_lens:
            self.lower_lens = [None for level in self.level_names]
        else:
            self.lower_lens = lower_lens

        # If overall clustering scores are not explicitly available, they can be computed from sample scores
        if self.sample_scores and self.metrics and not self.clustering_scores:
            self.clustering_scores = [[np.mean([scores[i] for scores in score_dict.values()])
                                       for i in range(len(self.metrics))] for score_dict in self.sample_scores]

        # If level_names are not explicitly provided, extract from sample name assuming naming conventions are followed
        if self.sample_scores and not self.level_names:
            self.level_names = [list(score_dict.keys())[0].rsplit('_', maxsplit=1)[0] for score_dict in self.sample_scores]

        # Check hierarchy level attribute lists are all == len
        if self.level_names and self.sample_scores and self.parents:
            if not (len(self.level_names) == len(self.sample_scores) == len(self.parents) == len(self.upper_lens) ==
                    len(self.lower_lens)):
                raise RuntimeError("Attribute lists for hierarchy levels must be of equal length!")


    def run_full_hc(
            self,
            top_min_seq_id,
            top_min_coverage,
            min_seq_ids,
            min_coverages,
            metrics,
            n_samples,
            linclust=False,
            top_name=None,
            level_names=None,
            sample_names=None,
            include_og_reps=False,
            parent_og_default=False,
            priority=None,
            in_fasta=None,
            upper_lens=None,
            lower_lens=None,
            max_threads=20,
            nice=None,
            quiet=True
    ):
        """ Perform hierarchical clustering at provided thresholds. If multi-sampling of cluster representatives is used
        for each cascade clustering, the first metric of self.metrics will be used to select the best clustering.
        Priority for certain sequences to be used as reps can be specified as a list of lists, where the order of the
        outer list dictates priority. Alternatively, upper and lower bounds can be set on representative sequences.

        NOTE: For now, assumes that only the original database is available (i.e., need to do the top level clustering
        and then all samples for all levels. """

        # Seems to want to re-use old clusterings at same %id - maybe clearing tmp prevents it?
        if os.path.exists('tmp'):
            shutil.rmtree('tmp')
        os.mkdir('tmp')

        if not top_name:  # Infer name from master seq DB, hc_name (if provided), and min seqID/coverage
            self.top_clustering = self.master_seq_db.split('DB')[0] + (f"_{self.hc_name}" if self.hc_name else "") \
                                  + f"_clust{int(top_min_seq_id*100)}" + f"_cov{int(top_min_coverage*100)}" + "DB"
        else:
            self.top_clustering = top_name.split("DB")[0] + "DB"

        if not level_names:
            self.level_names = \
                [self.master_seq_db.split('DB')[0] + (f"_{self.hc_name}" if self.hc_name else "") +
                 f"_clust{int(min_seq_ids[i]*100)}" + f"_cov{int(min_coverages[i]*100)}" for i in range(len(min_seq_ids))]
        else:
            self.level_names = level_names

        if isinstance(n_samples, int):  # Same number of samples for each level
            n_samples = [n_samples for i in range(len(min_seq_ids))]

        if not sample_names:  # Use generic sample naming if not provided
            sample_names = [[f"sample{j}" for j in range(n_samples[i])] for i in range(len(n_samples))]

        # Check if priority sets have been provided as internal (integer) IDs and otherwise map external to internal IDs
        if priority:
            try:
                int(priority[1][0])
                self.priority = priority   # If no error, already have internal IDs
            except ValueError:  # Need to convert to internal IDs
                id_map = db_id_map(self.master_seq_db, int2ext=False)
                self.priority = [[id_map[ext_id] for ext_id in level if ext_id in id_map.keys()] for level in priority]
        else:
            self.priority = []

        # Optional re-sampling of original mmseqs representatives (sample 0) if length is outside a defined range
        # Any length bounds are ignored if priority sets are provided

        if not self.priority:
            # Upper bounds
            if isinstance(upper_lens, list):
                self.upper_lens = upper_lens
            elif isinstance(upper_lens, int):   # All levels use same max
                self.upper_lens = [upper_lens for level in self.level_names]
            else:  # No upper limits
                self.upper_lens = [None for level in self.level_names]
            # Lower bounds
            if isinstance(lower_lens, list):
                self.lower_lens = lower_lens
            elif isinstance(lower_lens, int):   # All levels use same min
                self.lower_lens = [lower_lens for level in self.level_names]
            else:  # No lower limits
                self.lower_lens = [None for level in self.level_names]
        else:
            self.upper_lens = [None for level in self.level_names]
            self.lower_lens = [None for level in self.level_names]

        # If any length bounds are imposed, need to compute length dictionary of DB sequences

        if sum([True if value else False for value in self.upper_lens + self.lower_lens]) > 0:
            seq_len_dict = db_seq_lens(self.master_seq_db)
        else:
            seq_len_dict = None

        if not (len(min_seq_ids) == len(min_coverages) == len(self.level_names) == len(n_samples) ==
                len(self.upper_lens) == len(self.lower_lens)):
            raise RuntimeError("Parameter lists for hierarchy levels must be of equal length.")

        if not isinstance(metrics, list):
            self.metrics = [metrics]
        else:
            self.metrics = metrics

        if include_og_reps and parent_og_default:  # og reps are the first sample and should be used as parent clustering
            self.parent_default = True

        if self.hc_name:
            for file in glob.glob(f"{self.master_seq_db.split('DB')[0]}_{self.hc_name}*"):
                os.remove(file)

        # Run a single clustering at the top level
        cluster_custom(self.master_seq_db, self.top_clustering, min_seq_id=top_min_seq_id, min_coverage=top_min_coverage,
                       linclust=linclust, max_threads=max_threads, nice=nice, quiet=quiet)
        self.parents = [self.top_clustering]

        self.in_fasta = in_fasta
        self.clustering_scores = []
        self.sample_scores = []

        # Iteratively multi-sample cascade cluster from best clustering sample at each level
        for i in range(len(min_seq_ids)):

            # Apply sample names to new DB names for this level
            clust_db_names = [f"{self.level_names[i]}_{sample_names[i][j]}DB" for j in range(len(sample_names[i]))]

            # Run multi-sample cascade clustering
            cascade_cluster_multisample(
                self.master_seq_db,
                self.parents[-1],
                clust_db_names,
                len(clust_db_names),
                min_seq_id=min_seq_ids[i],
                min_coverage=min_coverages[i],
                linclust=linclust,
                include_og_reps=include_og_reps,
                priority=self.priority,
                upper_len=self.upper_lens[i],
                lower_len=self.lower_lens[i],
                seq_lens_dict=seq_len_dict,
                og_cluster_map=cluster_map(self.parents[-1]),
                max_threads=max_threads,
                nice=nice
            )

            # Need to get rep_maps
            # TODO: can avoid doing this again here by having multi-sample cascade clustering return computed rep_maps
            rep_maps = [rep_map_reconstruction(
                self.parents[-1],
                sub_cluster_db=clust_db
            ) for clust_db in clust_db_names]

            # Run cluster comparison metrics for this level

            # Metric score lists/dicts for this level (overall + individual samples)
            self.clustering_scores.append([])  # Will add overall metric means for each metric for this level
            # Will add individual metric scores for each sample and metric for this level
            self.sample_scores.append({name : [] for name in clust_db_names})

            for k in range(len(self.metrics)):    # Note: just using db names as sample names here
                comp_results = cluster_comparison_multi(clust_db_names, rep_maps, self.metrics[k], clust_db_names)
                self.clustering_scores[-1].append(comp_results[1])  # Average over all sample clusterings
                for j in range(len(clust_db_names)):
                    self.sample_scores[-1][clust_db_names[j]].append(comp_results[2][j][1])
                # Next parent is best performing under first listed metric unless original is default
                if k == 0 and not self.parent_default:
                    self.parents.append(comp_results[3])

            if self.parent_default:  # Next parent is the "sample" with original mmseqs cluster reps
                self.parents.append(clust_db_names[0])




    # TODO: Complete cluster mapping functions

    # def cluster_by_id(self, int_id, level=None, lowest=False, sub_db_name=None, fasta_name=None):
    #     """ Extract a cluster as a sub-database and/or fasta file at
    #     TODO: Need mapping from self.master_db - not working with sub-databases. Could read in entire ID mapping but going to be
    #     a problem with RAM for very large databases (e.g. Uniparc MBL domain hits is too big) """
    #
    #
    # def get_clustering_dict(self):
    #     """ Get a dictionary mapping each sequence to the list of representative sequences at each level of the
    #     clustering hierarchy. Sequences clustered at higher levels of the hierarchy are assumed to belong to the same
    #     clusters as their cluster's representatives at the next level. """
    #
    #     # Instantiate dictionary
    #     clust_dict = {ext_id : [] for ext_id in get_external_ids(self.master_seq_db)}
    #
    #     # Append parent clustering ID maps
    #     os.system(f"mmseqs createtsv {self.master_seq_db} {self.master_seq_db} {self.top_clustering} temp.tsv")
    #     with open("temp.tsv") as in_file:
    #         for line in in_file:
    #             data = line.strip().split("\t")
    #             clust_dict[data[1]] = data[0]
    #
    #     # Iteratively follow cascade clusterings to assign clusters for each sequence at each hierarchy level
    #     # for level in self.parents:   # TODO: Complete this function
    #
    #
    # def write_id_clustering(self, file_name, ext_ids=True):
    #     """ Write a hierarchy of cluster membership for each sequence in a tab-delimited format. External IDs (default)
    #      or ..... (need to see how easy it is to make this internal, although could just read in an ID mapping to do this,
    #       however that may pose problems for spikes in RAM requirements) """


    def write_hc_metrics(self, filename):
        """ Write hierarchical clustering metrics to file. """

        # TODO: Add upper/lower rep length params

        filename = filename.split('.hc')[0] + '.hc'   # Enforce hierarchical clustering extension

        with open(filename, 'w') as file:

            # Write global parameters
            file.write("Parameters\n")
            file.write(f"\tMaster Sequence DB:\t{self.master_seq_db}\n")
            metrics_str = "\t".join(self.metrics)
            file.write(f"\t{metrics_str}\n")
            file.write("//\n")

            # Write info on top clustering (currently only supports a single top-level clustering)
            file.write("Parent Clustering\n")
            file.write(f"\t{self.top_clustering}\t{getDBLen(self.top_clustering)}\n")
            file.write("//\n")

            # Write remaining (sub/cascade) clusters
            file.write("Sub Clusterings\n")
            for i in range(len(self.level_names)):
                scores_str = "\t".join([str(round(score, 4)) for score in self.clustering_scores[i]])
                file.write(f"\t{self.level_names[i]}\tParent Clustering: {self.parents[i]}\t{scores_str}\n")

                # Write results for each sampling
                for j in range(len(self.sample_scores[i])):
                    sample_name = list(self.sample_scores[i].keys())[j]
                    metric_str = "\t".join([str(round(score, 4)) for score in self.sample_scores[i][sample_name]])
                    # \t \t sample name \t db_len \t metric1 \t metric2 \t .... \t ?NEXT_PARENT?
                    file.write(f"\t\t{sample_name}\t{getDBLen(sample_name)}\t{metric_str}\t"
                               f"{'NEXT' if (self.parents[i+1] == sample_name and j < len(self.sample_scores[i])-1) else ''}\n")

            file.write("//\n")


    def parent_reps_to_fasta(
            self,
            levels=None,
            fasta_naming=None,
            priority_resample=True
    ):
        """ Extract representatives of each parent clustering for the specified levels (if None, extract for all levels). """

        priority = self.priority if priority_resample else None

        if levels:
            levels = [level.split('DB')[0] + 'DB' for level in levels]
        else:
            levels = self.parents

        if not fasta_naming:  # Function for transforming DB name to fasta output file name
            fasta_naming = lambda x: x.split('DB')[0] + '.fa'

        for parent in levels[:-1]:  # All except bottom level will already have any re-sampling applied
            extract_cluster_reps(
                self.master_seq_db,
                parent,
                fasta_name=fasta_naming(parent)
            )

        # May need priority re-sample of bottom level
        extract_cluster_reps(
            self.master_seq_db,
            levels[-1],
            fasta_name=fasta_naming(parent),
            priority_resample=priority
        )


    def extract_expanded_clusters(self, file_naming=None):
        """ For each level of the hierarchical clustering, extract all sequences clustered together at all higher
        levels. For example, an expanded cluster from a 70% clustering would contain the 80% clustering members of each
        70% member, and so on. If the bottom level of the hierarchy has been further re-sampled, provide the resampled
        representatives as a fasta file. Clusterings are written as .cluster files in the fasta-like cluster format
        used by mmseqs2. If representatives have been re-sampled using priority at each level, expanded clusters will
        use the re-sampled representative. """

        if not file_naming:
            file_naming = lambda x: x + '.cluster'

        # Get rep maps and cluster maps for each level
        cluster_maps = [cluster_map(self.parents[i]) for i in range(len(self.parents))]

        rep_maps = [rep_map_reconstruction(self.parents[i], sub_cluster_db=self.parents[i+1]) for i in range(len(self.parents)-1)]
        # Need to do priority resampling of lowest level clustering reps
        rep_maps.append(priority_resample_cluster_reps(self.parents[-1], 1, self.priority)[0])

        # Get seq DB ID map and map ext ID to seq string
        id_map = db_id_map(self.master_seq_db)
        ext_seq_map = {}
        for seq in SeqIO.parse(self.in_fasta, "fasta"):
            ext_seq_map[seq.name] = seq.seq


        for i in range(len(self.parents)):

            cluster_file = file_naming(self.parents[i])

            # Get rep map for this level
            level_rep_map = rep_maps[i]
            level_cluster_map = cluster_maps[i]

            if i == 0:  # Handle top level clustering separately

                # ext_rep_map = {}

                # Remove .cluster file if it exists
                if cluster_file in os.listdir():
                    os.remove(cluster_file)

                for old_rep, new_rep in level_rep_map.items():

                    members = set(level_cluster_map[old_rep])

                    # Convert expanded maps to external seq IDs
                    ext_rep = id_map[new_rep]
                    ext_members = [id_map[member] for member in members]
                    ext_members.remove(ext_rep)  # Rep will be handled separately

                    # Write .cluster file
                    with open(cluster_file, "a") as file:
                        file.write(f">{ext_rep}\n>{ext_rep}\n{ext_seq_map[ext_rep]}\n")
                        for member in ext_members:
                            file.write(f">{member}\n{ext_seq_map[member]}\n")

            else:  # Use previous level's cluster file

                # Read in previous level's .cluster file
                prev_file_name = file_naming(self.parents[i-1])
                prev_clust_ext = extract_cluster_ext(prev_file_name, clustering_map=True)

                # Remove .cluster file if it exists
                if cluster_file in os.listdir():
                    os.remove(cluster_file)

                for old_rep, new_rep in level_rep_map.items():

                    members = level_cluster_map[old_rep]
                    members_ext = [id_map[member] for member in members]
                    exp_members_ext = []
                    for ext_member in members_ext:
                        exp_members_ext.extend(prev_clust_ext[ext_member])
                    ext_rep = id_map[new_rep]
                    exp_members_ext.remove(ext_rep)

                    with open(cluster_file, 'a') as file:
                        file.write(f">{ext_rep}\n>{ext_rep}\n{ext_seq_map[ext_rep]}\n")
                        for member in exp_members_ext:
                            file.write(f">{member}\n{ext_seq_map[member]}\n")


                # for seq in SeqIO.parse(self.in_fasta, "fasta"):
                #     if seq.name == ext_rep:
                #         rep_str = seq.seq
                #     elif seq.name in ext_members:
                #         member_seq_strs[seq.name] = seq.seq
                #
                # with open(cluster_file, 'a') as file:
                #     file.write(f">{ext_rep}\n>{ext_rep}\n{rep_str}\n")
                #     for seq_id, seq_str in member_seq_strs.items():
                #         file.write(f">{seq_id}\n{seq_str}\n")

                #
                # # Create temporary sub-DB containing expanded cluster members and extract to fasta
                # with open("temp_list.txt", 'w') as temp_file:
                #     for member in expanded_members:
                #         temp_file.write(f"{member}\n")
                # os.system(f"mmseqs createsubdb temp_list.txt {self.master_seq_db} temp_listDB")
                # os.system(f"mmseqs convert2fasta temp_listDB temp_list.fa")
                #
                # # Start new cluster with >rep_name
                # with open(cluster_file, 'a') as file:
                #     file.write('>' + id_map[new_rep] + '\n')

                # # Append all cluster members and remove temp files
                # os.system(f"cat temp_list.fa >> {cluster_file}")
                # os.system(f"rm temp_list*")


# def reps_from_cluster_file(ext_cluster_file, out_file=None):
#     """
#     NOTE: this seems to be buggy. Have re-written below
#     Extract representative sequences from a .cluster file as a fasta file.
#      """
#
#     if not out_file:
#         out_file = ext_cluster_file.split('.cluster')[0] + '_reps.fa'
#
#     with open(ext_cluster_file) as in_file:
#         with open(out_file, 'w') as out_file:
#             header = None  # ID if previous line was a header
#             rep = None  # Representative of current cluster
#             for line in in_file:
#                 if header:
#                     if line.startswith('>'):  # Two header lines in a row means new cluster
#                         rep = header  # Will write rep sequence when it is encountered
#                         header = line.strip().strip('>')
#                     elif header == rep:
#                         out_file.write(f">{header}\n{line.strip()}\n")
#                         header = None
#                 else:
#                     if line.startswith('>'):
#                         header = line.strip().strip('>')
#                     else:
#                         raise RuntimeError(f"Cluster file {ext_cluster_file} is not correctly formatted.")


def reps_from_cluster_file(ext_cluster_file, rep_file=None):
    """ Extract representative sequences from a .cluster file as a fasta file. """

    with open(ext_cluster_file) as in_file:

        if not rep_file:
            rep_file = ext_cluster_file.splilt(".cluster")[0] + "_reps.fa"

        with open(rep_file, 'w') as out_file:

            rep = None
            prev_header = None
            found_rep = False
            for line in in_file:

                if found_rep:
                    if line.startswith(">"):
                        found_rep = False
                    else:
                        out_file.write(line)

                if prev_header and line.startswith(">"):
                    rep = prev_header
                    prev_header = None
                if line.startswith(">"):
                    if line.strip().split(">")[1] == rep:
                        found_rep = True
                        out_file.write(line)
                    prev_header = line.strip().split(">")[1]
                else:
                    prev_header = False


def extract_cluster_ext(
        ext_cluster_file,
        rep_fasta=None,
        clustering_map=False
):
    """ Parse fasta-like external cluster file into dictionary mapping {rep : [all_members]}, or write a fasta file
     containing representative sequences, or both.

     #TODO: cluster files are having whole header info written for members - need to fix"""

    if not (rep_fasta or clustering_map):
        raise RuntimeError("Either a clustering map must be returned, or a representative fasta file written, or both.")

    rep_fasta = rep_fasta if rep_fasta else "temp_reps.fa"

    clust_map = {}

    with open(ext_cluster_file) as in_file:
        with open(rep_fasta, 'w') as out_file:
            header = False  # Was there a header in the previous line?
            rep = None  # Representative of current cluster
            rep_seq = None
            for line in in_file:
                if header:
                    if line.startswith('>'):  # Two header lines in a row means new cluster
                        if rep:  # Record previous cluster if this isn't the first one
                            clust_map[rep] = members
                        rep = header.split()[0]
                        out_file.write(f">{rep}\n")
                        members = [line.strip().strip('>').split()[0]]
                        if members[0] == rep:
                            rep_seq = ""
                    else:  # This is just a continuation of current cluster
                        members.append(header)
                        rep_seq = None
                    header = False
                else:
                    if line.startswith('>'):
                        if isinstance(rep_seq, str):
                            out_file.write(rep_seq + '\n')
                            rep_seq = None
                        elif line.strip().strip('>').split()[0] == rep:
                            rep_seq = ''
                        header = line.strip().strip('>').split()[0]
                    elif isinstance(rep_seq, str):
                        rep_seq += line.strip()
        # Reached EOF and need to record the last cluster
        clust_map[rep] = members


    if rep_fasta == "temp_reps.fa":
        os.remove("temp_reps.fa")

    return clust_map


def get_hc_maps(cluster_files=None, name=None):
    """ Return a list of cluster maps of external IDs in the form {rep : [all_members]} at each level of hierarchical
    clustering. This function is only applied to fasta-like cluster files such that any re-sampling of cluster
    representative should be already reflected. Cluster file names can be provided as a list or a name can be provided
    if files follow the naming conventions of this module, in which case all .cluster files containing the name in the
    current directory are used. """

    if not (cluster_files or name):
        raise RuntimeError("File names must be provided either explicitly as a list or implicitly as a hc name.")

    if not cluster_files:  # Find all files following naming conventions with the clustering name and map to %
        file_list = [(file, int(file.split('_cov')[0].split('_clust')[1])) for file in os.listdir()
                     if (name in file and file.endswith('.cluster'))]
        file_list.sort(key=lambda x:x[1], reverse=True)
        cluster_files = [file[0] for file in file_list]

    return [extract_cluster_ext(file, clustering_map=True) for file in cluster_files]


def get_hc_level_info(cluster_files=None, name=None):
    """ Returns a list of tuples [(file_name, seq_id_threshold, seq_cov_threshold), ...] ordered high to low by
     sequence identity threshold. """

    if not cluster_files:
        cluster_files = [file for file in os.listdir() if (name in file and file.endswith(".cluster"))]

    hc_info = [(file, int(file.split('_cov')[0].split('_clust')[1])) for file in cluster_files]
    hc_info.sort(reverse=True, key=lambda x:x[1])

    return hc_info


def read_hc_file(filename):
    """ Create a HierarchicalClustering object from a .hc file. All relevant databases for the HC should be in the
    same directory. """


def hc_from_fasta(
        in_fasta,
        db_name=None,
        min_seq_ids=None,
        min_coverages=None,
        linclust=False,
        priority=None,
        max_threads=10,
        nice=10,
        hc_file_name=None,
        cluster_file_naming=None,
        skip_cluster_files=False
):
    """ Runs hierarchical clustering workflow directly from a Fasta file. Note that this workflow runs only a single
     clustering sample at each min seq identity and coverage and hence does not assess cluster confidence. The top level
      cluster parameters should be included as the first entry in parameter lists here. """

    # Create database
    if not db_name:
        db_name = in_fasta.split('.')[0] + 'DB'
    create_seq_db(in_fasta, db_name)

    # Priority can't be None
    if not priority:
        priority = []

    if not min_seq_ids:
        min_seq_ids = [0.9, 0.8, 0.7, 0.6, 0.5]
    if not min_coverages:
        min_coverages = [0.8, 0.8, 0.8, 0.7, 0.7]

    # Set up and run HC
    hc = HierarchicalClustering(db_name)
    hc.run_full_hc(
        top_min_seq_id=min_seq_ids[0],
        top_min_coverage=min_coverages[0],
        min_seq_ids=min_seq_ids[1:],
        min_coverages=min_coverages[1:],
        linclust=linclust, metrics=[],
        n_samples=1,
        priority=priority,
        include_og_reps=True,
        parent_og_default=True,
        max_threads=max_threads,
        nice=nice,
        in_fasta=in_fasta)

    # Write cluster files for each clustering level
    if not skip_cluster_files:
        if not cluster_file_naming:
            cluster_file_naming = lambda x: x.split("DB")[0].split("_sample")[0] + ".cluster"
        hc.extract_expanded_clusters(file_naming=cluster_file_naming)

    # Write hc metrics file
    if not hc_file_name:
        hc_file_name = db_name.split('DB')[0] + '.hc'
    hc.write_hc_metrics(hc_file_name)


def expand_cluster_reps(
        hc_maps,
        lower=None,
        upper=None,
        lower_reps=None
):
    """ Given select cluster representatives of a low-level clustering (the lowest by default, else an index can be
    provided corresponding to the order in the list), expand the selected clusters by climbing back up the hierarchy to
    the upper limit (the un-clustered level by default, else an index can be provided indicating the final level for
    which reps should be included in the expanded clusters. """

    if lower is None:
        lower = -1  # Starting point for cluster expansion is bottom level (i.e. reps from bottom level are expanded reps)

    if upper is None:
        upper = -1  # Expand all the way back up the hierarchy such that all sequences are included

    lower_clusters = {rep : members for rep, members in hc_maps[lower].items()}

    # Subset clustering
    sub_lower_clusters = {rep : members for rep, members in lower_clusters.items() if rep in lower_reps} if lower_reps \
        else lower_clusters

    # Expanded members in select lower level clusters
    sub_members = set()
    for clust in sub_lower_clusters.values():
        sub_members.update(clust)

    # Need to exclude sequences which AREN'T reps at upper level
    exclude = set()
    if upper >= 0:

        sub_upper_reps = set([rep for rep in hc_maps[upper].keys()])
        exclude = sub_members - sub_upper_reps

    if exclude:  # If we have an upper clustering level
        filt_lower_clusters = {}
        for rep, members in sub_lower_clusters.items():
            filt_lower_clusters[rep] = list(set(members) - exclude)  # Remove non-reps from upper clustering
    else:
        filt_lower_clusters = sub_lower_clusters

    return filt_lower_clusters


def cluster_aln(args):
    """ Function to be called by multi-processing for clustering multiple sub-alignments for clusters. The rep and list
     of members should be passed as a double (cluster_data). """

    rep, members, seq_fa, naming, remove_sub_fa, nice = args  # Unpack arguments

    if not naming:
        fa_name = f"temp_sub_fa_{rep}.fa"
        aln_name = f"{seq_fa.split('.')[0]}_sub_{rep}.aln"
    else:
        fa_name = naming(seq_fa, rep) + '.fa'
        aln_name = naming(seq_fa, rep) + '.aln'

    # Extract seqs to sub fasta and align (mafft-linsi)
    file_util.raw_extract_fasta(seq_fa, fa_name, members, id_format=None)

    if len(members) > 1:
        aln.linsi_aln(fa_name, aln_name, quiet=True, nice=nice)
    else:  # Just single sequence to this cluster so no need to align
        shutil.copy(fa_name, aln_name)

    if remove_sub_fa:
        os.remove(fa_name)


def cluster_sub_alns(
        seq_fa,
        cluster_map,
        naming=None,
        remove_sub_fa=True,
        max_processes=6,
        nice=None
):
    """ Generate alignments for each cluster from a provided cluster map. File naming can be specified as a lambda
     function acting on the sequence fasta name and the representative name of the cluster, else default naming will be
     be applied. Provided naming functions should output the desired file name excluding the extension. """

    # Create list of arguments to pass to cluster_aln in separate processes
    tasks = []
    for rep, members in cluster_map.items():
        tasks.append((rep, members, seq_fa, naming, remove_sub_fa, nice))

    # Create all cluster sub alignments using multi-processing
    with multiprocessing.Pool(max_processes) as pool:
        pool.map(cluster_aln, tasks)


def default_sub_aln_naming(seq_fa, rep):
    """ Helper function for automatically naming cluster sub-alignments. Lambda can't be used as it is passed as an
     argument to multi-processing. """

    return f"temp_sub_aln_{rep}"


def hc_merged_aln(
        full_seq_fa,
        rep_aln,
        out_file,
        cluster_files=None,
        hc_name=None,
        remove_sub_alns=True,
        custom_ext_clust_map=None,
        lower_rep_idx=None,
        upper_rep_idx=None,
        select_seqs=None,
        max_processes=6,
        nice=None
):
    """ Given an alignment of cluster representatives from a hierarchical clustering, produce alignments of all clusters
     independently and merge them as implemented in aln.merge_cluster_alns. Reps of the lower clustering level index
     should match the reps in the rep_aln. The upper clustering indicates the highest clustering level for which reps
      should be used in the merged alignment (i.e. some level of redundancy reduction can still be maintained). See
    the expand_cluster_reps function for further details. If the rep_aln is to be filtered, specific sequences names can
    be provided. """

    if not (cluster_files or hc_name):
        raise RuntimeError("File names must be provided either explicitly as a list or implicitly as a hc name.")

    if isinstance(rep_aln, str):  # Alignment needs to be read in from file
        rep_aln = AlignIO.read(rep_aln, "fasta")

    if select_seqs:  # Filter the alignment if required
        filt_seqs = [record for record in rep_aln if record.name in select_seqs]
        rep_aln = AlignIO.MultipleSeqAlignment(filt_seqs)

    reps = [seq.name for seq in rep_aln]

    if not cluster_files:  # Find all files following naming conventions with the clustering name and map to %
        file_list = [(file, int(file.split('_cov')[0].split('_clust')[1])) for file in os.listdir()
                     if (hc_name in file and file.endswith('.cluster'))]
        file_list.sort(key=lambda x:x[1], reverse=True)
        cluster_files = [file[0] for file in file_list]

    if custom_ext_clust_map:  # If re-allocations have been made relative to original clustering
        exp_clust_map = custom_ext_clust_map
    else:
        hc_maps = [extract_cluster_ext(file, clustering_map=True) for file in cluster_files]
        # Get expanded cluster map required for building sub-alignments
        exp_clust_map = expand_cluster_reps(hc_maps, lower=lower_rep_idx, upper=upper_rep_idx, lower_reps=reps)

    if remove_sub_alns:  # Cluster alignments will create many temp files
        clust_aln_names = [f"temp_sub_aln_{rep}.aln" for rep in reps]

        naming = default_sub_aln_naming

    else:  # Cluster alignment function will apply default names
        clust_aln_names = [f"{full_seq_fa.split('.')[0]}_sub_{rep}.aln" for rep in reps]
        naming = None

    # Generate all sub-alignments
    cluster_sub_alns(full_seq_fa, exp_clust_map, naming=naming, max_processes=max_processes, nice=nice)

    # Create merged alignment
    aln.merge_cluster_alns(rep_aln, clust_aln_names, out_file=out_file)

    if remove_sub_alns:
        for file in glob.glob(f"temp_sub_aln_*"):
            os.remove(file)


def dash_seed_cluster_aln(
        full_seq_fa,
        merged_aln_name,
        cluster_files=None,
        hc_name=None,
        rep_aln_name=None,
        remove_rep_aln=False,
        remove_sub_alns=True,
        lower_rep_idx=None,
        upper_rep_idx=None,
        include_seqs=None,
        max_processes=6,
        nice=None,
        quiet=False
):
    """ For a given hierarchical clustering, specified either by a list of cluster files, or the name of the .hc file,
     align the representatives of a lower level clustering using MAFFT-DASH (with local-pair). This is then used as the
      representative alignment for hc_merged_aln, with the upper level of the clustering for which sequences are to be
    included in the cluster sub-alignments also specified here. """

    if not cluster_files:  # Find all files following naming conventions with the clustering name and map to %
        file_list = [(file, int(file.split('_cov')[0].split('_clust')[1])) for file in os.listdir()
                     if (hc_name in file and file.endswith('.cluster'))]
        file_list.sort(key=lambda x: x[1], reverse=True)
        cluster_files = [file[0] for file in file_list]
    else:
        cluster_files.sort(reverse=True, key=lambda x:int(x.split('_clust')[1].split("_cov")[0]))

    # Use lowest possible level clustering (last in list) if none provided
    lower_rep_idx = lower_rep_idx if lower_rep_idx else -1

    cluster_f_name = cluster_files[lower_rep_idx]

    if not rep_aln_name:
        lower_cov = cluster_f_name.split("cov")[1].split("_")[0]
        rep_aln_name = cluster_f_name.split("_cov")[0] + f"_cov{lower_cov}.aln"

    rep_fa_name = rep_aln_name.split(".aln")[0] + ".fa"

    if include_seqs:   # Get full expanded cluster map - we may need to switch cluster reps should a lower level rep be excluded

        hc_maps = [extract_cluster_ext(file, clustering_map=True) for file in cluster_files]

        # # TODO: For debugging - REMOVE
        # with open("test_log.txt", 'w') as log:
        #     log.write("\n".join(str(hc_map) for hc_map in hc_maps))

        exp_clust_map = expand_cluster_reps(hc_maps, lower=lower_rep_idx, upper=upper_rep_idx)

        # # TODO: For debugging - REMOVE
        # with open("test_log.txt", 'a') as log:
        #     log.write(str(exp_clust_map))

        revised_exp_clust_map = {}
        for og_rep, og_members in exp_clust_map.items():
            new_members = [member for member in og_members if member in include_seqs]
            if new_members:
                if og_rep in new_members:
                    new_rep = og_rep
                else:  # Take the longest sequence if we need to replace rep
                    member_lens = [(seq.name, len(seq)) for seq in SeqIO.parse(full_seq_fa, "fasta") if seq.name in new_members]
                    new_rep = max(member_lens, key=lambda x:x[1])[0]
                revised_exp_clust_map[new_rep] = new_members

        # TODO: For debugging - REMOVE
        # with open("test_log.txt", 'w') as log:
        #     log.write(str(revised_exp_clust_map))

        # Get updated rep sequences
        rep_seqs = [seq for seq in SeqIO.parse(full_seq_fa, "fasta") if seq.name in revised_exp_clust_map.keys()]
        SeqIO.write(rep_seqs, rep_fa_name, "fasta")

    else:  # Using all original clusters
        reps_from_cluster_file(cluster_f_name, rep_fa_name)
        revised_exp_clust_map = None


    # # Extract rep sequences from desired lowest clustering level
    # if filt_seqs:  # These sequences are to be excluded
    #     reps_from_cluster_file(cluster_f_name, "temp_reps_abc.fa")
    #     unfilt_aln = list(SeqIO.parse("temp_reps_abc.fa", "fasta"))
    #     filt_aln = [seq for seq in unfilt_aln if seq.name not in filt_seqs]
    #     SeqIO.write(filt_aln, rep_fa_name, "fasta")
    #     os.remove("temp_reps_abc.fa")
    # if include_seqs:  # ONLY these sequences should be included
    #     reps_from_cluster_file(cluster_f_name, "temp_reps_abc.fa")
    #     unfilt_aln = list(SeqIO.parse("temp_reps_abc.fa", "fasta"))
    #     # print(f"UNFILT ALN LEN: {len(unfilt_aln)}")
    #     filt_aln = [seq for seq in unfilt_aln if seq.name in include_seqs]
    #     # print(f"FILT ALN LEN: {len(filt_aln)}")
    #     SeqIO.write(filt_aln, rep_fa_name, "fasta")
    #     os.remove("temp_reps_abc.fa")
    # else:
    #     reps_from_cluster_file(cluster_f_name, rep_fa_name)

    nice = str(int(nice)) if nice else '0'

    # Perform MAFFT-DASH alignment for these reps

    args = [
        "nice",
        "-n", nice,
        "mafft",
        "--dash",
        "--reorder",
        "--thread", str(max_processes),
        "--threadit", str(max([max_processes,30])),
        "--originalseqonly"
    ]

    if quiet:
        args.append("--quiet")

    args.append(rep_fa_name)

    with open(rep_aln_name, 'w') as out_file:
        subprocess.run(args, stdout=out_file)

    # os.system(f"nice -n {nice} mafft --dash --reorder --localpair --thread {max_processes} --threadit "
    #           f"{max([max_processes,30])} --originalseqonly {'--quiet' if quiet else ''}{rep_fa_name} > {rep_aln_name}")

    # Construct merged alignment from sub-alignments of clusters
    hc_merged_aln(
        full_seq_fa,
        rep_aln_name,
        merged_aln_name,
        cluster_files=cluster_files,
        custom_ext_clust_map=revised_exp_clust_map,
        remove_sub_alns=remove_sub_alns,
        lower_rep_idx=lower_rep_idx,
        upper_rep_idx=upper_rep_idx,
        max_processes=max_processes,
        nice=nice
    )
