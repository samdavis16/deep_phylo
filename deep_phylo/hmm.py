""" Profile HMM operations and workflows, mostly utilising HMMER3 and other modules in this package. Many functions are
specific to construction and searching of pHMMs from monophyletic groups in phylogenetic trees. """

import concurrent.futures
import glob
import multiprocessing as mp
import os
import random
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from pathlib import Path

import ete3
import pandas as pd
import treeswift
from Bio import SeqIO

from . import aln
from . import annots
from . import file_util
from . import tree
from . import cluster


def build_hmm(aln_file, out_file, quiet=True):
    """ Construct a profile hidden markov model over all extant sequences in a tree using corresponding multiple
    sequence alignment. """

    hmm_name = out_file.split('.')[0]

    args = ["hmmbuild", "-n", hmm_name, out_file, aln_file]

    subprocess.run(args,
                   check=True,
                   stdout=subprocess.DEVNULL if quiet else None)


def count_profiles(hmm_file):
    """ Count the number of individual profile HMMs in a given .hmm file. """

    # Check if files have already been concatenated
    with open(hmm_file) as in_file:
        cnt = 0
        for line in in_file:
            if line.strip().startswith("//"):
                cnt += 1
    return cnt


################## TODO: Sort these out (just need downstream to properly call the new function)

def build_clade_hmm_old(tree_file, aln_file, in_nodes=None, out_node=None, excl_nodes=None, realign=True, clade_name=None,
                        out_path='', stored_aln=False):
    """
    TODO: More control over file naming so this can be called from more locations easily
    Build a profile HMM over all sequences under the common ancestor of extant labels in in_nodes.
    A single outgroup extant ID can initially be used root the input tree. If no out_node is provided,
    assumes that the clade of interest is already monophyletic. If less than two in_nodes are provided,
    construct a pHMM for the entire tree. If exclusion nodes are provided, remove these sequences before HMM building.
    If a subtree is used, sequences can be optionally realigned. If re-alignnment is already available, set stored_aln
    to True.
    TODO: Currently not handling case where clade name is not provided"""

    t = tree.load_tree(tree_file) if not isinstance(tree_file, ete3.Tree) else tree_file

    # Define name for hmm file if not provided
    if not clade_name:
        if isinstance(in_nodes, list):
            clade_name = aln_file.split('.')[0] + '_sub_' + '_'.join(in_nodes)
        else:  # Specific clade not specified - will build pHMM on entire alignment
            clade_name = aln_file.split('.')[0] + '_full'


    if len(in_nodes) < 2 or not isinstance(in_nodes, list):  # Target clade not implied, therefore use whole tree
        build_hmm(aln_file, clade_name + '.hmm')

    else:  # Need to make a subtree

        if out_node:  # Initially root the tree on out_node
            t.set_outgroup(t & out_node)

        t = tree.get_subtree(t, clade_bounds=in_nodes)
        leaf_names = [node.name for node in t.iter_descendants() if node.is_leaf()]
        # Remove any nodes which are to be excluded from profile
        if excl_nodes:
            for node in excl_nodes:
                leaf_names.remove(node)

        # Extract relevant seqs from original aln
        sub_aln_file = out_path + aln_file.split('.')[0].split('/')[-1] + '_sub.aln' if not clade_name else clade_name + '.aln'
        file_util.extract_fasta(aln_file, sub_aln_file, leaf_names)


        if realign:  # structural re-aligment of sequences in subtree

            ### May want to store structural alignments for clades if performing analyses with the same clades ###
            ### multiple times. ###
            # aln_name = sub_aln_file.split('.aln')[0]
            # out_file = '../../struct_alns/' + aln_name + '_re.aln'

            out_file = sub_aln_file.split('.aln')[0] + '_re.aln'
            if stored_aln:  # Avoids structural alignment bottleneck if several analyses are performed with same clades
                pass
            else:
                aln.struct_aln(sub_aln_file, out_file=out_file)

        else:  # Just using the extracted sequences w/out realignment
            out_file = sub_aln_file

        # Construct profile HMM from relevant aln
        build_hmm(out_file, clade_name + '.hmm')


def build_clade_hmm(
        full_aln,
        sub_aln,
        out_hmm,
        full_tree,
        target_clade_node=None,
        target_clade_bounds=None,
        out_leaf=None,
        reroot=False,
        include=None,
        exclude=None,
        realign=False,
        structural=False,
        quiet=True,
        missing_error=True
):
    """ Build a profile HMM over all sequences in a target clade, specified either by a set of extant nodes on either
    side of the LCA, or an internal node label. The root can optionally be re-rooted to ensure the target clade is
    monophyletic. If only a subset of clade sequences should be used to construct the HMM, labels can be provided for
    either explicit inclusion or exclusion. The sub-sequences can be optionally realigned. """

    if not include:
        include = []
    if not exclude:
        exclude = []

    aln.clade_sub_aln(
        full_aln,
        sub_aln,
        full_tree,
        target_clade_node=target_clade_node,
        target_clade_bounds=target_clade_bounds,
        out_leaf=out_leaf,
        reroot=reroot,
        include=include,
        exclude=exclude,
        realign=realign,
        structural=structural,
        quiet=quiet,
        missing_error=missing_error)

    build_hmm(sub_aln, out_hmm)


# TODO: Below is a re-written version of build_clade_hmm_old - code calling build_clade_hmm_old should be redesigned to call below function instead
# TODO: OR, if it's decided that this functionality is unecessary for workflows, remove and have functions calling build_clade_hmm_old do the clade extraction followed by hmmbuild

# def build_clade_hmm_NEW(tree_file, aln_file, in_nodes=None, out_node=None, excl_nodes=None, realign=True, clade_name=None, stored_aln=False):
#     """ TODO: More control over file naming so this can be called from more locations easily
#     Build a profile HMM over all sequences under the common ancestor of extant labels in in_nodes.
#     A single outgroup extant ID can initially be used root the input tree. If no out_node is provided,
#     assumes that the clade of interest is already monophyletic. If less than two in_nodes are provided,
#     construct a pHMM for the entire tree. If exclusion nodes are provided, remove these sequences before HMM building.
#     If a subtree is used, sequences can be optionally realigned. If re-alignnment is already available, set stored_aln
#     to True.
#     TODO: Currently not handling case where clade name is not provided"""
#
#     t = tree.load_tree(tree_file) if not isinstance(tree_file, ete3.Tree) else tree_file
#
#     # Define name for hmm file if not provided
#     if not clade_name:
#         if isinstance(in_nodes, list):
#             clade_name = aln_file.split('.')[0] + '_sub_' + '_'.join(in_nodes)
#         else:  # Specific clade not specified - will build pHMM on entire alignment
#             clade_name = aln_file.split('.')[0] + '_full'
#
#
#     if len(in_nodes) < 2 or not isinstance(in_nodes, list):  # Target clade not implied, therefore use whole tree
#         build_hmm(aln_file, clade_name + '.hmm', out_path=out_path)
#
#     else:  # Need to make a subtree
#
#         if out_node:  # Initially root the tree on out_node
#             t.set_outgroup(t & out_node)
#
#         t = tree.get_subtree(t, ext_nodes=in_nodes)
#         leaf_names = [node.name for node in t.iter_descendants() if node.is_leaf()]
#         # Remove any nodes which are to be excluded from profile
#         if excl_nodes:
#             for node in excl_nodes:
#                 leaf_names.remove(node)
#
#         # Extract relevant seqs from original aln
#         sub_aln_file = out_path + aln_file.split('.')[0].split('/')[-1] + '_sub.aln' if not clade_name else clade_name + '.aln'
#         parsers.extract_fasta(aln_file, sub_aln_file, leaf_names, id_format=None)
#
#
#         if realign:  # structural re-aligment of sequences in subtree
#
#             ### May want to store structural alignments for clades if performing analyses with the same clades ###
#             ### multiple times. ###
#             # aln_name = sub_aln_file.split('.aln')[0]
#             # out_file = '../../struct_alns/' + aln_name + '_re.aln'
#
#             out_file = sub_aln_file.split('.aln')[0] + '_re.aln'
#             if stored_aln:  # Avoids structural alignment bottleneck if several analyses are performed with same clades
#                 pass
#             else:
#                 aln.struct_aln(sub_aln_file, out_file=out_file)
#
#         else:  # Just using the extracted sequences w/out realignment
#             out_file = sub_aln_file
#
#         # Construct profile HMM from relevant aln
#         build_hmm(out_file, clade_name + '.hmm', out_path=out_path)



def hmm_press(db_name, hmm_names, out_path=None):
    # TODO: Moving towards allowing greater flexibility with where temp files are stored
    # TODO: hmmsearch and hmmscan seem to have similar run times for equivalent profile/sequence sets. Maybe just use search in future
    """ Press a collection of HMMs to an indexed db to be used with hmmscan.
    Assumes HMM files are stored in'{hmm_name}/{hmm_name}.hmm' relative to current directory """

    out_path = 'out' if not out_path else out_path

    # Concatenate all relevant HMM profiles
    hmm_files = ""
    for name in hmm_names:
        hmm_files += f"{out_path}/{name}/{name}.hmm "
    db_path = f"{out_path}/{db_name}.hmm"
    os.system(f"cat {hmm_files}> {db_path}")

    # Press concatenated HMMs to db
    os.system(f"hmmpress {db_path}")


def hmm_concat(in_files, out_file):
    """ Concatenate multiple HMM files for use by hmmsearch. Note: this process is not the same as hmmpress, which
    prepares various files used to run hmmscan. """
    hmm_files = ""
    for hmm in in_files:
        hmm_files += f"{hmm} "
    os.system(f"cat {hmm_files}> {out_file}")


def parse_tblout(results_file, id_format=None, best_domain=True):
    """ Parse tbl-formatted hmmsearch output into df (for multiple profiles) or dict (for a single profile). The score
    of the best domain is parsed by default, however the score for the whole sequence can be used if desired. """
    with open(results_file, 'r') as file:

        # Get non-redundant sets of pHMMs and sequences and construct df
        data = [line for line in file.readlines() if not line.startswith('#')]
        profiles = set([line.split()[2] for line in data])

        # Convert sequence ID as required (only uniref50 currently supported)
        if id_format == 'uniref50':
            seqs = set([line.split()[0].split('UniRef50_')[1] for line in data])
        else:   # Just use ID as is
            seqs = set([line.split()[0] for line in data])

        # Df for multiple profiles, or else dict if searching with a single one
        scores = pd.DataFrame(index=list(seqs), columns=list(profiles), dtype=object) if len(profiles) > 1 else {}

        # Collect scores
        for line in data:
            seq_id = line.split()[0]

            # Adjust seq_id for id format of database
            if id_format == 'uniref50':
                seq_id = seq_id.split('UniRef50_')[1]

            if len(profiles) > 1:   # scores is a dataframe for many profiles
                scores[line.split()[2]][seq_id] = float(line.split()[8 if best_domain else 5])
            else:  # scores is a dict for one profile
                scores[seq_id] = float(line.split()[8 if best_domain else 5])

        return scores


def parse_dombtblout(results_file, threshold=0, envelope=False):
    """ Parse domtbl-formatted hmmsearch output into a dictionary {seq_id : (score, (boundary_start, boundary_end))} for
    single profile, or a nested dictionary
    {seq_id : {profile1 : (score, (boundary_start, boundary_end)), profile2 : ( .... } ... } for multiple profiles.
    Boundaries can be the domain (default) or envelope boundaries. A domain score threshold can be applied if desired.

    TODO: Only currently reports best single hit
    """

    with open(results_file) as file:

        data = [line for line in file.readlines() if not line.startswith('#')]

        # Get set of pHMMs for which search results are present
        profiles = set([line.split()[3] for line in data])

        # Get set of sequence IDs for which search results are present
        seqs = set([line.split()[0] for line in data])

        # Initialise dictionary
        if len(profiles) > 1:   # Multiple profiles -> nested dict
            result_dict = {seq : {profile : None for profile in profiles} for seq in seqs}
        else:  # Just one profile
            result_dict = {seq : None for seq in seqs}

        # Get scores and boundaries for each hit
        for line in data[::-1]:   # If multiple domain hits for seq, best will be the last encountered
            fields = line.split()
            if len(profiles) > 1 and float(fields[13]) > threshold:  # {seq : {profile : (score, (boundaries)) }}
                result_dict[fields[0]][fields[3]] = (float(fields[13]), ((int(fields[19]), int(fields[20]))
                                                                         if envelope else (int(fields[17]), int(fields[18]))))
            elif float(fields[13]) > threshold:
                result_dict[fields[0]] = (float(fields[13]), ((int(fields[19]), int(fields[20])) if envelope else
                                                              (int(fields[17]), int(fields[18]))))

    return result_dict


def parse_domtblout_full(results_file):
    """ Parse domtbl-formatted hmmsearch output including profile and hit alignment coordinates for all domain hits for
     each profile-query pair. Returns a nested dictionary:
     {profile1 : {seq1 : [(domain_score_hit1, (hmm_coords_hit1), (ali_coords_hit1), (env_coords_hit1)), ...], ...}, ...}
     """
    with open(results_file) as file:

        data = [line for line in file.readlines() if not line.startswith('#')]

        # Get set of pHMMs for which search results are present
        profiles = set([line.split()[3] for line in data])

        # Get set of sequence IDs for which search results are present
        seqs = set([line.split()[0] for line in data])

        result_dict = {profile : {seq : [] for seq in seqs} for profile in profiles}

        for hit in data:
            fields = hit.split()
            # (score, (hmm_start, hmm_end), (ali_start, ali_end), (env_start, env_end))
            hit_info = (float(fields[13]),
                        (int(fields[15]), int(fields[16])),
                        (int(fields[17]), int(fields[18])),
                        (int(fields[19]), int(fields[20]))
                        )
            result_dict[fields[3]][fields[0]].append(hit_info)

    return result_dict


def chain_profile_coverage(chain):
    """ For a chain of hits of the form [score, (hmm_start, hmm_end), ... + optionally target coordinates], return a
    list of hmm consensus positions covered in at least one hit. """

    positions = set()

    for hit in chain:
        positions.update(range(hit[1][0], hit[1][1]+1))

    positions = list(positions)
    positions.sort()

    return positions


def get_profile_chains(
        dom_results,
        region_type="env",
        overlap_tol_prop=0.1,
        hmm_files=None
):
    """ For a set of domain hits for a given target sequence against a profile HMM, group hits which represent
    individual instances of paths through both target sequence coordinate space and profile HMM coordinate space. The implementation
    scans domain hits from left to right by the target start coordinate (ali or env as specified) and extends the
    current chain if hmm coordinates suggest progression through the same domain instance, allowing for insertions
    (discontinuity in target sequence coordinates) and deletions (discontinuity in hmm coordinates) relative to the hmm.
     A new chain is initialised if hmm coordinates of an encountered hit are to the left of the previous hit, with some
    specified level of overlap tolerance for fuzzy boundaries. """

    if not hmm_files:  # hmm file names should be consistent with what's in results file
        hmm_files = {profile : f"{profile}.hmm" for profile in dom_results}

    # Get lengths of profiles in order to define overlap tolerance as number of positions
    hmm_lens = {}
    for profile, file_name in hmm_files.items():
        with open(file_name) as file:
            for line in file:
                if line.startswith("LENG"):
                    hmm_lens[profile] = int(line.split("LENG")[1].strip())
                    break

    chain_dict = {profile : {} for profile in dom_results}

    for profile in dom_results:

        tol = overlap_tol_prop * hmm_lens[profile]  # Get profile-specific hmm overlap tolerance

        for target_seq, hits in dom_results[profile].items():

            if not hits:  # If not hits for this profile-target pair, hits will be empty
                continue

            hits_ordered = hits
            hits_ordered.sort(key=lambda x:x[3 if region_type == "env" else 2][0])

            chains = []  # Store chains once finalised
            active_chain = [hits_ordered[0]]  # Store the current active chain
            prospective_chain = []  # An apparent new chain to be confirmed

            for hit in hits_ordered[1:]:

                if prospective_chain:
                    # Check if this is a possible continuation to the hit of the prospective chain
                    if hit[1][0] >= prospective_chain[0][1][1] - tol:
                        prospective_chain.append(hit)

                # Now check for continuation of main active chain
                if hit[1][0] >= active_chain[-1][1][1] - tol:
                    active_chain.append(hit)

                    if len(prospective_chain) > 1:  # Need to decide to maintain either active or prospective chain

                        if len(chain_profile_coverage(prospective_chain)) > len(chain_profile_coverage(active_chain)):
                            # Adopt prospective chain as new active chain
                            chains.append(active_chain)
                            active_chain = prospective_chain
                            prospective_chain = []

                        else:
                            chains.append(prospective_chain)
                            prospective_chain = []

                    elif len(prospective_chain) == 1:  # Prospective chain wasn't extended
                        chains.append(prospective_chain)
                        prospective_chain = []

                else:

                    if len(prospective_chain) > 0:  # This must be second hit in a row not to extend active chain
                        chains.append(active_chain)
                        active_chain = prospective_chain
                        if len(prospective_chain) == 1:  # Current hit didn't extend prospective chain
                            prospective_chain = [hit]
                        else:  # Current hit already in prospective chain - start a new one
                            prospective_chain = []

                    else:  # No current prospective chain - maintain main chain for now
                        prospective_chain = [hit]

            if active_chain:
                chains.append(active_chain)
            if prospective_chain:
                chains.append(prospective_chain)
            chain_dict[profile][target_seq] = chains

    return chain_dict


def chain_coords_per_target(chains, region_type_idx):
    """ Function to be called per-thread when multi-threading target sequences and pre-complied chains for a given
    profile, for obtaining per-chain coordinates and spans. """

    coords = set()
    spans = []

    for chain in chains:

        # Coords - These are for subsequence extraction and hmmalign - don't need to separate by individual hit
        for hit in chain:
            coords.add(hit[region_type_idx])

        # Spans per chain
        chain_start = min([hit[region_type_idx][0] for hit in chain])
        chain_end = max([hit[region_type_idx][1] for hit in chain])
        spans.append((chain_start, chain_end))

    return coords, spans


def hmmalign_subseqs_per_profile(
        profile_name,
        subseqs,
        temp_tag,
        hmm_file
):
    """ Function to be called per process (per-profile) when multiprocessing for target subsequence extraction and
    profile alignment with hmmalign. """

    # Write subsequence file
    SeqIO.write(subseqs, f"temp_{temp_tag}_{profile_name}_subseqs.fa", "fasta")

    # Run hmmalign for this profile against all relevant target sub-sequences
    args = ["hmmalign", "--outformat", "A2M", hmm_file, f"temp_{temp_tag}_{profile_name}_subseqs.fa"]
    with open(f"temp_{temp_tag}_{profile_name}_alns.aln", 'w') as out_file:
        subprocess.run(args, stdout=out_file)


def pcpt_init(
        profile_name,
        profile_alns,
        hmm_params,
        seq_lens,
        profile_alpha,
        region_type_idx,
        target_spans
):
    """ Instantiate global variables available to each weighted chain processing process. """

    global _profile_name
    _profile_name = profile_name

    global _profile_alns
    _profile_alns = profile_alns

    global _hmm_params
    _hmm_params = hmm_params

    global _seq_lens
    _seq_lens = seq_lens

    global _profile_alpha
    _profile_alpha = profile_alpha

    global _region_type_idx
    _region_type_idx = region_type_idx

    global _target_spans
    _target_spans = target_spans


def process_chains_per_target(target_name, chains):
    """ Function to be called per-process when multi-processing over target sequences and pre-compiled chains for a
    given profile, for processing weighted chain scores. """

    weighted_scores = []

    for k in range(len(chains)):

        chain_alns = []

        for hit in chains[k]:
            hit_target_coords = hit[_region_type_idx]
            subseq_name = f"{target_name}_{hit_target_coords[0]}_{hit_target_coords[1]}"
            aln_str = str([seq for seq in _profile_alns if seq.name == subseq_name][0].seq)
            chain_alns.append(aln_str)

        # Count per-hmm-position and per-target-position depths
        # calculated below on basis of re-computed alignments
        profile_depths = [0 for _ in range(_hmm_params[_profile_name][0])]
        target_depths = [0 for _ in range(_seq_lens[target_name])]

        # Approximate position-wise score contributions using match state emission log-odds

        per_hit_pos_conts = []  # To store all per-hit, per-aligned-position score contributions

        for i in range(len(chains[k])):  # Per hit in chain

            # Get raw (non-normalised) match state emission log-odds for each aligned position
            raw_contributions = []  # Store contribution at each position as (raw_sub_score, hmm_pos, target_pos)

            # Initialise coordinates at beginning of hit
            hmm_pos = 1  # Alignments span over full profile coordinates regardless of where hit sits
            target_pos = chains[k][i][_region_type_idx][0]

            try: # TODO: For debugging - REMOVE try/except

                for j in range(len(chain_alns[i])):

                    # Scores in hmmer file are s = -ln(p) for both background frequencies and match state probs
                    # Trivial to show that log-odds is = transformed_background_score - transformed_match_score
                    char = chain_alns[i][j]
                    if char in _profile_alpha:  # Either a match state or hmm gap
                        if char != "-":  # Only attribute weight to match states, assume gaps have 0 weight
                            profile_depths[hmm_pos - 1] += 1
                            target_depths[target_pos - 1] += 1
                            sub_score = _hmm_params[_profile_name][1][char] - _hmm_params[_profile_name][2][hmm_pos][char]
                            raw_contributions.append((sub_score, hmm_pos, target_pos))
                            target_pos += 1  # Only progress in target coordinate space if it's a match state
                        hmm_pos += 1  # Gaps still progress profile coordinate
                    else:
                        target_pos += 1

                # Need to normalise to ensure decomposition sums to original domain score
                try:
                    norm_factor = sum([pos_data[0] for pos_data in raw_contributions]) / chains[k][i][0]
                    log_odds_norm = [(score / norm_factor, hmm_pos, target_pos)
                                     for score, hmm_pos, target_pos in raw_contributions]
                except ZeroDivisionError:  # Edge case where a domain hit score is 0.0
                    log_odds_norm = [(0, hmm_pos, target_pos) for score, hmm_pos, target_pos in raw_contributions]
                per_hit_pos_conts.extend(log_odds_norm)

            except IndexError:
                with open(f"failure.log", "a") as fail_log:
                    fail_log.write(_profile_name+"\n\n")
                    fail_log.write(target_name+"\n\n")


        # Sum contributions over hmm positions weighted by position depth
        tot_weighted_hmm = 0
        for hmm_pos in range(1, _hmm_params[_profile_name][0] + 1):
            for hit_pos in [hit_pos for hit_pos in per_hit_pos_conts if hit_pos[1] == hmm_pos]:
                tot_weighted_hmm += hit_pos[0] / profile_depths[hmm_pos - 1]

        # Do the same over target positions
        tot_weighted_target = 0
        for target_pos in range(1, _seq_lens[target_name] + 1):
            for hit_pos in [hit_pos for hit_pos in per_hit_pos_conts if hit_pos[2] == target_pos]:
                tot_weighted_target += hit_pos[0] / target_depths[target_pos - 1]

        # Get total weighted score accounting for both hmm and target overlaps
        # Also record the overall start and end coordinates (span) of the chain
        tot_weighted = 0.5 * (tot_weighted_hmm + tot_weighted_target)
        weighted_scores.append((tot_weighted, _target_spans[_profile_name][target_name][k]))

    return weighted_scores


def get_weighted_chain_scores(
        chain_dict,
        full_seq_fa,
        region_type="env",
        temp_tag=None,
        hmm_files=None,
        max_workers=6
):
    """ For a dictionary of profile-target chains of the form returned by get_profile_chains, calculate a cumulative,
    weighted domain score accounting for overlaps in both the profile and target coordinates. Per-position contributions
    to a domain bit-score are estimated by ratios of match state emission scores across a computed alignment between the
    profile and target interval (either ali or env). """

    if not temp_tag:
        temp_tag = random.randint(int(1e7), int(1e8)-1)  # Avoid identically named temp files if threading

    if not hmm_files:  # hmm file names should be consistent with what's in results file
        hmm_files = {profile : f"{profile}.hmm" for profile in chain_dict}

    region_type_idx = 3 if region_type == "env" else 2

    # Get profile length and equilibrium and match state probabilities for use in position-specific score attribution
    # { profile : (profile_len, {eq_probs}, {pos : {ms_probs} ) }
    hmm_params = {profile : [None for _ in range(3)] for profile in chain_dict}
    for profile, file_name in hmm_files.items():
        with open(file_name) as file:
            for line in file:
                if line.strip().startswith("LENG"):
                    this_len = int(line.split("LENG")[1].strip())
                    hmm_params[profile][0] = this_len
                    hmm_params[profile][2] = {i+1 : None for i in range(this_len)}
                elif line.strip().startswith("HMM "):
                    alphabet = [char for char in line.strip().split()[1:]]
                elif line.strip().startswith("COMPO "):
                    data = line.strip().split()[1:]
                    eq_probs = {}
                    for i in range(len(alphabet)):
                        eq_probs[alphabet[i]] = float(data[i])
                    hmm_params[profile][1] = eq_probs
                else:
                    try:  # Check if this is a match state emission prob line
                        col = line.strip().split()[0]
                        if float(col) == int(col):  # Whole number suggests that this is indeed a column number
                            data = line.strip().split()[1:21]
                            ms_probs = {alphabet[i] : float(data[i]) for i in range(len(alphabet))}
                            hmm_params[profile][2][int(col)] = ms_probs
                    except ValueError:
                        continue

    # For each profile, collect sub-target-sequences for which hmmalign will later need to be run
    # Also for each profile-target-chain, pair, get the complete span of the chain for later subsequence extraction
    target_coords = {}
    target_spans = {}
    for profile, target_chains in chain_dict.items():

        target_coords[profile] = {}
        target_spans[profile] = {}

        with ThreadPoolExecutor(max_workers=max_workers) as executor:

            future_store = {executor.submit(chain_coords_per_target, chains, region_type_idx) : target
                              for target, chains in target_chains.items()}

            for future in as_completed(future_store):
                target = future_store[future]
                coords, spans = future.result()
                target_coords[profile][target] = coords
                target_spans[profile][target] = spans

        # TODO: Now re-implemented above with threading - DELETE BELOW ONCE TESTED
        # for target, chains in target_chains.items():
        #     for chain in chains:
        #         for hit in chain:
        #
        #             # Sub-sequence coordinates
        #             try:
        #                 target_coords[profile][target].add(hit[region_type_idx])
        #             except KeyError:
        #                 target_coords[profile][target] = {hit[region_type_idx]}
        #
        #         # Chain spans
        #         chain_start = min([hit[region_type_idx][0] for hit in chain])
        #         chain_end = max([hit[region_type_idx][1] for hit in chain])
        #         try:
        #             target_spans[profile][target].append((chain_start, chain_end))
        #         except KeyError:
        #             target_spans[profile][target] = [(chain_start, chain_end)]

    # Extract all required sub-sequences and record sequence lengths
    seq_lens = {}
    subseqs = {profile : [] for profile in target_coords}
    for seq in SeqIO.parse(full_seq_fa, "fasta"):
        seq_lens[seq.name] = len(seq)
        for profile in target_coords:
            if seq.name in target_coords[profile]:
                for coords in target_coords[profile][seq.name]:
                    subseq = seq[coords[0] - 1:coords[1]]
                    subseq.id = f"{seq.name}_{coords[0]}_{coords[1]}"
                    subseqs[profile].append(subseq)

    # for profile in target_coords:
    #     SeqIO.write(subseqs[profile], f"temp_{temp_tag}_{profile}_subseqs.fa", "fasta")

    # Write subseqs and run hmmalign per-profile in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:

        future_store = [executor.submit(hmmalign_subseqs_per_profile,
                                        profile,
                                        subseqs[profile],
                                        temp_tag,
                                        hmm_files[profile])
                        for profile in hmm_files]

        for future in as_completed(future_store):
            future.result()

    # Full results for return
    all_weighted_scores = {}

    # Compute weighted scores for each chain using per-position score contribution approximations
    # (based on per-pos emission state probs and computed hmm alignments)
    for profile, target_chains in chain_dict.items():

        all_weighted_scores[profile] = {}
        profile_alpha = list(hmm_params[profile][1].keys()) + ["-"]

        # # Run hmmalign for this profile against all relevant target sub-sequences
        # args = ["hmmalign", "--outformat", "A2M", hmm_files[profile], f"temp_{temp_tag}_{profile}_subseqs.fa"]
        # with open(f"temp_{temp_tag}_{profile}_alns.aln", 'w') as out_file:
        #     subprocess.run(args, stdout=out_file)

        profile_alns = list(SeqIO.parse(f"temp_{temp_tag}_{profile}_alns.aln", "fasta"))

        with ProcessPoolExecutor(
                max_workers=max_workers,
                initializer=pcpt_init,
                initargs=(
                        profile,
                        profile_alns,
                        hmm_params,
                        seq_lens,
                        profile_alpha,
                        region_type_idx,
                        target_spans)
                             ) as executor:

            future_target = {executor.submit(process_chains_per_target, target, chains) : target
                             for target, chains in target_chains.items()}

            for future in as_completed(future_target):
                all_weighted_scores[profile][future_target[future]] = future.result()


        # for target, chains in target_chains.items():
        #
        #     all_weighted_scores[profile][target] = []  # (weighted score, (chain_start, chain_end)) for each chain
        #
        #     for k in range(len(chains)):
        #
        #         chain_alns = []
        #         for hit in chains[k]:
        #             hit_target_coords = hit[region_type_idx]
        #             subseq_name = f"{target}_{hit_target_coords[0]}_{hit_target_coords[1]}"
        #             aln_str = str([seq for seq in profile_alns if seq.name == subseq_name][0].seq)
        #             chain_alns.append(aln_str)
        #
        #         # Count per-hmm-position and per-target-position depths
        #         # calculated below on basis of re-computed alignments
        #         profile_depths = [0 for _ in range(hmm_params[profile][0])]
        #         target_depths = [0 for _ in range(seq_lens[target])]
        #
        #         # Approximate position-wise score contributions using match state emission log-odds
        #
        #         per_hit_pos_conts = []  # To store all per-hit, per-aligned-position score contributions
        #
        #         for i in range(len(chains[k])):   # Per hit in chain
        #
        #             # Get raw (non-normalised) match state emission log-odds for each aligned position
        #             raw_contributions = []  # Store contribution at each position as (raw_sub_score, hmm_pos, target_pos)
        #
        #             # Initialise coordinates at beginning of hit
        #             hmm_pos = 1  # Alignments span over full profile coordinates regardless of where hit sits
        #             target_pos = chains[k][i][region_type_idx][0]
        #
        #             for j in range(len(chain_alns[i])):
        #
        #                 # Scores in hmmer file are s = -ln(p) for both background frequencies and match state probs
        #                 # Trivial to show that log-odds is = transformed_background_score - transformed_match_score
        #                 char = chain_alns[i][j]
        #                 if char in profile_alpha:  # Either a match state or hmm gap
        #                     if char != "-":  # Only attribute weight to match states, assume gaps have 0 weight
        #                         profile_depths[hmm_pos-1] += 1
        #                         target_depths[target_pos-1] += 1
        #                         sub_score = hmm_params[profile][1][char] - hmm_params[profile][2][hmm_pos][char]
        #                         raw_contributions.append((sub_score, hmm_pos, target_pos))
        #                         target_pos += 1  # Only progress in target coordinate space if it's a match state
        #                     hmm_pos += 1  # Gaps still progress profile coordinate
        #                 else:
        #                     target_pos += 1
        #
        #             # Need to normalise to ensure decomposition sums to original domain score
        #             try:
        #                 norm_factor = sum([pos_data[0] for pos_data in raw_contributions]) / chains[k][i][0]
        #                 log_odds_norm = [(score / norm_factor, hmm_pos, target_pos)
        #                                 for score, hmm_pos, target_pos in raw_contributions]
        #             except ZeroDivisionError:  # Edge case where a domain hit score is 0.0
        #                 log_odds_norm = [(0, hmm_pos, target_pos) for score, hmm_pos, target_pos in raw_contributions]
        #             per_hit_pos_conts.extend(log_odds_norm)
        #
        #         # Sum contributions over hmm positions weighted by position depth
        #         tot_weighted_hmm = 0
        #         for hmm_pos in range(1, hmm_params[profile][0]+1):
        #             for hit_pos in [hit_pos for hit_pos in per_hit_pos_conts if hit_pos[1] == hmm_pos]:
        #                 tot_weighted_hmm += hit_pos[0] / profile_depths[hmm_pos-1]
        #
        #         # Do the same over target positions
        #         tot_weighted_target = 0
        #         for target_pos in range(1, seq_lens[target]+1):
        #             for hit_pos in [hit_pos for hit_pos in per_hit_pos_conts if hit_pos[2] == target_pos]:
        #                 tot_weighted_target += hit_pos[0] / target_depths[target_pos-1]
        #
        #         # Get total weighted score accounting for both hmm and target overlaps
        #         # Also record the overall start and end coordinates (span) of the chain
        #         tot_weighted = 0.5 * (tot_weighted_hmm + tot_weighted_target)
        #         all_weighted_scores[profile][target].append((tot_weighted,
        #                                                      (target_spans[profile][target][k])))

    for temp_file in glob.glob(f"temp_{temp_tag}*"):
        os.remove(temp_file)

    return all_weighted_scores


def hmm_search_chain(
        hmm_files,
        target_seqs,
        results_file=None,
        annot_file=None,
        region_type="env",
        overlap_tol_prop=0.1,
        quiet=True,
        max_on=True,
        max_threads=2,
        nice=None
):
    """ Search one or multiple pHMMs against a sequence database and return the score of the best scoring chain as
      determined by profile chain assembly and chain weighting functions of this module. """

    rand_temp_tag = random.randint(int(1e7), int(1e8)-1)  # Avoid identically named temp files if threading

    if not isinstance(hmm_files, list):

        # Check if files have already been concatenated
        with open(hmm_files) as in_file:
            cnt = 0
            for line in in_file:
                if line.strip().startswith("//"):
                    cnt+=1
        if cnt <= 1:  # Single hmm query
            hmm_files = [hmm_files]
            # else leave the file as is - will note that it's pre-concatenated at search time

    if isinstance(hmm_files, list):  # Need to concatenate
        if len(hmm_files) > 1:
            query_file = f"temp_{rand_temp_tag}_all.hmm"
            hmm_concat(hmm_files, query_file)
        else:
            query_file = hmm_files[0]

    else:  # Pre-concatenated file
        query_file = hmm_files

    if not results_file:
        results_file = f"temp_{rand_temp_tag}_results.txt"

    # Run hmmsearch
    args = []
    if nice:
        args.extend(["nice", "-n", str(nice)])
    args.extend(["hmmsearch", "--cpu", str(max_threads), "--domtblout", results_file])
    if max_on:
        args.append("--max")
    if quiet:
        args.extend(["-o", f"temp_{rand_temp_tag}_output.txt"])
    args.extend([query_file, target_seqs])
    subprocess.run(args)

    # Partition hits into domain chains for each query-target pair
    full_dom_results = parse_domtblout_full(results_file)
    dom_chains = get_profile_chains(full_dom_results, region_type=region_type, overlap_tol_prop=overlap_tol_prop)

    # Get weighted chain scores and store the best score per hmm-target pair
    # Chains scores are returned alongside the absolute (start, end) coordinates for the chain
    all_chain_scores = get_weighted_chain_scores(dom_chains, target_seqs, region_type="env",
                                                 temp_tag=rand_temp_tag, max_workers=max_threads)
    best_chain_scores = {profile : {} for profile in all_chain_scores}
    for profile, target_chain_scores in all_chain_scores.items():
        for target, chain_scores in target_chain_scores.items():
            best_chain_scores[profile][target] = max(chain_scores, key=lambda x:x[0])

    if annot_file:
        annot_dict = {}
        for profile, target_chains in best_chain_scores.items():
            for target, chain_data in target_chains.items():
                try:
                    annot_dict[target][profile] = chain_data[0]
                except KeyError:
                    annot_dict[target] = {profile : chain_data[0]}
        annots.annot_file_from_dict(annot_file, annot_dict)
        annots.create_itol_metadata(annot_file.split(".")[0] + ".itol", annot_file)

    for temp_file in glob.glob(f"temp_{rand_temp_tag}*"):
        os.remove(temp_file)

    return best_chain_scores


def get_domain_hit_sequence():
    """ For a one or more of domain hits for a given target sequence against a profile HMM, """
    return


def hmm_search(
        profile_path,
        seq_db,
        no_return=False,
        results_file=None,
        annot_file=None,
        region_type="full",
        id_format=None,
        best_domain=True,
        quiet=True,
        max_on=False,
        max_threads=2,
        nice=None
):
    """ Search a single pHMM against a sequence database. Either return a dictionary of {seq_id : score}, optionally,
    saving the hmmsearch results file by specifying a results file name, OR use no_return flag to create results file
    only. A tab-delimited file (FigTree annotation format) can also optionally be written. """

    if no_return and not results_file:
        raise RuntimeError("Output file must be specified if no_return flag is used.")

    out_file = results_file if results_file else f"{profile_path.split('.')[0]}_results.txt"

    nice_option = f"nice -n {nice} " if nice else ""

    # Region type (full sequence or domain alignmnet/envelope) affects required results format
    if region_type == "full":
        output_option = f"--tblout {out_file}"
    elif region_type in ["ali", "env"]:
        output_option = f"--domtblout {out_file}"
    else:
        raise RuntimeError(f"Invalid region type: {region_type}.")

    # Run hmmsearch
    os.system(f"{nice_option}hmmsearch --cpu {max_threads} {'--max ' if max_on else ''}{'-o temp.txt ' if quiet else ''}"
              f"{output_option} {profile_path} {seq_db}")

    if annot_file or not no_return:
        # Scores from results file
        if region_type == "full":
            scores = parse_tblout(out_file, id_format, best_domain=best_domain)
        else:
            scores = parse_dombtblout(out_file, envelope=True if region_type == "env" else False)

    if annot_file:  # TODO: Not adapted for cases where domtblout used
        annot_dict = {seq : {profile_path.split('.hmm')[0]+"_score" : score} for seq, score in scores.items()}
        annots.annot_file_from_dict(annot_file, annot_dict)

    if no_return:  # Only want the results file
        return

    if not results_file:
        os.remove(out_file)

    return scores


def hmm_threshold_search(
        profile,
        seq_db,
        in_seqs,
        annot_name=None,
        max_threads=2,
        nice=None
):
    """ Generate annotation dictionary containing hmmsearch scores for a given profile against sequences (generally
    across a tree, differentiating between 'in' sequences, presumably with which the profile was constructed, and all
    other ('threshold') sequences. In-sequences can be provided as a list of IDs or as an alignment file from which
    IDs will be obtained. """

    if isinstance(in_seqs, str):  # Must be a file
        in_seqs = [seq.name for seq in SeqIO.parse(in_seqs, 'fasta')]

    if not annot_name:  # Use default name
        annot_name = profile.split('.hmm')[0] + "_scores"

    # Get scores for all seqs
    all_scores = hmm_search(profile, seq_db, max_on=True, max_threads=max_threads, nice=nice)

    annot_dict = {seq : {annot_name : score} for seq, score in all_scores.items()}  # Annotation field for all sequences

    # Annotation fields for in vs threshold sequences
    for seq in annot_dict:
        score = annot_dict[seq][annot_name]
        if seq in in_seqs:
            annot_dict[seq][annot_name+'_in'] = score
        else:
            annot_dict[seq][annot_name+'_out'] = score

    return annot_dict


def hmm_search_dom(
        profiles,
        target_seqs,
        no_return=False,
        results_file=None,
        quiet=True,
        thresholds=None,
        region_type="env",
        max_on=False,
        max_cpu=2,
        nice=None
):
    """ Run hmmsearch inlcuding domain information in results. Results can be written to a file or returned as a
    dictionary mapping {seq_id : (domain_score : (hit_start, hit_end))}. Boundaries can be the domain (default) or
     envelope boundaries. Where multiple domain hits are found, the best (determined by score) is incldued only. """

    temp_tag = random.randint(1e8, 1e9-1)

    # Determine if profiles need to be pressed
    if isinstance(profiles, list):  # Individual profiles provided, so concatenate
        hmm_file = f"all_{temp_tag}.hmm"
        hmm_concat(profiles, hmm_file)
        rm_tmp_profiles = True
    else:  # Profiles already concatenated
        hmm_file = profiles
        rm_tmp_profiles=False

    n_profiles = count_profiles(hmm_file)

    if not thresholds:  # Apply 0 thresh to all profs
        thresholds = [0 for _ in range(n_profiles)]
    elif not isinstance(thresholds, list):  # Apply single thresh to all profs
        thresholds = [thresholds for _ in range(n_profiles)]
    elif len(thresholds) != n_profiles:
        raise RuntimeError(f"Number of provided thresholds is not equal to number of profiles.")

    if region_type == "ali":
        region_type_idx = 2
    elif region_type == "env":
        region_type_idx = 3
    else:
        raise RuntimeError(f"Region type {region_type} is not implemented.")

    if no_return and not results_file:
        raise RuntimeError("Output file must be specified if no_return flag is used.")

    out_file = results_file if results_file else f"{temp_tag}_results.txt"

    args = []

    if nice:
        args.extend(["nice", "-n", str(nice)])

    args.extend(["hmmsearch", "--cpu", str(max_cpu), "--domtblout", out_file])

    if quiet:
        args.extend(["-o", f"{temp_tag}_out.txt"])

    if max_on:
        args.append("--max")

    args.extend([hmm_file, target_seqs])

    subprocess.run(args, check=True)

    if no_return:  # Only want the results file

        if not results_file:
            os.remove(out_file)
        if rm_tmp_profiles:
            os.remove(hmm_file)
        if quiet:
            os.remove(f"{temp_tag}_out.txt")

        return

    else:

        # Get full results (including multiple hits and all region types)
        full_results = parse_domtblout_full(out_file)

        best_results = {}
        for prof, target_results in full_results.items():
            best_results[prof] = {}
            for target, results in target_results.items():
                if results:
                    best_hit = max(results, key=lambda x:x[0])
                    best_results[prof][target] = (best_hit[0], best_hit[region_type_idx])

        if not results_file:
            os.remove(out_file)
        if rm_tmp_profiles:
            os.remove(hmm_file)
        if quiet:
            os.remove(f"{temp_tag}_out.txt")

        return best_results


# def hmm_weighted_domain_search(profiles, seqs, results_file=None, quiet=True, threshold=0, envelope=False,
#                       max_on=False, max_threads=2, nice=None):
#     """ Run hmmsearch for a single pHMM with the score being the maximum of {the best single domain score; the sum of
#     domain scores weighted by per-position profile coverage depth}.
#
#      TODO: Currently just implemented for a single pHMM. """
#
#     if isinstance(profiles, str):
#         profiles = [profiles]
#
#     if len(profiles) > 1:
#         hmm_concat(profiles, "temp_all.hmm")
#         profile_option = "temp_all.hmm"
#     else:
#         profile_option = profiles[0]
#
#     out_file = results_file if results_file else "temp_results.txt"
#
#
#     os.system(
#         f"{nice_option}hmmsearch --cpu {max_threads} {'--max ' if max_on else ''}{'-o temp.txt ' if quiet else ''}"
#         f"--domtblout {out_file} {profile} {seqs}")
#
#     hmmer_args = []
#
#     if nice:
#         hmmer_args.extend(["nice", "-n", nice])
#
#     hmmer_args.extend(["hmmsearch", "--cpu", str(max_threads), "--domtblout", out_file, "--max"])
#
#     if quiet:
#         hmmer_args.extend(["-o", "temp_output.txt"])
#
#     hmmer_args.extend([profile_option, seqs])
#
#     subprocess.run(hmmer_args)
#
#
#
#
#
#
#
#
#
#     # Run hmmsearch
#     hmmer_args = ["hmsearch", ]
#     subprocess.run(())
#
#
#     if not results_file:
#         os.remove("temp_results.txt")
#
#     if len(profiles) > 1:
#         hmm_concat(profiles, "temp_all.hmm")


def hmm_search_multi(
        profiles,
        seq_db,
        no_return=False,
        results_file=None,
        id_format=None,
        best_domain=True,
        region_type="full",
        name=None,
        quiet=True,
        max_on=False,
        max_threads=2,
        nice=None
):
    """ Perform hmmsearch over a sequence database for several previously defined profile HMMs. pHMM input can either be
     in the form of a list of file names of single HMMs to be concatenated, or the file name of pre-concatenated HMMs.
     Either returns None and outputs tblout- or dombtblout-formatted results file, or returns a data frame profiles as
     columns, indexed by sequence ID. """

    if no_return and not results_file:
        raise RuntimeError("Output file must be specified if no_return flag is used.")

    out_file = results_file if results_file else "results.txt"

    nice_option = f"nice -n {nice} " if nice else ""

    # Determine if profiles need to be pressed
    if isinstance(profiles, list):  # Individual profiles provided, so concatenate
        hmm_file = f"{name}.hmm" if name else "all.hmm"
        hmm_concat(profiles, hmm_file)
    else:  # Profiles already concatenated
        hmm_file = profiles

    rand_temp_tag = random.randint(int(1e7), int(1e8)-1)  # Avoid identical output file names if threading

    # Region type (full sequence or domain aligmnet/envelope) affects required results format
    if region_type == "full":
        output_option = f"--tblout {out_file}"
    elif region_type in ["ali", "env"]:
        output_option = f"--domtblout {out_file}"
    else:
        raise RuntimeError(f"Invalid region type: {region_type}.")

    # Run hmmsearch
    os.system(f"{nice_option}hmmsearch --cpu {max_threads} {'--max ' if max_on else ''}{f'-o temp{rand_temp_tag}.txt ' if quiet else ''}"
              f"{output_option} {hmm_file} {seq_db}")

    if no_return:  # Only want the results file
        return

    # Get score df from results file
    if region_type == "full":
        scores = parse_tblout(out_file, id_format, best_domain=best_domain)
    else:
        scores = parse_dombtblout(out_file, envelope=True if region_type == "env" else False)

    if not results_file:
        os.remove(out_file)

    if quiet:
        os.remove(f'temp{rand_temp_tag}.txt')

    return scores


def hmm_search_custom(
        profiles,
        target_seqs,
        mode="dom",
        region_type="env",
        chain_overlap_tol=0.1,
        max_on=False,
        annot_file=None,
        annot_prefix=None,
        cpu_per_search=2,
        nice=None
):
    """ In future, to be used as a single call point for hmmsearch functionality (all modes) outside of this module.
      For all modes, this function returns a dictionary: {profile : {target_seq : (score, (start, end . """

    if mode == "dom":
        results = hmm_search_dom(
            profiles,
            target_seqs,
            region_type=region_type,
            max_on=max_on,
            max_cpu=cpu_per_search,
            nice=nice
        )


    elif mode == "chain":
        results = hmm_search_chain(
            profiles,
            target_seqs,
            region_type=region_type,
            overlap_tol_prop=chain_overlap_tol,
            max_on=max_on,
            max_threads=cpu_per_search,
            nice=nice
        )

    else:
        raise RuntimeError(f"hmm_search mode {mode} is not implemented.")

    if annot_file:

        score_annots = {}

        for prof in results:
            if annot_prefix:
                annot_prefix = annot_prefix.rsplit('_', maxsplit=1)[0] + '_'
            else:
                annot_prefix = ""
            annot_label = annot_prefix + prof

            for seq, hit_data in results[prof].items():
                try:
                    score_annots[seq][annot_label] = hit_data[0]
                except KeyError:
                    score_annots[seq] = {annot_label : hit_data[0]}

        # Check if annot file exists already and score annots can be merged
        if Path(annot_file).exists():
            in_file_arg = [annot_file]
        else:
            in_file_arg = None

        annots.merge_annots(
            annot_dicts=[score_annots],
            annot_files=in_file_arg,
            out_file=annot_file
        )

        annots.create_itol_metadata(annot_file.rsplit('.')[0] + ".itol",
                                    annot_file)

    return results


def profile_db(
        profiles,
        thresholds,
        seq_db,
        results_file=None,
        hits_file=None,
        best_domain=True,
        extract_as="full",
        max_on=False,
        full_hit_file=None,
        id_format=None,
        name=None,
        no_return=True,
        max_threads=2,
        nice=None
):

    # TODO: Currently outputs hits for all profiles to the same fasta. Implement separate files if desired
    """ Run hmmsearch for a set of given profiles and thresholds. Inputs equivalent to hmm_search_multi.
    Output sequences either as a dictionary: {seq_id : {profile_1 : score}, ... {profile_n : score}, .... }  for all
    sequence/profile combinations above the respective thresholds, or write all hits to a fasta file, or both. The full
    sequence containing the hit, or the best domain's alignment or envelope only, can be written to the hits file. """

    if no_return and not hits_file:
        raise RuntimeError("If no_return is selected, a hits file name must be provided.")

    if not isinstance(profiles, list):  # Just one profile
        profiles = [profiles]
        if not isinstance(thresholds, list):
            thresholds = [thresholds]

    if extract_as not in ["full", "ali", "env"]:
        raise RuntimeError("Hits may be extracted as full sequences, or the domain alignment/envelope only.")

    # Need to provide results file explicitly to hmm_search as using no_return option
    remove_temp = False if results_file else True
    results_file = results_file if results_file else "results.txt"

    if len(profiles) > 1:  # Multiple profiles

        # Run search (not applying thresholds yet)
        scores = hmm_search_multi(profiles, seq_db, results_file=results_file, best_domain=best_domain,
                                  id_format=id_format, region_type=extract_as, name=name, max_on=max_on,
                                  max_threads=max_threads, nice=nice)

        # Profile names in tblout should correspond to hmm file names
        profile_names = [profile.split('.')[0] for profile in profiles]

        # Map profile name to respective threshold
        thresh_map = {profile_names[i] : thresholds[i] for i in range(len(profiles))}  # {Profile name : threshold}

        if extract_as == "full":  # Results will be in DataFrame - TODO: standardise the output of hmmsearch_multi to be nested dict

            try:  # TODO: For debugging - remove this outer exception block

                # Collect all hits over respective profile thresholds
                hits = {}  # {Seq_ID : {profile1:score1, ... }, ... } for seq-profile combos over provided thresholds
                for profile in scores.columns:
                    for seq_id in scores.index:
                        if scores[profile][seq_id] >= thresh_map[profile]:  # Apply threshold
                            try:
                                hits[seq_id][profile] = scores[profile][seq_id]
                            except KeyError:
                                hits[seq_id] = {profile : scores[profile][seq_id]}

            except Exception as e:
                print(f"Error for segment {seq_db}:\n{e}")

            boundary_map = None

        else:  # Results will be in nested dict with domain boundary information

            hits = {}
            boundary_map = {}  # NOTE: uses boundary of highest scoring hit for extraction
            for seq_id, score_data in scores.items():
                for profile, hit_info in score_data.items():
                    if hit_info[0] >= thresh_map[profile]:
                        try:
                            hits[seq_id][profile] = hit_info[0]
                        except KeyError:
                            hits[seq_id] = {profile : hit_info[0]}
                        try:
                            if hit_info[1] > boundary_map[seq_id]:
                                boundary_map[seq_id] = hit_info[1]
                        except KeyError:
                            boundary_map[seq_id] = hit_info[1]

    else:   # Single profile

        scores = hmm_search(profiles[0], seq_db, results_file=results_file, best_domain=best_domain,
                            region_type=extract_as, id_format=id_format, max_threads=max_threads, nice=nice)
        hits = {}

        if extract_as == "full":  # We just have dictionary mapping {sequence : score}

            for name, score in scores.items():
                if score >= thresholds[0]:
                    hits[name] = score

            boundary_map = None

        else:  # Have dictionary mapping {sequence : (score, (boundary_start, boundary_end))}
            boundary_map = {}
            for name, hit_info in scores.items():
                if hit_info[0] > thresholds[0]:
                    hits[name] = hit_info[0]
                    boundary_map[name] = hit_info[1]

    # Remove full results file if not required
    if remove_temp:
        os.remove(results_file)

    # Extract hits from seq_db if required
    if hits_file:
        # TODO: Currently just putting all hits in the same file

        target_seqs = list(hits.keys())

        if boundary_map:
            boundaries = [boundary_map[target] for target in target_seqs]
        else:
            boundaries = None

        file_util.extract_subset_fasta(seq_db, hits_file, target_seqs, boundaries)
        # parsers.raw_extract_fasta(seq_db, hits_file, target_seqs, boundaries, id_format)

        # TODO: Add option to additionally save full sequences for hits when boundaries are used

    if not no_return:
        return hits


def get_hit_names(
        profiles,
        thresholds,
        seq_db,
        id_format=None
):
    """ Return a non-redundant list of hits for a sequence database against a collection of profiles and thresholds. """

    hits = profile_db(profiles, thresholds, seq_db, id_format=id_format, no_return=False)
    return list(hits.keys())


def partial_clade_strap(
        tree_file,
        aln_file,
        parent_boundaries,
        in_boundaries,
        profile_names,
        out_boundaries,
        excl_nodes=None,
        realign=True,
        out_path=None
):
    """ Given confident internal clade(s) of a parent clade, identify sequences outside confident clades for removal on
       basis of thresholds defined by sequences known to be exterior to the clade of interest. All clades are defined
       by a list of sequences for which the LCA is the boundary. Optionally, exclude certain clade sequences.
       Assumes the tree is rooted. Currently, requires an input tree and alignment containing both the parent and any
       external threshold clades. """

    if len(profile_names) != len(in_boundaries):
        raise RuntimeError("Number of profile names and boundaries do not match.")

    if excl_nodes:
        if len(excl_nodes) != len(in_boundaries):
            raise RuntimeError("Incorrect exclusion list length. Use empty sub-lists for clades without exclusions.")
    else:
        excl_nodes = [[] for _ in in_boundaries]

    if out_path:
        shutil.rmtree(out_path, ignore_errors=True)
        os.mkdir(out_path)
        os.chdir(out_path)

    for i in range(len(in_boundaries)):
        build_clade_hmm_old(tree_file, aln_file, in_nodes=in_boundaries[i], excl_nodes=excl_nodes[i], realign=realign,
                            clade_name=profile_names[i])

    t = tree.load_tree(tree_file)

    # Determine sequences internal to parent clade to be interrogated
    parent_names = tree.get_subtree_leaves(t, ext_nodes=parent_boundaries)
    in_nodes = []   # Nodes making up confident profile(s) - don't need to check scores of these
    for boundary in in_boundaries:  # Add nodes from all internal clades
        in_nodes.extend(tree.get_subtree_leaves(t, ext_nodes=boundary))
    test_nodes = [name for name in parent_names if name not in in_nodes]
    # Add exclusion nodes which are currently internal to defined in_clades
    for excl_seqs in excl_nodes:
        test_nodes.extend(excl_seqs)
    # Extract seqs to fasta
    file_util.extract_fasta(aln_file, "profile_strap_test.fa", test_nodes, ungap=True)

    # Get out_clade sequences to define threshold
    out_nodes = []
    for boundary in out_boundaries:
        out_nodes.extend(tree.get_subtree_leaves(t, ext_nodes=boundary))
    file_util.extract_fasta(aln_file, "profile_strap_out.fa", out_nodes, ungap=True)

    # Maps seqs in test set to number of in_profiles for which it scores above threshold
    test_hits = {name : 0 for name in test_nodes}   # TODO: For now just trying this method of evaluating/reporting.

    for profile in profile_names:

        # score out seqs to determine threshold
        out_scores = list(hmm_search(f"{profile}.hmm", 'profile_strap_out.fa', id_format=None).values())
        out_scores.sort(reverse=True)
        threshold = out_scores[0]  # Take threshold as highest out-node score

        # Score test seqs against threshold
        test_scores = hmm_search(f"{profile}.hmm", 'profile_strap_test.fa', id_format=None)
        for seq_name, score in test_scores.items():
            if score >= threshold:   # If test seq is passing for this profile, add to number of hits
                test_hits[seq_name] += 1

    if out_path:
        os.chdir("../..")
        shutil.rmtree(out_path, ignore_errors=True)

    return test_hits


def cs_worker_init(
        out_path,
        subtree_path,
        target_leaves,
        full_fa,
        thresh_fa,
        target_aln,
        holdout_sets,
        hmm_mode,
        max_on,
        chain_overlap_tol,
        full_results,
        lock,
        max_clade_proportion,
        threads_per_search,
        nice
):
    """ Ensure availability of subtree and objects defined in the parent process for each worker. """

    # # TODO: This logging is just for debugging - delete this once it's working
    # pid = mp.current_process().pid
    # logfile = f"worker_{pid}.log"
    # sys.stdout = open(logfile, "w", buffering=1)  # line-buffered
    # sys.stderr = sys.stdout

    global _out_path,_subtree_path,_target_leaves,_full_fa,_thresh_fa,_target_aln,_hmm_mode,_chain_overlap_tol,_max_on
    global _holdout_sets,_full_results,_lock,_max_clade_proportion,_threads_per_search,_nice

    _out_path = out_path
    _subtree_path = subtree_path
    _target_leaves = target_leaves
    _full_fa = full_fa
    _thresh_fa = thresh_fa
    _target_aln = target_aln
    _hmm_mode = hmm_mode
    _chain_overlap_tol = chain_overlap_tol
    _max_on = max_on
    _holdout_sets = holdout_sets
    _full_results = full_results
    _lock = lock
    _max_clade_proportion = max_clade_proportion
    _threads_per_search = threads_per_search
    _nice = nice


def assess_sub_clade(node_id):
    """ Holdout evaluation process for each eligible sub-clade in the broader target clade. """

    # Read in target clade from subtree file



    subtree = tree.load_tree(_subtree_path, nwk_format=1)

    node = subtree & node_id


    if len(node) <= len(subtree) * _max_clade_proportion:  # Make sure clade isn't > max proportion

        assessment_id = random.randint(1000000000000, 9999999999999)  # Tag for temp files

        # Define holdout seqs under this ancestral node
        if node.is_leaf():
            holdout_seqs = [node.name]
        else:
            holdout_seqs = node.get_leaf_names()

        # holdout_seqs_fa = f"{_out_path}/{assessment_id}_holdout.fa"
        # retain_sub_aln = f"{_out_path}/{assessment_id}_retain.aln"
        # sub_hmm = f"{_out_path}/{assessment_id}.hmm"

        holdout_seqs_fa = f"{_out_path}_{assessment_id}_holdout.fa"
        retain_sub_aln = f"{_out_path}_{assessment_id}_retain.aln"
        sub_hmm = f"{_out_path}_{assessment_id}.hmm"

        # These are seqs to exclude from profile building and test against
        file_util.extract_fasta(_full_fa, holdout_seqs_fa, holdout_seqs)

        # Remaining seqs under target clade can be retained for pHMM building - extract from target_aln
        retain_seqs = [name for name in _target_leaves if name not in holdout_seqs]
        file_util.extract_fasta(_target_aln, retain_sub_aln, retain_seqs)

        # Build pHMM from retained sequences
        build_hmm(retain_sub_aln, sub_hmm)

        # Score threshold (outer) seqs
        if _hmm_mode == "chain":  # TODO: Fix this after modularising all hmmsearch modes
            # TODO: Expose option for dom vs env region window (currently uses default - env)
            chain_results = hmm_search_chain(sub_hmm,
                                             _thresh_fa,
                                             overlap_tol_prop=_chain_overlap_tol,
                                             max_on=_max_on,
                                             max_threads=_threads_per_search,
                                             nice=_nice)

            chain_results = list(chain_results.values())[0]  # Only one profile in results
            thresh_scores = {seq_id : chain_info[0] for seq_id, chain_info in chain_results.items()}

        else:  # Uses best domain score if not chain
            thresh_scores = hmm_search(
                sub_hmm,
                _thresh_fa,
                id_format=None,
                max_on=_max_on,
                max_threads=_threads_per_search,
                nice=_nice)

        thresh_scores = list(thresh_scores.items())
        thresh_scores.sort(reverse=True, key=lambda x: x[1])  # Sort high to low

        # Score heldout seqs
        if _hmm_mode == "chain":

            chain_results = hmm_search_chain(
                sub_hmm,
                holdout_seqs_fa,
                overlap_tol_prop=_chain_overlap_tol,
                max_on=_max_on,
                max_threads=_threads_per_search,
                nice=_nice)

            chain_results = list(chain_results.values())[0]  # Only one profile in results
            test_scores = {seq_id: chain_info[0] for seq_id, chain_info in chain_results.items()}

        else:

            test_scores = hmm_search(
                sub_hmm,
                holdout_seqs_fa,
                id_format=None,
                max_on=_max_on,
                max_threads=_threads_per_search,
                nice=_nice)

        test_scores = list(test_scores.items())
        test_scores.sort(key=lambda x: x[1])

        # Check if any test seqs score lower than threshold seqs
        sub_results = {}  # Results for this iteration / set of excluded sequences
        for test_name, test_score in test_scores:
            thresh_above = []
            if test_score >= thresh_scores[0][1]:  # All remaining test seqs score above max threshold score
                break
            for thresh_name, thresh_score in thresh_scores:
                if test_score >= thresh_score:  # No more thresh_scores will be above this test_score
                    break
                else:  # This sequence scores higher than the test_seq score
                    thresh_above.append(thresh_name)
            if thresh_above:  # Threshold sequences scoring better than this test seq
                sub_results[test_name] = thresh_above

        with _lock:  # Controlled simultaneous addition of holdout sets and corresponding results to master lists
            _holdout_sets.append(holdout_seqs)
            _full_results.append(sub_results)

        # Remove temp files
        os.remove(holdout_seqs_fa)
        os.remove(retain_sub_aln)
        os.remove(sub_hmm)


def full_clade_strap(
        in_tree,
        aln_file,
        clade_name=None,
        use_int_labels=False,
        target_clade_node=None,
        target_clade_boundary=None,
        out_leaf=None,
        relabeled_tree_file=None,
        external_out_seqs=None,
        hmm_mode="dom",
        chain_overlap_tol=0.1,
        max_on=False,
        max_clade_proportion=0.5,
        realign=False,
        thresh_cluster_id=None,
        thresh_cluster_cov=None,
        linclust=False,
        max_processes=2,
        threads_per_search=6,
        nice=None):
    """ Given an apparent clade in a rooted phylogenetic tree, hold out all sub-clades with less than a defined maximum
    proportion of the larger clade's extant sequences. Iteratively hold out sequences, and build a new profile. Report
    any withheld sequences which score below a threshold defined by profile scores against all sequences which fall
    outside the clade of interest (or cluster representatives of said sequences at a provided threshold), plus optional
    further external sequences provided as SeqRecords. Return a list of dictionaries mapping failing holdout sequences
    to threshold seqs which score above them. This list is indexed by another list containing the set of all sequences
    heldout for each iteration.

    TODO: External out seqs (seqs to threshold against that aren't in current tree) currently only taken as SeqRecords
    """

    # Analysis creates a lot of temporary files
    out_path = clade_name + "_strap"
    shutil.rmtree(out_path, ignore_errors=True)
    os.mkdir(out_path)
    # os.chdir(out_path)

    # Get full tree and target clade
    full_tree = tree.load_tree(in_tree)

    if not external_out_seqs:
        external_out_seqs = []

    # If internal nodes already exist, assume re-rooting is not required
    if not use_int_labels:  # If not, re-root to ensure target clade is monophyletic and label internal nodes

        # NOTE: tree is not deep copied during rooting so target_clade_node is preserved if used
        full_tree = tree.root_tree(full_tree, mode="outgroup", og_bounds=target_clade_boundary, in_leaf=out_leaf,
                                   og_node=target_clade_node, no_copy=True)

        # Label / re-label internal nodes
        full_tree = tree.relabel_internal_nodes(full_tree, out_tree_file=relabeled_tree_file, no_return=False)

    # Get subtree corresponding to clade of interest
    target_clade = tree.get_subtree(full_tree, int_node=target_clade_node, clade_bounds=target_clade_boundary)

    # Name clade as internal node label if name is not explicitly provided
    if not clade_name:
        if len(target_clade.get_children()) == 1:  # Root is preserved during pruning - true LCA is root's child
            clade_name = target_clade.get_children()[0].name
        else:  # If root of original tree is still the LCA of the clade
            clade_name = target_clade.name

    # Subtree needs to be available for individual subprocesses to load and use
    subtree_path = f"{out_path}/{clade_name}.nwk"
    target_clade.write(outfile=subtree_path, format=1)  # TODO: Check this is the correct format

    target_leaves = target_clade.get_leaf_names()

    # Write ungapped aln to file (full seq db)
    full_fa = f"{out_path}/{clade_name}.fa"
    file_util.extract_fasta(aln_file, full_fa, full_tree.get_leaf_names(), ungap=True)

    # Write thresholding seqs to file
    thresh_seqs = [name for name in full_tree.get_leaf_names() if name not in target_clade.get_leaf_names()]
    thresh_fa = f"{out_path}/{clade_name}_thresh.fa"
    if thresh_cluster_id and thresh_cluster_cov:  # Threshold sequences are to be clustered
        file_util.extract_fasta(full_fa, f"{out_path}/temp_thresh_full.fa", thresh_seqs, ungap=True)
        cluster.rr_seqs(f"{out_path}/temp_thresh_full.fa", thresh_fa, thresh_cluster_id, thresh_cluster_cov,
                        linclust=linclust, max_threads=threads_per_search, nice=nice)

    else:  # Use all sequences to threshold
        file_util.extract_fasta(full_fa, thresh_fa, thresh_seqs, ungap=True)

    # Add external out seqs to thresholding set if provided
    # TODO: Probably buggy based on format of external_out_seqs. Currently don't need them but need to fix
    if external_out_seqs:
        if isinstance(external_out_seqs, str):  # Need to read these sequences in
            external_out_seqs = list(SeqIO.parse(external_out_seqs, "fasta"))
        file_util.merge_seqs(in_fastas=[thresh_fa], in_records=[external_out_seqs], out_file=thresh_fa, no_return=True)

    # TODO: Need per-holdout file names when multi-processing
    # Temp file names for each iteration
    # temp_train_fa = "cladestrap_temp_train.fa"
    temp_test_fa = "cladestrap_temp_test.fa"
    temp_aln = "cladestrap_temp.aln"
    temp_hmm = "cladestrap_temp.hmm"

    # If realigning, perform structure-guided alignment once only, and extract relevant sequences in each clade-strap
    # cycle, else just extract aligned targets seqs from the provided alignment and remove all-gap columns
    target_fa = f"{out_path}/cladestrap_target.fa"
    target_aln = f"{out_path}/cladestrap_target.aln"
    if realign:
        file_util.extract_fasta(full_fa, target_fa, target_leaves)
        aln.struct_aln(target_fa, target_aln, quiet=True)
    else:
        file_util.extract_fasta(aln_file, target_aln, target_leaves)
        aln.remove_gap_only_cols(target_aln, target_aln)

    # Implementing multi-processing for per-clade assessments

    manager = mp.Manager()
    # Index for holdout sets. Order should match results output
    holdout_sets = manager.list()
    full_results = manager.list()  # Indexed by holdout_sets
    lock = mp.Lock()  # Ensure the list of holdout sequences and corresponding results are appended simultaneously

    # Set up multi-processing for robustness assessment over all sub-clades
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=max_processes,
            initializer=cs_worker_init,
            initargs=(
                    out_path,
                    subtree_path,
                    target_leaves,
                    full_fa,
                    thresh_fa,
                    target_aln,
                    holdout_sets,
                    hmm_mode,
                    max_on,
                    chain_overlap_tol,
                    full_results,
                    lock,
                    max_clade_proportion,
                    threads_per_search,
                    nice)
            ) as executor:

        future_list = [executor.submit(assess_sub_clade, node.name)
                       for node in target_clade.iter_descendants()]
        for future in as_completed(future_list):
            future.result()

    # os.chdir("..")
    shutil.rmtree(out_path)

    # os.remove(f"{clade_name}.fa")

    return list(holdout_sets), list(full_results)


def clade_association_strength(t, holdout_sets, full_results):
    """ Test the association of each leaf internal to a target clade on the basis of results from a full_clade_strap.
     Association strength is assessed using a framework inspired by item response theory. Each internal leaf node (i) is
      assessed by its ability to 'beat' all external 'threshold' leaves (t) when scored against profiles built during
       iterations where that i is held out.

       TODO: Could compute mean of internals' losses to each threshold across their respective folds (Gi) to address topology effects somewhat differently  """

    # # TODO: For testing, DELETE
    # cas_log = open(f"{str(holdout_sets)[:20]}_{str(full_results)[:20]}", 'a')
    # cas_log.write(t.name + "\n\n")
    # cas_log.write(str(holdout_sets) + "\n\n")
    # cas_log.write(str(full_results) + "\n\n")
    # cas_log.close()

    if len(holdout_sets) != len(full_results):
        raise RuntimeError("Length of holdout sets list and threshold results must be equal and order must correspond to "
                           "one another")

    if not isinstance(t, ete3.Tree):
        t = tree.load_tree(t)

    # Define internals as set of all sequences in holdout set - all should be held out at least once
    internals = set()
    for holdouts in holdout_sets:
        internals.update(holdouts)

    # Define thresholds as any sequence appearing outside that clade - get from the tree
    # TODO: Note that if the holdout sets and results aren't consistent with the tree, results here will not be valid
    thresholds = [leaf for leaf in t.get_leaf_names() if leaf not in internals]

    # Params for strength association
    # n_folds = len(holdout_sets)  # Number of different holdouts
    # n_thresholds = len(thresholds)

    # Determine win rate of each threshold t against heldout internals in each fold (fold-wise difficulty of t)
    diff_g_t = []  # [{t1 : t1_diff for fold 0, t2 : t2_dif.....}, {t1 : t1_diff for fold 1, ....}, .... ]
    for i in range(len(holdout_sets)):

        # Record instances of a threshold sequence beating an internal holdout in this fold
        thresh_win_cnts = {t : 0 for t in thresholds}
        for thresh_seqs in full_results[i].values():
            for thresh_seq in thresh_seqs:
                thresh_win_cnts[thresh_seq] += 1

        diff_g_t.append({})
        for thresh_seq, win_cnt in thresh_win_cnts.items():
            diff_g_t[-1][thresh_seq] = thresh_win_cnts[thresh_seq] / len(holdout_sets[i])   # Per-fold win rate of each threshold seq

    # Calculate global difficulty of each t (per-fold diffs averaged over all folds)
    diff_t = {}
    for t in thresholds:
        diff_t[t] = sum([fold_diff[t] for fold_diff in diff_g_t]) / len(diff_g_t)

    # Fold-wise association weakness of an internal based on losses to threshold sequences weighted by their difficulty
    weak_g_i = []  # [{int1:weak_g1(i_1), int2 : weak_g1(i_2) ...}, {int1 : weak_g2(i_1)},...] (Ordered by holdout_sets)
    for i in range(len(holdout_sets)):

        weak_g_i.append({})
        for int_heldout in holdout_sets[i]:

            if int_heldout in full_results[i].keys():
                weighted_losses = 0  # Weight losses to threshold sequences based on difficulty of threshold
                for high_thresh in full_results[i][int_heldout]:
                    weighted_losses += (1-diff_t[high_thresh])  # Penalise losses to more difficult thresholds less
                # Normalise association weakness by total number of threshold sequences
                weak_g_i[-1][int_heldout] = weighted_losses / len(thresholds)

            else:  # No losses of this internal against threshold seqs
                weak_g_i[-1][int_heldout] = 0

    # Compute final association strength for each internal by averaging over folds which it is held out in
    ass_scores = {}
    for int_seq in internals:
        n_folds_i = 0
        total_weakness = 0
        for fold in weak_g_i:
            if int_seq in fold.keys():
                n_folds_i += 1
                total_weakness += fold[int_seq]
        ass_scores[int_seq] = 1 - total_weakness / n_folds_i

    return ass_scores


def full_clade_assessment(
        tree_file,
        aln_file,
        use_int_labels=False,
        target_clade_node=None,
        target_clade_boundary=None,
        clade_name=None,
        out_leaf=None,
        external_out_seqs=None,
        hmm_mode="dom",
        max_on=False,
        chain_overlap_tol=0.1,
        max_clade_proportion=0.5,
        realign=False,
        thresh_cluster_id=None,
        thresh_cluster_cov=None,
        linclust=False,
        cas_pass_score=0.9,
        min_cas_pass=0.75,
        max_processes=2,
        threads_per_search=6,
        nice=None
):
    """
    Assess the representation robustness of a clade and suitability for use in further sequence curation by full
    cladestrapping and subsequent clade association screening.

    :param tree_file:
    :param aln_file:
    :param use_int_labels:
    :param target_clade_node:
    :param target_clade_boundary:
    :param clade_name:
    :param out_leaf:
    :param external_out_seqs:
    :param hmm_mode:
    :param max_on:
    :param chain_overlap_tol:
    :param max_clade_proportion:
    :param realign:
    :param thresh_cluster_id:
    :param thresh_cluster_cov:
    :param max_processes:
    :param threads_per_search:
    :param nice:
    :return:
    """

    if not external_out_seqs:
        external_out_seqs = []

    # Run full cladestrapping
    holdout_sets, full_results = full_clade_strap(
        tree_file,
        aln_file,
        use_int_labels=use_int_labels,
        target_clade_node=target_clade_node,
        target_clade_boundary=target_clade_boundary,
        clade_name=clade_name,
        out_leaf=out_leaf,
        external_out_seqs=external_out_seqs,
        hmm_mode=hmm_mode,
        max_on=max_on,
        chain_overlap_tol=chain_overlap_tol,
        max_clade_proportion=max_clade_proportion,
        realign=realign,
        thresh_cluster_id=thresh_cluster_id,
        thresh_cluster_cov=thresh_cluster_cov,
        linclust=linclust,
        max_processes=max_processes,
        threads_per_search=threads_per_search,
        nice=nice)

    # Compute per-leaf Clade Association Scores (CASs) on basis of cladestrap results
    cas = clade_association_strength(tree_file, holdout_sets, full_results)

    # Check how many sequences have CAS < min pass threshold
    passed_seqs = []
    for seq, score in cas.items():
        if score >= cas_pass_score:
            passed_seqs.append(seq)

    # Return list of passing sequences if the proportion of passing sequences is above threshold, else return None
    if len(passed_seqs) / len(cas) >= min_cas_pass:
        return passed_seqs
    else:
        return None


def score_dist_compare(
        t,
        seqs_fasta,
        profile,
        in_nodes,
        out_node=None,
        score_metric="chain",
        max_threads=8,
        nice=None
):
    """ For a given sub-clade and associated pHMM, compute the score and distance from the LCA of the sub-clade for each
    external (threshold) sequence within the full tree. Returns two lists corresponding to LCA distances and profile
    scores with equivalent order for each external sequence. """

    t = tree.load_tree(t)
    if out_node:
        t = tree.outgroup_root(t, in_nodes, out_node)
    in_names = tree.get_subtree_leaves(t, ext_nodes=in_nodes)

    out_names = [leaf for leaf in t.get_leaf_names() if leaf not in in_names]

    # Score profile over all external sequences
    if score_metric == "chain":  # Best weighted chain score
        out_dict = hmm_search_chain(
            profile,
            seqs_fasta,
            max_threads=max_threads,
            max_on=True
        )
        prof_name = list(out_dict.keys())[0]
        full_scores = {seq : best_chain[0] for seq, best_chain in out_dict[prof_name].items()}

    elif score_metric == "domain":  # Best domain directly from hmmsearch (HMMER)
        full_scores = hmm_search(
            profile,
            seqs_fasta,
            best_domain=True,
            max_threads=max_threads,
            max_on=True
        )

    else:  # Full domain score
        full_scores = hmm_search(
            profile,
            seqs_fasta,
            best_domain=False,
            max_threads=max_threads,
            max_on=True
        )

    out_scores = {seq : score for seq, score in full_scores.items() if seq in out_names}

    # Fast pairwise patristic distances using treeswift
    nwk_str = t.write()
    ts_tree = treeswift.read_tree_newick(nwk_str)
    all_dists = ts_tree.distance_matrix(leaf_labels=True)

    # Get LCA node for profile clade
    ext1_node = t & in_nodes[0]
    ext2_node = t & in_nodes[1]
    lca = ext1_node.get_common_ancestor(ext2_node)

    # Use a representative of the profile clade to derive all clade LCA to external leaf distances
    ext1_lca_dist = lca.get_distance(ext1_node)  # Distance from LCA to clade rep
    ext1_dists = all_dists[in_nodes[0]]
    lca_ext_dists = {out_name : ext1_dists[out_name]-ext1_lca_dist for out_name in out_names}

    return lca_ext_dists, out_scores


######### TODO: Experimenting with automatic threshold determination below #########################


# def plot_cumulative_correlation(x, y, bin_size=100):
#
#     paired_data = list(zip(x, y))
#     paired_data.sort(key=lambda x:x[1], reverse=True)
#     n = []
#     corrs = []
#     for i in range(0, len(paired_data), bin_size):
#         n.append(i+bin_size)
#         corrs.append(spearmanr([point[0] for point in paired_data[i:i+bin_size]],
#                                    [point[1] for point in paired_data[i:i+bin_size]])[0])
#
#     fig, ax = plt.subplots()
#     ax.plot(n, corrs)
#     ax.invert_yaxis()
#     plt.show()
#
#
# def plot_cumulative_correlation_score(x, y, n_score_bins=20):
#
#     paired_data = list(zip(x, y))
#     paired_data.sort(key=lambda x:x[1], reverse=True)
#     min_score = min(y)
#     max_score = max(y)
#     bin_size = (max_score - min_score) / n_score_bins
#     score_centroids = []
#     corrs = []
#     lower_bound = min_score
#     upper_bound = lower_bound + bin_size
#
#     for i in range(n_score_bins):
#
#         score_centroids.append((lower_bound+upper_bound)/2)
#
#         point_subset = [point for point in paired_data if lower_bound <= point[1] < upper_bound]
#         corrs.append(spearmanr([point[0] for point in point_subset],
#                                    [point[1] for point in point_subset])[0])
#
#     fig, ax = plt.subplots()
#     ax.plot(score_centroids, corrs)
#     ax.invert_yaxis()
#     plt.show()
#
#
# def plot_corr_threshold_differential(x, y, n_thresholds=50):
#
#     max_score = max(y)
#     min_score = min(y)
#     step = (max_score - min_score) / (n_thresholds + 1)
#     thresholds = [min_score + (i+1)*step for i in range(n_thresholds)]
#
#     paired_data = list(zip(x, y))
#     corr_diffs = []
#     for thresh in thresholds:
#         sub_above = [point for point in paired_data if point[1] >= thresh]
#         sub_below = [point for point in paired_data if point[1] < thresh]
#
#         corr_above = spearmanr([point[0] for point in sub_above],
#                                [point[1] for point in sub_above])[0]
#
#         corr_below = spearmanr([point[0] for point in sub_below],
#                                [point[1] for point in sub_below])[0]
#
#         corr_diffs.append(corr_above - corr_below)
#
#     fig, ax = plt.subplots()
#     ax.plot(thresholds, corr_diffs)
#     ax.invert_xaxis()
#     plt.show()
#
#
# def bootstrap_cumulative_correlation(x, y, step_n=100, bootstrap_m=20, bootstrap_k=100):
#
#     paired_data = list(zip(x, y))
#     paired_data.sort(key=lambda x: x[1], reverse=True)
#     n = []
#     corrs = []
#     for i in range(step_n, len(paired_data), step_n):
#
#         n.append(i)
#
#         pool = paired_data[:i]
#         bootstrap_corrs = []
#
#         for j in range(bootstrap_k):
#             sample = random.choices(pool, k=bootstrap_m)
#
#             bootstrap_corrs.append(spearmanr([point[0] for point in sample],
#                                [point[1] for point in sample])[0])
#
#         corrs.append(np.mean(bootstrap_corrs))
#
#
#     fig, ax = plt.subplots()
#     ax.plot(n, corrs)
#     ax.invert_yaxis()
#     plt.show()
#
#
# def sliding_window_correlation(x, y, corr_window=50, kernel_window=50):
#
#     paired_data = list(zip(x, y))
#     paired_data.sort(key=lambda x: x[1], reverse=True)
#     centroids = []
#     corrs = []
#     for i in range(int(len(paired_data) - corr_window)):
#
#         window_data = paired_data[i:i + corr_window]
#         centroids.append(i + corr_window / 2)
#         corrs.append(spearmanr([point[0] for point in window_data],
#                                [point[1] for point in window_data])[0])
#
#     # Compute kernel density
#     kernel_centroid = []
#     kernel_corrs = []
#     for i in range(len(centroids) - kernel_window):
#         kernel_centroid.append(centroids[i] + kernel_window / 2)
#         kernel_corrs.append(np.mean(corrs[i:i+kernel_window]))
#
#     fig, ax = plt.subplots()
#     ax.plot(kernel_centroid, kernel_corrs)
#     ax.invert_yaxis()
#     plt.show()
#
#
# def local_regression(x, y, slope_window=100, smoothing_window=):
#
#     paired_data = list(zip(x, y))
#     paired_data.sort(key=lambda x: x[1], reverse=True)
#     n = []
#     slopes = []
#     for i in range(0, len(paired_data) - slope_window):
#
#         window_data = paired_data[i : i + slope_window]
#         n.append(i + slope_window / 2)
#         slopes.append(linregress([point[0] for point in window_data],
#                                [point[1] for point in window_data])[0])
#
#     fig, ax = plt.subplots()
#     ax.plot(n, slopes)
#     ax.invert_yaxis()
#     plt.show()
#
#
# def partial_spearmans_corr(x, y, prop)
#
#
# def dist_corr():
#     """ Compute the Spearman rank correlation between scores for a clade-specific profile HMM on sequences outside
#      that clade with patristic distance of those sequences from the last common ancestor of the clade on which the
#       profile was built. """
#
#     return


def profile_cousin_clades(
        t,
        seq_fa,
        in_nodes,
        profile,
        root_on=None,
        include_target=False,
        best_domain=False,
        score_metric="chain",
        max_threads=8,
        nice=None
):
    """ Get scores for a clade pHMM binned by cousin clades as defined for tree.cousin_clade_sets. Returns a dictionary
    mapping {cousin_clade_level : [list_of_scores]} beginning at level 0 (target clade itself) if specified with
     include_target_clade, else at 1 (sister clade to target). """

    # Run hmm_search
    if score_metric == "chain":
        out_dict = hmm_search_chain(
            profile,
            seq_fa,
            max_threads=max_threads,
            max_on=True)
        prof_name = list(out_dict.keys())[0]
        all_scores = {seq : best_chain[0] for seq, best_chain in out_dict[prof_name].items()}

    elif score_metric == "domain":
        all_scores = hmm_search(
            profile,
            seq_fa,
            best_domain=True,
            max_on=True,
            max_threads=max_threads,
            nice=nice)

    else:
        all_scores = hmm_search(
            profile,
            seq_fa,
            best_domain=False,
            max_on=True,
            max_threads=max_threads,
            nice=nice)

    # Get cousin clade sets
    cousin_clades = tree.cousin_clade_sets(
        t, in_nodes,
        root_on=root_on,
        include_target=include_target)

    cousin_clade_scores = {}
    i = 0 if include_target else 1
    for cousin_leaves in cousin_clades:
        clade_scores = []
        for cousin in cousin_leaves:
            clade_scores.append(all_scores[cousin])
        cousin_clade_scores[i] = clade_scores
        i+=1

    # TODO: Currently just returning a dictionary mapping level to list of scores for that level
    return cousin_clade_scores


def threshold_clade_hmm(
        profile,
        in_tree,
        in_fasta,
        beat_seq_prop=0.2,
        target_clade_node=None,
        target_clade_bounds=None,
        hmm_mode="dom",
        max_on=False,
        chain_overlap_tol=0.1,
        max_cpu=6,
        nice=None):
    """ Determine an appropriate threshold for a profile HMM representing a clade in a phylogeny using the distribution
    of scores across leaf sequences outside the clade of interest.
    NOTE: MSA can be provided as in_fasta
     TODO: Currently very simplistic. Working on the maths to implement some more sophisticated ideas. """

    temp_tag = random.randint(int(1e7), int(1e8)-1)

    in_tree = tree.load_tree(in_tree)  # Assumes a rooted tree

    # Get leaf labels within the target clade
    if target_clade_node:
        if not isinstance(target_clade_node, ete3.TreeNode):  # Assume node label has been provided
            target_clade_node = in_tree & target_clade_node
        in_seqs = tree.get_subtree_leaves(in_tree, anc=target_clade_node)
    else:
        in_seqs = tree.get_subtree_leaves(in_tree, ext_nodes=target_clade_bounds)

    out_seqs = [seq.name for seq in SeqIO.parse(in_fasta, "fasta") if seq.name not in in_seqs]

    # Get temp file containing only sequences outside of clade
    out_seqs_fa = f"out_{temp_tag}.fa"
    file_util.extract_fasta(in_fasta, out_seqs_fa, out_seqs, ungap=True)

    # Score profile over all external sequences

    if hmm_mode == "domain":  # Best domain directly from hmmsearch
        scores = hmm_search(
            profile,
            out_seqs_fa,
            best_domain=True,
            max_threads=max_cpu,
            max_on=max_on,
            nice=nice)

    elif hmm_mode == "chain":  # Best weighted chain score
        out_dict = hmm_search_chain(
            profile,
            out_seqs_fa,
            overlap_tol_prop=chain_overlap_tol,
            max_threads=max_cpu,
            nice=nice,
            max_on=max_on)
        prof_name = list(out_dict.keys())[0]
        scores = {seq: best_chain[0] for seq, best_chain in out_dict[prof_name].items()}

    else:  # Full seq score
        scores = hmm_search(
            profile,
            out_seqs_fa,
            best_domain=False,
            max_threads=max_cpu,
            nice=nice,
            max_on=max_on)

    # Order seqs by score
    scores = list(scores.items())
    scores.sort(key=lambda x:x[1], reverse=True)

    # Determine index of sequence to determine threshold from
    thresh_idx = round(len(scores) * beat_seq_prop) - 1

    os.remove(out_seqs_fa)

    return scores[thresh_idx][1]





#### TODO: This is the old version before multi-processing implementation - retain for now
# def full_clade_strap_parallel(tree_file, aln_file, target_clade_boundary, clade_name, external_out_seqs=[],
#                      max_clade_proportion=0.5, realign=False):
#     """ Given an apparent clade in a rooted phylogenetic tree, hold out all sub-clades with less than a defined maximum
#     proportion of the larger clade's extant sequences. Iteratively hold out sequences, and build a new profile. Report
#     any withheld sequences which score below a threshold defined by profile scores against all sequences which fall
#     outside the clade of interest, plus optional further external sequences provided as SeqRecords. Return a list of
#     dictionaries mapping failing holdout sequences to threshold seqs which score above them. This list is indexed by
#     another list containing the set of all sequences heldout for each iteration.
#
#     TODO: External out_seqs (seqs to threshold against that aren't in current tree) currently only taken as Biopython SeqRecords
#     """
#
#     # Analysis creates a lot of temporary files
#     out_path = clade_name + "_strap"
#     shutil.rmtree(out_path, ignore_errors=True)
#     os.mkdir(out_path)
#     os.chdir(out_path)
#
#     # Get full tree and target clade
#     if isinstance(tree_file, ete3.Tree):  # Check if an ete3 Tree object has been provided
#         full_tree = tree_file
#     else:    # Must be a file
#         full_tree = tree.load_tree("../"+tree_file)
#     target_clade = tree.get_subtree(full_tree, ext_nodes=target_clade_boundary)
#     target_leaves = target_clade.get_leaf_names()
#
#     # Write ungapped aln to file (full seq db)
#     full_fa = f"{clade_name}.fa"
#     parsers.extract_fasta("../"+aln_file, full_fa, full_tree.get_leaf_names(), id_format=None, ungap=True)
#
#     # Write thresholding seqs to file
#     thresh_seqs = [name for name in full_tree.get_leaf_names() if name not in target_clade.get_leaf_names()]
#     thresh_fa = f"{clade_name}_thresh.fa"
#     parsers.extract_fasta(full_fa, thresh_fa, thresh_seqs, id_format=None, ungap=True)
#
#     # Add external out seqs to thresholding set if provided if provided
#     if external_out_seqs:
#         if isinstance(external_out_seqs, list):  # External seqs provided as Bio Seq Records
#             parsers.merge_seqs(in_fastas=[thresh_fa], in_records=[external_out_seqs], out_file=thresh_fa, no_return=True)
#         else:  # Must instead be a fasta file
#             parsers.merge(in_fastas=[thresh_fa, external_out_seqs], out_file=thresh_fa, no_return=True)
#
#     # Temp file names for each iteration
#     temp_train_fa = "cladestrap_temp_train.fa"
#     temp_test_fa = "cladestrap_temp_test.fa"
#     temp_aln = "cladestrap_temp.aln"
#     temp_hmm = "cladestrap_temp.hmm"
#
#     # Index for holdout sets. Order should match results output
#     idx = []
#
#     full_results = []  # Indexed by idx (set of sequences heldout)
#
#     # If realigning, perform structure-guided alignment once only, and extract relevant sequences in each clade-strap
#     # cycle, else just extract aligned targets seqs from the provided alignment and remove all-gap columns
#     target_fa = "cladestrap_target.fa"
#     target_aln = "cladestrap_target.aln"
#     if realign:
#         parsers.extract_fasta(full_fa, target_fa, target_leaves, id_format=None)
#         aln.struct_aln(target_fa, target_aln, quiet=True)
#     else:
#         parsers.extract_fasta('../'+aln_file, target_aln, target_leaves, id_format=None)
#         aln.remove_gap_only_cols(target_aln, target_aln)
#
#
#     # Visit all nodes, including leaves and excluding the root, check the clade under it should be clade-strapped
#     for node in target_clade.iter_descendants():
#
#         # Does this sub-clade comprise > max_clade_proportion
#         if len(node) > len(target_clade) * max_clade_proportion:
#             continue
#
#         # Define holdout seqs under this ancestral node
#         if node.is_leaf():
#             holdout_seqs = [node.name]
#         else:
#             holdout_seqs = node.get_leaf_names()
#
#         # These are seqs to exclude from profile building and test against
#         parsers.extract_fasta(full_fa, temp_test_fa, holdout_seqs, id_format=None)
#
#         # Remaining seqs under target clade can be retained for pHMM building - extract from target_aln
#         retain_seqs = [name for name in target_leaves if name not in holdout_seqs]
#         parsers.extract_fasta(target_aln, temp_aln, retain_seqs, id_format=None)
#
#         # Build pHMM
#         build_hmm(temp_aln, temp_hmm)
#
#         # Score threshold (outer) seqs
#         thresh_scores = hmm_search(temp_hmm, thresh_fa, id_format=None)
#         thresh_scores = list(thresh_scores.items())
#         thresh_scores.sort(reverse=True, key=lambda x:x[1])  # Sort high to low
#
#         # Score test seqs
#         test_scores = hmm_search(temp_hmm, temp_test_fa, id_format=None)
#         test_scores = list(test_scores.items())
#         test_scores.sort(key=lambda x:x[1])
#
#         # Check if any test seqs score lower than threshold seqs
#         sub_results = {}   # Results for this iteration / set of excluded sequences
#         for test_name, test_score in test_scores:
#             thresh_above = []
#             if test_score >= thresh_scores[0][1]:   # All remaining test seqs score above max threshold score
#                 break
#             for thresh_name, thresh_score in thresh_scores:
#                 if test_score >= thresh_score:  # No more thresh_scores will be above this test_score
#                     break
#                 else:  # This sequence scores higher than the test_seq score
#                     thresh_above.append(thresh_name)
#             if thresh_above:  # Threshold sequences scoring better than this test seq
#                 sub_results[test_name] = thresh_above
#
#         if sub_results:   # At least one leaf in this clade is not supported
#             idx.append(holdout_seqs)
#             full_results.append(sub_results)
#
#         # Remove temp files
#         os.remove(temp_test_fa)
#         os.remove(temp_aln)
#         os.remove(temp_hmm)
#
#     os.chdir("..")
#     shutil.rmtree(out_path)
#
#     return idx, full_results





def clade_support(
        tree_file,
        aln_file,
        target_clade_boundary,
        clade_name,
        external_out_seqs=None,
        max_clade_proportion=0.5
):
    """ For a given clade in a phylogenetic tree, return True if full_clade_strap returns no non-supported leaves,
     False otherwise.
     # TODO: Now deprecated """

    if not external_out_seqs:
        external_out_seqs = []

    non_support_leaves = full_clade_strap(
            tree_file,
            aln_file,
            target_clade_boundary,
            clade_name,
            external_out_seqs,
            max_clade_proportion)[1]


    if len(non_support_leaves) > 0:
        return True
    else:
        return False


def get_domain_boundaries(
        seqs,
        profile,
        results_file=None,
        envelope=False,
        threshold=0
):
    """ Return a dictionary mapping sequence ID to profile boundaries as (start_position, end_position). The domain
     (default) or envelope boundary can be used. A bit score threshold can be optionally applied.

     TODO: Currently supports case where profile file contains single pHMM only"""

    # Run hmmsearch and extract domain/envelope score and boundary information
    data = hmm_search_dom(
        profile,
        seqs,
        no_return=False,
        results_file=results_file,
        quiet=True,
        threshold=threshold,
        envelope=envelope)

    # Get dict {seq_id : (start, end)}
    bounds = {seq_id : info[1] for seq_id, info in data.items()}

    return bounds


def subseqs_by_profile(
        full_seq_file,
        subseq_file,
        profile,
        hmm_mode="dom",
        chain_overlap_tol=0.1,
        threshold=0,
        region_type="env",
        filt_full_seq_fa=None,
        max_cpu=2,
        nice=None
):
    """ Extract sub-sequences from sequences in fasta format based with boundaries informed by a pHMM. A bit score
    threshold can be applied for inclusion if desired. """

    full_results = hmm_search_custom(
        profile,
        full_seq_file,
        mode=hmm_mode,
        region_type=region_type,
        chain_overlap_tol=chain_overlap_tol,
        max_on=False,
        cpu_per_search=max_cpu,
        nice=nice
    )

    # Only a single profile - just need {seq : hit_info}
    full_results = list(full_results.items())[0][1]

    filt_seqs = []
    filt_coords = []

    for seq, hit_info in full_results.items():
        if hit_info[0] >= threshold:
            filt_seqs.append(seq)
            filt_coords.append(hit_info[1])

    file_util.extract_fasta(
        full_seq_file,
        subseq_file,
        filt_seqs,
        coords=filt_coords
    )

    if filt_full_seq_fa:
        file_util.extract_fasta(
            full_seq_file,
            filt_full_seq_fa,
            filt_seqs
        )
