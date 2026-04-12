import os
import random
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO

from . import annots
from . import tree


""" Functionality for sequence-sequence searching with either BLAST or MMSeqs2. """


def top_easy_search(
        query_fa,
        target_fa,
        results_file=None
):
    """ Using the MMSeqs easy-search module, map each sequence in query_fa to the best hit (on basis of e-value) in
    target_fa. Returns a dictionary mapping each query sequence to [best_hit_name, e-val, prop_id, prop_cov] where the
    coverage is the minimum of query and target coverage. """

    if results_file:
        out_file = results_file
    else:
        out_file = f"tmp_{random.randint(int(1e8), int(1e9)-1)}.tsv"

    tmp_dir = f"tmp_{random.randint(int(1e8), int(1e9)-1)}"

    args = [
        "mmseqs",
        "easy-search",
        query_fa,
        target_fa,
        out_file,
        tmp_dir,
        "-v", '0',
        "--format-output",
        "query,target,pident,qcov,tcov,evalue"
    ]

    subprocess.run(args, check=True)

    # Process output file
    best_hits = {}
    with open(out_file) as result_f:
        for line in result_f:
            data = line.strip().split('\t')
            if data[0] not in best_hits:
                best_hits[data[0]] = [
                    data[1],
                    float(data[5]),
                    round(float(data[2]), 1),
                    round(min(float(data[3]), float(data[4])), 1)
                ]

    Path(out_file).unlink()
    shutil.rmtree(tmp_dir, ignore_errors=True)

    return best_hits


def blastp_single(
        query_fa,
        target_fa,
        query_name=None,
        threads=8,
        nice=10,
        fast=False,
        max_targets="all"
):

    if len(list(SeqIO.parse(query_fa, "fasta"))) > 1:
        if not query_name:
            raise RuntimeError("A target sequence name must be provided if more than one sequence is present in the "
                               "query file.")
        try:
            query_seq = [[seq for seq in SeqIO.parse(query_fa, "fasta") if seq.name == query_name][0]]
            SeqIO.write(query_seq, f"temp_blast_query.fa", "fasta")
            query_fa = f"temp_blast_query.fa"
        except IndexError:
            raise RuntimeError(f"Sequence {query_name} not found in {query_fa}.")

    if max_targets == "all":
        max_targets = 0
        for seq in SeqIO.parse(target_fa, "fasta"):
            max_targets += 1

    subprocess.run(["makeblastdb", "-in", target_fa, "-dbtype", "prot", "-out", "temp_blast"])

    # Perform full BLAST search (no score/e-value filters)
    format_str = "6 qseqid sseqid qlen pident length evalue"
    args = ["nice", "-n", str(nice), "blastp", "-query", query_fa, "-db", "temp_blast", "-out",
            "temp_blast_results.txt", "-max_target_seqs", str(max_targets),
            "-outfmt", format_str, "-evalue", "10000", "-num_threads", str(threads)]
    if fast:
        args.extend(['-task', 'blastp-fast'])
    subprocess.run(args)

    result_dict = {}  # Store results for each query/subject pair
    with open("temp_blast_results.txt") as result_file:
        for line in result_file:
            data = line.strip().split('\t')
            result_list = [float(data[5]), round(float(data[3])), round(float(data[4]) / float(data[2]) * 100)]
            result_dict[data[1]] = result_list

    # Delete all temporary files
    os.system("rm temp_blast*")

    return result_dict


def blast_best_hits(query_fa, subject_fa, threads=8, nice=10, fast=False):
    """ Find the best hit (by BLAST E-value) to each query sequence in the subject sequences. Output a dictionary
     mapping {query : [best_subject_id, e_val, %id, %coverage]}. """

    # Create temporary BLAST database TODO: Maybe want an option to retain these?
    subprocess.run(["makeblastdb", "-in", subject_fa, "-dbtype", "prot", "-out", "temp_blast"])
    # os.system(f"makeblastdb -in {subject_fa} -dbtype prot -out temp_blast")

    # Perform full BLAST search (no score/e-value filters)
    format_str = "6 qseqid sseqid qlen pident length evalue"
    args = ["nice", "-n", str(nice), "blastp", "-query", query_fa, "-db", "temp_blast", "-out", "temp_blast_results.txt",
            "-outfmt", format_str, "-evalue", "10000", "-num_threads", str(threads)]
    if fast:
        args.extend(['-task', 'blastp-fast'])
    subprocess.run(args)
    # os.system(f"blastp -query {query_fa} -db temp_blast -out temp_blast_results.txt -outfmt {format_str} -evalue 10000")

    full_dict = {}  # Store results for each query/subject pair
    with open("temp_blast_results.txt") as result_file:
        for line in result_file:
            data = line.strip().split('\t')
            result_list = [data[1], float(data[5]), round(float(data[3])), round(float(data[4])/float(data[2])*100)]
            try:
                full_dict[data[0]].append(result_list)
            except KeyError:
                full_dict[data[0]] = [result_list]

    # Delete all temporary files
    os.system("rm temp_blast*")

    # Get best hit for each query sequence (by minimum E-value)
    return {query : min(results, key=lambda x:x[1]) for query, results in full_dict.items()}


def pw_tree_mapping(
        seq_file_1,
        seq_file_2,
        tree_1=None,
        tree_2=None,
        clade_names_1=None,
        clade_names_2=None,
        clade_bounds_1=None,
        clade_bounds_2=None,
        annot_file_1=None,
        annot_file_2=None,
        mapping_annot_label="tree_mapping",
        clade_map_prefix="map",
        clade_label_prefix="",
        threads=8,
        nice=10,
        min_id=80,
        min_cov=80,
        no_return=True,
        fast=False):
    """ Get BLAST-based mapping between sequences in two trees. Optionally, annotate mappings for specifically selected
    clades in each tree to the closest hits in the other. """

    if not (annot_file_1 or annot_file_2 or no_return):
        raise RuntimeError("If no return is specified, annotation file name(s) must be provided.")

    if not clade_names_1:
        clade_names_1 = []
    if not clade_names_2:
        clade_names_2 = []
    if not clade_bounds_1:
        clade_bounds_1 = []
    if not clade_bounds_2:
        clade_bounds_2 = []

    if tree_1 and tree_2:
        tree_1 = tree.load_tree(tree_1)
        tree_2 = tree.load_tree(tree_2)

    if ((clade_names_1 and clade_bounds_1) and
        len(clade_names_1) != len(clade_bounds_1)) or ((clade_names_2 and clade_bounds_2) and
                                                       len(clade_names_2) != len(clade_bounds_2)):
        raise RuntimeError("Each provided clade name must correspond to a pair of extant labels representing its "
                           "boundaries.")

    # Define these sets and extract relevant results later for annotations
    clades_1 = [tree.get_subtree_leaves(tree_1, ext_nodes=extant_pair) for extant_pair in clade_bounds_1]
    clades_2 = [tree.get_subtree_leaves(tree_2, ext_nodes=extant_pair) for extant_pair in clade_bounds_2]

    # # Get full BLAST hits in both directions
    # t1_t2_results = blast_best_hits(seq_file_1, seq_file_2, threads=threads, nice=nice, fast=fast)  # t1 query, t2 target
    # t2_t1_results = blast_best_hits(seq_file_2, seq_file_1, threads=threads, nice=nice, fast=fast)  # t2 query, t1 target

    # Get best hit in both directions with MMSeqs easy-search
    t1_t2_results = top_easy_search(seq_file_1, seq_file_2)
    t2_t1_results = top_easy_search(seq_file_2, seq_file_1)

    # Extract clade-specific results
    t1_t2_mappings = [{} for clade in clade_names_1]
    t2_t1_mappings = [{} for clade in clade_names_2]

    # Dicts containing mappings and clade definitions - NOTE: annot_dict_1 will actually be annot file for tree2 and v.v
    annot_dict_1 = {}
    annot_dict_2 = {}

    for i in range(len(clades_1)):
        for leaf in clades_1[i]:
            if t1_t2_results[leaf][2] >= min_id and t1_t2_results[leaf][3] >= min_cov:
                t1_t2_mappings[i][leaf] = t1_t2_results[leaf][0]
            try:
                annot_dict_2[leaf][f"{clade_label_prefix + '_' if clade_label_prefix else ''}{clade_names_1[i]}"] = True
            except KeyError:
                annot_dict_2[leaf] = {f"{clade_label_prefix + '_' if clade_label_prefix else ''}{clade_names_1[i]}":True}

    for i in range(len(clades_2)):
        for leaf in clades_2[i]:
            if t2_t1_results[leaf][2] >= min_id and t2_t1_results[leaf][3] >= min_cov:
                t2_t1_mappings[i][leaf] = t2_t1_results[leaf][0]
            try:
                annot_dict_1[leaf][f"{clade_label_prefix + '_' if clade_label_prefix else ''}{clade_names_2[i]}"] = True
            except KeyError:
                annot_dict_1[leaf] = {f"{clade_label_prefix + '_' if clade_label_prefix else ''}{clade_names_2[i]}":True}

    # Create annotation files if provided
    if annot_file_1:  # Mapping of tree 1 leaves/clades to tree2 (annot file for tree2)

        # First, iterate over all results (regardless of presence in defined clade or not)
        for t1_seq, mapping_data in t1_t2_results.items():
            if mapping_data[2] >= min_id and mapping_data[3] >= min_cov:
                try:
                    annot_dict_1[mapping_data[0]][mapping_annot_label] = t1_seq
                except KeyError:
                    annot_dict_1[mapping_data[0]] = {mapping_annot_label : t1_seq}

        for i in range(len(clades_1)):  # For each defined clade in t1 map each hit to it's query in other tree
            for t1_seq, mapping in t1_t2_mappings[i].items():
                try:
                    annot_dict_1[mapping][f"{clade_map_prefix}_{clade_names_1[i]}"] = t1_seq
                except KeyError:
                    annot_dict_1[mapping] = {f"{clade_map_prefix}_{clade_names_1[i]}" : t1_seq}


    if annot_file_2:  # Mapping of tree 2 leaves/clades to tree1 (annot file for tree1)

        # First, iterate over all results (regardless of presence in defined clade or not)
        for t2_seq, mapping_data in t2_t1_results.items():
            if mapping_data[2] >= min_id and mapping_data[3] >= min_cov:
                try:
                    annot_dict_2[mapping_data[0]][mapping_annot_label] = t2_seq
                except KeyError:
                    annot_dict_2[mapping_data[0]] = {mapping_annot_label : t2_seq}

        for i in range(len(clades_2)):  # For each defined clade in t1 map each hit to it's query in other tree
            for t2_seq, mapping in t2_t1_mappings[i].items():
                try:
                    annot_dict_2[mapping][f"{clade_map_prefix}_{clade_names_2[i]}"] = t2_seq
                except KeyError:
                    annot_dict_2[mapping] = {f"{clade_map_prefix}_{clade_names_2[i]}" : t2_seq}


    if annot_file_1:
        annots.annot_file_from_dict(annot_file_1, annot_dict_1)

    if annot_file_2:
        annots.annot_file_from_dict(annot_file_2, annot_dict_2)

    if not no_return:
        return t1_t2_mappings, t2_t1_mappings




