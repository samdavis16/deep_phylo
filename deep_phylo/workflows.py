import glob
import math
import os
import subprocess
from pathlib import Path

from Bio import SeqIO

from . import hmm
from . import phylo_partition
from . import aln
from . import curation
from . import annots
from . import aln_processing
from . import file_util
from . import tree
from . import cluster

""" Workflows for automating complex tasks making use of various PhyloKit modules. """


#################################################  SEQUENCE ANNOTATIONS   ##############################################

def generate_all_annots(
        in_fasta,
        out_annot_file=None,
        raw_uparc_file=None,
        raw_uprot_file=None,
        tax_file=None,
        full_sprot_map_file=None,
        max_threads=16
):
    """ Workflow for generation of all available annotation files from an input collection of Uniparc IDs. """

    # Naming for output .annot and .itol files
    if not out_annot_file:
        out_annot_file = in_fasta.split(".")[0] + ".annot"
    out_itol_file = out_annot_file.split(".")[0] + ".itol"

    # Naming for Uniparc and Uniprot raw annotation files
    # NOTE: if others databases are added, they should be put here as well
    if not raw_uparc_file:
        raw_uparc_file = in_fasta.split(".")[0] + "_uparc.txt"
    if not raw_uprot_file:
        raw_uprot_file = in_fasta.split(".")[0] + "_uprot.txt"

    # Naming for taxonomy files
    if not tax_file:
        tax_file = in_fasta.split(".")[0] + "_tax.txt"

    # Retrieve raw Uniparc and Uniprot annotation files
    seq_ids = [seq.name for seq in SeqIO.parse(in_fasta, "fasta")]
    annots.retrieve_all_annots(
        seq_ids,
        raw_uparc_file,
        raw_uprot_file,
        max_threads=max_threads
    )

    # Map Uniparc IDs to annotations from associated Uniprot entries
    parc2prot_map_file = raw_uparc_file.split(".txt")[0] + "_parc2prot.tsv"
    annots.map_and_extract_annots(
        parc2prot_map_file,
        raw_uprot_file,
        annots.UP_ANNOTS,
        out_file=out_annot_file,
        full_sprot_map_file=full_sprot_map_file,
        sprot_prefix="sp",
        no_return=True,
        max_threads=max_threads
    )

    # Fetch NCBI taxonomy annotations for Uniparc sequences, extract rank annotations and merge to .annot file
    annots.fetch_uparc_taxonomy(raw_uparc_file, tax_file)
    annots.extract_tax_lineage(
        tax_file,
        out_file=out_annot_file,
        annot_merge=True
    )

    # # Get GTDB annotations mapped to Uniparc IDs (where entry is associated with UProt proteome and NCBI assembly)
    # gtdb_annots = annots.uparc_to_gtdb_annots(raw_uparc_file)
    # annots.merge_annots(annot_dicts=[gtdb_annots],  # Merge to existing .annot file
    #                     annot_files=out_annot_file,
    #                     out_file=out_annot_file)

    # Create equivalent itol metadata file
    annots.create_itol_metadata(out_itol_file, out_annot_file)



#################################################  ALIGNMENT QUALITY   #################################################

def aln_quality_filter(in_aln,
                       out_aln,
                       valid_chars="ACDEFGHIKLMNPQRSTVWY-",
                       trim_retention_gt=0.1,
                       trim_retention_pc_cutoff=0.01,
                       trim_gap_prop_gt=0.1,
                       trim_gap_prop_pc_cutoff=0.01
                       ):
    """
    Default alignment quality control pipeline for use during iterations of large-phylogenetic-scope evolutionary
    analyses.
    NOTE: Not all functionality of each filtering function are exposed here.

    :param in_aln:
    :param out_aln:
    :param valid_chars: Iterable of allowable characters (should include gap character)
    :param trim_retention_gt: Gap threshold for column trimming in trim retention filtering
    :param trim_retention_pc_cutoff: Proportion of worst scoring sequences to exclude during trim retention filtering
    :param trim_gap_prop_gt: Gap threshold for column trimming in trim gap proportion filtering
    :param trim_gap_prop_pc_cutoff: Proportion of worst scoring sequences to exclude during gap proportion filtering

    :return: None
    """

    to_remove = set()

    # Non-standard character filtering
    for seq in SeqIO.parse(in_aln, "fasta"):
        for char in seq:
            if char not in valid_chars:
                to_remove.add(seq.name)
                break

    # Trim retention filtering
    # Currently only filters sequences by worst x% of sequences. To implement more sophisticated cutoff determination
    tr_scores = aln_processing.trim_retention(in_aln, gt=trim_retention_gt)
    tr_scores = [(seq, score) for seq, score in tr_scores.items()]
    tr_scores.sort(key=lambda x:x[1])   # Order worst to best
    tr_rank_cutoff = int(trim_retention_pc_cutoff * len(tr_scores))
    for i in range(len(tr_scores)):
        if i >= tr_rank_cutoff:
            break
        to_remove.add(tr_scores[i][0])

    # Trim gap proportion filtering
    tgp_scores = aln_processing.trim_gap_pc(in_aln, gt=trim_gap_prop_gt)
    tgp_scores = [(seq, score) for seq, score in tgp_scores.items()]
    tgp_scores.sort(reverse=True, key=lambda x:x[1])   # Order worst to best
    tgp_rank_cutoff = int(trim_gap_prop_pc_cutoff * len(tgp_scores))
    for i in range(len(tr_scores)):
        if i >= tgp_rank_cutoff:
            break
        to_remove.add(tgp_scores[i][0])

    # Filter sequences and remove gap only columns
    retain_seqs = [seq.name for seq in SeqIO.parse(in_aln, "fasta") if seq.name not in to_remove]
    aln.sub_aln(in_aln, out_aln, target_seqs=retain_seqs)


#################################################  CLADE EXPANSION   ##################################################

def profiles_to_tree(
        profiles,
        thresholds,
        full_hit_file,
        dom_hit_file,
        seg_files,
        iter_log_file=None,
        uparc_file_path=None,
        segment_idx_suffix=None,
        hmm_mode="dom",
        region_type="env",
        chain_overlap_tol=0.1,
        skip_curation=False,
        max_on=False,
        annotations=None,
        full_sprot_map_file=None,
        full_uprot_map_file=None,
        clust_min_seq_ids=None,
        clust_min_seq_covs=None,
        linclust=True,
        max_rep_seqs=2500,
        rep_level=None,
        max_seed_seqs=750,
        seed_aln_level=None,
        init_aln_mode="default",  # TODO: Test feasibility of less heuristic aligners,
        re_aln_mode="dash_cluster",
        valid_chars="ACDEFGHIKLMNPQRSTVWY-",
        trim_retention_gt=0.1,
        trim_retention_pc_cutoff=0.01,
        trim_gap_prop_gt=0.1,
        trim_gap_prop_pc_cutoff=0.01,
        aln_trim_mode="trimal",
        aln_trim_gt=0.1,
        tree_inference_mode="fasttree",
        max_download_threads=8, max_work_threads=48, nice=None
):

    """
    Beginning with one or more profile HMMs and associated thresholds, run curation and all steps of the representation
    expansion workflow up to and including tree inference.

    :param profiles:
    :param thresholds:
    :param full_hit_file:
    :param dom_hit_file:
    :param seg_files:
    :param iter_log_file:
    :param uparc_file_path:
    :param hmm_mode:
    :param region_type:
    :param chain_overlap_tol:
    :param skip_curation:
    :param max_on:
    :param annotations:
    :param full_sprot_map_file:
    :param full_uprot_map_file:
    :param clust_min_seq_ids:
    :param clust_min_seq_covs:
    :param linclust:
    :param max_rep_seqs:
    :param rep_level: HC level from which to extract representatives for subsequent analysis (overrides max_rep_seqs)
    :param max_seed_seqs:
    :param seed_aln_level: HC level at which to perform the structure-guided seed alignment (overrides max_seed_seqs)
    :param init_aln_mode:
    :param valid_chars: Iterable of allowable characters (should include gap character)
    :param trim_retention_gt: Gap threshold for column trimming in trim retention filtering
    :param trim_retention_pc_cutoff: Proportion of worst scoring sequences to exclude during trim retention filtering
    :param trim_gap_prop_gt: Gap threshold for column trimming in trim gap proportion filtering
    :param trim_gap_prop_pc_cutoff: Proportion of worst scoring sequences to exclude during gap proportion filtering
    :param aln_trim_mode:
    :param aln_trim_gt: Gap threshold for trimming refined alignment when TrimAl is selected as trimming mode
    :param tree_inference_mode:
    :param max_download_threads:
    :param max_work_threads:
    :param nice:
    :return: None

    """

    if annotations not in ["full", "reps", None]:
        raise ValueError(f"Annotations mode {annotations} is not implemented. "
                         f"Use option 'full' for to generate for all hits, or "
                         f"'reps' to generate only for sequences in the final tree.")

    if isinstance(profiles, str):
        profiles = [profiles]
    if isinstance(thresholds, str):
        thresholds = [thresholds]

    # Default HC level thresholds
    if not (clust_min_seq_ids and clust_min_seq_covs):
        clust_min_seq_ids = [0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3],
        clust_min_seq_covs = [0.95, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7]

    # Sequence curation
    if skip_curation:  # Ensure domain hit file is provided - if only a full hit file, set dom_hit_file as full hit file
        dir_files = os.listdir()
        if dom_hit_file not in dir_files:
            if full_hit_file in dir_files:
                dom_hit_file = full_hit_file
            else:
                raise RuntimeError("Existing sequence hits must be provided if the skip_curation flag is set.")

    else:

        # TODO: Quick fix - Need to refactor hmm_mode options in cluster_phylo_test + hmm and other downstream modules
        if hmm_mode == "dom":
            best_domain = True
            best_chain = False
        elif hmm_mode == "chain":
            best_domain = False
            best_chain = True
        else:
            raise RuntimeError(f"hmm search mode {hmm_mode} is not implemented.")


        # TODO: In process of switching over to single call point for segmented db curation

        cpu_per_search = max_work_threads // max_download_threads

        if iter_log_file:
            with open(iter_log_file, 'a') as log_f:
                log_f.write(f"Beginning curation with profiles and thresholds:"
                            f"\n{str(profiles)}\n{str(thresholds)}\n\n")

        curation.profile_search_segmented_db(
            profiles,
            thresholds,
            seg_files,
            full_hit_file,
            idx_suffix=segment_idx_suffix,
            dom_hit_file=dom_hit_file,
            hmm_mode=hmm_mode,
            region_type=region_type,
            chain_overlap_tol=chain_overlap_tol,
            remove_seg_files=False,
            max_on=max_on,
            cpu_per_search=cpu_per_search,
            max_search_processes=max_download_threads,
            nice=nice
        )

    hit_names = [seq.name for seq in
                 SeqIO.parse(full_hit_file, "fasta")]

        # # Curate with specified profiles and thresholds
        # cluster_phylo_test.profile_uniparc_search_new(profiles,
        #                                               full_hit_file,
        #                                               thresholds,
        #                                               dom_hit_file=dom_hit_file,
        #                                               uparc_file_path=uparc_file_path,
        #                                               best_domain=best_domain,
        #                                               max_download_threads=max_download_threads,
        #                                               max_hmm_threads=max_work_threads,
        #                                               nice=nice,
        #                                               best_chain=best_chain)

    # # TODO: This is a temporary fix to non-chains searches not producing a domain file
    # TODO: This was causing re-expansion of domain hits back to full hits. Unsure if removing will break something else
    # if not best_chain:
    #     dom_hit_file = full_hit_file

    if annotations == "full":

        # Generate all relevant annotation files
        generate_all_annots(
            full_hit_file,
            full_sprot_map_file=full_sprot_map_file
        )

        # Get priority sets for hierarchical clustering
        priority = annots.map_up_priority(
            uprot_annot_file=full_hit_file.split(".")[0] + "_uprot.txt",
            prot2parc_map=full_hit_file.split(".")[0] + "_uparc_prot2parc.tsv",
            full_sprot_map_file=full_sprot_map_file,
            full_uprot_map_file=full_uprot_map_file,
            include=hit_names
        )

    elif (full_sprot_map_file and full_uprot_map_file):
        # Can still generate priority sets for clustering
        priority = annots.map_up_priority(
            full_sprot_map_file=full_sprot_map_file,
            full_uprot_map_file=full_uprot_map_file,
            include=hit_names
        )

    else:
        # No information provided for priority clustering sets
        priority = None  # Annotations not available to determine priority sequence sets

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning hierarchical clustering with domain hits: {dom_hit_file}\n\n")

    # Run hierarchical clustering
    cluster.hc_from_fasta(dom_hit_file,
                          min_seq_ids=clust_min_seq_ids,
                          min_coverages=clust_min_seq_covs,
                          linclust=linclust,
                          priority=priority,
                          max_threads=max_work_threads,
                          nice=nice)

    # Determine the sequence identity threshold to extract representatives at
    hc_name = dom_hit_file.split(".")[0]
    hc_level_info = cluster.get_hc_level_info(name=hc_name)
    hc_maps = cluster.get_hc_maps(name=hc_name)

    if rep_level:
        rep_fa = hc_level_info[rep_level][0].split('.')[0] + "_reps.fa"
    else:
        rep_fa = None
        clust_sizes = [len(cluster_reps) for cluster_reps in hc_maps]
        for i in range(len(clust_sizes)):
            if clust_sizes[i] <= max_rep_seqs:
                rep_level = i
                rep_fa = hc_level_info[i][0].split('.')[0] + "_reps.fa"
                break
        if not rep_fa:  # All clusterings larger than max_rep_seqs --> use lowest level clustering
            rep_level = len(hc_maps) - 1
            rep_fa = hc_level_info[-1][0].split('.')[0] + "_reps.fa"

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Extracting representatives at level: {rep_fa.split('.')[0]}\n\n")

    # Extract cluster representative sequences
    clust_db_file = glob.glob(f"{rep_fa.split('_reps.fa')[0]}*DB")[0]  # DB file name may also contain sample number
    cluster.extract_cluster_reps(dom_hit_file.split('.')[0] + "DB",
                                 clust_db_file,
                                 fasta_name=rep_fa,
                                 priority_resample=priority)

    # Run initial MSA  # TODO: Currently just runs MAFFT with default params
    init_aln_name = rep_fa.split('.')[0] + "_init.aln"

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning initial alignment with mode {init_aln_mode} --> {init_aln_name}\n\n")

    if init_aln_mode == "default":
        aln.default_aln(rep_fa,
                        init_aln_name,
                        quiet=True,
                        nice=nice,
                        max_threads=max_work_threads)

    else:
        raise RuntimeError(f"Alignment mode {init_aln_mode} is not implemented.")

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning QC on initial alignment.\n\n")

    # Run automated alignment QC
    filt_aln_name = init_aln_name.split(".aln")[0] + "_filt.aln"
    filt_fa_name = filt_aln_name.split(".aln")[0] + ".fa"
    aln_quality_filter(init_aln_name,
                       filt_aln_name,
                       valid_chars=valid_chars,
                       trim_retention_gt=trim_retention_gt,
                       trim_retention_pc_cutoff=trim_retention_pc_cutoff,
                       trim_gap_prop_gt=trim_gap_prop_gt,
                       trim_gap_prop_pc_cutoff=trim_gap_prop_pc_cutoff)
    file_util.ungap_fasta(filt_aln_name, filt_fa_name)  # Ungap


    # Re-align
    re_aln_name = init_aln_name.split("_init")[0] + "_filt_realn.aln"
    seed_aln_name = re_aln_name.split("_realn")[0] + "_seed.aln"

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning re-alignment of filtered representatives --> {re_aln_name}\n\n")

    if re_aln_mode == "dash_cluster":

        clust_sizes = [len(cluster_reps) for cluster_reps in hc_maps]

        # Determine the sequence identity threshold to build seed aln at if not specified
        if not seed_aln_level:
            for i in range(rep_level, len(clust_sizes)):
                if clust_sizes[i] <= max_seed_seqs:
                    seed_aln_level = i
                    break
            if not seed_aln_level:
                seed_aln_level = len(hc_maps) - 1

        # If all clusters are too large for MAFF-DASH (>3000 reps) --> just use default aln
        if clust_sizes[seed_aln_level] > 3000:
            aln.default_aln(
                filt_fa_name,
                re_aln_name,
                quiet=True,
                nice=nice,
                max_threads=max_work_threads)

        else:  # Otherwise good to go with dash_cluster alignment

            # List of names for seqs in previously filtered alignment
            retained_seqs = [seq.name for seq in SeqIO.parse(filt_fa_name, "fasta")]

            cluster.dash_seed_cluster_aln(
                filt_fa_name,
                re_aln_name,
                hc_name=hc_name,
                rep_aln_name=seed_aln_name,
                lower_rep_idx=seed_aln_level,
                upper_rep_idx=rep_level,
                include_seqs=retained_seqs,
                max_processes=max_work_threads,
                nice=nice,
                quiet=True)

    elif re_aln_mode == "default":
        aln.default_aln(
            filt_fa_name,
            re_aln_name,
            quiet=True,
            nice=nice,
            max_threads=max_work_threads)

    else:
        raise RuntimeError(f"Alignment mode {re_aln_mode} is not implemented.")

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning alignment trimming with mode {aln_trim_mode}.")

    # Alignment trimming
    if aln_trim_mode == "trimal":
        trimmed_aln_name = re_aln_name.split(".")[0] + f"_t{int(100*aln_trim_gt)}.aln"
        aln.run_trimal(re_aln_name, trimmed_aln_name, gt=aln_trim_gt)

    else:
        raise RuntimeError(f"Alignment trimming mode {aln_trim_mode} is not implemented.")

    # Construct phylogeny
    tree_name = re_aln_name.split('.')[0] + ".nwk"

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning tree inference with mode {tree_inference_mode} --> {tree_name}.\n\n")

    if tree_inference_mode == "fasttree":
        tree.run_fasttree(re_aln_name,
                          tree_name,
                          nice=nice,
                          log_file=re_aln_name.split('.')[0] + ".log")

    else:
        raise RuntimeError(f"Tree inference mode {tree_inference_mode} is not implemented.")

    if annotations == "reps":
        generate_all_annots(
            filt_fa_name,
            full_sprot_map_file=full_sprot_map_file
        )

    # For passing to next module if a full iteration is run
    return tree_name, re_aln_name


def tree_to_profiles(
        in_tree_file,
        aln_file,
        iter_log_file=None,
        out_dir=None,
        rooted_tree_file=None,
        branch_support_file=None,
        sup_tree_file=None,
        reroot_strategy="midpoint",
        outgroup_leaves=None,
        ingroup_leaf=None,
        hmm_mode_cs="dom",
        max_on=False,
        chain_overlap_tol=0.1,
        min_split_support=90,
        min_clade_prop=0,
        max_clade_prop=0.5,
        max_holdout_prop=0.5,
        realign_target_clade=False,
        threshold_cluster_id=None,
        threshold_cluster_cov=None,
        linclust=False,
        cas_pass_score=0.9,
        min_cas_pass=0.75,
        profile_prefix="",
        beat_seq_prop=0.2,
        hmm_mode_thresh="dom",
        max_processes=8,
        threads_per_search=6,
        nice=None):
    """
    From a phylogeny (possibly generated by the complementary profiles_to_tree sub-module), identify robust
    clade representations and determine threshold for corresponding profiles for use in subsequent curation iterations.

    :param in_tree_file:
    :param aln_file:
    :param iter_log_file:
    :param out_dir: All output files can optionally be placed in a sub-directory
    :param rooted_tree_file:
    :param sup_tree_file: Path for writing the scaffold subtree containing only nodes deemed 'supported' by pipeline
    :param reroot_strategy: Strategy for re-rooting input tree or None if not required
    :param outgroup_leaves: Pair of leaves from either side of outgroup for use if outgroup re-rooting is specified
    :param ingroup_leaf: Leaf from ingroup for use if outgroup re-rooting is specified
    :param hmm_mode_cs:
    :param max_on:
    :param chain_overlap_tol:
    :param min_split_support: Minimum (percentage) support for consideration as a clade for constructing representation
    :param min_clade_prop: Minimum proportion of total tree for consideration as a clade for constructing representation
    :param max_clade_prop: Maximum proportion of total tree for consideration as a clade for constructing representation
    :param max_holdout_prop:
    :param realign_target_clade:
    :param threshold_cluster_id:
    :param threshold_cluster_cov:
    :param linclust:
    :param cas_pass_score:
    :param min_cas_pass:
    :param profile_prefix:
    :param beat_seq_prop:
    :param max_processes:
    :param threads_per_search:
    :param nice:

    :return:
    """

    # Load and process tree for partitioning (re-rooting, construction of supported sub-tree)
    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Loading tree {in_tree_file} for robust clade assessment.\n\n")
    in_tree = tree.load_tree(in_tree_file)

    # TODO: File convention changes will require significant refactor - to do later
    out_dir = Path(out_dir) if out_dir else Path.cwd()

    if not rooted_tree_file:
        rooted_tree_file = in_tree_file.split(".")[0] + "_rooted.nwk"
        # rooted_tree_file = out_dir / rooted_tree_file
    if not branch_support_file:
        branch_support_file = rooted_tree_file.split(".")[0] + "_supports.tsv"
        # branch_support_file = out_dir / branch_support_file
    if not sup_tree_file:
        sup_tree_file = rooted_tree_file.split(".")[0] + "_sup.nwk"
        # sup_tree_file = out_dir / sup_tree_file

    if not profile_prefix:
        profile_prefix = in_tree_file.split(".nwk")[0]

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning preliminary identification of supported clade candidates.\n\n")

    sup_tree, rooted_tree = phylo_partition.partition_internal_nodes(
        in_tree,
        min_support=min_split_support,
        reroot_strategy=reroot_strategy,
        out_nodes=outgroup_leaves,
        in_node=ingroup_leaf,
        min_prop=min_clade_prop,
        max_prop=max_clade_prop,
        full_out_file=rooted_tree_file,
        sup_tree_file=sup_tree_file,
        branch_support_file=branch_support_file)

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning robustness assessment for clades in supported subtree.\n\n")

    final_clades = phylo_partition.sup_tree_cladestrap(
        sup_tree,
        rooted_tree,
        aln_file,
        max_clade_proportion=max_holdout_prop,
        realign_target_clade=realign_target_clade,
        hmm_mode=hmm_mode_cs,
        max_on=max_on,
        chain_overlap_tol=chain_overlap_tol,
        thresh_cluster_id=threshold_cluster_id,
        thresh_cluster_cov=threshold_cluster_cov,
        linclust=linclust,
        cas_pass_score=cas_pass_score,
        min_cas_pass=min_cas_pass,
        max_processes=max_processes,
        threads_per_search=threads_per_search,
        nice=nice)

    # Construct pHMMs for supported clades which also pass CAS screening and determine thresholds

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Beginning threshold determination for final clades:\n{str(list(final_clades.keys()))}\n\n")

    profile_thresholds = {}
    final_clade_annots = {}
    for int_node, leaves in final_clades.items():

        sub_aln_file_name = f"{profile_prefix}_{int_node}.aln"
        hmm_file_name = f"{profile_prefix}_{int_node}.hmm"

        hmm.build_clade_hmm(
            aln_file,
            sub_aln_file_name,
            hmm_file_name,
            rooted_tree,
            target_clade_node=int_node,
            include=leaves)

        profile_thresholds[hmm_file_name] = hmm.threshold_clade_hmm(
            hmm_file_name,
            rooted_tree,
            in_fasta=aln_file,
            beat_seq_prop=beat_seq_prop,
            target_clade_node=int_node,
            hmm_mode=hmm_mode_thresh,
            max_on=max_on,
            chain_overlap_tol=chain_overlap_tol,
            max_cpu=threads_per_search,
            nice=nice)

        # Add leaves for this clade to annot dict
        for leaf in leaves:
            final_clade_annots[leaf] = {"pHMM_clades" : int_node}

    # Write annotation files for final clades
    annot_file_name = rooted_tree_file.rsplit('.', maxsplit=1)[0] + ".annot"
    itol_file_name = annot_file_name.rsplit('.', maxsplit=1)[0] + ".itol"
    annots.annot_file_from_dict(annot_file_name, final_clade_annots)
    annots.create_itol_metadata(itol_file_name, annot_file_name)

    if iter_log_file:
        with open(iter_log_file, 'a') as log_f:
            log_f.write(f"Iteration complete with final profiles/thresholds:\n{str(profile_thresholds)}\n")

    return profile_thresholds


def full_expansion_iter(
        profiles,
        thresholds,
        full_hit_file,
        dom_hit_file,
        seg_files,
        iter_log_file=None,
        out_dir=None,
        uparc_file_path=None,
        segment_idx_suffix=None,
        hmm_mode_search="dom",
        hmm_mode_cs="dom",
        region_type="env",
        max_on=False,
        annotations=None,
        full_sprot_map_file=None,
        full_uprot_map_file=None,
        clust_min_seq_ids=None,
        clust_min_seq_covs=None,
        linclust=True,
        max_rep_seqs=2500,
        rep_level=None,
        max_seed_seqs=750,
        seed_aln_level=None,
        init_aln_mode="default",
        re_aln_mode="dash_cluster",
        valid_chars="ACDEFGHIKLMNPQRSTVWY-",
        trim_retention_gt=0.1,
        trim_retention_pc_cutoff=0.01,
        trim_gap_prop_gt=0.1,
        trim_gap_prop_pc_cutoff=0.01,
        aln_trim_mode="trimal",
        aln_trim_gt=0.1,
        tree_inference_mode="fasttree",
        chain_overlap_tol=0.1,
        min_split_support=90,
        min_clade_prop=0,
        max_clade_prop=0.5,
        max_holdout_prop=0.5,
        realign_target_clade=False,
        threshold_cluster_id=None,
        threshold_cluster_cov=None,
        cas_pass_score=0.9,
        min_cas_pass=0.75,
        profile_prefix="",
        beat_seq_prop=0.2,
        max_download_threads=8,
        max_work_threads=48,
        nice=None
):
    """
    Runs a single, full iteration of the clade expansion workflow, optionally skipping annotation generation.

    :param profiles:
    :param thresholds:
    :param full_hit_file:
    :param dom_hit_file:
    :param out_dir:
    :param uparc_file_path:
    :param segment_idx_suffix:
    :param seg_files:
    :param iter_log_file:
    :param hmm_mode_search:
    :param hmm_mode_cs:
    :param region_type:
    :param max_on:
    :param annotations:
    :param full_sprot_map_file:
    :param full_uprot_map_file:
    :param clust_min_seq_ids:
    :param clust_min_seq_covs:
    :param linclust:
    :param max_rep_seqs:
    :param rep_level:
    :param max_seed_seqs:
    :param seed_aln_level:
    :param init_aln_mode:
    :param re_aln_mode:
    :param valid_chars:
    :param trim_retention_gt:
    :param trim_retention_pc_cutoff:
    :param trim_gap_prop_gt:
    :param trim_gap_prop_pc_cutoff:
    :param aln_trim_mode:
    :param aln_trim_gt:
    :param tree_inference_mode:
    :param chain_overlap_tol:
    :param min_split_support:
    :param min_clade_prop:
    :param max_clade_prop:
    :param max_holdout_prop:
    :param realign_target_clade:
    :param threshold_cluster_id:
    :param threshold_cluster_cov:
    :param cas_pass_score:
    :param min_cas_pass:
    :param profile_prefix:
    :param beat_seq_prop:
    :param max_download_threads:
    :param max_work_threads:
    :param nice:
    :return:
    """

    # Default HC level thresholds
    if not (clust_min_seq_ids and clust_min_seq_covs):
        clust_min_seq_ids = [0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3],
        clust_min_seq_covs = [0.95, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7]

    # Run profiles --> tree segment
    unrooted_tree, re_aln_file = profiles_to_tree(
        profiles=profiles,
        thresholds=thresholds,
        full_hit_file=full_hit_file,
        dom_hit_file=dom_hit_file,
        uparc_file_path=uparc_file_path,
        segment_idx_suffix=segment_idx_suffix,
        iter_log_file=iter_log_file,
        seg_files=seg_files,
        hmm_mode=hmm_mode_search,
        region_type=region_type,
        max_on=max_on,
        chain_overlap_tol=chain_overlap_tol,
        skip_curation=False,
        annotations=annotations,
        full_sprot_map_file=full_sprot_map_file,
        full_uprot_map_file=full_uprot_map_file,
        clust_min_seq_ids=clust_min_seq_ids,
        clust_min_seq_covs=clust_min_seq_covs,
        linclust=linclust,
        max_rep_seqs=max_rep_seqs,
        rep_level=rep_level,
        max_seed_seqs=max_seed_seqs,
        seed_aln_level=seed_aln_level,
        init_aln_mode=init_aln_mode,
        re_aln_mode=re_aln_mode,
        valid_chars=valid_chars,
        trim_retention_gt=trim_retention_gt,
        trim_retention_pc_cutoff=trim_retention_pc_cutoff,
        trim_gap_prop_gt=trim_gap_prop_gt,
        trim_gap_prop_pc_cutoff=trim_gap_prop_pc_cutoff,
        aln_trim_mode=aln_trim_mode,
        aln_trim_gt=aln_trim_gt,
        tree_inference_mode=tree_inference_mode,
        max_download_threads=max_download_threads,
        max_work_threads=max_work_threads,
        nice=nice)

    threads_per_search = max_work_threads // max_download_threads

    # Pass unrooted tree and alignment to tree --> profiles segment
    profile_thresholds = tree_to_profiles(
        in_tree_file=unrooted_tree,
        aln_file=re_aln_file,
        iter_log_file=iter_log_file,
        out_dir=out_dir,
        rooted_tree_file=None,
        branch_support_file=None,
        sup_tree_file=None,
        reroot_strategy="midpoint",
        hmm_mode_cs=hmm_mode_cs,
        max_on=max_on,
        chain_overlap_tol=chain_overlap_tol,
        min_split_support=min_split_support,
        min_clade_prop=min_clade_prop,
        max_clade_prop=max_clade_prop,
        max_holdout_prop=max_holdout_prop,
        realign_target_clade=realign_target_clade,
        threshold_cluster_id=threshold_cluster_id,
        threshold_cluster_cov=threshold_cluster_cov,
        linclust=linclust,
        cas_pass_score=cas_pass_score,
        min_cas_pass=min_cas_pass,
        profile_prefix=profile_prefix,
        beat_seq_prop=beat_seq_prop,
        hmm_mode_thresh=hmm_mode_search,
        max_processes=max_download_threads,
        threads_per_search=threads_per_search,
        nice=nice)

    return profile_thresholds


def initialise_search_space(
        profiles,
        init_hit_file,
        uparc_segs_path=None,
        max_search_processes=6,
        cpu_per_search=6,
        max_on=True,
        nice=None

):
    """ Generate an initial search space
     TODO: Currently only uses UniParc with domain search (all thresholds = 0, extracting full sequence) """

    max_hmm_threads = max_search_processes * cpu_per_search

    curation.profile_uniparc_search_new(
        profiles,
        init_hit_file,
        [0 for _ in range(len(profiles))],
        dom_hit_file=None,
        uparc_file_path=uparc_segs_path,
        best_chain=False,
        best_domain=True,
        max_download_threads=max_search_processes,
        max_hmm_threads=max_hmm_threads,
        nice=nice,
        max_on=max_on,
        remove_segments=True,
        region_type="env"
    )


def run_n_iters(
        init_profiles,
        init_thresholds,
        run_name,
        n_full_iters=5,
        convergence_prop=None,
        max_iters=None,
        run_db_init=False,
        init_db_segments=None,
        init_db_prefix=None,
        max_on_init=False,
        max_on_iters=False,
        n_init_segments=8,
        full_db_dir=None,
        segment_idx_suffix=None,
        annotations=None,
        full_sprot_map_file=None,
        full_uprot_map_file=None,
        hmm_mode_search="dom",
        hmm_mode_cs="dom",
        region_type="env",
        clust_min_seq_ids=None,
        clust_min_seq_covs=None,
        linclust=True,
        max_rep_seqs=2500,
        rep_level=None,
        max_seed_seqs=750,
        seed_aln_level=None,
        init_aln_mode="default",
        re_aln_mode="default",
        valid_chars="ACDEFGHIKLMNPQRSTVWY-",
        trim_retention_gt=0.1,
        trim_retention_pc_cutoff=0.01,
        trim_gap_prop_gt=0.1,
        trim_gap_prop_pc_cutoff=0.01,
        aln_trim_mode="trimal",
        aln_trim_gt=0.1,
        tree_inference_mode="fasttree",
        chain_overlap_tol=0.1,
        min_split_support=90,
        min_clade_prop=0.02,
        max_clade_prop=0.5,
        max_holdout_prop=0.5,
        realign_target_clade=False,
        threshold_cluster_id=0.4,
        threshold_cluster_cov=0.7,
        cas_pass_score=0.9,
        min_cas_pass=0.75,
        profile_prefix="",
        beat_seq_prop=0.2,
        max_search_processes=6,
        cpu_per_search=6,
        nice=None
):
    """ Interim main entry point for full clade expansion tool.
    TODO: Various downstream refactoring planned and will tidy this up after that. """

    run_log = run_name + ".runlog"

    with open(run_log, 'w') as run_log_f:
        run_log_f.write(f"Beginning run {run_name}\n\n")

    # Check if convergence params supplied
    if convergence_prop:

        if convergence_prop > 100 or convergence_prop <= 0:
            raise ValueError(f"Invalid convergence proportion: {convergence_prop}. "
                             f"Must be provided as a proportion or percentage.")
        elif convergence_prop > 1:  # Need to convert to proportion
            convergence_prop = convergence_prop / 100

        if max_iters:
            if max_iters >= 3:
                n_full_iters = max_iters  # Override default n_iters
                with open(run_log, 'a') as run_log_f:
                    run_log_f.write(f"Convergence proportion: {convergence_prop}\n")
                    run_log_f.write(f"Maximum iterations: {max_iters}\n\n")
            else:
                raise ValueError("Current implementation of convergence assessment is only compatible with "
                                 "max_iterations >= 3.")
        else:
            raise ValueError("If a convergence proportion criteria is set, a maximum number of iterations must also be "
                             f"elected.")

    else:
        with open(run_log, 'a') as run_log_f:
            run_log_f.write(f"# iterations: {n_full_iters}\n\n")

    this_iter = 0  # TODO: In future implement update to init_db after m iterations

    max_work_threads = max_search_processes * cpu_per_search

    if not init_db_prefix:
        init_db_prefix = run_name + "_initdb"
    full_init_hit_file = init_db_prefix + f"_it{this_iter}_all.fa"

    # Default clustering identity/coverage thresholds
    if not (clust_min_seq_ids and clust_min_seq_covs):
        clust_min_seq_ids = [0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
        clust_min_seq_covs = [0.95, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 0.6]

    # Skip search space initialisation if segments already available
    if init_db_segments:
        seg_files = init_db_segments

    # If explicitly stated not to run initdb
    elif not run_db_init:

        # TODO: Hard-coded for Uniparc only
        if full_db_dir:
            seg_files = []
            for i in range(1,201):
                file_name = f"uniparc_active_p{i}.fasta"
                if file_name not in os.listdir():
                    subprocess.run(["ln", "-s", f"{full_db_dir.rstrip('/')}/{file_name}", '.'])
                if segment_idx_suffix:  # TODO: Quick fix - to be tidied up
                    segment_idx_suffix = segment_idx_suffix.lstrip('.')
                    seg_idx_name = file_name.split('.')[0] + '.' + segment_idx_suffix
                    if seg_idx_name in os.listdir(full_db_dir.rstrip('/')):
                        subprocess.run(["ln", "-s", f"{full_db_dir.rstrip('/')}/{seg_idx_name}", '.'])

                seg_files.append(file_name)

        else:   # Need to download each segment
            download_threads = min(16, max_search_processes * cpu_per_search)
            curation.fetch_uniparc_fasta(
                unzip=True,
                max_threads=download_threads)
            seg_files = [f"uniparc_active_p{i}.fasta" for i in range(1,201)]
            # raise RuntimeError(f"Only Uniparc database currently implemented for full search per-iteraiton.")

    else:  # Run search space initialisation

        # NOTE: Not currently used

        # Run initialiation of search space searching over Uniparc
        initialise_search_space(
            init_profiles,
            full_init_hit_file,
            uparc_segs_path=full_db_dir,
            max_search_processes=max_search_processes,
            cpu_per_search=cpu_per_search,
            max_on=max_on_init,
            nice=nice
        )

        # Split segments
        seg_files = []
        if n_init_segments > 1:

            n_init_hits = 0
            for _ in SeqIO.parse(full_init_hit_file, "fasta"):
                n_init_hits += 1
            hits_per_seg = math.ceil(n_init_hits / n_init_segments)
            this_seg = 1
            this_hit = 0
            init_seg_file = full_init_hit_file.split("_all")[0] + f"_{this_seg}.fa"
            seg_f = open(init_seg_file, 'w')
            seg_files.append(init_seg_file)

            for seq in SeqIO.parse(full_init_hit_file, "fasta"):

                # Check if we need new segment
                if this_hit and this_hit % hits_per_seg == 0:
                    seg_f.close()
                    this_seg += 1
                    init_seg_file = full_init_hit_file.split("_all")[0] + f"_{this_seg}.fa"
                    seg_files.append(init_seg_file)
                    seg_f = open(init_seg_file, 'w')


                seg_f.write(f">{seq.id}\n{str(seq.seq)}\n")
                this_hit += 1

            seg_f.close()

        else:
            seg_files = [full_init_hit_file]

    per_iter_hits = {}

    # Set variables for first iteration
    current_profiles = init_profiles
    current_thresholds = init_thresholds

    # Run full expansion iterations
    for iter_i in range(n_full_iters):

        # Start logging file for this iter
        iter_log_file = f"{run_name}_it{iter_i+1}.iterlog"
        with open(iter_log_file, 'w') as log_f:
            pass

        with open(run_log, 'a') as run_log_f:
            run_log_f.write(f"Beginning full iteration {iter_i + 1}:\n")
            run_log_f.write(f"Profiles:\n{str(current_profiles)}\n")
            run_log_f.write(f"Thresholds:\n{str(current_thresholds)}\n")

        full_hit_file = run_name + f"_it{iter_i+1}_hits.fa"
        dom_hits_file = run_name + f"_it{iter_i+1}_dom_hits.fa"

        iter_results = full_expansion_iter(
            profiles=current_profiles,
            thresholds=current_thresholds,
            full_hit_file=full_hit_file,
            dom_hit_file=dom_hits_file,
            seg_files=seg_files,
            iter_log_file=iter_log_file,
            out_dir=None,
            uparc_file_path=None,
            segment_idx_suffix=segment_idx_suffix,
            hmm_mode_search=hmm_mode_search,
            hmm_mode_cs=hmm_mode_cs,
            region_type=region_type,
            max_on=max_on_iters,
            chain_overlap_tol=chain_overlap_tol,
            annotations=annotations,
            full_sprot_map_file=full_sprot_map_file,
            full_uprot_map_file=full_uprot_map_file,
            clust_min_seq_ids=clust_min_seq_ids,
            clust_min_seq_covs=clust_min_seq_covs,
            linclust=linclust,
            max_rep_seqs=max_rep_seqs,
            rep_level=rep_level,
            max_seed_seqs=max_seed_seqs,
            seed_aln_level=seed_aln_level,
            init_aln_mode=init_aln_mode,
            re_aln_mode=re_aln_mode,
            valid_chars=valid_chars,
            trim_retention_gt=trim_retention_gt,
            trim_retention_pc_cutoff=trim_retention_pc_cutoff,
            trim_gap_prop_gt=trim_gap_prop_gt,
            trim_gap_prop_pc_cutoff=trim_gap_prop_pc_cutoff,
            aln_trim_mode=aln_trim_mode,
            aln_trim_gt=aln_trim_gt,
            tree_inference_mode=tree_inference_mode,
            min_split_support=min_split_support,
            min_clade_prop=min_clade_prop,
            max_clade_prop=max_clade_prop,
            max_holdout_prop=max_holdout_prop,
            realign_target_clade=realign_target_clade,
            threshold_cluster_id=threshold_cluster_id,
            threshold_cluster_cov=threshold_cluster_cov,
            cas_pass_score=cas_pass_score,
            min_cas_pass=min_cas_pass,
            profile_prefix=profile_prefix,
            beat_seq_prop=beat_seq_prop,
            max_download_threads=max_search_processes,
            max_work_threads=max_work_threads,
            nice=nice
        )

        # Record hit accessions per iteration
        per_iter_hits[iter_i + 1] = set(
            (seq.name for seq in SeqIO.parse(
                dom_hits_file,
                "fasta")
             )
        )

        with open(run_log, 'a') as run_log_f:
            run_log_f.write(f"Total hits: {len(per_iter_hits[iter_i + 1])}\n")

        # Check convergence criteria if specified
        if convergence_prop and iter_i >= 2:

            # Collect union of hits from previous three iterations
            # Keys are counts from 1, not idxs
            prev_hits = set()
            for j in range(iter_i, iter_i-2, -1):
                prev_hit_file = run_name + f"_it{j}_dom_hits.fa"
                prev_hits.update(
                    set((seq.name for seq in
                         SeqIO.parse(prev_hit_file, "fasta")
                         ))
                )

            # Check if proportion of sequences not seen in previous
            # two iterations is < convergence proportion
            this_iter_hits = per_iter_hits[iter_i + 1]
            new_hits = this_iter_hits - prev_hits
            prop_new = len(new_hits) / len(this_iter_hits)
            with open(run_log, 'a') as run_log_f:
                run_log_f.write(f"New hit proportion: {prop_new}\n")

            if prop_new < convergence_prop:
                # Converged!
                with open(run_log, 'a') as run_log_f:
                    run_log_f.write(f"Completed iteration {iter_i+1}.\n\n")
                    run_log_f.write(f"Sequence space exploration converged after {iter_i+1} iterations.\n")

                return iter_results

        with open(run_log, 'a') as run_log_f:
            run_log_f.write(f"Completed iteration {iter_i+1}\n\n")

        # Profiles and thresholds for next iteration
        current_profiles, current_thresholds = zip(*iter_results.items())
        current_profiles = list(current_profiles)
        current_thresholds = list(current_thresholds)

    with open(run_log, 'a') as run_log_f:
        if convergence_prop:
            run_log_f.write(f"Run completed - exploration did not converge.\n")
        else:
            run_log_f.write(f"Run completed.")

    return iter_results
