from copy import deepcopy

import ete3

from . import hmm
from . import annots
from . import tree

""" A collection of methods and accessory functions for the (possibly partial) partitioning of phylogenetic trees. A 
 motivating use case is the curation/expansion workflow, implemented in this library. """


def longest_branch_root(t):
    """ Set the root on the longest branch of an input unrooted tree and return the rooted version of the tree.
     NOTE: this is not bisecting the branch to the selected outgroup - one of the branches from the new node will be
     0-length. """

    if isinstance(t, ete3.Tree):
        t = deepcopy(t)
    else:
        t = tree.load_tree(t)

    # Identify the longest branch
    longest_node = None
    max_len = 0
    for node in t.traverse("preorder"):
        if node.dist > max_len and not node.is_leaf():
            max_len = node.dist
            longest_node = node

    # Note: no bisection of this branch.
    t.set_outgroup(longest_node)
    return t


def partition_internal_nodes(
        t,
        min_support=90,
        reroot_strategy=None,
        return_as="nodes",
        out_nodes=None,
        in_node=None,
        min_prop=0,
        max_prop=1,
        sup_tree_file=None,
        full_out_file=None,
        branch_support_file=None
):
    """ Hierarchically partition a tree depth-first simply based on various criteria. Assumes that the input tree has
    branch supports annotated in the internal node label fields.
    """

    t = tree.load_tree(t)

    if reroot_strategy:  # TODO: Need to write a single rooting function which handles all modes in tree module

        if reroot_strategy == "longest":
            t = longest_branch_root(t)

        elif reroot_strategy == "lca":  # Identify outgroup by LCA of two extant labels
            if not (out_nodes and in_node):
                raise RuntimeError(f"An outgroup seq label for each side of the outgroup LCA and an ingroup seq are "
                                   f"required when rooting using the LCA method.")
            t = tree.outgroup_root(t, out_nodes, in_node)

        elif reroot_strategy == "midpoint":
            t = tree.midpoint_root(t, no_return=False)

        else:
            raise RuntimeError(f"Re-rooting strategy '{reroot_strategy}' is not implemented.")

    # Re-label internal nodes
    t = tree.relabel_internal_nodes(t,
                                    full_out_file,
                                    no_return=False)

    # Check if branch supports are in pct or decimal format, convert if required and store supports for writing
    pct_supports = tree.support_value_type(t)
    if branch_support_file or not pct_supports:
        int_supports = []
        for node in t.traverse(strategy="preorder"):
            if not node.is_leaf():
                if not pct_supports:
                    node.support = round(100 * node.support, 1)
                int_supports.append((node.name, node.support))

    # Write branch supports to file if specified
    if branch_support_file:
        int_supports.sort(key=lambda x: int(x[0][1:]))
        with open(branch_support_file, 'w') as bs_out_file:
            for label, support in int_supports:
                bs_out_file.write(f"{label}\t{support}\n")

    # Maintain hierarchy of supported nodes as Tree
    supported_subtree = ete3.Tree()

    for node in t.traverse(strategy="preorder"):

        # Check each internal node against each criteria and do not consider if all are not satisfied

        if node.support < min_support:  # Clade support criteria
            continue

        if min_prop != 0 or max_prop != 1:  # Clade's proportion of total sequences criteria
            if not min_prop <= len(node.get_leaves()) / len(t) <= max_prop:
                continue

        # Add this node to supported node tree if it passes all criteria
        found_anc = False
        for anc in node.get_ancestors():
            if anc.name in [n.name for n in supported_subtree.traverse()]:
                (supported_subtree & anc.name).add_child(name=node.name)
                found_anc = True
                break
        if not found_anc:  # New descendant directly off root
            supported_subtree.add_child(name=node.name)

        if sup_tree_file:
            supported_subtree.write(outfile=sup_tree_file, format=1)

    # Returns the supported subtree alongside the full tree with equivalent internal node names
    return supported_subtree, t


def annotate_supported_clades(
        supported_subtree: ete3.Tree,
        full_tree: ete3.Tree,
        out_format="annot",
        annot_prefix=None
):
    """ Output sets of extants sequences associated with given supported internal nodes. """

    clades = {}

    for sup_node in supported_subtree.traverse(strategy="preorder"):
        # Get node in full tree
        if sup_node.name:
            full_node = full_tree & sup_node.name
            clades[sup_node.name] = full_node.get_leaf_names()

    if out_format == "annot":  # Annot dictionary format i.e. annotating leaves by internal nodes
        annot_dict = {}
        for node, leaves in clades.items():
            for leaf in leaves:
                try:
                    annot_dict[leaf][f"{annot_prefix + '_' if annot_prefix else ''}{node}"] = "True"
                except KeyError:
                    annot_dict[leaf] = {f"{annot_prefix + '_' if annot_prefix else ''}{node}" : "True"}
        return annot_dict

    elif out_format == "dict":
        return clades

    else:
        raise RuntimeError(f"Output format {out_format} is not implemented.")


def annotate_supported_clades_max(
        supported_subtree: ete3.Tree,
        full_tree: ete3.Tree,
        out_format="annot",
        annot_label="max_sup_clades"
):
    """ Output the set of largest monophyletic clades in a given supported tree. """

    clades = {}

    for root_sup_node in supported_subtree.get_children():
        full_node = full_tree & root_sup_node.name
        clades[root_sup_node.name] = full_node.get_leaf_names()

    if out_format == "annot":
        annot_dict = {}
        for node, leaves in clades.items():
            for leaf in leaves:
                annot_dict[leaf] = {annot_label : node}
        return annot_dict

    elif out_format == "dict":
        return clades

    else:
        raise RuntimeError(f"Output format {out_format} is not implemented.")


def partition_and_annotate(
        t,
        annot_file_name,
        min_support=0.9,
        re_root=None,
        return_as="nodes",
        out_nodes=None,
        in_node=None,
        min_prop=0,
        max_prop=1,
        sup_tree_file=None,
        full_tree_labeled_file=None,
        sup_annot_prefix=None,
        max_annot_label="max_sup_clades"
):

    # Run partitioning and write internal node-labelled tree and supported trees if specified
    sup_tree, full_tree = partition_internal_nodes(
        t,
        min_support=min_support,
        reroot_strategy=re_root,
        return_as=return_as,
        out_nodes=out_nodes,
        in_node=in_node,
        min_prop=min_prop,
        max_prop=max_prop,
        sup_tree_file=sup_tree_file,
        full_out_file=full_tree_labeled_file)

    # Create annotations for all supported clades and the set of largest supported clades
    full_sup_annots = annotate_supported_clades(sup_tree, full_tree, out_format="annot", annot_prefix = sup_annot_prefix)
    max_sup_annots = annotate_supported_clades_max(sup_tree, full_tree, out_format="annot", annot_label=max_annot_label)
    # return full_sup_annots, max_sup_annots
    annots.merge_annots(annot_dicts=[full_sup_annots, max_sup_annots], out_file=annot_file_name)


def sup_tree_cladestrap(
        sup_tree,
        full_tree,
        aln_file,
        max_clade_proportion=0.5,
        realign_target_clade=False,
        hmm_mode="dom",
        max_on=False,
        chain_overlap_tol=0.1,
        thresh_cluster_id=None,
        thresh_cluster_cov=None,
        linclust=False,
        cas_pass_score=0.9,
        min_cas_pass=0.75,
        max_processes=8,
        threads_per_search=6,
        nice=None
):
    """
    From a sup_tree containing only nodes passing criteria for cladestrapping, cladestrap nodes depth-first. Descendant
    nodes (i.e. representing sub-clades) of those which pass clade association screening are not considered for further
    cladestrapping.

    :param sup_tree:
    :param full_tree: Rooted and labeled tree (such is returned after supported clade partitioning)
    :param aln_file:
    :param max_clade_proportion:
    :param realign_target_clade:
    :param hmm_mode:
    :param max_on:
    :param chain_overlap_tol:
    :param thresh_cluster_id:
    :param thresh_cluster_cov:
    :param cas_pass_score:
    :param min_cas_pass:
    :param max_processes:
    :param threads_per_search:
    :param nice:

    :return:
    """

    sup_tree = tree.load_tree(sup_tree, nwk_format=1)
    full_tree = tree.load_tree(full_tree)
    final_clades = {}  # Map node name to list of passing seqs for all clades passing CAS screening

    def assess_sup_clade(node : ete3.TreeNode):
        """ For recursive assessment of supported nodes. Where a parent in the supported sub-tree passes CAS screening,
        its descendants are no longer considered. If a node fails, it's children are each assessed, and so on until tips
        are reached for a given lineage. """

        # Run CAS for this node/clade. If it fails screening, None will be returned
        node_name = node.name
        passing_leaves = hmm.full_clade_assessment(
            full_tree,
            aln_file,
            use_int_labels=True,
            target_clade_node=node_name,
            clade_name=node_name,
            max_clade_proportion=max_clade_proportion,
            realign=realign_target_clade,
            hmm_mode=hmm_mode,
            max_on=max_on,
            chain_overlap_tol=chain_overlap_tol,
            thresh_cluster_id=thresh_cluster_id,
            thresh_cluster_cov=thresh_cluster_cov,
            linclust=linclust,
            cas_pass_score=cas_pass_score,
            min_cas_pass=min_cas_pass,
            max_processes=max_processes,
            threads_per_search=threads_per_search,
            nice=nice)

        if passing_leaves:  # Record this node and its passing seqs and terminate assessment of this lineage
            final_clades[node_name] = passing_leaves
        else:
            if not node.is_leaf():  # Recursively assess descendants if this node failed CAS screening
                for child_node in node.get_children():
                    assess_sup_clade(child_node)

    top_nodes = sup_tree.get_children()  # Nodes at top level of supported node hierarchy
    for node in top_nodes:  # Recursively assess all nodes underneath
        assess_sup_clade(node)

    return final_clades
