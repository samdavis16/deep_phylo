import subprocess
from copy import deepcopy

import ete3
from ete3.parser.newick import NewickError

""" Collection of processes for dealing with trees and alignments. Many of these functions combine ete3 functions and 
 handle the problem of ete3 arbitrarily selecting a root on unrooted input trees. """


def load_tree(t, nwk_format=0):
    """ Load an ete3 Tree object from a Newick file. If input is already an ete3 Tree object, return the provided tree.
     """
    if isinstance(t, ete3.Tree):
        return t
    else:
        if nwk_format == 1:
            return ete3.Tree(t, format=1)  # format 1 allows internal node labels
        else:
            try:
                return ete3.Tree(t, format=0)  # format 0 assumes support values
            except NewickError:
                return ete3.Tree(t, format=1)  # User may not know tree has internal node labels


def has_internal_labels(t):
    """ Checks if a tree has an internal labels annotated. If a ete3.TreeNode object is provided, checks if the .name
    attribute is populated for any internal node. If a file is provided, checks if internal nodes are labelled in the
    Newick string (possibly with support values). """

    t = load_tree(t)  # This will add internal node labels in Newick string to .name attribute if t is a file

    labelled_int_nodes = sum([(True if node.name else False) for node in t.traverse(strategy="preorder")
                              if not node.is_leaf()])
    if labelled_int_nodes > 0:
        return True
    else:
        return False


def check_nn_nodes(t):
    """ Checks if a tree contains any non-numeric internal node labels. Returns True if any identified, False otherwise.
     """

    t = load_tree(t)

    # Need to ensure supports are attached to .support attribute rather than .name so that its safe to re-root
    non_numeric_labels = False  # If any internal labels are non-numeric, assume support labels aren't in .name
    for node in t.traverse(strategy="preorder"):
        if node.name and not node.is_leaf():
            try:
                float(node.name)
            except TypeError:
                non_numeric_labels = True
                break

    if non_numeric_labels:
        return True
    else:
        return False


def support_value_type(t):
    """ Checks the form of support values annotated to internal nodes. Returns 0 for decimal, 1 for percentage. """

    t = load_tree(t)

    # Check if branch supports are in pct or decimal format - if latter, convert to former
    pct_supports = False
    for node in t.traverse():
        if not node.is_leaf():
            if node.support:
                if float(node.support) > 1.0:
                    pct_supports = True
                    break

    if pct_supports:
        return 1
    else:
        return 0


def relabel_internal_nodes(in_tree, out_tree_file=None, no_return=True):
    """ For a tree which may have internal nodes annotated with split supports, transfer splits to the .support field
     of the ete3.TreeNode if necessary and optionally write the file (with internal labels but without supports) to a
     new Newick file. Internal nodes are then re-labeled N{i} with i ordered by preorder traversal. If a branch support
     file is specified, supports are written to a tsv mapped to the corresponding internal node name. """

    in_tree = load_tree(in_tree)
    # if branch_support_file:
    #     bs_out_file = open(branch_support_file, 'w')

    # Check if we need to transfer split supports to .support attribute
    int_labels = has_internal_labels(in_tree)
    non_numeric_labels = check_nn_nodes(in_tree)

    if int_labels and not non_numeric_labels:
        in_tree = transfer_supports(in_tree)

    i = 0
    for node in in_tree.traverse(strategy="preorder"):
        if node.is_leaf():  # TODO: This may be necessary for tree_partitioning
            node.support = 0
        else:
            node.name = f"N{i}"
            # if node.name == "":  # This may be a 0 branch length node or a node added artificially during of re-rooting
            #     node.support = 100
            # else:
            #     if pct_supports:
            #         node.support = float(node.name)
            #     else:
            #         node.support = float(node.name) * 100
            # node.name = f"N{i}"
            # if branch_support_file:
            #     bs_out_file.write(f"{node.name}\t{node.support}\n")
            i+=1

    # if branch_support_file:
    #     bs_out_file.close()

    if out_tree_file:
        in_tree.write(outfile=out_tree_file, format=1)
    if not no_return:
        return in_tree


def transfer_supports(t: ete3.TreeNode):
    """ Transfers stored values from .name attribute to .support attribute. Accepts only an ete3 TreeNode
    object (i.e. does not read in from tree files). """

    for node in t.traverse(strategy="preorder"):
        if not node.is_leaf():
            node.support = float(node.name)
            node.name = ""

    return t


"""
def get_node(label)
    ''' TODO: node-level equivalent of load_tree '''

"""


def root_tree(in_tree, out_file=None, mode="midpoint", og_bounds=None, in_leaf=None, og_node=None,
              out_leaf=None, in_format=0, out_format=0, no_copy=False):
    """ Return and optionally write a rooted tree to a file. The input tree is deep-copied before rooting unless no_copy
    is specified. """

    in_tree = load_tree(in_tree, nwk_format=in_format)
    if no_copy:
        rooted_tree = in_tree
    else:
        rooted_tree = deepcopy(in_tree)

    if mode == "midpoint":
        outgroup = rooted_tree.get_midpoint_outgroup()
        rooted_tree.set_outgroup(outgroup)

    elif mode == "outgroup":
        if og_node:  # Node or node name for outgroup is explicitly provided
            if not isinstance(og_node, ete3.TreeNode):  # Assume this is a node label
                og_node = rooted_tree & og_node
                rooted_tree.set_outgroup(og_node)
        elif og_bounds:
            if in_leaf:  # Optionally set preliminary root on in-group leaf to ensure intended outgroup is monophyletic
                if not isinstance(in_leaf, ete3.TreeNode):
                    in_leaf = rooted_tree & in_leaf
                rooted_tree.set_outgroup(in_leaf)
            if not isinstance(og_bounds[0], ete3.TreeNode):
                og_bounds = [rooted_tree & og_bounds[i] for i in range(len(og_bounds))]
            lca = og_bounds[0].get_common_ancestor(*og_bounds[1:])  # Get LCA of all outgroup bounds
            rooted_tree.set_outgroup(lca)
        else:
            raise RuntimeError(f"Either an outgroup nodes or outgroup bounds (as two extant leaves on either side of the "
                               f"split must be provided for outgroup rooting.")

    elif mode == "single":
        if not isinstance(out_leaf, ete3.TreeNode):
            out_leaf = rooted_tree & out_leaf
        rooted_tree.set_outgroup(out_leaf)

    else:
        raise RuntimeError(f"Rooting mode {mode} is not implemented.")

    if out_file:
        rooted_tree.write(format=out_format)

    return rooted_tree


def single_root(tree, out_node: str):
    """ Return a rooted tree with a single taxa as outgroup. Node provided as a label (str)

    TODO: now redundant but still some usages - switch to using the universal rooting function """

    # Check if tree is already loaded, else assume filename for nwk tree is provided
    t = load_tree(tree) if not isinstance(tree, ete3.Tree) else tree

    t.set_outgroup(t & out_node)
    return t


def outgroup_root(tree, out_nodes, in_node=None):
    """ Return a rooted tree object, with outgroup corresponding to the last common ancestor of out_nodes.
    If an in_node is provided, it can be set as the initial outgroup to prevent unexpected behaviour by ete3

    TODO: now redundant but still some usages - switch to using the universal rooting function"""

    # Check if tree is already loaded, else assume filename for nwk tree is provided
    t = load_tree(tree) if not isinstance(tree, ete3.Tree) else tree

    # If in_node provided, set as initial outgroup so actual outgroup is monophyletic in re-rooted tree
    if in_node:
        t.set_outgroup(t & in_node)

    out_nodes = [out_nodes] if not isinstance(out_nodes, list) else out_nodes

    # Need to get node instances using node label
    out_nodes = [t & label for label in out_nodes]

    lca = out_nodes[0].get_common_ancestor(*out_nodes[1:])
    t.set_outgroup(lca)

    return t


def midpoint_root(t, out_file=None, no_return=True):
    """ Performs midpoint rooting and returns the rooted tree or writes it to a Newick file. The input tree is deep
     copied before re-rooting.

     TODO: now redundant but still some usages - switch to using the universal rooting function"""

    unrooted = load_tree(t)
    rooted = deepcopy(unrooted)

    # Need to ensure supports are attached to .support attribute rather than .name so that its safe to re-root
    # If any internal labels are non-numeric, assume support labels aren't in .name
    non_numeric_labels = check_nn_nodes(rooted)
    int_labels = has_internal_labels(rooted)

    if int_labels and not non_numeric_labels:
        # Assume we need to transfer support vals on internal .name attributes to .support
        rooted = transfer_supports(rooted)

    outgroup = rooted.get_midpoint_outgroup()
    rooted.set_outgroup(outgroup)

    if out_file:
        rooted.write(outfile=out_file, format=0)  # TODO: May need to change output format here
    if not no_return:
        return rooted


def get_subtree(tree, int_node=None, clade_bounds=None):
    """ Return a pruned tree containing only nodes under a given ancestral node. If the ancestral node is not
    explicitly provided, infer it as the LCA of a collection of ext_nodes. """

    t = load_tree(tree) if not isinstance(tree, ete3.Tree) else tree

    t2 = deepcopy(t)  # Prevent pruning of original tree if part of iterative process

    if int_node:
        int_node = t2 & int_node if not isinstance(int_node, ete3.Tree) else t2 & int_node.name  # Need equivalent anc in t2
    else:
        if not clade_bounds:
            raise RuntimeError("Either a clade ancestor or a set of representative extants must be provided.")
        else:
            clade_bounds = [clade_bounds] if not isinstance(clade_bounds, list) else clade_bounds
            clade_bounds = [t2 & clade_bounds[i] for i in range(len(clade_bounds))]
            int_node = clade_bounds[0].get_common_ancestor(*clade_bounds[1:])
    t2.prune([node for node in int_node.iter_descendants()], preserve_branch_length=True)

    return t2


def get_subtree_leaves(tree, anc=None, ext_nodes=None):
    """ Return a list of extant node names for a subtree defined by a last common ancestor """
    subtree = get_subtree(tree, anc, ext_nodes)
    return subtree.get_leaf_names()


def get_subtree_boundaries(tree, outgroup_bounds=[], in_node=None, singleton_clades=True):
    """ For a given phylogenetic tree, get boundaries for all subtrees in the form of two extant nodes whose common
      ancestor represents the ancestral clade boundary. For single leaves representing a sister to a multi-leaf clade
      ('singleton' leaves), optionally record said singletons as a clade indicated by its ID (default=True). Optionally,
      re-root the tree using outgroup and ingroup node params. """

    t = load_tree(tree)

    if outgroup_bounds:  # Re-root if required
        t = outgroup_root(t, outgroup_bounds, in_node)

    clade_bounds = []

    # Traverse tree, recording clade boundaries if criteria are met
    for node in t.traverse():

        if not node.is_leaf():  # Internal node, so record boundary as one leaf from either side of bp
            clade_bounds.append([child.get_leaf_names()[0] for child in node.get_children()])

        else:
            if singleton_clades:
                if len(node.get_sisters()[0]) > 1:  # Nested singleton leaf - define as own clade
                    clade_bounds.append(node)

    return clade_bounds


# def annotate_clade_bounds(tree, outgroup_bounds=[], in_nodes=None, singleton_clades)


def cousin_clade_sets(t, in_nodes, root_on=None, include_target=False):
    """ For a monophyletic target clade in a phylogenetic tree, the first cousin clade is defined as the sister clade
     to the target clade. The second cousin clade is defined as the sister of (target + cousin clade 1), and so on
      until all sequences up to the outgroup are incorporate as cousin clades.  """

    t = load_tree(t)

    # Root if root_on is specified
    if root_on:
        if isinstance(root_on, str):
            t = single_root(t, root_on)
        elif isinstance(root_on, list):
            t = outgroup_root(t, root_on, in_node=in_nodes[0])

    cousin_clades = []

    # Get LCA of target clade
    in_node_1 = t & in_nodes[0]
    in_node_2 = t & in_nodes[1]
    lca = in_node_1.get_common_ancestor(in_node_2)
    if include_target:
        cousin_clades.append(lca.get_leaf_names())

    found_labels = set(lca.get_leaf_names())  # Keep track of names to differentiate between sister clades each time

    while True:
        lca = lca.up
        cousin_clades.append([leaf for leaf in lca.get_leaf_names() if leaf not in found_labels])
        found_labels.update(cousin_clades[-1])
        if lca.is_root():
            break

    return cousin_clades


def run_fasttree(in_aln, out_tree, nice=None, log_file=None, quiet=False):
    """

    :param in_aln:
    :param out_tree:
    :param nice:
    :param log_file:
    :return:
    """

    ft_args = []

    if nice:
        ft_args.extend(["nice", "-n", str(int(nice))])

    ft_args.append("FastTree")

    if log_file:
        ft_args.extend(["-log", log_file])

    if quiet:
        ft_args.append("-quiet")

    ft_args.append(in_aln)

    with open(out_tree, 'w') as out_file:
        subprocess.run(ft_args, stdout=out_file)


# def phylo_annot_correlation(t, annots, dist_dict=None, measure="standard", top_pct=100, sample=1):
#
#     """ Prototype for testing phylogenetic conformity of some continuous-valued annotation to each node.
#     Note: currently only implemented for extant sequences. """
#
#     t = load_tree(t)
#     if not isinstance(annots, dict):
#         annots = parsers.tab_del_file_to_dict(annots)
#
#     # TODO: Currently using non-transformed values for annotation - may need max/min scaling or similar
#     annots_raw = [(label, float(val)) for label,val in annots.items()]
#
#     tree_dists = []
#     annot_diffs = []
#
#     for i in range(len(annots_raw)):
#         for j in range(i, len(annots_raw)):
#
#             annot_diffs.append(np.abs(annots_raw[i][1] - annots_raw[j][1]))
#
#             if tree_dists:
#                 tree_dists.append(dist_dict[(annots_raw[i][0], annots_raw[j][0])])
#             else:
#                 node1 = t.get_leaves_by_name(annots_raw[i][0])[0]
#                 node2 = t.get_leaves_by_name(annots_raw[j][0])[0]
#                 tree_dists.append(t.get_distance(node1, node2))
#
#     plt.plot(tree_dists, annot_diffs, '.')
#     plt.show()
#
#
# def place_hyrbid_seqs(t, s1, s2, node_positions=None, node_names=None, n_nodes=None, out_nwk=None):
#
#     t = load_tree(t)  # Read in tree if not already ete3.Tree
#
#     # Get s1 and s2 nodes
#     s1 = t & s1
#     s2 = t & s2
#
#     # Find most recent common ancestor (MRCA)
#     mrca = s1.get_common_ancestor(s2)
#     mrca_dist = t.get_distance(s1, mrca)  # Get distance to MRCA along trajectory
#
#     if node_positions:
#         int_dists = node_positions
#     else:  # Naively place n nodes at equal distance between s1 and s2
#         # Step size = total_path / (n-1)
#         total_dist = t.get_distance(s1, s2)
#         step = total_dist / (n_nodes+1)
#         int_dists = np.arange(step, total_dist-0.5*step, step)
#
#     if not node_names:
#         node_names = [f"int_{i}" for i in range(1, len(int_dists)+1)]
#
#     # As tree will change after adding new nodes, need to evaluate the path between s1 and s2 each time
#     for int_node, int_name in zip(int_dists, node_names):
#
#         # Get current path (node list) between s1 and s2
#         path1 = []
#         current = s1
#         while current != mrca:
#             path1.append(current)
#             current = current.up
#
#         # Path from node2 to MRCA (excluding MRCA), then reverse
#         path2 = []
#         current = s2
#         while current != mrca:
#             path2.append(current)
#             current = current.up
#         path2.reverse()
#
#         # Combine: path1 → MRCA → path2
#         full_path = path1 + [mrca] + path2
#
#         # Get list of boundaries between existing nodes between s1 and s2
#         bound_dists = [0]
#         current_dist = 0
#         for node in path1:
#             current_dist += node.dist
#             bound_dists.append(current_dist)
#         for node in path2:
#             current_dist += node.dist
#             bound_dists.append(current_dist)
#         mrca_idx = len(path1)
#
#         # Work out which two nodes we're in between
#         next_idx = None
#         for i in range(len(bound_dists)):
#             if bound_dists[i] > int_node:
#                 next_idx = i  # Index of next character along trajectory
#                 break
#         if not next_idx:
#             next_idx = len(bound_dists)-1
#
#         if int_node >= mrca_dist:  # We're already past the MRCA
#
#             branch_len = bound_dists[next_idx] - bound_dists[next_idx-1]
#
#             to_detach = full_path[next_idx]  # This is the node that will need to be detached
#             # Add intermediate node to parent of detached node
#             this_int = to_detach.up.add_child(dist=int_node - bound_dists[next_idx - 1])
#
#             # Detach and adjust branch length
#             detached = to_detach.detach()
#             detached.dist = branch_len - this_int.dist
#
#             # Add detached branch and a zero length branch to the intermediate
#             this_int.add_child(dist=0, name=int_name)
#             this_int.add_child(detached)
#
#         else:
#
#             branch_len = bound_dists[next_idx] - bound_dists[next_idx-1]
#
#             to_detach = full_path[next_idx-1]
#             # Add intermediate node to parent of detached node
#             this_int = to_detach.up.add_child(dist=bound_dists[next_idx]-int_node)
#
#             # Detach and adjust branch length
#             detached = to_detach.detach()
#             detached.dist = branch_len - this_int.dist
#
#             # Add detached branch and a zero length branch to the intermediate
#             this_int.add_child(dist=0, name=int_name)
#             this_int.add_child(detached)
#
#         if out_nwk:
#             t.write(outfile=out_nwk, format=1)
#
#
# # def (in_file):
# #
# #     with json.
#
# def parse_distance_json(json_path):
#     with open(json_path, 'r') as f:
#         data = json.load(f)
#     # keep original structure
#     distance_dict = data
#     all_sequences = {}
#     method = "equidistant"
#     for id_pair, info in data[method].items():
#         for key, value in info.items():
#             if key == "distance":
#                 continue
#             if key not in all_sequences:
#                 all_sequences[key] = value
#     return distance_dict, all_sequences
#
#
# def apply_trajectory(in_aln, in_tree, out_aln, out_tree, int_data, dist_type_idx=0):
#
#     """ TODO: Currently just takes a single distance type and outputs alignment and tree for it. Need to call multiple
#      times if you want all distance types. """
#
#     # Read in intermediate node data
#     data, seqs = parse_distance_json(int_data)
#
#     # Check which distance types are present
#     dist_types = list(data.keys())
#
#     # Get the start and end of the trajectory (s1+s2)  # NOTE: This was when start and end nodes were explicitly stated in json
#     # s1 = list(list(data[dist_types[dist_type_idx]].items())[0][1].items())[1][0]
#     # s2 = list(list(data[dist_types[dist_type_idx]].items())[-1][1].items())[2][0]
#
#     # TODO: Decide on how start and end nodes should be specified in the file type
#     # Currently just pulling from the json file name
#     int_data = int_data.split("/")[-1]  # If file isn't in current directory
#     s1 = int_data.split("_")[0]
#     s2 = int_data.split("_")[1]
#
#     # Get cumulative distances from s1 for all intermediate nodes
#     type_data = list(data.items())[dist_type_idx][1]
#     dist_from_s1 = 0
#     cum_dists = []  # Cumulative distance of each intermediate node from s1
#     int_names = []  # Names of intermediate nodes to place
#
#     for branch_id, branch_data in list(type_data.items())[:-1]:
#
#         int_names.append(list(branch_data.keys())[2])
#         dist_from_s1 += branch_data["distance"]
#         cum_dists.append(dist_from_s1)
#
#     # Append intermediate sequences to alignment
#     with open(in_aln, 'r') as in_file:
#         with open(out_aln, 'w') as out_file:
#             for line in in_file:
#                 out_file.write(line)
#             for name, seq_str in list(seqs.items())[1:-1]:
#                 out_file.write(f">{name}\n{seq_str}\n")
#
#     # Add intermediate nodes
#     place_hyrbid_seqs(in_tree, s1, s2, node_positions=cum_dists, node_names=int_names, out_nwk=out_tree)
#
#


# def place_hybrid_seqs_og(t, s1, s2, n_nodes):
#
#     t = load_tree(t)  # Read in tree if not already ete3.Tree
#
#     # Get s1 and s2 nodes
#     s1 = t&s1
#     s2 = t&s2
#
#     # Find most recent common ancestor (MRCA)
#     mrca = s1.get_common_ancestor(s2)
#
#     path1 = []
#     current = s1
#     while current != mrca:
#         path1.append(current)
#         current = current.up
#
#     # Path from node2 to MRCA (excluding MRCA), then reverse
#     path2 = []
#     current = s2
#     while current != mrca:
#         path2.append(current)
#         current = current.up
#     path2.reverse()
#
#     # Combine: path1 → MRCA → path2
#     full_path = path1 + [mrca] + path2
#
#     # Step size = total_path / (n-1)
#     total_dist = t.get_distance(s1, s2)
#     step = total_dist / (n_nodes -1)
#
#     bound_dists = []
#     current_dist = 0
#     for node in path1:
#         current_dist += node.dist
#         bound_dists.append(current_dist)
#     for node in path2:
#         current_dist += node.dist
#         bound_dists.append(current_dist)
#     mrca_idx = len(path1)
#
#     current_int_node = 0
#     current_anc_node = 0  # Next boundary point to be encountered
#     current_dist = 0
#     while current_int_node < n_nodes:
#         crossed_nodes = []
#         current_dist += step
#         for bound in bound_dists[current_anc_node:]:
#             if current_dist > bound:
#                 crossed_nodes.append(bound)
#             else:
#                 break
#
#         # Update node which is next to be encountered
#         current_anc_node += len(crossed_nodes)
#
#         # Need to insert the next intermediate node
#         # If moving down path, node to insert under is PREVIOUS node to be encountered
#         if current_anc_node > mrca_idx:
#
#             to_detach = path2[current_anc_node - mrca_idx - 1]  # This is the node that will need to be detached
#             # Add intermediate node to parent of detached node
#             this_int = to_detach.up.add_child(dist=current_dist-bound_dists[current_anc_node-1])
#
#             # Detach and adjust branch length
#             detached = to_detach.detach()
#             detached.dist = bound_dists[current_anc_node] - this_int.dist
#
#             # Add detached branch and a zero length branch to the intermediate
#             this_int.add_child(dist=0)
#             this_int.add(detached)
#
#         else:
#
#             to_detach = path1[current_anc_node]  # Index of path1 is shifted by 1 relative to bound dists
#
#             # Work out what the distance split will be between the two ancestral nodes
#             detached_dist = current_dist - bound_dists[current_anc_node-1] if current_anc_node > 0 else current_dist
#             int_dist = bound_dists[current_anc_node] - detached_dist
#             this_int = to_detach.up.add_child(dist=int_dist)
#
#
#             this_int = to_detach.up.add_child(dist=current_dist-bound_dists[current_anc_node-1])
#
#
#
#
#             # Insert node between root and A
#             a = t & "A"
#             intermediate = a.up.add_child(dist=0.0)
#             intermediate.add_child(a.detach())
#
#
#             # Add new node under ancestral node corresponding to PREVIOUS boundary point
#             # int_dist = current_dist - prev_boundary
#             # add node to prev_boundary at length=int_dist
#             # detach node corresponding to NEXT boundary point
#             # Add 0 length branch (fake) to new node, add detached node with length=NEXT_node_dist-int_dist
#
#         else:
#             # Add new node under NEXT boundary point
#             # int_dist = current_dist - PREVIOUS_node_dist
#             # add node to next_boundary at length = int_dist
#             # detach node corresponding to PREVIOUS boundary point
#             # Add 0 length branch to new node, add detached node with length = PREV_node-int_dist
#
#
#         insert_nodes(crossed_nodes)
#         if current_dist + step >
#
#
#     def insert_node(crossed_nodes)





