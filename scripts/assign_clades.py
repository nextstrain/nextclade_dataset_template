from Bio import Phylo
import pandas as pd
import argparse
from collections import defaultdict
import json

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Assign clades to a tree')
    parser.add_argument('--tree', type=str, help='Newick tree file')
    parser.add_argument('--clade-key', type=str, default='clade', help='Key to use for clade lookup in metadata')
    parser.add_argument('--metadata', type=str, help='metadata file with subtype information')
    parser.add_argument('--output', type=str, help='output file')
    args = parser.parse_args()

    tree = Phylo.read(args.tree, 'newick')
    metadata = pd.read_csv(args.metadata, sep='\t', index_col=0)

    node_to_subtype = {k: x[args.clade_key] for k, x in metadata.iterrows()}
    nodes_by_subtype = defaultdict(list)

    for node in tree.get_terminals():
        nodes_by_subtype[node_to_subtype[node.name]].append(node)

    clade_demarcations = {}
    for subtype, nodes in nodes_by_subtype.items():
        clade = tree.is_monophyletic(nodes)
        if clade:
            clade_demarcations[subtype] = clade
        else:
            clade_demarcations[subtype] = tree.common_ancestor(nodes)
            print("Subtype {} is not monophyletic".format(subtype))
        clade_demarcations[subtype].clade_label = subtype

    for node in tree.find_clades(order = 'preorder'):
        node.clade = 'unassigned'

    for node in tree.find_clades(order = 'preorder'):
        if hasattr(node, "clade_label"):
            node.clade = node.clade_label
            for child in node.find_clades():
                child.clade = node.clade_label

    node_data = {'nodes': {}}
    for node in tree.find_clades(order = 'preorder'):
        if node.is_terminal():
            if node.clade != node_to_subtype[node.name]:
                print("Node {} has clade {} but subtype {}".format(node.name, node.clade, node_to_subtype[node.name]))
                node.clade = node_to_subtype[node.name]
        datum = {'clade_membership': node.clade}
        if hasattr(node, 'clade_label'):
            datum['clade_annotation'] = node.clade_label
        node_data['nodes'][node.name] = datum


    with open(args.output, 'w') as f:
        json.dump(node_data, f, indent=2)