"""
Script for retrieving molecule from PubChem and drawing it using networkX
"""

breakpoint()

import networkx as nx
import matplotlib.pyplot as plt

# This will be for coloring nodes in the visualization
color_map = {"C": "gray", "O": "red", "N": "blue", "H": "#D3D3D3"}


def create_graph(nodes, edges):
    """
    Create a graph from a list of molecular nodes (atoms) and edges (bonds).
    
    Parameters
    ----------
    nodes: list[tuple]
        A list of elements in the molecule.

    edges: list[tuple]
        A list of bonds in the molecule. The first two elements should indicate
        the atom indices, with the third indicating the bond order.

    Returns
    -------
    nx.Graph
        A networkX graph.
    """
    G = nx.Graph()

    for i, node in enumerate(nodes):
        atom = (0, node)
        G.add_node(atom)

    for edge in edges:
        # Have to have nodes which have the same
        # names as those added above.
        node_1 = (edge[0], nodes[edge[0]])
        node_2 = (edge[1], nodes[edge[1]])

        G.add_edge(node_1, node_2, weight=edge[2])

    return G


def retrieve_sdf(molecule_name: str) -> str:
    """
    Retrieve molecule structure from pubchem based on name.

    Parameters
    ----------
    molecule_name: str
        The name of the molecule to retrieve

    Returns
    -------
    data : str
        The text for the sdf file for the molecule
    """

    from urllib import request
  
    sdf_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molecule_name}/SDF"
    )

    with request.urlopen(sdf_url) as response:
        text = response.read()
        data = text.split("\n")

    return data


def process_sdf(sdf_text: str, include_hydrogen: bool = False):

    num_atoms = int(sdf_text[3].split()[0])
    num_bonds = int(sdf_text[3].split()[1])

    atom_end = 4 + num_atoms
    bond_end = 4 + num_atoms + num_bonds

    atoms = sdf_text[4:atom_end]
    bonds = sdf_text[atom_end:bond_end]

    elements = [x.split()[3] for x in atoms]

    bond_pairs = [
        (int(x.split()[0]), int(x.split()[1]), int(x.split()[2])) for x in bonds
    ]

    if include_hydrogen is not True:
        bond_pairs = [
            x
            for x in bond_pairs
            if elements[x[0] - 1] != "H" and elements[x[1] - 1] != "H"
        ]
        elements = [x for x in elements if x != "H"]

    # Fix so that bond indices use python indexing
    bond_pairs = [(x[0] - 1, x[1] - 1, x[2]) for x in bond_pairs]

    return elements, bond_pairs


def create_image(G, colors, widths, labels, molecule_name):

    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx(
        G, node_color=colors, width=widths, with_labels=True, labels=labels
    )

    plt.savefig(f"{molecule_name}.png", dpi=300)


if __name__ == "__main__":

    import sys

    molecule_name = sys.argv[1]
  
    sdf = retrieve_sdf(molecule_name)

    nodes, edges = process_sdf(sdf)

    G = create_graph(nodes, edges)

    colors = [ color_map[x] for x in nodes ]
    widths = [ 2 * bond[2] for bond in edges ]
    node_labels = {node: node[1] for node in G.nodes}

    create_image(G, colors, widths, node_labels, molecule_name)
