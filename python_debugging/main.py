"""Script for retrieving a molecule from PubChem and drawing it using NetworkX."""
import networkx as nx
import matplotlib.pyplot as plt

color_map = {"C": "gray", "O": "red", "N": "blue", "H": "#D3D3D3"}


def create_graph(nodes, edges):
    """Create a graph from element symbols and (atom index, atom index, order) bonds."""
    G = nx.Graph()
    for i, node in enumerate(nodes):
        atom = (0, node)
        G.add_node(atom)
    for edge in edges:
        node_1 = (edge[0], nodes[edge[0]])
        node_2 = (edge[1], nodes[edge[1]])
        G.add_edge(node_1, node_2, weight=edge[2])
    return G


def retrieve_sdf(molecule_name: str, local_path=None) -> list[str]:
    """Retrieve SDF lines by molecule name, or read an optional local SDF."""
    from urllib import request
    sdf_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{molecule_name}/SDF"
    source = open(local_path, "rb") if local_path else request.urlopen(sdf_url)
    with source as response:
        text = response.read()
        data = text.split("\n")
    return data


def process_sdf(sdf_text: list[str], include_hydrogen: bool = False):
    """Parse the first V2000 SDF record into elements and zero-based bond indices."""
    num_atoms = int(sdf_text[3][:3])
    num_bonds = int(sdf_text[3][3:6])
    atom_end = 4 + num_atoms
    bond_end = atom_end + num_bonds
    atoms = sdf_text[4:atom_end]
    bonds = sdf_text[atom_end:bond_end]
    elements = [x.split()[3] for x in atoms]
    bond_pairs = [(int(x[:3])-1, int(x[3:6])-1, int(x[6:9])) for x in bonds]
    # Preserve atom identity when removing hydrogen, even if H atoms are interleaved.
    keep = [i for i, element in enumerate(elements) if include_hydrogen or element != "H"]
    old_to_new = {old: new for new, old in enumerate(keep)}
    bond_pairs = [(old_to_new[a], old_to_new[b], order)
                  for a, b, order in bond_pairs if a in old_to_new and b in old_to_new]
    return [elements[i] for i in keep], bond_pairs


def create_image(G, colors, widths, labels, molecule_name):
    pos = nx.kamada_kawai_layout(G, weight=None)
    nx.draw_networkx(G, node_color=colors, width=widths, with_labels=True, labels=labels)
    plt.savefig(f"{molecule_name}.png", dpi=300)


if __name__ == "__main__":
    import sys
    molecule_name = sys.argv[1]
    local_path = sys.argv[2] if len(sys.argv) > 2 else None
    sdf = retrieve_sdf(molecule_name, local_path)
    nodes, edges = process_sdf(sdf)
    G = create_graph(nodes, edges)
    colors = [color_map[x] for x in nodes]
    widths = [2 * bond[2] for bond in edges]
    node_labels = {node: node[1] for node in G.nodes}
    create_image(G, colors, widths, node_labels, molecule_name)
