import functools

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import lib.files
import lib.map_names



def read_network(network_name):
    network_filename = lib.files.make_network_filename(network_name)
    network_filepath = lib.files.make_filepath_to_networks(network_filename)
    return read_weighted_edgelist(network_filepath)


def read_STRING():
    return nx.read_weighted_edgelist(lib.files.make_filepath_to_networks('4932.protein.links.v11.5.txt'))


def read_weighted_edgelist(filepath, nodetype=str):
    return nx.read_weighted_edgelist(filepath, nodetype=nodetype)


def read_inviable_proteins(with_organism_name=False):
    lines = lib.files.read_filelines(lib.files.make_filepath_to_data('inviable_proteins.csv'))
    lines = list(filter(None, lines))  # Remove empty lines
    lines = [line.split(',') for line in lines]  # Split into SGD, systemic name, common name, _, _
    systemic_names = [line[1].replace('"', '') for line in lines]
    if with_organism_name:
        systemic_names = ['4932.' + name for name in systemic_names] # Add in the organism name when needed.
    return systemic_names  # Return systemic name


def write_weighted_edgelist(network, filepath):
    nx.write_weighted_edgelist(network, path=filepath, delimiter=' ')


def make_from_nodes(nodes):
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    return graph


def make_from_union(graphs):
    graph_new = nx.Graph()
    for graph_old in graphs:
        graph_new.update(graph_old.edges(), graph_old.nodes())
    return graph_new


def make_from_intersection(graphs):
    return functools.reduce(
        lambda G, H: nx.intersection(G, H), graphs)


def get_edge_weight(graph, edge):
    return list(graph.get_edge_data(edge[0], edge[1]).values())[0]


def remove_diff(parent, subgraph):
    nodes_to_remove = [node for node in subgraph.nodes() if node not in parent.nodes()]
    subgraph.remove_nodes_from(nodes_to_remove)
    return subgraph


def remove_edges_below_threshold(graph, threshold):
    for edge in graph.edges:
        if get_edge_weight(graph, edge) <= threshold:
            graph.remove_edge(edge[0], edge[1])
    return graph


def remove_nodes_with_degree_lte(graph, k):
    nodes_to_remove = [node for node in graph.nodes if get_degree(graph, node) <= k]
    graph.remove_nodes_from(nodes_to_remove)
    return graph


def remove_nodes_from_subgraph(graph, subgraph):
    graph.remove_nodes_from(subgraph.nodes())
    return graph


def remove_nodes_from_list(graph, nodes):
    graph.remove_nodes_from(nodes)
    return graph


def get_largest_cc_with_node(graph, source):
    found = False
    # Largest to smallest.
    connected_components = sorted(nx.connected_components(graph), key=len, reverse=True)
    for cc in connected_components:
        if source in cc:
            found = True
            break
    if found:
        return graph.subgraph(cc)
    else:
        raise AssertionError("Node not in any connected component?")


def get_shortest_path_length_from_source_to_targets(graph, source, targets):
    shortest_paths = nx.single_source_shortest_path_length(graph, source)
    shortest_paths = [(node, length) for (node, length) in shortest_paths.items() if node in targets]
    return shortest_paths


def get_number_of_nodes(graph):
    return graph.number_of_nodes()


def get_number_of_edges(graph):
    return graph.number_of_edges()


def get_degree(graph, node):
    return graph.degree(node)


def get_neighbours(graph, node):
    return list(graph.neighbours(node))


def is_connected(graph):
    return nx.connected.is_connected(graph)


def get_number_connected_components(graph):
    return nx.connected.number_connected_components(graph)


def get_largest_connected_component(graph):
    return graph.subgraph(
        max(nx.connected_components(graph), key=len)
    )


def get_neighbourhood(graph, node, path_length=1, frozen=True):
    # TODO: Change this to use the DataFrame data to avoid calculating everytime.
    nodes = [node, *get_nodes_m_or_less_away(graph, node, path_length)]
    return graph.subgraph(nodes)


def get_giant(graph):
    return get_largest_connected_component(graph)


def get_degree_sequence(graph):
    degrees = [graph.degree()[node] for node in list(graph.nodes())]
    degrees = sorted(degrees)
    return np.array(degrees)


def get_mean_degree(graph):
    return np.average(get_degree_sequence(graph))


def get_nodes_of_degree(graph, degree):
    return [node for node in graph.nodes() if graph.degree()[node] == degree]


def get_nodes_m_or_less_away(graph, node, m):
    shortest_paths = nx.single_source_shortest_path_length(graph, node)
    return [node for (node, length) in shortest_paths.items() if length <= m]


def plot_degrees_vs_threshold(graph):
    # 4 thresholds
    # Increase threshold per plot to remove excess nodes.
    thresholds = range(600, 1000, 100)

    fig, axs = plt.subplots(2, 2)

    for i, threshold in enumerate(thresholds):
        graph = remove_edges_below_threshold(graph, threshold)
        degrees = get_degree_sequence(graph)

        row = 0 if i < 2 else 1
        col = i % 2
        axs[row][col].plot(degrees, color='blue', marker='o', markersize=4, linewidth=2)

    plt.show()


def get_node_degrees_by_threshold(graph, node, start=700, stop=1000, step=100):
    degrees = []
    for threshold in range(start, stop, step):
        graph = lib.graph.remove_edges_below_threshold(graph, threshold)
        degrees.append(graph.degree(node))
    return degrees

def rename_with_gene_names(network):
    network_systematic_names = list(network.nodes())
    network_gene_names = lib.map_names.map_names_descriptions(network_systematic_names, 'systematic_name_nonumber', 'gene_name')
    mapping = dict(zip(network_systematic_names, network_gene_names))
    network = nx.relabel_nodes(network, mapping) # Idaelly better to update in place but can't use frozen graphs then.
    return network

#Turn a path represented as a list of nodes into a list of tuples which are edges
def path_nodes_to_edges(path_nodes):
    
    path_edges = []

    for i in range(len(path_nodes) - 1):
    
        edge = (path_nodes[i], path_nodes[i + 1])
        path_edges.append(edge)
        
    return(path_edges)

#Get edge weights from a path represented as a list of nodes
#Returns 0 if no edge
def path_nodes_to_edgeweights(path_nodes, network):

    path_edges = lib.graph.path_nodes_to_edges(path_nodes)

    weight_dict = {}

    for edge in path_edges:
        
        if edge in network.edges():
            weight_dict[edge] = network[edge[0]][edge[1]]["weight"]
        else:
            weight_dict[edge] = 0

    return(weight_dict)
