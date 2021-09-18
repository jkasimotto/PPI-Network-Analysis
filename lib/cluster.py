import itertools

import markov_clustering as mcl
import networkx as nx

import lib.files


def read_yhtp2008():
    """
    Read the csv file into a list of sets representing protein complexes verified experimentally.
    """
    lines = lib.files.read_filelines(lib.files.make_filepath_to_data("yhtp2008_complex.csv"))
    complexes = []
    for line in lines:
        cid, orf, _ = line.split(',')
        # The appearance of a new number in the cid column indicates a new complex.
        if cid:
            complexes.append(set())
        complexes[-1].add(orf)
    return complexes


def get_complexes_containing_protein(complexes, protein):
    return [complex for complex in complexes if protein in complex]


def get_complexes_size_sequence(complexes):
    return sorted([len(complex) for complex in complexes])


def count_complexes_containing_protein(complexes, orf):
    return len(get_complexes_containing_protein(complexes, orf))


def get_proteins_in_complexes(complexes):
    return list(itertools.chain.from_iterable(complexes))


def non_overlapping(old_clusters):
    """Very inefficient"""
    # Keep a list of all nodes and the index of their largest cluster
    old_clusters = sorted(old_clusters, key=len, reverse=True)  # Largest to smallest
    new_clusters = []
    nodes = {k: None for k in itertools.chain.from_iterable(old_clusters)}
    for i, cluster in enumerate(old_clusters):
        for node in cluster:
            if nodes[node] is None:
                nodes[node] = i
    for i in range(len(old_clusters)):
        new_clusters.append([])
        for node, index in nodes.items():
            if index == i:
                new_clusters[-1].append(node)
    return new_clusters


def remove_nodes_far_from_source(network, cluster, source, path_length):
    shortest_paths = nx.single_source_shortest_path_length(network, source)
    return [node for node in cluster if shortest_paths[node] <= path_length]


def run_mcl(graph, inflation=2):
    # Convert to sparse matrix to run the algorithm.
    matrix = nx.to_scipy_sparse_matrix(graph)
    result = mcl.run_mcl(matrix, inflation=inflation)
    clusters = mcl.get_clusters(result)
    clusters_semantic = mcl_semantic_clusters(graph, clusters)
    modularity = mcl.modularity(matrix=result, clusters=clusters)
    return MCLData(matrix, result, clusters, clusters_semantic, modularity)


def mcl_semantic_clusters(graph, clusters):
    # Convert clusters of indexes back to clusters of systemic names.
    nodes = list(graph.nodes())
    semantic_clusters = [[nodes[index] for index in cluster] for cluster in clusters]
    return semantic_clusters


def neighbourhood_clusters(clusters, shortest_path_lengths, path_length=1):
    local_clusters = [[node for node in cluster if shortest_path_lengths[node] <= path_length] for cluster in clusters]
    local_clusters = lib.cluster.non_overlapping(local_clusters)
    local_clusters = sorted(local_clusters, key=len)
    return list(filter(None, local_clusters))

def write_to_file(filepath, clusters):
    lines = []
    for i, cluster in enumerate(clusters):
        for node in cluster:
            lines.append(f"{i},{node}\n")
    with open(filepath, "w") as f:
        f.writelines(lines)


class MCLData:
    def __init__(self, matrix, result, sparse_clusters, semantic_clusters, modularity=None):
        self.matrix = matrix
        self.result = result
        self.sparse_clusters = sparse_clusters  # Obtained from mcl.get_clusters(result)
        self.clusters = semantic_clusters  # Obtained from mcl_semantic_clusters(network, sparse_clusters)
        self.modularity = None

