import itertools

import markov_clustering as mcl
import networkx as nx

import lib.files


def read_csv(filepath):
    lines = lib.files.read_filelines(filepath)

    clusters = []
    for line in lines:
        cluster_id, protein = line.split(',')
        cluster_id = int(cluster_id)
        if len(clusters) == cluster_id:
            clusters.append(set())
        clusters[-1].add(protein)
    return clusters


def read_yhtp2008():
    """
    Read the csv file into a list of sets representing protein clusters verified experimentally.
    """
    lines = lib.files.read_filelines(lib.files.make_filepath_to_clusters("yhtp2008_cluster.csv"))
    clusters = []
    for line in lines:
        cid, orf, _ = line.split(',')
        # The appearance of a new number in the cid column indicates a new cluster.
        if cid:
            clusters.append(set())
        clusters[-1].add(orf)
    return clusters


def clusters_with_protein(clusters, protein):
    return [cluster for cluster in clusters if protein in cluster]


def lengths(clusters):
    return sorted([len(cluster) for cluster in clusters])


def number_clusters_with_protein(clusters, orf):
    return len(clusters_with_protein(clusters, orf))


def proteins(clusters):
    return list(itertools.chain.from_iterable(clusters))


def neighbourhood_clusters(clusters, shortest_path_lengths, path_length=1):
    # Remove node far away
    local_clusters = [[node for node in cluster if shortest_path_lengths[node] <= path_length] for cluster in clusters]
    # Make non-overlapping
    local_clusters = lib.cluster.non_overlapping(local_clusters)
    # Sort from smallest to largest
    local_clusters = sorted(local_clusters, key=len)
    # Remove empty clusters
    return list(filter(None, local_clusters))


def non_overlapping(old_clusters):
    """Very inefficient. TODO: Write 1-3 tests."""
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


def remove_clusters_of_size_lte(clusters, size=3):
    return [cluster for cluster in clusters if len(cluster) > size]


def run_mcl(graph, inflation=2):
    # Convert to sparse matrix to run the algorithm.
    matrix = nx.to_scipy_sparse_matrix(graph)
    result = mcl.run_mcl(matrix, inflation=inflation)
    clusters = mcl.get_clusters(result)
    clusters_semantic = mcl_semantic_clusters(graph, clusters)
    modularity = mcl.modularity(matrix=result, clusters=clusters)
    return MCLData(matrix, result, clusters, clusters_semantic, modularity)


def run_mcl_and_write_to_file(graph, filepath, inflation=2):
    mcl_data = run_mcl(graph, inflation)
    write_to_file(filepath, mcl_data.clusters)


def mcl_semantic_clusters(graph, clusters):
    # Convert clusters of indexes back to clusters of systemic names.
    nodes = list(graph.nodes())
    semantic_clusters = [[nodes[index] for index in cluster] for cluster in clusters]
    return semantic_clusters


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
