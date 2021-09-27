import itertools
import math

import markov_clustering as mcl
import networkx as nx
import pandas as pd

import lib.files


# PROPERTIES

def intersection(cluster1, cluster2):
    return list(set(cluster1).intersection(set(cluster2)))


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


# CENTRALITIES

def compute_centralities(cluster_filepath, network, clusters, centrality_function, centrality_name):
    """
    This function computes the centralities for each cluster and appends it as a column to the cluster file.

    :param cluster_filepath:  The filepath from which the clusters were read
    :param network: The network from which the clusters were computed
    :param clusters: # The clusters
    :param centrality_function: A function accepting a graph as an argument and returning a dict of node:value items.
    :param centrality_name: The name of the centrality measure which will be stored as a column header in cluster_filepath
    :return:
    """
    # Cast clusters as subgraphs to perform analysis on them.
    clusters = [network.subgraph(cluster) for cluster in clusters]
    # Perform analysis
    clusters_centralities = [centrality_function(cluster) for cluster in clusters]
    # Write to the cluster file
    cids = []
    orfs = []
    centralities = []
    for cid, cluster in enumerate(clusters_centralities):
        for orf, centrality in cluster.items():
            cids.append(cid)
            orfs.append(orf)
            centralities.append(centrality)
    append_column(cluster_filepath, centrality_name, cids, orfs, centralities)


# VALIDATION MEASURES

def accuracy(sensitivity, ppv):
    return math.sqrt(sensitivity * ppv)


def sensitivity(complexes, contingency_table):
    """
    :param contingency_table: A table with entries t(i,j) representing the number of shared proteins between complex i and cluster j
    """
    rows = range(len(contingency_table))
    cols = range(len(contingency_table[0]))
    # The sum of maximum match for each complex
    numerator = sum(max(contingency_table[i]) for i in rows)
    # Divided by the length of each complex
    denominator = (sum(len(complexes[i]) for i in rows))
    return numerator / denominator


def positive_predictive_value(contingency_table):
    rows = range(len(contingency_table))
    cols = range(len(contingency_table[0]))
    # The sum of maximum match for each cluster
    numerator = sum(max(contingency_table[i][j] for i in rows) for j in cols)
    # The sum of the numbers of all matches
    denominator = sum(sum(contingency_table[i][j] for i in rows) for j in cols)
    return numerator / denominator


def contingency_table(complexes, clusters):
    """
    :param complexes:  Our validation clusters. These are only protein complexes however and don't have functional modules.
    :param clusters: Our clusters.
    :return: A matrix with entries t(i, j) representing the number of nodes in complex i and cluster j
    """
    return [
        [len(intersection(cluster, complex)) for cluster in clusters]
        for complex in complexes
    ]


# MCL
class MCLData:
    def __init__(self, matrix, result, sparse_clusters, semantic_clusters, modularity=None):
        self.matrix = matrix
        self.result = result
        self.sparse_clusters = sparse_clusters  # Obtained from mcl.get_clusters(result)
        self.clusters = semantic_clusters  # Obtained from mcl_semantic_clusters(network, sparse_clusters)
        self.modularity = None


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


# FILE HANDLING

def read_csv(filepath, as_df=False):
    """


    :param filepath: The filepath to the data.
    :param as_df: Boolean. If true return pandas dataframe else return list of (cluster_id, protein_name).
    :return:
    """
    df = pd.read_csv(filepath, header=0, index_col=0)
    if as_df:
        return df

    clusters = []
    for cid, orf in list(zip(df['cid'], df['orf'])):
        if cid == len(clusters):
            clusters.append([])
        clusters[-1].append(orf)
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


def write_to_file(filepath, clusters):
    lines = []
    for i, cluster in enumerate(clusters):
        for node in cluster:
            lines.append(f"{i},{node}\n")
    with open(filepath, "w") as f:
        f.writelines(lines)


def append_column(filepath, column_name, cids, orfs, column_values):
    df1 = read_csv(filepath, as_df=True)
    df2 = pd.DataFrame.from_records(zip(cids, orfs, column_values), columns=['cid', 'orf', column_name])
    df = df1.merge(df2, on=['cid', 'orf'])
    df.to_csv(filepath)
