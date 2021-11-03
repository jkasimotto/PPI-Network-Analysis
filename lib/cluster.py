import itertools
import math

import markov_clustering as mcl
import networkx as nx
import pandas as pd

import PC2P.PC2P_ParallelMultiprocess
import PC2P.PC2P_Sequential
import lib.cluster
import lib.constants
import lib.files
import lib.graph


# PROPERTIES

def intersection(cluster1, cluster2):
    return list(set(cluster1).intersection(set(cluster2)))


def cluster_idxs_with_protein(clusters, protein):
    return [i for (i, cluster) in enumerate(clusters) if protein in cluster]


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


def get_connected_clusters(clusters):
    return [cluster for cluster in clusters if nx.is_connected(cluster)]


def percentage_connected(clusters):
    return len(get_connected_clusters(clusters)) / len(clusters)


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
    append_dataframe_columns(cluster_filepath, centrality_name, cids, orfs, centralities)


def append_average_clustering_coefficient(cluster_df_filepath, network, clusters):
    cluster_ids = []
    cluster_avg_clustering = []
    for cluster_id, cluster in enumerate(clusters):
        # Change to a subgraph to apply to networkx algorithm.
        cluster = network.subgraph(cluster)
        # Add data
        cluster_ids.append(cluster_id)
        cluster_avg_clustering.append(nx.average_clustering(cluster))  # Change cluster to a subgraph
    # Append it to our dataframe.
    append_dataframe_column(cluster_df_filepath, cluster_ids, 'avg_clust_coeff', cluster_avg_clustering)


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
    def __init__(self, matrix, result, sparse_clusters, systematic_clusters, modularity=None):
        self.matrix = matrix
        self.result = result
        self.sparse_clusters = sparse_clusters  # Obtained from mcl.get_clusters(result)
        self.clusters = systematic_clusters  # Obtained from mcl_systematic_clusters(network, sparse_clusters)
        self.modularity = None


def run_mcl(graph, inflation=2, weight=None):
    # Set all the edge weights
    # Convert to sparse matrix to run the algorithm.
    matrix = nx.to_scipy_sparse_matrix(graph, weight=weight)
    result = mcl.run_mcl(matrix, inflation=inflation)
    # These clusters use the sparse representation of a node i.e. an index
    # We want to convert them back to the systematic names of the proteins.
    clusters = mcl.get_clusters(result)
    clusters_systematic = mcl_systematic_clusters(graph, clusters)
    modularity = mcl.modularity(matrix=result, clusters=clusters)
    return MCLData(matrix, result, clusters, clusters_systematic, modularity)


def run_mcl_and_write_to_file(graph, filepath, inflation=2, weight=None):
    mcl_data = run_mcl(graph, inflation, weight)
    write_to_file(filepath, mcl_data.clusters)


def mcl_systematic_clusters(graph, clusters):
    # Convert clusters of indexes back to clusters of systemic names.
    nodes = list(graph.nodes())
    systematic_clusters = [[nodes[index] for index in cluster] for cluster in clusters]
    return systematic_clusters


# PC2P
def run_pc2p(network):
    G = network.copy()
    edge_cut = PC2P.PC2P_Sequential.Find_CNP(G)
    # PC2P clusters nodes into connected components by removing edges
    G_copy = G.copy()
    G_copy.remove_edges_from(edge_cut)
    # Save predicted clusters in
    G_cnp_components = list(nx.connected_components(G_copy))
    G_cnp_components.sort(key=len, reverse=True)
    return G_cnp_components


def run_pc2p_parallel(network):
    G = network.copy()
    edge_cut = PC2P.PC2P_ParallelMultiprocess.Find_CNP(G)
    # PC2P clusters nodes into connected components by removing edges
    G_copy = G.copy()
    G_copy.remove_edges_from(edge_cut)
    # Save predicted clusters in
    G_cnp_components = list(nx.connected_components(G_copy))
    G_cnp_components.sort(key=len, reverse=True)
    return G_cnp_components


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
    for cid, orf in list(zip(df['cluster'], df['protein'])):
        if cid == len(clusters):
            clusters.append([])
        clusters[-1].append(str(orf))
    return clusters


def read_yhtp2008():
    """
    Read the csv file into a list of sets representing protein clusters verified experimentally.
    """
    lines = lib.files.read_filelines(lib.files.make_filepath_to_clusters("validation/yhtp2008_complex.csv"))
    clusters = []
    for line in lines:
        cid, orf, _ = line.split(',')
        # The appearance of a new number in the cid column indicates a new cluster.
        if cid:
            clusters.append(set())
        clusters[-1].add(orf)
    return list(filter(lambda x: len(x) >= 3, clusters))


def read_cyc2008():
    lines = lib.files.read_filelines(lib.files.make_filepath_to_clusters("validation/CYC2008_complexes.txt"))
    clusters = [line.split(' ') for line in lines]
    return clusters


def read_sgd():
    lines = lib.files.read_filelines(lib.files.make_filepath_to_clusters("validation/SGD_complexes.txt"))
    clusters = [line.split(' ') for line in lines]
    clusters = list(map(set, clusters))
    return list(filter(lambda x: len(x) >= 3, clusters))


def write_to_file(filepath, clusters):
    cluster_ids = []
    proteins = []
    for i, cluster in enumerate(clusters):
        for protein in cluster:
            cluster_ids.append(i)
            proteins.append(protein)
    df = pd.DataFrame.from_records(list(zip(cluster_ids, proteins)), columns=['cluster', 'protein'])
    df.to_csv(filepath)


def read_clusters(network_name, clusters_name):
    clusters_filename = lib.files.make_clusters_filename(network_name, clusters_name)
    clusters_filepath = lib.files.make_filepath_to_clusters(clusters_filename)
    return read_csv(clusters_filepath)


# DATAFRAMES

def append_dataframe_column(df_filepath, cluster_ids, column_name, column_values):
    df1 = pd.read_csv(df_filepath, header=0, index_col=0)
    assert len(cluster_ids) == len(df1)
    df2 = pd.DataFrame.from_records(zip(cluster_ids, column_values), columns=['cluster', column_name])
    df = df1.merge(df2, on=['cluster'])
    df.to_csv(df_filepath)


# PLUG AND PLAY

def generate_dataframe(network, clusters):
    """
    This function is a plug and play method to take a network + clustering and creating a dataframe with 
    interesting variables for analysis.
    :param network: Networkx network - should be connected.
    :param clusters: List of nodes belonging to each cluster.
    :return: pandas Dataframe
    """

    # billy_metric = []

    # The features to be generated.
    ids = []
    size = []
    icp55_shell = []
    pim1_shell = []
    avg_clust_coeff = []
    is_connected = []

    # These are used several times and stored as a variable.
    icp55_shortest_paths = nx.shortest_path_length(network, lib.constants.ICP55)
    pim1_shortest_paths = nx.shortest_path_length(network, lib.constants.PIM1)

    # Add the data
    for id, cluster in enumerate(clusters):
        # Save the cluster as a subgraph to apply networkx algorithms.
        cluster = network.subgraph(cluster)
        # Save values as lists to reuse below.
        cluster_icp55_shortest_paths = list((icp55_shortest_paths[node] for node in cluster))
        cluster_pim1_shortest_paths = list((pim1_shortest_paths[node] for node in cluster))
        # Create records
        ids.append(id)
        size.append(len(cluster))
        icp55_shell.append(min(cluster_icp55_shortest_paths))
        pim1_shell.append(min(cluster_pim1_shortest_paths))
        avg_clust_coeff.append(nx.average_clustering(cluster))
        is_connected.append(nx.is_connected(cluster))

    # Return a dataframe of the records.
    # Make sure the names match the data in order.
    return pd.DataFrame.from_records(
        data=list(zip(
            # billy_metric,
            ids,
            size,
            icp55_shell,
            pim1_shell,
            avg_clust_coeff,
            is_connected
        )),
        columns=[
            'cluster',
            'size',
            'icp55_shell',
            'pim1_shell',
            'avg_clust_coeff',
            'is_connected'
        ]
    )


def generate_and_save_dataframe(network, clusters, filepath):
    df = generate_dataframe(network, clusters)
    df.to_csv(filepath)
