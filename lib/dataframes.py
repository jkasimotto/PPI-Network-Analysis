from os.path import exists

import networkx as nx
import pandas as pd

import PC2P.Analysis.PredictedClusters_Analysis as pc2p_analysis
import lib.centrality
import lib.cluster
import lib.constants
import lib.files
import lib.graph


def node_dataframe(
        network_name,
        clusters_name,
        _degree=False,
        _inviable=False,
        _icp55_shell=False,
        _pim1_shell=False,
        _betweenness=False,
        _closeness=False,
        _bridging=False,
        _eigenvector=False,
        _cluster_id=False,
        _cluster_size=False,
        _cluster_degree=False,
        _cluster_eigenvector=False,
        _cluster_closeness=False,
        _cluster_betweenness=False
):
    """

    :param network: The filename (without the extension) of the network. E.g. icp55-cc-900-inv
    :param clusters: The filename (without the path or extension) of the clusters. E.g. mcl-2.5
    :param _columns: The desired computations.
    :return:
    """

    # Read in the network.
    network = lib.graph.read_weighted_edgelist(lib.files.make_filepath_to_networks(f"{network_name}.txt"))

    # Read in the clusters.
    clusters_filename = lib.files.make_clusters_filename(network_name, clusters_name)
    clusters_filepath = lib.files.make_filepath_to_clusters(clusters_filename)
    clusters = lib.cluster.read_csv(clusters_filepath)

    # Make sure the clusters partition the network.
    clusters_proteins = lib.cluster.proteins(clusters)
    assert len(clusters_proteins) == len(network.nodes())
    assert set(clusters_proteins) == set(network.nodes())

    # Cast clusters to a subgraph datatype for algorithms.
    clusters = list(map(network.subgraph, clusters))

    # Define the columns of the dataframe.
    columns = {
        'protein': [],
        'degree': [],
        'inviable': [],
        'icp55_shell': [],
        'pim1_shell': [],
        'betweenness': [],
        'eigenvector': [],
        'closeness': [],
        'bridging': [],
        'cluster_id': [],
        'cluster_size': [],
        'cluster_degree': [],
        'cluster_eigenvector': [],
        'cluster_closeness': [],
        'cluster_betweenness': []
    }

    # Prepare data for algorithms.
    print("Preparing data")
    inviable_proteins = lib.graph.read_inviable_proteins()
    if _icp55_shell or _pim1_shell:
        shortest_path_lengths = dict(nx.shortest_path_length(network))
    if _icp55_shell:
        icp55_shell = shortest_path_lengths[lib.constants.ICP55]
    if _pim1_shell:
        pim1_shell = shortest_path_lengths[lib.constants.PIM1]

    # Run flagged algorithms.
    if _degree:
        degree = network.degree()
    if _inviable:
        inviable = {node: 1 if node in inviable_proteins else 0 for node in network}
    if _betweenness:
        print("Computing betweenness")
        betweenness = nx.betweenness_centrality(network)
    if _closeness:
        print("Computing closeness")
        closeness = nx.closeness_centrality(network)
    if _bridging:
        print("Computing bridging")
        if _betweenness:
            bridging = lib.centrality.bridging_centrality(network,
                                                          betweeness=betweenness)  # Pass betweeness to avoid recalc.
        else:
            bridging = lib.centrality.bridging_centrality(network)  # We haven't computed betweenness.
    if _eigenvector:
        print("Computing eigenvector")
        eigenvector = nx.eigenvector_centrality(network)

    # Run flagged algorithms in clusters.
    if _cluster_closeness:
        print("Computing closeness in cluster")
        cluster_closeness = [nx.closeness_centrality(cluster) for cluster in clusters]
    if _cluster_eigenvector:
        print("Computing eigenvector in cluster")
        cluster_eigenvector = [nx.eigenvector_centrality(cluster) for cluster in clusters]
    if _cluster_betweenness:
        print("Computing betweenness in cluster")
        cluster_betweenness = [nx.betweenness_centrality(cluster) for cluster in clusters]

    # Iterate over cluster_ids (it's index) and clusters
    #   Since clusters form a partition of the network we therefore iterate over every node exactly once.
    print("Arranging columns")
    for cluster_id, cluster in enumerate(clusters):
        for node in cluster:
            columns['protein'].append(node)
            if _degree:
                columns['degree'].append(degree[node])
            if _inviable:
                columns['inviable'].append(inviable[node])
            if _icp55_shell:
                columns['icp55_shell'].append(icp55_shell[node])
            if _pim1_shell:
                columns['pim1_shell'].append(pim1_shell[node])
            if _betweenness:
                columns['betweenness'].append(betweenness[node])
            if _closeness:
                columns['closeness'].append(closeness[node])
            if _bridging:
                columns['bridging'].append(bridging[node])
            if _eigenvector:
                columns['eigenvector'].append(eigenvector[node])
            if _cluster_id:
                columns['cluster_id'].append(cluster_id)
            if _cluster_size:
                columns['cluster_size'].append(len(cluster))
            if _cluster_degree:
                columns['cluster_degree'].append(cluster.degree()[node])
            if _cluster_closeness:
                columns['cluster_closeness'].append(cluster_closeness[cluster_id][node])
            if _cluster_eigenvector:
                columns['cluster_eigenvector'].append(cluster_eigenvector[cluster_id][node])
            if _cluster_betweenness:
                columns['cluster_betweenness'].append(cluster_betweenness[cluster_id][node])

    # Remove empty columns
    columns = {k: v for (k, v) in columns.items() if v != []}

    # Write the dataframe to a file
    print("Adding dataframe to file")
    filename = f"{network_name}.{clusters_name}.nodes.dataframe.csv"
    filepath = lib.files.make_path_to_dataframes(filename)
    add_dataframe_columns(filepath, columns, on=['protein'])


def cluster_dataframe(network_name, clusters_name):
    # These are all the columns
    columns = {
        'id': [],
        'size': [],
        'icp55_shell': [],
        'pim1_shell': [],
        'percent_connected': [],
        'modularity': [],
        'sensitivity_yhtp': [],
        'positive_predicted_value_yhtp': [],
        'accuracy_yhtp': [],
        'fraction_matched_yhtp': [],
        'separation_yhtp': [],
        'precision_yhtp': [],
        'recall_yhtp': [],
        'f_measure_yhtp': [],
        'sensitivity_sgd': [],
        'positive_predicted_value_sgd': [],
        'accuracy_sgd': [],
        'fraction_matched_sgd': [],
        'separation_sgd': [],
        'precision_sgd': [],
        'recall_sgd': [],
        'f_measure_sgd': []
    }

    network = lib.graph.read_network(network_name)

    clusters = lib.cluster.read_clusters(network_name, clusters_name)

    # Validated complexes
    yhtp2008 = lib.cluster.read_yhtp2008()
    sgd = lib.cluster.read_sgd()

    columns['percent_connected'].append(len(
        [cluster for cluster in clusters if nx.is_connected(network.subgraph(cluster))]
    ) / len(clusters))

    columns['sensitivity_yhtp'].append(pc2p_analysis.clusteringwise_sensitivity(yhtp2008, clusters))
    columns['positive_predicted_value_yhtp'].append(pc2p_analysis.positive_predictive_value(yhtp2008, clusters))
    columns['accuracy_yhtp'].append(pc2p_analysis.accuracy(yhtp2008, clusters))
    columns['fraction_matched_yhtp'].append(pc2p_analysis.fraction_matched(yhtp2008, clusters))
    columns['separation_yhtp'].append(pc2p_analysis.clusteringwise_separation(yhtp2008, clusters))
    columns['precision_yhtp'].append(pc2p_analysis.precision_Jaccard(yhtp2008, clusters))
    columns['recall_yhtp'].append(pc2p_analysis.recall_Jaccard(yhtp2008, clusters))
    columns['f_measure_yhtp'].append(pc2p_analysis.F_measure_Jaccard(yhtp2008, clusters))

    columns['sensitivity_sgd'].append(pc2p_analysis.clusteringwise_sensitivity(sgd, clusters))
    columns['positive_predicted_value_sgd'].append(pc2p_analysis.positive_predictive_value(sgd, clusters))
    columns['accuracy_sgd'].append(pc2p_analysis.accuracy(sgd, clusters))
    columns['fraction_matched_sgd'].append(pc2p_analysis.fraction_matched(sgd, clusters))
    columns['separation_sgd'].append(pc2p_analysis.clusteringwise_separation(sgd, clusters))
    columns['precision_sgd'].append(pc2p_analysis.precision_Jaccard(sgd, clusters))
    columns['recall_sgd'].append(pc2p_analysis.recall_Jaccard(sgd, clusters))
    columns['f_measure_sgd'].append(pc2p_analysis.F_measure_Jaccard(sgd, clusters))

    filename = lib.files.make_clusters_dataframe_filename(network_name)
    filepath = lib.files.make_path_to_dataframes(filename)
    add_dataframe_rows(filepath, columns)


def add_dataframe_columns(filepath, columns, on=None):
    """
    This function appends columns to an existing dataframe file or creates a brand new file is none exists.
    :param filepath:
    :param columns: A dict of {col_name: rows}
    :param on: The columns to merge old and new data on (same as df.merge)
    :return:
    """
    # Create the dataframe.
    df = pd.DataFrame.from_dict(columns)

    # If no file exists, write to the filepath.
    if not exists(filepath):
        df.to_csv(filepath)
        return

    # If a file exists, we need to know what columns to merge the new data on.
    assert on is not None

    # Then read the existing file.
    df2 = pd.read_csv(filepath, index_col=0, header=0)
    df2 = df2.astype({'protein': str})  # Networks such as karate-club and gnp-100-0.5 use integer node names.

    # Assert the same number of rows exist (i.e. all proteins accounted for)
    assert len(df) == len(df2)

    # Remove overlapping columns from the old dataframe.
    overlap = [col for col in df.columns if col in df2.columns and col not in on]
    df2 = df2.drop(columns=overlap)

    # Merge the dataframes on the index
    df3 = df2.merge(df, on=on)

    # Write to the filepath
    df3.to_csv(filepath)


def add_dataframe_rows(filepath, columns):
    # Create the dataframe.
    df = pd.DataFrame.from_dict(columns)

    # If no file exists, write to the filepath.
    if not exists(filepath):
        df.to_csv(filepath)
        return

    # Then read the existing file.
    df2 = pd.read_csv(filepath, index_col=0, header=0)
    df2 = df2.astype({'protein': str})  # Networks such as karate-club and gnp-100-0.5 use integer node names.

    # Assert we have the same columns
    assert df.columns == df2.columns

    # Concatenate the new rows to the old rows
    df3 = pd.concat([df2, df], ignore_index=True)

    df3.to_csv(filepath)



def add_min_shell(filepath):
    """df must have icp55_shell and pim1_shell"""
    df = pd.read_csv(filepath, header=0, index_col=0)
    assert "icp55_shell" in df.columns and "pim1_shell" in df.columns
    df['min-shell'] = df['icp55_shell'].combine(df['pim1_shell'], min, 0)
    df.to_csv(filepath)


def add_min_shell_weighted_centralities(filepath):
    df = pd.read_csv(filepath, header=0, index_col=0)
    df['betweenness/min-shell'] = df['betweenness'] / df['min-shell']
    df['closeness/min-shell'] = df['closeness'] / df['min-shell']
    df['eigenvector/min-shell'] = df['eigenvector'] / df['min-shell']
    df['bridging/min-shell'] = df['bridging'] / df['min-shell']
    df.to_csv(filepath)
