
import time
import numpy as np
import pandas as pd
import networkx as nx

import lib.graph
import lib.cluster
import PC2P.Analysis.PredictedClusters_Analysis as pc2p_analysis
if __name__ == "__main__":
    network_name = 'icp55-cc-900-inv'
    network_filepath = lib.files.make_filepath_to_networks(f'{network_name}.txt')
    network = lib.graph.read_weighted_edgelist(network_filepath)

    # These are the two gold standard validation sets to use
    # We need to convert each cluster from a list to a set for pc2p_analysis
    yhtp2008 = lib.cluster.read_yhtp2008()
    yhtp2008 = list(map(set, yhtp2008))
    sgd = lib.cluster.read_sgd()
    sgd = list(map(set, sgd))

    # These are the metrics to calculate.
    # Anything with a suffix _yhtp will use the yhtp2008 validation set.
    # Anything with a suffix _sgd will use the sgd validation set.
    # These are the independant metrics
    percentage_connected = []
    modularity = []
    num_clusters = []
    # YHTP metrics
    sensitivity_yhtp = []
    positive_predicted_value_yhtp = []
    accuracy_yhtp = []
    fraction_matched_yhtp = []
    separation_yhtp = []
    precision_jaccard_yhtp = []
    recall_jaccard_yhtp = []
    f_measure_jaccard_yhtp = []
    # SGD metrics
    sensitivity_sgd = []
    positive_predicted_value_sgd = []
    accuracy_sgd = []
    fraction_matched_sgd = []
    separation_sgd = []
    precision_jaccard_sgd = []
    recall_jaccard_sgd = []
    f_measure_jaccard_sgd = []

    # These are the different inflation values we will use.
    inflations = np.linspace(3.1, 4, num=10, endpoint=True)

    # These are the times taken for different inflation values
    times = []

    # These are filepaths to each dataset
    filepaths = []

    # Calculate total time taken
    time_0 = time.time()

    # for inflation in inflations:
    for i in range(1):
        start = time.time()

        # Run markov clustering
        # mcl_data = lib.cluster.run_mcl(network, inflation)
        # clusters = mcl_data.clusters

        # Run PC2P clustering
        clusters = lib.cluster.run_pc2p(network)

        # We need to convert each cluster from a list to a set for pc2p_analysis.
        clusters = list(map(set, clusters))

        # Save the clusters to a file
        # filepath = lib.files.make_filepath_to_mcl_clusters(f"{network_name}.mcl-{inflation}.csv")
        filepath = lib.files.make_filepath_to_clusters(f"{network_name}.pc2p.csv")
        filepaths.append(filepath)
        lib.cluster.write_to_file(filepath, clusters)

        # # Compute the metrics for each clustering.
        # modularity.append(mcl_data.modularity)


        num_clusters.append(len(clusters))
        percentage_connected.append(len(
            [cluster for cluster in clusters if nx.is_connected(network.subgraph(cluster))]
        ) / len(clusters))

        # Compute yhtp metrics
        sensitivity_yhtp.append(pc2p_analysis.clusteringwise_sensitivity(yhtp2008, clusters))
        positive_predicted_value_yhtp.append(pc2p_analysis.positive_predictive_value(yhtp2008, clusters))
        accuracy_yhtp.append(pc2p_analysis.accuracy(yhtp2008, clusters))
        fraction_matched_yhtp.append(pc2p_analysis.fraction_matched(yhtp2008, clusters))
        separation_yhtp.append(pc2p_analysis.clusteringwise_separation(yhtp2008, clusters))
        precision_jaccard_yhtp.append(pc2p_analysis.precision_Jaccard(yhtp2008, clusters))
        recall_jaccard_yhtp.append(pc2p_analysis.recall_Jaccard(yhtp2008, clusters))
        f_measure_jaccard_yhtp.append(pc2p_analysis.F_measure_Jaccard(yhtp2008, clusters))

        # Compute sgd metrics
        sensitivity_sgd.append(pc2p_analysis.clusteringwise_sensitivity(sgd, clusters))
        positive_predicted_value_sgd.append(pc2p_analysis.positive_predictive_value(sgd, clusters))
        accuracy_sgd.append(pc2p_analysis.accuracy(sgd, clusters))
        fraction_matched_sgd.append(pc2p_analysis.fraction_matched(sgd, clusters))
        separation_sgd.append(pc2p_analysis.clusteringwise_separation(sgd, clusters))
        precision_jaccard_sgd.append(pc2p_analysis.precision_Jaccard(sgd, clusters))
        recall_jaccard_sgd.append(pc2p_analysis.recall_Jaccard(sgd, clusters))
        f_measure_jaccard_sgd.append(pc2p_analysis.F_measure_Jaccard(sgd, clusters))
        #
        end = time.time()
        seconds = end - start
        times.append(seconds)

    # Make a dataframe with each row representing a clustering.
    # Save the dataframe to a file.
    df2 = pd.DataFrame.from_records(
        list(zip(
            # inflations,
            # modularity,
            percentage_connected,
            sensitivity_yhtp,
            positive_predicted_value_yhtp,
            accuracy_yhtp,
            fraction_matched_yhtp,
            separation_yhtp,
            precision_jaccard_yhtp,
            recall_jaccard_yhtp,
            f_measure_jaccard_yhtp,
            sensitivity_sgd,
            positive_predicted_value_sgd,
            accuracy_sgd,
            fraction_matched_sgd,
            separation_sgd,
            precision_jaccard_sgd,
            recall_jaccard_sgd,
            f_measure_jaccard_sgd,
            filepaths,
            num_clusters
        )),
        columns=[
            # "inflation",
            # "modularity",
            "percentage_connected",
            "sensitivity_yhtp",
            "positive_predicted_value_yhtp",
            "accuracy_yhtp",
            "fraction_matched_yhtp",
            "separation_yhtp",
            "precision_jaccard_yhtp",
            "recall_jaccard_yhtp",
            "f_measure_jaccard_yhtp",
            "sensitivity_sgd",
            "positive_predicted_value_sgd",
            "accuracy_sgd",
            "fraction_matched_sgd",
            "separation_sgd",
            "precision_jaccard_sgd",
            "recall_jaccard_sgd",
            "f_measure_jaccard_sgd",
            "filepath",
            "num_clusters"
        ]
    )

    print(f"Seconds: {(time.time()- time_0) / 60}")

    df_filepath = lib.files.make_filepath_to_clusters(f"{network_name}.pc2p.global.csv")
    # df = pd.read_csv(df_filepath)
    df.to_csv(df_filepath)
    # df3 = pd.concat([df, df2], ignore_index=True)
    #
    # df3_filepath = lib.files.make_filepath_to_clusters(f"{network_name}.mcl.global.concat.csv")
    # df3.to_csv(df3_filepath)
