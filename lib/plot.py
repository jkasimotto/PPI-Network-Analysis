import colorsys
import itertools
import math

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

import lib.cluster
import lib.constants
import lib.files
import lib.graph
import lib.map_names


def generate_n_colours(n):
    # Using HSV to RGB
    hues = np.linspace(0, 1, n + 2) # Use 2 extra to skip first and last.
    return [colorsys.hsv_to_rgb(hue, 1, 1) for hue in hues][1:-1]
    # Generating permutations of RGB colours.
    # num = math.ceil(math.pow(n, (1 / 3)))
    # colours = [x for x in itertools.product(np.linspace(0, 1, num), repeat=3)]
    # return colours[:n]


def plot_histograms_n_by_2(xs, num_bin, xscale, yscale, titles=[], xlims=None, ylims=None):
    n = len(xs)
    fig, axs = plt.subplots(math.ceil(n / 2), 2)

    hb = np.arange(1, 10)
    hb = np.append(hb, np.logspace(np.log10(10), np.log10(10000), num_bin))

    row = 0
    for i, x in enumerate(xs):
        if n == 1:
            ax = axs[i % 2]
        else:
            ax = axs[row][i % 2]
        ax.hist(x, density=True, log=True, bins=hb)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if xlims is not None:
            ax.set_xlim(1, xlims[i])
        if ylims is not None:
            ax.set_ylim(ylims[i])
        ax.set_title(str(titles[i]))
        if (i + 1) % 2 == 0:
            row += 1

    plt.show()


def hist(x, title, xlabel, ylabel, bins=None, density=False, title_size=20, xlabel_size=20,
         ylabel_size=20, xticks=None, yticks=None, xticks_size=10, yticks_size=10, xlim=None, ylim=None):
    plt.hist(x, density=density, bins=bins)
    plt.title(title, size=title_size)
    plt.xlabel(xlabel, size=xlabel_size)
    plt.ylabel(ylabel, size=ylabel_size)
    plt.xticks(xticks, size=xticks_size)
    plt.yticks(yticks, size=yticks_size)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.show()


def bar(x, y, title, xlabel, ylabel, width=None, density=False, title_size=20, xlabel_size=20, ylabel_size=20,
        xticks=None, yticks=None, xticks_size=10, yticks_size=10):
    plt.bar(x, y, width=width)
    plt.title(title, size=title_size)
    plt.xlabel(xlabel, size=xlabel_size)
    plt.ylabel(ylabel, size=ylabel_size)
    plt.xticks(xticks, size=xticks_size)
    plt.yticks(yticks, size=yticks_size)
    plt.show()


def network_layers(network, subgraph_kwargs, base=None, ax=None):
    # TODO: Introduce node_size into subgraph_KWargs to make icp55 clusters larger.
    """

    :param subgraph_kwargs: A list of [(subgraph, colour)] from background to foreground.
    :param base: The base graph such that only subgraphs of will be plotted.
    :return:
    """

    # If no base is provided we take the subgraph of all desired nodes AND a neighbourhood of 2 around them so we can
    # see their connections.
    if base is None:
        subgraphs = [kwarg['graph'] for kwarg in subgraph_kwargs]
        base = network.subgraph(itertools.chain.from_iterable(
            subgraph.nodes() for subgraph in subgraphs
        ))

    pos = nx.spring_layout(base)
    # pos = nx.random_layout(base)
    # pos = nx.shell_layout(base)
    # pos = nx.spiral_layout(base)
    # pos = nx.planar_layout(base)
    # pos = nx.bipartite_layout(base)

    nx.draw(base, pos=pos, node_color="black", node_size=1, ax=ax)
    for kwarg in subgraph_kwargs:
        graph = kwarg['graph']
        kwarg = {k: v for (k, v) in kwarg.items() if k != 'graph'}  # Remove the graph kwarg before pasing to nx.draw.
        nx.draw(
            graph,
            pos=pos,
            **kwarg
        )


def targets_with_clusters(network_name, clusters_name, targets, all_shorps=False, top_size=200, target_size=200,
                          base_size=200, ax=None, top_colour='pink', base_colour='yellow'):
    """
    This function plots a target protein within its cluster in relation to ICP55 and PIM1
    :param network_name: The name of the network we are in
    :param clusters_name: The name of the clusters we are in
    :param targets: A list of proteins we are targeting
    :param all_shorps: If true, plot all shortest paths between all targets. If False plot shortest paths only from ICP55/PIM1 to targets.
    :param top_size: Size of ICP55 / PIM1 nodes
    :param target_size: Size of target proteins.
    :param base_size: Base of proteins on the shortest paths. (These appear as squares).
    :param top_colour: Colour of ICP55 / PIM1.
    :param base_colour: Colour of proteins on the shortest paths.
    :param ax: An axis to plot on if you wish.
    :return:
    """
    network_filename = lib.files.make_network_filename(network_name)
    network_filepath = lib.files.make_filepath_to_networks(network_filename)
    network = lib.graph.read_weighted_edgelist(network_filepath)

    clusters_filename = lib.files.make_clusters_filename(network_name, clusters_name)
    clusters_filepath = lib.files.make_filepath_to_clusters(clusters_filename)
    clusters = lib.cluster.read_csv(clusters_filepath)

    # This is our subgraph layer of icp55 and pim1
    top_layer = network.subgraph([lib.constants.ICP55, lib.constants.PIM1])

    # This is our layer of target proteins.
    # I create a subgraph for each one because we want the target layers to have the same colour as their cluster
    target_layers = [network.subgraph([target]) for target in targets]

    # This gets a list of all the "cluster-ids" aka indexes we will need.
    all_clusters = [lib.cluster.cluster_idxs_with_protein(clusters, target)[0] for target in targets]
    # There may be targets in the same cluster so I create a set of unique cluster indexes here.
    # This is used to assign each target the same colour as it's cluster.
    unique_clusters = list(set(all_clusters))
    # This creates a subgraph layer for each cluster we will plot.
    cluster_layers = [network.subgraph(clusters[idx]) for idx in unique_clusters]

    # This creates a base layer of a neighbourhood around icp55 and pim1 to show the links

    # This creates a base layer of nodes in the shortest paths between ICP55, PIM1 and targets.
    base_layers = []
    base_layers += list(nx.shortest_path(network, lib.constants.ICP55, lib.constants.PIM1))
    base_layers += list(
        itertools.chain.from_iterable([nx.shortest_path(network, lib.constants.ICP55, target) for target in targets]))
    base_layers += list(
        itertools.chain.from_iterable([nx.shortest_path(network, lib.constants.PIM1, target) for target in targets]))
    if all_shorps:
        base_layers += list(itertools.chain.from_iterable(itertools.chain.from_iterable(
            [nx.shortest_path(network, target1, target2) for target1 in targets] for target2 in targets
        )))
    base_layer = network.subgraph(base_layers)

    # This creates a base layer of nodes in the 2-shell around ICP55/PIM1
    # base_layer = network.subgraph(
    #     list(lib.graph.get_neighbourhood(network, lib.constants.ICP55, neighbourhood_dist).nodes())
    #     # list(lib.graph.get_neighbourhood(network, lib.constants.PIM1, neighbourhood_dist).nodes())
    # )

    # We will colour icp55 and pim1 pink with large nodes and labels
    top_layer_kwargs = {
        'graph': lib.graph.rename_with_gene_names(top_layer),
        'node_color': top_colour,
        'node_size': top_size,
        'with_labels': True
    }

    # We will give each cluster a visibly separate colour.
    colours = generate_n_colours(len(cluster_layers))
    cluster_layer_kwargs = [{
        'graph': lib.graph.rename_with_gene_names(cluster),
        'node_color': colours[i],
        'node_size': 10
    } for i, cluster in enumerate(cluster_layers)]

    # This lambda function get's the index of our target's cluster to get the same colour.
    get_target_colour = lambda i: unique_clusters.index(all_clusters[i])
    target_layer_kwargs = [{
        'graph': lib.graph.rename_with_gene_names(target),
        'node_color': colours[get_target_colour(i)],
        'node_size': target_size,
        'with_labels': True
    } for i, target in enumerate(target_layers)]

    base_layer_kwargs = {
        'graph': lib.graph.rename_with_gene_names(base_layer),
        'node_color': base_colour,
        'node_size': base_size,
        'with_labels': True,
        'node_shape': 's'
    }

    # Because we renamed all the subgraphs we must rename the original network too for the plot function.
    network_renamed = lib.graph.rename_with_gene_names(network)

    # Now we want to plot the corresponding layers in the following order:
    # Base
    # Clusters
    # Targets
    # Top layer
    network_layers(network_renamed,
                   [base_layer_kwargs,
                    *cluster_layer_kwargs,
                    *target_layer_kwargs,
                    top_layer_kwargs],
                   ax=ax)
    
    
#Plot a given cluster and its shortest paths to icp55/pim1, if within a threshold path length
def closest_clusters_with_paths_vis(network_name, clusters_name, master_df_name, cluster_id, shorp_threshold, all_shorps=False, top_size=200, target_size=200,
                          base_size=200, ax=None, top_colour='pink', base_colour='yellow',
                                   save = False,
                                   save_directory = "..\\data/cluster_validation/clusters_and_shorps/"):
    """
    This function plots a target protein within its cluster in relation to ICP55 and PIM1
    :param network_name: The name of the network we are in
    :param clusters_name: The name of the clusters we are in
    :param targets: A list of proteins we are targeting
    :param all_shorps: If true, plot all shortest paths between all targets. If False plot shortest paths only from ICP55/PIM1 to targets.
    :param top_size: Size of ICP55 / PIM1 nodes
    :param target_size: Size of target proteins.
    :param base_size: Base of proteins on the shortest paths. (These appear as squares).
    :param top_colour: Colour of ICP55 / PIM1.
    :param base_colour: Colour of proteins on the shortest paths.
    :param ax: An axis to plot on if you wish.
    :return:
    """
    network_filename = lib.files.make_network_filename(network_name)
    network_filepath = lib.files.make_filepath_to_networks(network_filename)
    network = lib.graph.read_weighted_edgelist(network_filepath)

    clusters_filename = lib.files.make_clusters_filename(network_name, clusters_name)
    clusters_filepath = lib.files.make_filepath_to_clusters(clusters_filename) 
    clusters = lib.cluster.read_csv(clusters_filepath)
    
    master_df = pd.read_csv(master_df_name)
    
    
    ##Get proteins in cluster than are within shorp_threshold
    
    #Set ip top_layer_nodes for later
    top_layer_nodes = []
    
    #Icp55
    icp55_connectors = master_df.loc[(master_df["icp55_shell"] <= shorp_threshold) &
                                     (master_df["cluster_id"] == cluster_id), "protein"]
    icp55_connector_paths = []
    if len(icp55_connectors) == 0:
        icp55_connector_path_nodes = []
    else:
        top_layer_nodes.append(lib.constants.ICP55)
        
        for connector in icp55_connectors:
            icp55_connector_paths.extend(nx.all_shortest_paths(network,
                                                              lib.constants.ICP55,
                                                              connector))
        icp55_connector_path_nodes = list(set(itertools.chain.from_iterable(icp55_connector_paths)))
    
    
    #Pim1
    pim1_connectors = master_df.loc[(master_df["pim1_shell"] <= shorp_threshold) &
                                     (master_df["cluster_id"] == cluster_id), "protein"]
    pim1_connector_paths = []
    if len(pim1_connectors) == 0:
        pim1_connector_path_nodes = []
    else:
        top_layer_nodes.append(lib.constants.PIM1)
        
        for connector in pim1_connectors:
            pim1_connector_paths.extend(nx.all_shortest_paths(network,
                                                              lib.constants.PIM1,
                                                              connector))
        pim1_connector_path_nodes = list(set(itertools.chain.from_iterable(pim1_connector_paths)))
        
    #Combine both
    connector_path_nodes = list(itertools.chain.from_iterable([icp55_connector_path_nodes,
                           pim1_connector_path_nodes]))
    
    
    # This is our subgraph layer of icp55 and pim1
    
    top_layer = network.subgraph(top_layer_nodes)
    
    #This is our layer of clusters
    cluster_layers = [network.subgraph(clusters[cluster_id])]
    
    #This is our base layer of connector proteins
    base_layer = network.subgraph(connector_path_nodes)
    
    #This is a network from combined all
    all_nodes = list(top_layer.nodes())
    all_nodes.extend(cluster_layers[0].nodes())
    all_nodes.extend(base_layer.nodes())
    all_node_network = network.subgraph(all_nodes)
    
    # We will colour icp55 and pim1 pink with large nodes and labels
    top_layer_kwargs = {
        'graph': lib.graph.rename_with_gene_names(top_layer),
        'node_color': top_colour,
        'node_size': top_size,
        'with_labels': True
    }

    # We will give each cluster a visibly separate colour.
    colours = generate_n_colours(len(cluster_layers))
    cluster_layer_kwargs = [{
        'graph': lib.graph.rename_with_gene_names(cluster),
        'node_color': colours[i],
        'node_size': 10
    } for i, cluster in enumerate(cluster_layers)]

    base_layer_kwargs = {
        'graph': lib.graph.rename_with_gene_names(base_layer),
        'node_color': base_colour,
        'node_size': base_size,
        'with_labels': True,
        'node_shape': 's'
    }

    # Because we renamed all the subgraphs we must rename the original network too for the plot function.
    network_renamed = lib.graph.rename_with_gene_names(network)

    # Now we want to plot the corresponding layers in the following order:
    # Base
    # Clusters
    # Targets
    # Top layer
    f = plt.figure()
    
    network_layers(network_renamed,
                   [base_layer_kwargs,
                    *cluster_layer_kwargs,
                    top_layer_kwargs],
                   ax=ax)
    
    #Save subgraph for input into STRING
    if save:
        
        #Save image
        f.savefig(save_directory + "cluster" + str(cluster_id) + 'shorps_graph.png')
        
        #Save nodes for STRING
        with open(save_directory + "cluster" + str(cluster_id) + 'shorps_nodes.txt', 'w') as f:
            for node in list(all_node_network.nodes()):
                f.write("%s\n" % node)
 
    
                   