import colorsys
import itertools
import math

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt

import lib.cluster
import lib.constants
import lib.files
import lib.graph


def generate_n_colours(n):
    # Using HSV to RGB
    hues = np.linspace(0, 1, n + 1)
    return [colorsys.hsv_to_rgb(hue, 1, 1) for hue in hues][:-1]
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
        kwarg = {k:v for (k,v) in kwarg.items() if k != 'graph'} # Remove the graph kwarg before pasing to nx.draw.
        nx.draw(
            graph,
            pos=pos,
            **kwarg
        )


def targets_with_clusters(network_name, clusters_name, targets, ax=None, top_size=200, target_size=200):
    """
    This function plots a target protein within it's cluster in relation to ICP55 and PIM1
    :param network_name: The name of the network we are in
    :param clusters_name: The name of the clusters we are in
    :param targets: A list of proteins we are targeting
    :param ax: An axis to plot on if you wish.
    :param top_size: Size of ICP55 / PIM1 nodes
    :param target_size: Size of target proteins.
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
    base_layer_1 = list(itertools.chain.from_iterable([nx.shortest_path(network, lib.constants.ICP55, target) for target in targets]))
    base_layer_2 = list(itertools.chain.from_iterable([nx.shortest_path(network, lib.constants.PIM1, target) for target in targets]))
    base_layer_3 = list(itertools.chain.from_iterable(itertools.chain.from_iterable(
        [nx.shortest_path(network, target1, target2) for target1 in targets] for target2 in targets
    )))
    base_layer = network.subgraph(base_layer_1 + base_layer_2 + base_layer_3)

    # This creates a base layer of nodes in the 2-shell around ICP55/PIM1
    # base_layer = network.subgraph(
    #     list(lib.graph.get_neighbourhood(network, lib.constants.ICP55, neighbourhood_dist).nodes())
    #     # list(lib.graph.get_neighbourhood(network, lib.constants.PIM1, neighbourhood_dist).nodes())
    # )

    # We will colour icp55 and pim1 pink with large nodes and labels
    top_layer_kwargs = {
        'graph': top_layer,
        'node_color': 'pink',
        'node_size': top_size,
        'with_labels': True
    }

    # We will give each cluster a visibly separate colour.
    colours = generate_n_colours(len(cluster_layers))
    cluster_layer_kwargs = [{
       'graph': cluster,
       'node_color': colours[i],
       'node_size': 10
    } for i, cluster in enumerate(cluster_layers)]

    # This lambda function get's the index of our target's cluster to get the same colour.
    get_target_colour = lambda i: unique_clusters.index(all_clusters[i])
    target_layer_kwargs = [{
        'graph': target,
        'node_color': colours[get_target_colour(i)],
        'node_size': target_size,
        'with_labels': True
    } for i, target in enumerate(target_layers)]

    base_layer_kwargs = {
        'graph': base_layer,
        'node_color': 'black',
        'node_size': 10,
    }

    # Now we want to plot the corresponding layers in the following order:
    # Base
    # Clusters
    # Targets
    # Top layer
    network_layers(network,
                   [base_layer_kwargs,
                    *cluster_layer_kwargs,
                    *target_layer_kwargs,
                    top_layer_kwargs],
                   ax=ax)




