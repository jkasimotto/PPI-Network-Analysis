import colorsys
import itertools
import math
import random

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt

import lib.graph


def generate_n_colours(n):
    # Using HSV to RGB
    hues = np.linspace(0, 1, n)
    return [colorsys.hsv_to_rgb(hue, 1, 1) for hue in hues]
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


def network_layers(network, subgraph_args, base=None, ax=None):
    # TODO: Introduce node_size into subgraph_KWargs to make icp55 clusters larger.
    """

    :param subgraph_kwargs: A list of [(subgraph, colour)] from background to foreground.
    :param base: The base graph such that only subgraphs of will be plotted.
    :return:
    """
    # So far the only arguments catered for are node_colors
    subgraphs, colours = list(zip(*subgraph_args))

    # We only plot nodes on top of the base.
    if base is None:
        base = network.subgraph(itertools.chain.from_iterable(
            subgraph.nodes() for subgraph in subgraphs
        ))
    else:
        subgraphs = [lib.graph.make_from_intersection([base, subgraph]) for subgraph in subgraphs]

    pos = nx.spring_layout(base)
    # pos = nx.random_layout(base)
    # pos = nx.shell_layout(base)
    # pos = nx.spiral_layout(base)
    # pos = nx.planar_layout(base)
    # pos = nx.bipartite_layout(base)

    nx.draw(base, pos=pos, node_color="black", ax=ax)
    for i in range(len(subgraphs)):
        subgraph = subgraphs[i]
        colour = colours[i]
        nx.draw(subgraph, pos=pos, node_color=colour, ax=ax)