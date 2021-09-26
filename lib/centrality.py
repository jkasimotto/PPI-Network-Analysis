import networkx as nx


def bridging_centrality(graph,
                        betweenness_approx_nodes=False):
    """

    :param graph: The graph on which bridging centrality will be computed. Needs to be connected.
    :param betweenness_approx_nodes: False specifies that betweenness centrality should be computed exactly (computationally costly). Otherise, an integer value needs to be entered, which is the number of nodes used to approximate betweenness centrality
    :return: A dictionary where keys are node names, and values are bridging centrality values for those nodes
    """
    # Calculate betweenness centrality
    if betweenness_approx_nodes == False:
        betweeness = nx.betweenness_centrality(graph)
    else:
        betweeness = nx.betweenness_centrality(graph,
                                               k=betweenness_approx_nodes)

    # Calculate bridging coefficient (Hwang et al. 2006)
    bridging_coeff = dict()

    for node in list(graph.nodes()):
        bridging_coeff[node] = graph.degree(node) ** (-1) / sum(
            [degree ** (-1) for (node, degree) in graph.degree(list(graph.neighbors(node)))])

    # Calculate brdiging centrality
    bridging_cent = dict()
    for node in bridging_coeff.keys():
        bridging_cent[node] = betweeness[node] * bridging_coeff[node]

    # Output
    return bridging_cent
