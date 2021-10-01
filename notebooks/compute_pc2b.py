import lib.graph
import lib.cluster
import lib.files

def go():
    import os
    import sys
    print(os.getcwd())
    print(sys.path)
    network_name = 'icp55-cc-900-inv'
    network_filepath = lib.files.make_filepath_to_networks(f'{network_name}.txt')
    network = lib.graph.read_weighted_edgelist(network_filepath)

    # Because of the issue of hub nodes, I'm going to remove all nodes with high degree
    # nodes_to_remove = [node for node in network if network.degree()[node] >= 50]
    # network.remove_nodes_from(nodes_to_remove)
    # network = lib.graph.get_largest_cc_with_node(network, lib.constants.ICP55)
    # print(f"Numbewr of nodes {len(network)}")

    clusters = lib.cluster.run_pc2p(network)
    # Save the clusters to a file
    filepath = lib.files.make_filepath_to_clusters(f"{network_name}.pc2p.csv")
    lib.cluster.write_to_file(filepath, clusters)
if __name__ == "__main__":
    go()
