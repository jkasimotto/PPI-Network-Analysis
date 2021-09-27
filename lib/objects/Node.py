class Node:
    """
    Single data structure to keep track of all node data we might need.
    """

    def __init__(self, name, betweeness_centrality=None):
        self.name = name
        self.betweeness_centrality = betweeness_centrality

    def __hash__(self):
        """
        networkx uses the node name as a hash for dictionaries.
        """
        return hash(self.name)
