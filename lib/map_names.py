import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import lib.constants
import lib.files


def map_names_descriptions(names,
                          name_type,
                          output):
    """
        This function maps yeast protein gene names to  formal names, and also fetches descriptions of  proteins
        :param names: A list of names that you will input (accepts a list of length one e.g. ["ICP55"]
        :param name_type: the type of the names that you will input. This can be "gene_name" (e.g. ICP55), "systematic_name_nonumber" (e.g. YER078C), or "systematic_name_wnumber" (e.g. 4932.YER078C)
        :param output: what you want to map the inputted names to. Accepts the same values as for name_type. "description" will give a description of the protein, and "all" will return everything as a pd array
        :return:
        """
    
    ###Get name matching df
    name_mapping_df = pd.read_csv(lib.files.make_path_to_dataframes('4932.genenames.systematicnames.descriptions.csv'), header=0)
    
    ###Get df by looping. It is ordered
    for i in range(len(names)):

        if i == 0:
            matched_df = name_mapping_df[name_mapping_df[name_type] == names[i]]

        else:
            matched_df = matched_df.append(name_mapping_df[name_mapping_df[name_type] == names[i]])

    ###Return
    if output == "all":
        return(matched_df)
    else:
        return(list(matched_df[output]))