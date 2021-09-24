import gzip

PATH_TO_PROJECT = "..\\"
PATH_TO_DATA = f"{PATH_TO_PROJECT}/data"
PATH_TO_GRAPHS = f"{PATH_TO_DATA}/graphs"
PATH_TO_NETWORKS = f"{PATH_TO_DATA}/networks"
PATH_TO_CLUSTERS = f"{PATH_TO_DATA}/clusters"
PATH_TO_MCL_CLUSTERS = f"{PATH_TO_CLUSTERS}/mcl"


def make_filepath_to_data(file_name):
    return f"{PATH_TO_DATA}/{file_name}"


def make_filepath_to_graphs(file_name):
    return f"{PATH_TO_GRAPHS}/{file_name}"


def make_filepath_to_networks(file_name):
    return f"{PATH_TO_NETWORKS}/{file_name}"


def make_filepath_to_clusters(file_name):
    return f"{PATH_TO_CLUSTERS}/{file_name}"


def make_filepath_to_mcl_clusters(file_name):
    return f"{PATH_TO_MCL_CLUSTERS}/{file_name}"


def read_zipped_filelines(filepath):
    with gzip.open(filepath, 'rt') as f:
        lines = f.readlines()
    return lines


def read_filelines(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    return [line.strip() for line in lines]  # Strip newline characters


def write_lines_to_file(filepath, lines):
    with open(filepath, 'w+') as f:
        f.write('\n'.join(lines))


def preprocess_txt_gz_file(filepath):
    pass
    # lines = read_zipped_filelines(filepath)
    # lines = lines[1:]  # Remove header
    # write_lines_to_file(make_unzipped_filepath(file_name), lines)


def read_data(file_name):
    return read_filelines(make_filepath_for_nb(file_name))
