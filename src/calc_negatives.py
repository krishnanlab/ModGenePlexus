import numpy as np
import json
import os.path as osp
from typing import Dict
from typing import Any
from typing import Literal
from scipy.stats import hypergeom
from itertools import compress

def _load_json_file(file_loc: str, file_name: str) -> Dict[str, Any]:
    """Load JSON into dictionary.

    Args:
        file_loc: Location of data files.
        file_name: Name of the file.

    """
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)
    return json.load(open(file_path, "rb"))


def check_file(path: str):
    """Check existence of a file.

    Args:
        path: Path to the file.

    Raises:
        FileNotFoundError: if file not exist.

    """
    if not osp.isfile(path):
        raise FileNotFoundError(path)




# def load_gsc(
#     file_loc: str,
#     gsc: config.GSC_TYPE,
#     net_type: config.NET_TYPE,
# ) -> config.GSC_DATA_TYPE:
#     """Load gene set collection dictionary.

#     Args:
#         file_loc: Location of data files.
#         target_set: Target gene set collection.
#         net_type: Network used.

#     """
#     file_name = f"GSC_{gsc}_{net_type}_GoodSets.json"
#     return _load_json_file(file_loc, file_name)

def load_gsc(
    file_loc: str,
    gsc: str,
    net_type: str,
):
    """Load gene set collection dictionary.

    Args:
        file_loc: Location of data files.
        target_set: Target gene set collection.
        net_type: Network used.

    """
    file_name = f"GSC_{gsc}_{net_type}_GoodSets.json"
    return _load_json_file(file_loc, file_name)



def _load_np_file(
    file_loc: str,
    file_name: str,
    load_method: Literal["npy", "txt"],
) -> np.ndarray:
    """Check np file existence and load.

    Args:
        file_loc: Location of data files.
        file_name: Name of the file.
        load_method: How to load the file ('npy' or 'txt').

    """
    file_path = osp.join(file_loc, file_name)
    check_file(file_path)

    if load_method == "npy":
        return np.load(file_path)
    elif load_method == "txt":
        return np.loadtxt(file_path, dtype=str)
    else:
        raise ValueError(f"Unknwon load method: {load_method!r}")


# def load_genes_universe(
#     file_loc: str,
#     gsc: config.GSC_TYPE,
#     net_type: config.NET_TYPE,
# ) -> np.ndarray:
#     """Load gene universe a given network and GSC.

#     Args:
#         file_loc: Location of data files.
#         gsc: Gene set collection.
#         net_type: Network used.

#     """
#     file_name = f"GSC_{gsc}_{net_type}_universe.txt"
#     return _load_np_file(file_loc, file_name, load_method="txt")

def load_genes_universe(
    file_loc: str,
    gsc: str,
    net_type: str,
):
    """Load gene universe a given network and GSC.

    Args:
        file_loc: Location of data files.
        gsc: Gene set collection.
        net_type: Network used.

    """
    file_name = f"GSC_{gsc}_{net_type}_universe.txt"
    return _load_np_file(file_loc, file_name, load_method="txt")


#train_genes are high/med seen genes. train_sub says to use that to subset Remy's files or not
def get_negatives(file_loc, net_type, gsc, pos_genes_in_net, train_genes,train_sub: bool):
    uni_genes = load_genes_universe(file_loc, gsc, net_type)
    good_sets = load_gsc(file_loc, gsc, net_type)

    #subset the above to only have
    if train_sub == True:
        #an array
        #https://numpy.org/doc/stable/reference/generated/numpy.isin.html
        mask=np.isin(uni_genes,train_genes)
        uni_genes=uni_genes[mask]
        #dictionary, needs to subset the values of each key
        #multiple values, Name and Genes
        #the values are stored as lists, so loop for each key, do the above
        for k in good_sets:
            newmask=np.isin(good_sets[k]['Genes'],train_genes)
        	#test_dict[k]['Genes']=good_sets[k][newmask]
            good_sets[k]['Genes']=list(compress(good_sets[k]['Genes'],newmask))
        #good_sets
        #test_dict={key:val for key, val in good_sets.items() if val in train_genes}

    M = len(uni_genes)
    N = len(pos_genes_in_net)
    genes_to_remove = pos_genes_in_net
    for akey in good_sets:
        n = len(good_sets[akey]["Genes"])
        k = len(np.intersect1d(pos_genes_in_net, good_sets[akey]["Genes"]))
        pval = hypergeom.sf(k - 1, M, n, N)
        if pval < 0.05:
            genes_to_remove = np.union1d(genes_to_remove, good_sets[akey]["Genes"])
    negative_genes = np.setdiff1d(uni_genes, genes_to_remove)
    return negative_genes

