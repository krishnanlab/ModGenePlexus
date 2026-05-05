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
        raise ValueError(f"Unknown load method: {load_method!r}")


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


def get_negatives(file_loc, net_type, gsc, pos_genes_in_net, train_genes, train_sub: bool):
    """Compute negative genes by excluding positives and related gene sets.

    Args:
        file_loc: Location of data files.
        net_type: Network type.
        gsc: Gene set collection name.
        pos_genes_in_net: Positive genes present in the network.
        train_genes: High/medium study-bias genes used to subset the universe.
        train_sub: If True, subset universe and GSC to train_genes only.

    """
    uni_genes = load_genes_universe(file_loc, gsc, net_type)
    good_sets = load_gsc(file_loc, gsc, net_type)

    if train_sub:
        mask = np.isin(uni_genes, train_genes)
        uni_genes = uni_genes[mask]
        for k in good_sets:
            newmask = np.isin(good_sets[k]['Genes'], train_genes)
            good_sets[k]['Genes'] = list(compress(good_sets[k]['Genes'], newmask))

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

