import numpy as np
import pandas as pd
from itertools import combinations
import networkx as nx
import csv
import argparse

def parse_GSC_file(GSC, min_set_size, max_set_size, net):
    # This section will parse the input file and build two
    # dictionaries to be used later

    # data_dict = The key is the ontology ID and the value is a list of genes assoicated with that ID
    # names_dict = The key is the ontology ID and the value is the more understandable name
    # It is in this section that the max and min fo the set size is applied
    if 'GO' in GSC:
        FN = "../data/propogated/annotations_" + net + "_GO.tsv"
    elif 'mondo' == GSC:
        FN = '../data/propogated/propogated_annotations_mondo_' + net + '.tsv'
    elif 'diff' == GSC:
        FN='../data/differential/input_gmts/creeds_'+net+'_disease.tsv'
    elif 'creeds_drug' == GSC:
        FN='../data/differential/input_gmts/creeds_'+net+'_drug.tsv'
    elif 'creeds_gene' == GSC:
        FN='../data/differential/input_gmts/creeds_'+net+'_gene.tsv'
    elif 'magma_1e2' == GSC:
        FN='../data/gwas/input_gmts/GWAS_1e2.tsv'
    elif 'magma_1e5' == GSC:
        FN='../data/gwas/input_gmts/GWAS_1e5.tsv'
    elif 'magma_1e8' == GSC:
        FN='../data/gwas/input_gmts/GWAS_1e8.tsv'
    elif 'magma_1e1' == GSC:
        FN='../data/gwas/input_gmts/GWAS_1e1.tsv'
    elif 'magma_5e2' == GSC:
        FN='../data/gwas/input_gmts/GWAS_5e2.tsv'
    else:
        print('Not a valid geneset collection')
    data_dict = {}
    names_dict = {}
    all_annot_genes = []
    if FN.strip().split('.')[-1] == 'tsv':
        with open(FN, 'r') as f:
            for idx, line in enumerate(f):
                if idx == 0:
                    continue
                line = line.strip().split('\t')
                ID = line[0]
                genes = line[3].strip().split(', ')
                all_annot_genes = list(set(all_annot_genes + genes))
                if (len(genes) >= min_set_size) and (len(genes) <= max_set_size):
                    data_dict[ID] = genes
                    names_dict[ID] = line[1]
    elif FN.strip().split('.')[-1] == 'gmt':
        with open(FN, 'r') as f:
            for idx, line in enumerate(f):
                line = line.strip().split('\t')
                ID = line[0]
                genes = line[2:]
                all_annot_genes = list(set(all_annot_genes + genes))
                if (len(genes) >= min_set_size) and (len(genes) <= max_set_size):
                    data_dict[ID] = genes
                    names_dict[ID] = line[1].replace(' ', '_')
    print('The number of IDs in the original file is',idx-1)
    print('The number of IDs after min and max threshold is',len(data_dict))
    print('The total number of genes annotted to the GSC is',len(all_annot_genes))
    return data_dict, names_dict, all_annot_genes
    

def make_edgelist(data_dict, jac_thresh, over_thresh):
    keys = list(data_dict)
    all_key_combos = combinations(keys, 2)

    edges = []
    for idx, akey_combo in enumerate(all_key_combos):
        list0 = data_dict[akey_combo[0]]
        list1 = data_dict[akey_combo[1]]
        inter = len(np.intersect1d(list0,list1,assume_unique=True))
        union = len(np.union1d(list0,list1))
        min_set_size = np.min(np.array([len(list0),len(list1)]))
        jac = inter/union
        over = inter/min_set_size
        if (jac >= jac_thresh) and (over >= over_thresh):
            edges.append([akey_combo[0],akey_combo[1]])
    df_edges = pd.DataFrame(edges, columns=['Node1', 'Node2'])
    print('The number of possible edges tried is', idx - 1)
    return keys, df_edges


def get_CCs(df_edges, keys):
    G = nx.convert_matrix.from_pandas_edgelist(df_edges, 'Node1', 'Node2')
    nodelist = list(G.nodes())
    good_ids = np.setdiff1d(keys, nodelist)
    print('The number of nodes in the graph is', len(nodelist))
    print('The number of edges in the graph is', df_edges.shape[0])
    print('The number of IDs with no edges is', len(good_ids))
    CCs = sorted(nx.connected_components(G), key=len, reverse=True)
    print('The number of connected components is', len(CCs))
    return G, CCs, good_ids

def find_sets_from_CCs(CCs, G, data_dict):
    good_from_CCs = []
    for aCC in CCs:
        aCC_tmp = np.array(list(aCC))
        aCC_tmp.sort()
        print(aCC_tmp)
        aCC_scores = []
        for idx1, anode1 in enumerate(aCC_tmp):
            score_tmp = 0
            neighs_tmp = list(G.neighbors(anode1))
            for anode2 in neighs_tmp:
                score_tmp = score_tmp + len(np.intersect1d(data_dict[anode1], data_dict[anode2])) / len(data_dict[anode2])
            aCC_scores.append(score_tmp)
        aCC_scores = np.array(aCC_scores)

        sorted_args = np.flip(np.argsort(aCC_scores))
        aCC_tmp_sorted = aCC_tmp[sorted_args]
        while len(aCC_tmp_sorted) > 0:
            good_from_CCs.append(aCC_tmp_sorted[0]) # add highest score node to useable from CC ID list
            neighs_tmp = np.array(list(G.neighbors(aCC_tmp_sorted[0])))
            aCC_tmp_sorted = np.setdiff1d(aCC_tmp_sorted,aCC_tmp_sorted[0]) # remove highest score node from the list
            aCC_tmp_sorted = np.setdiff1d(aCC_tmp_sorted,neighs_tmp) # remove all neighbors of highest score node from the list
    print('The number of genesets from the CCs is',len(good_from_CCs))
    return good_from_CCs

def get_sets_df(good_ids, good_from_CCs, names_dict, data_dict):
    good_ids = np.union1d(good_ids, good_from_CCs)
    print('The final number of IDs to use is',len(good_ids))

    final_results = []
    final_results2 = []
    for aid in good_ids:
        name_tmp = names_dict[aid]
        num_genes_tmp = len(data_dict[aid])
        genes_tmp = '\t'.join(data_dict[aid])
        genes_tmp2 = ', '.join(data_dict[aid])
        final_results.append([aid, name_tmp, num_genes_tmp, genes_tmp])
        final_results2.append([aid, name_tmp, num_genes_tmp, genes_tmp2])
    colnames = ['ID', 'Name', 'Gene_Count', 'Gene_IDs']
    df_sets = pd.DataFrame(final_results, columns=colnames)
    df_sets2 = pd.DataFrame(final_results2, columns=colnames)
    df_sets = df_sets.sort_values(by=['Gene_Count'], ascending=False)
    df_sets2 = df_sets2.sort_values(by=['Gene_Count'], ascending=False)
    print('The shape of the df_final is',df_sets.shape)

    subset_annot_genes = []
    data_dict = {}
    names_dict = {}
    for idx, line in df_sets.iterrows():
        ID = line[0]
        genes = line[3].strip().split('\t')
        subset_annot_genes = list(set(subset_annot_genes + genes))
    df_sets = df_sets2
    return df_sets, subset_annot_genes
    

def get_gene_bins(GSC, all_annot_genes):
    df_counts = pd.read_csv('../data/20200923_gene2pubmed.txt', sep='\t')
    if 'GO' in GSC:
        df_counts = df_counts[df_counts['#tax_id'] == names_pubmed_dict['hs']]
    elif GSC == 'mondo':
        df_counts = df_counts[df_counts['#tax_id']==names_pubmed_dict['hs']]
    df_counts = df_counts.groupby(['GeneID']).size().reset_index(name='counts')
    df_counts = df_counts.sort_values(by=['counts'],ascending=False)
    df_counts['GeneID'] = df_counts['GeneID'].apply(str)
    print('The number of genes that have a count is',df_counts.shape[0])
    df_counts = df_counts[df_counts['GeneID'].isin(all_annot_genes)]
    in_counts = df_counts['GeneID'].to_list()
    not_in_counts = np.setdiff1d(all_annot_genes, in_counts)
    print('The number of genes that have an annotation and are in the count file is', len(in_counts))
    print('The number of genes that have an annotation and are not in the count file is', len(not_in_counts))
    final_count_genes = np.array(list(in_counts) + list(not_in_counts))
    print('The number of total number of genes before intersection with networks genes is',len(final_count_genes))
    num_per_bin = int(len(final_count_genes)/3)
    top_genes = final_count_genes[0:num_per_bin]
    med_genes = final_count_genes[num_per_bin:2*num_per_bin]
    low_genes = final_count_genes[2*num_per_bin:3*num_per_bin]
    gene_bins_dict = {'top_genes': top_genes, 'med_genes': med_genes, 'low_genes': low_genes}
    print('The number of genes in the high, medium and low bin is', len(top_genes), len(med_genes), len(low_genes))
    return gene_bins_dict


def make_label_dict(df_sets, gene_bins_dict):
    label_dict = {}
    for idx in range(df_sets.shape[0]):
        df_tmp = df_sets.iloc[[idx]]
        genes_tmp = df_tmp['Gene_IDs'].tolist()
        genes_tmp = genes_tmp[0].strip().split(', ')
        genes_tmp = np.array([str(item) for item in genes_tmp])
        train_pos_tmp = np.intersect1d(genes_tmp, gene_bins_dict['top_genes'])
        train_pos_tmp = np.append(train_pos_tmp, np.intersect1d(genes_tmp, gene_bins_dict['med_genes']))

        test_pos_tmp = np.intersect1d(genes_tmp, gene_bins_dict['low_genes'])
        if (len(train_pos_tmp) >= 10) and (len(test_pos_tmp)>=10):
            ID_tmp = df_tmp['ID'].tolist()[0]
            label_dict[ID_tmp] = {}
            name_tmp = df_tmp['Name'].tolist()[0]
            label_dict[ID_tmp]['Name'] = name_tmp
            label_dict[ID_tmp]['Train-Positive'] = train_pos_tmp
            label_dict[ID_tmp]['Test-Positive'] = test_pos_tmp
            train_neg_tmp = np.setdiff1d(gene_bins_dict['top_genes'], train_pos_tmp)
            test_neg_tmp = np.setdiff1d(gene_bins_dict['low_genes'], test_pos_tmp)
            label_dict[ID_tmp]['Train-Negative'] = train_neg_tmp
            label_dict[ID_tmp]['Test-Negative'] = test_neg_tmp
    print('The number of sets with at least 10 positive in train and test is',len(label_dict))
    return label_dict



  
if __name__ == "__main__":

    names_pubmed_dict = {'hs': 9606, 'mm': 10090, 'ce': 6239,
                         'dm': 7227, 'dr': 7955, 'sc': 559292}

    parser = argparse.ArgumentParser()
    parser.add_argument('--network', type=str, required=True)
    parser.add_argument('--gs', type=str, required=True)
    args = parser.parse_args()
    net = args.network
    GSC = args.gs

    if GSC == "GO":
        min_set_size = 20
        max_set_size = 100
    elif GSC == "mondo":
        min_set_size = 50
        max_set_size = 1000
    else:
        min_set_size = 100
        max_set_size = 50000

    # Jaccard and overlap thresholds for non-redundancy filtering
    jac_thresh = 0.3
    over_thresh = 0.5

    print('Getting non-redundant genesets')
    data_dict, names_dict, all_annot_genes = parse_GSC_file(GSC, min_set_size, max_set_size, net)
    keys, df_edges = make_edgelist(data_dict, jac_thresh, over_thresh)
    G, CCs, good_ids = get_CCs(df_edges, keys)
    good_from_CCs = find_sets_from_CCs(CCs, G, data_dict)
    df_sets, subset_annot_genes = get_sets_df(good_ids, good_from_CCs, names_dict, data_dict)

    gene_bins_dict = get_gene_bins(GSC, all_annot_genes)

    fr = pd.DataFrame.from_dict(gene_bins_dict)
    dictowrite = "../data/study_bins/" + net + "_" + GSC + "_bins.tsv"
    fr.to_csv(dictowrite, sep="\t")

    label_dict = make_label_dict(df_sets, gene_bins_dict)

    results_t = []
    for idx, akey in enumerate(label_dict):
        trainpositives = ', '.join([*label_dict[akey]['Train-Positive']])
        testpos = ', '.join([*label_dict[akey]['Test-Positive']])
        results_t.append([akey, label_dict[akey]['Name'], trainpositives, testpos])
    col_names = ["ID", "Name", "Train", "Test"]
    dff = pd.DataFrame(results_t, columns=col_names)

    dirtowrite = "../data/splits/" + GSC + "_splits_" + net + ".tsv"
    dff.to_csv(dirtowrite, sep="\t", header=True, index=False)
