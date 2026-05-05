import geneplexus

def run_plex(genelist, gsc, curdisease, neggenes, neg_bool=False):
    gp = geneplexus.GenePlexus(net_type="STRING", features="Adjacency", gsc=gsc, file_loc="../data/data_2022_05_17")
    gp.load_genes(genelist, disease=curdisease, neg_genes=neggenes, orig_calc_negs=neg_bool)
    mdl_weights, df_probs, avgps, cvframe = gp.fit_and_predict(cross_validate=False)
    return df_probs, mdl_weights

