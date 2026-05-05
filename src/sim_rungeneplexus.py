import geneplexus

def run_plex(genelist,gsc,curdisease,neggenes,neg_bool=True):
    gp = geneplexus.GenePlexus(net_type="STRING", features="Adjacency", gsc=gsc,file_loc="/mnt/home/mckimale/play/plexpy/data_2022_05_17/")
    gp.load_genes(genelist,disease=curdisease,neg_genes="",orig_calc_negs=neg_bool)
    mdl_weights, df_probs, avgps, cvframe = gp.fit_and_predict(cross_validate=True)
    return df_probs, avgps
