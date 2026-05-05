library(tidyverse)
library(reticulate)
use_condaenv("base")
library(parallel)
source("moduleFunctions.R")
args = commandArgs(trailingOnly=TRUE)
disease=as.character(args[1])
nettoget=as.character(args[2])
typenet=as.character(args[3])
typediff=as.character(args[4])
pyfile="runGenePlexus.py"
source_python(pyfile)

run_rest<-function(genelist){
  terms=unique(clusters$Term)
  terms=terms[!grepl("_Neutral",terms)]

  all_out=tibble()
  for(term in terms){
    cur=genelist %>% filter(Term==term)
    postrain=unlist(strsplit(filter(cur,Type=="Postrain")$Gene,split=", "))
    negtrain=unlist(strsplit(filter(cur,Type=="Negtrain")$Gene,split=", "))
    out=run_plex(postrain,"DisGeNet",disease,neggenes=negtrain,neg_bool=F)

    coef=out[[2]] %>% as_tibble %>% bind_cols(nodeorder) %>% select(Entrez,value)
    colnames(coef)[2]="Beta"
    out=out[[1]] %>% as_tibble
    out$Disease=disease
    out$Cluster=term
    out = out %>% left_join(coef,by="Entrez")
    all_out=bind_rows(all_out,out)
  }
  return(all_out)
}

getRuns<-function(clusters, typediff){
  restout=run_rest(clusters)
  folder=paste0("../plexout/",nettoget,"_",typenet,"/",disease,"_out.txt")
  write_tsv(restout,folder)
  return("Done")
}

clusters = read_tsv(paste0("../data/pos_neg_combinations/",typenet,"_",typediff,"_neutral_diff_string_allclustype.tsv"),
                    col_types="cccc") %>% filter(Trait==disease)

nodeorder=read_tsv("../data/NodeOrder_STRING.txt",col_types="c",col_names=F)
colnames(nodeorder)[1]="Entrez"

starttime=Sys.time()
done=getRuns(clusters,typediff)
print(Sys.time()-starttime)
