#Running genepleus with mondo data
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

#Running for the types in the passed split information
run_rest<-function(genelist){
  #This will give me "AllAssign" and clusters from 1:n
  terms=unique(clusters$Term)
  terms=terms[!grepl("_Neutral",terms)]
  
  all_out=tibble()
  all_coef=tibble()
  for(term in terms){
    cur=genelist %>% filter(Term==term)
    postrain=unlist(strsplit(filter(cur,Type=="Postrain")$Gene,split=", "))
    negtrain=unlist(strsplit(filter(cur,Type=="Negtrain")$Gene,split=", "))
    if(term=="AllAssign" | grepl("_AllClus", term)){
      out=run_plex(postrain,"DisGeNet",disease,neggenes=negtrain,neg_bool=F)
    }else{
      out=run_plex(postrain,"DisGeNet",disease,neggenes=negtrain,neg_bool=F)
    }
    
    #Get stuff, out is predictions and coef are model coefficients
    #Model coefficient gene order is based on node order file
    coef=out[[2]] %>% as_tibble %>% bind_cols(nodeorder) %>% select(Entrez,value)
    colnames(coef)[2]="Beta"
    out=out[[1]] %>% as_tibble
    
    
    out$Disease=disease
    out$Cluster=term
    
    #Join out with coef to add beta column
    out = out %>% left_join(coef,by="Entrez")
    all_out=bind_rows(all_out,out)
  }
  return(all_out)
}

getRuns<-function(clusters, typediff){
  #GenePlexus with all genes. Will have to load stuff since this split information in 
  #different location

  restout=run_rest(clusters)
  

  towrite=restout

  
  #Hacky solution because of defining creeds_disease as "diff"
  folder=paste0("../plexout/",nettoget,"_",typenet,"/",disease,"_out.txt")

  write_tsv(towrite,folder)

  return("Done")
}

#clusters = read_tsv(paste0("/mnt/research/compbio/krishnanlab/projects/module_expansion/data/pos_neg_combinations/",typenet,"_neutral_diff_string_allclustype.tsv"),
#                    col_types="cccc") %>% filter(Trait==disease)

clusters = read_tsv(paste0("../data/pos_neg_combinations/",typenet,"_",typediff,"_neutral_diff_string_allclustype.tsv"),
                    col_types="cccc") %>% filter(Trait==disease)

nodeorder=read_tsv("../data/NodeOrder_STRING.txt",col_types="c",col_names=F)
colnames(nodeorder)[1]="Entrez"

starttime=Sys.time()
done=getRuns(clusters,typediff)
print(Sys.time()-starttime)
