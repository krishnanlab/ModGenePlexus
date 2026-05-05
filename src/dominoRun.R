library(tidyverse)
library(parallel)

formatoutput<-function(domout,disease){
  colnames(domout)[1]="Gene"
  domout = domout %>% mutate(Cluster = as.character(1:nrow(domout)))
  domout=domout %>% separate_rows(Gene,convert=F)
  domout = domout %>% filter(Gene != "")
  domout = domout %>% filter(Gene!="e")
  
  cluscounts = domout %>% group_by(Cluster) %>% count %>% arrange(n%>%desc) %>% ungroup
  cluscounts$NewCluster=as.character(1:nrow(cluscounts))
  clusassign = cluscounts %>% select(-n) 
  domout = domout %>% left_join(clusassign) %>% select(-Cluster)
  colnames(domout)=c("Gene","Cluster")
  domout$Disease=disease
  
  domout=domout %>% mutate(Cluster = paste("domino_", Cluster, sep = ""))
  
  return(domout %>% select(Disease, Cluster, Gene))
}

runDomino<-function(disease,typediff,typenet,dissplits){
  
  network = paste0("part0slices/",typenet,".sif")
  slices = paste0("part0slices/",typenet,"_slices.txt")
  outputfolder = paste0("./domino/output/strict_normal/",typediff,"/")
  
  agfilename = paste0("./domino/input/strict_normal/",typediff,"/",disease,".txt")
  ag = filter(dissplits,Name==disease)$Train %>% strsplit(", ") %>% unlist %>% paste0("e_",.) %>% as_tibble
  
  write_tsv(ag,agfilename,col_names=F)
  
  domcall = paste0("domino --active_genes_files ",agfilename," --network_file ",network,
                   " --slices_file ",slices," --output_folder ",outputfolder, " --visualization False")
  system(domcall)
  domout = read_tsv(paste0(outputfolder,"/",disease,"/modules.out"),col_names=F)
  if(nrow(domout)==0){
    return(tibble())
  }
  newdomout=formatoutput(domout,disease)
  
  return(newdomout)
}

make_dom_minusprop<-function(disease,output,dissplits){
  cursplit=dissplits %>% filter(Name==disease)
  curtrain=cursplit$Train %>% strsplit(", ") %>% unlist
  curtest=cursplit$Test %>% strsplit(", ") %>% unlist
  curgenes=c(curtrain,curtest)
  domminusprop=filter(output,Disease==disease) %>% filter(.,Gene %in% curgenes)
  domminusprop = domminusprop %>% mutate(Cluster = str_replace(Cluster,"domino_","noprop_"))
  return(domminusprop)
}
make_dom_minustest<-function(disease,output,dissplits){
  cursplit=dissplits %>% filter(Name==disease)
  curtrain=cursplit$Train %>% strsplit(", ") %>% unlist
  curgenes=curtrain
  domminusprop=filter(output,Disease==disease) %>% filter(.,Gene %in% curgenes)
  domminusprop = domminusprop %>% mutate(Cluster = str_replace(Cluster,"domino_","notest_"))
  return(domminusprop)
  
  
}

typenet="strict_normal"
typediffs=c("diff","creeds_drug","creeds_gene")

for(typediff in typediffs){
  
  dissplits=read_tsv(paste0("../data/splits/",typediff,"_splits_string.tsv"))
  dissplits=dissplits %>% arrange(Name)
  diseases=unique(dissplits$Name)
  
  out=mclapply(diseases,runDomino,typediff,typenet,dissplits,mc.cores=detectCores()) %>% bind_rows
  write_tsv(out,paste0("../data/clusters/","strict_normal_",typediff,"_domino.tsv"))
  dis=out$Disease %>% unique
  nopropout=mclapply(dis,make_dom_minusprop,out,dissplits,mc.cores=detectCores()) %>% bind_rows
  write_tsv(nopropout,paste0("./../data/clusters/","strict_normal_",typediff,"_noprop.tsv"))

  
}


