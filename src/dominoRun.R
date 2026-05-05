#This facilitates running domino using the command line

library(tidyverse)
library(parallel)
#Domino output is one line each module. Need to edit this. Remove the e_
formatoutput<-function(domout,disease){
  #Change to Gene
  colnames(domout)[1]="Gene"
  #cluster 1 - nrow
  domout = domout %>% mutate(Cluster = as.character(1:nrow(domout)))
  #separate the rows on the comma
  domout=domout %>% separate_rows(Gene,convert=F)
  #The brackets get turned into "". Get rid of these
  domout = domout %>% filter(Gene != "")
  #separate_rows above puts the genes into their own rows. Get rid of these
  domout = domout %>% filter(Gene!="e")
  
  #Renumber clusters based on how big they are
  cluscounts = domout %>% group_by(Cluster) %>% count %>% arrange(n%>%desc) %>% ungroup
  cluscounts$NewCluster=as.character(1:nrow(cluscounts))
  clusassign = cluscounts %>% select(-n) 
  #Remove old cluster
  domout = domout %>% left_join(clusassign) %>% select(-Cluster)
  #Change NewCluser to Cluster
  colnames(domout)=c("Gene","Cluster")
  #Add disease
  domout$Disease=disease
  
  #In Cluster column, modify it to have "domino_" in front
  domout=domout %>% mutate(Cluster = paste("domino_", Cluster, sep = ""))
  
  return(domout %>% select(Disease, Cluster, Gene))
}

#Actual running of the stuff and loading/writing files based on traits
runDomino<-function(disease,typediff,typenet,dissplits){
  
  #Setting up needed files
  network = paste0("part0slices/",typenet,".sif")
  slices = paste0("part0slices/",typenet,"_slices.txt")
  #Domino creates a subfolder in here and the modules inside of it
  outputfolder = paste0("./domino/output/strict_normal/",typediff,"/")
  
  #Create active gene file
  agfilename = paste0("./domino/input/strict_normal/",typediff,"/",disease,".txt")
  ag = filter(dissplits,Name==disease)$Train %>% strsplit(", ") %>% unlist %>% paste0("e_",.) %>% as_tibble
  
  write_tsv(ag,agfilename,col_names=F)
  
  domcall = paste0("domino --active_genes_files ",agfilename," --network_file ",network,
                   " --slices_file ",slices," --output_folder ",outputfolder, " --visualization False")
  #call domino
  system(domcall)
  #Load output
  domout = read_tsv(paste0(outputfolder,"/",disease,"/modules.out"),col_names=F)
  #Check if blank. If so just return now
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
  #Filter for training and test
  domminusprop=filter(output,Disease==disease) %>% filter(.,Gene %in% curgenes)
  #Replace "domino_" with "noprop_"
  domminusprop = domminusprop %>% mutate(Cluster = str_replace(Cluster,"domino_","noprop_"))
  return(domminusprop)
}
make_dom_minustest<-function(disease,output,dissplits){
  cursplit=dissplits %>% filter(Name==disease)
  curtrain=cursplit$Train %>% strsplit(", ") %>% unlist
  #curtest=cursplit$Test %>% strsplit(", ") %>% unlist
  #curgenes=c(curtrain,curtest)
  curgenes=curtrain
  #Filter for training and test
  domminusprop=filter(output,Disease==disease) %>% filter(.,Gene %in% curgenes)
  #Replace "domino_" with "noprop_"
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
  #writing output
  write_tsv(out,paste0("../data/clusters/","strict_normal_",typediff,"_domino.tsv"))
  #Only keep genes that are in original clusters
  dis=out$Disease %>% unique
  nopropout=mclapply(dis,make_dom_minusprop,out,dissplits,mc.cores=detectCores()) %>% bind_rows
  write_tsv(nopropout,paste0("./../data/clusters/","strict_normal_",typediff,"_noprop.tsv"))

  
}


