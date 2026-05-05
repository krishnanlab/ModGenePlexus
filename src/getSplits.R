library(tidyverse)
source("moduleFunctions.R")
library(parallel)
library(reticulate)
use_condaenv("base")
args <- commandArgs(TRUE)
pyfile="calc_negatives.py"
source_python(pyfile)

typediff=as.character(args[1])


call_negative_calc<-function(posgenes,gsctouse,net,traingenes,trainsub){
  file_loc="../plexpy/data_2022_05_17/"
  if(net=="string"){
    net_type="STRING"
  }else if(net=="string_exp"){
    net_type="STRING-EXP"
  }else if(net=="consensus"){
    net_type="ConsensusPathDB"
  }else if(net=="giant"){
    net_type="GIANT-TN"
  }
  
  if(gsctouse=="mondo"){
    gsc="DisGeNet"
  }else{
    gsc=gsctouse
  }
  pos_genes_in_net=posgenes
  out=get_negatives(file_loc,net_type,gsc,pos_genes_in_net,traingenes,trainsub) %>% as.vector
  return(out)
}

get_assignment<-function(postrain,postest,neguni){
  postest = postest %>% strsplit(", ") %>% unlist
  pos_train_percent= length(postrain) / length(c(postrain,postest))
  pos_test_percent=1-pos_train_percent
  #the indexing at the end is for an edge case
  negbins=(c(rep("Negative_Train", floor((pos_train_percent)*length(neguni))),
             rep("Negative_Test",ceiling((pos_test_percent)*length(neguni)))) %>% sample)[1:length(neguni)]
  negs=tibble(Neg_gene=neguni,Bin=negbins)
  negtrain=filter(negs,Bin=="Negative_Train")$Neg_gene
  negtest=filter(negs,Bin=="Negative_Test")$Neg_gene
  
  return(negs)
}

bin_gene_function<-function(pos_train_percent,pos_test_percent,neguni,reqsize,type){
  negbins=(c(rep("Negative_Train", floor((pos_train_percent)*reqsize)),
             rep("Negative_Test",ceiling((pos_test_percent)*reqsize))) %>% sample)[1:reqsize]
  negs=tibble(Neg_gene=neguni,Bin=negbins,TypeNeg=type)
  
  return(negs)
}

get_assignment_multi<-function(disease,diffpos,alltypenegs,typecluss){
  curtypenegs=alltypenegs %>% filter(Disease==disease) %>% filter(!grepl("_Neutral",NegType))
  
  neggenetib=tibble()
  commongenes=c()
  #Filtering on typeclus will get AllAssign and typeclus_AllClus
  for(typeclus in typecluss){
    curnegs=curtypenegs %>% filter(Typeclus==typeclus)
    curneggenes=lapply(curnegs$Genes, strsplit,", ") %>% unlist %>% unique

    neggenetib=bind_rows(neggenetib,tibble(Typeclus=typeclus,Genes=curneggenes))
    
  }
  
  allgenes=neggenetib$Genes %>% unique
  for(typeclus in typecluss){
    curgenes=filter(neggenetib,Typeclus==typeclus)$Genes
    commongenes=intersect(allgenes,curgenes)
  }
  untib = neggenetib %>% filter(!(Genes %in% commongenes))

  curpos=diffpos%>%filter(Name==disease)
  postrain=curpos$Postrain %>% strsplit(", ") %>% unlist
  postest = curpos$Postest %>% strsplit(", ") %>% unlist

  pos_train_percent= length(postrain) / length(c(postrain,postest))
  pos_test_percent=1-pos_train_percent
  
  totalnegsplit=bin_gene_function(pos_train_percent,pos_test_percent,commongenes,length(commongenes),"All")

  if(nrow(untib)>0){
    print(disease)
  for(typeclus in typecluss){
    curtib=untib %>% filter(Typeclus==typeclus)
    genes=curtib$Genes
    totalnegsplit=bind_rows(totalnegsplit,bin_gene_function(pos_train_percent,pos_test_percent,genes,length(genes),typeclus))
  }
  }
  totalnegsplit$Disease=disease
  print(totalnegsplit)
  return(totalnegsplit)
}

get_output<-function(disease,term,postrain,postest,negassigns,termneguni,aaneguni){
  termnegs=negassigns %>% filter(Neg_gene %in% termneguni)
  negtrain=filter(termnegs,Bin=="Negative_Train")$Neg_gene
  negtest = filter(negassigns,Bin=="Negative_Test")$Neg_gene
  
  postrain=paste0(postrain,collapse=", ")
  postest=paste0(postest,collapse=", ")
  negtrain=paste0(negtrain,collapse=", ")
  negtest=paste0(negtest,collapse=", ")

  out=tibble(Trait=disease,Term=term,Gene=c(postrain,postest,negtrain,negtest),
             Type=c("Postrain","Postest","Negtrain","Negtest"))
  
  return(out)
}

get_negs<-function(disease,diffpos,highmedgenes,typenet,typeclus,clusters){

  print(disease)
  curclusters=clusters %>% filter(Disease==disease)
  curpos=diffpos%>%filter(Name==disease)
  aapostrain=curpos$Postrain %>% strsplit(", ") %>% unlist %>% as.double
  acpostrain = curclusters$Gene
  postest = curpos$Postest %>% strsplit(", ") %>% unlist
  aaneguni=call_negative_calc(aapostrain,"mondo","string",highmedgenes,trainsub=T)
  acneguni=call_negative_calc(acpostrain,"mondo","string",highmedgenes,trainsub=T)

  negout = tibble(Disease=disease,Typenet=typenet,NegType="AllAssign",Genes=aaneguni%>%paste0(collapse=", "))
  neutralgenes=c()
  for(cluster in unique(curclusters$Cluster)){
    curcluspos=filter(curclusters,Cluster==cluster)$Gene
    if(length(curcluspos)<CLUSTER_SIZE){
      next
    }
    clusneguni=call_negative_calc(curcluspos,"mondo","string",highmedgenes,trainsub=T)
    neutralgenes=c(neutralgenes,setdiff(acneguni,clusneguni)) %>% unique
  }
  subacneguni=acneguni[!(acneguni %in% neutralgenes)]
  neutralgenes=paste0(neutralgenes,collapse=", ")

  negout = negout %>% bind_rows(tibble(Disease=disease,Typenet=typenet,NegType=paste0(typeclus,"_AllClus"),Genes=subacneguni %>% paste0(collapse=", ")))
  negout = negout %>% bind_rows(tibble(Disease=disease,Typenet=typenet,NegType=paste0(typeclus,"_Neutral"),Genes=neutralgenes))
}


get_split<-function(disease,diffpos,allnegs,typenets,typeclus,assignments,typediff){
  curpos=diffpos%>%filter(Name==disease)
  aapostrain=curpos$Postrain %>% strsplit(", ") %>% unlist
  postest = curpos$Postest
  disneg=filter(allnegs,Disease==disease,NegType!=paste0(typeclus,"_Neutral"))$Genes %>% strsplit(", ") %>% unlist %>% unique
  negassigns = assignments%>%filter(Disease==disease) %>% filter(TypeNeg=="All" | TypeNeg==typeclus) %>% select(-TypeNeg)
  print(any(duplicated(negassigns$Neg_gene)))
  aaneguni=filter(allnegs,Disease==disease & NegType=="AllAssign")$Genes[1] %>% strsplit(", ") %>% unlist %>% unique
  if(typeclus=="domino"){
    outtib=get_output(disease,term="AllAssign",aapostrain,postest,negassigns,aaneguni) %>% mutate(Net="All")
  }else{
    outtib=tibble()
  }
  for(typenet in typenets){
    clusters=read_tsv(paste0("../data/clusters/",typenet,"_",typediff,"_",typeclus,".tsv"),col_types="ccc")
    curclusters=clusters %>% filter(Disease==disease)
    if(nrow(curclusters)==0){
      print("Trait does not have big enough cluster in this typenet")
      next
    }
    acpostrain = curclusters$Gene
    acneguni=filter(allnegs,Disease==disease & Typenet==typenet & NegType==paste0(typeclus,"_AllClus"))$Genes[1] %>% strsplit(", ") %>% unlist %>% unique
    outtib = bind_rows(outtib,get_output(disease,term=paste0(typeclus,"_AllClus"),acpostrain,postest,negassigns,acneguni,aaneguni) %>% mutate(Net=typenet))
    clusnegtrain=filter(outtib,Term==paste0(typeclus,"_AllClus") & Net==typenet & Type=="Negtrain")$Gene
    clusnegtest=filter(outtib,Term==paste0(typeclus,"_AllClus") & Net==typenet &Type=="Negtest")$Gene
    for(cluster in unique(curclusters$Cluster)){
      curcluspos=filter(curclusters,Cluster==cluster)$Gene
      if(length(curcluspos)<CLUSTER_SIZE){
        next
      }
      curpostrain=paste0(curcluspos,collapse=", ")
      cl=as.character(cluster)
      outtib = outtib %>% bind_rows(tibble(Trait=disease,Term=cl,
                                           Gene=c(curpostrain,postest,clusnegtrain,clusnegtest),
                                           Type=c("Postrain","Postest","Negtrain","Negtest"),
                                           Net=typenet))
    }
    outtib = outtib %>% bind_rows(tibble(Trait=disease,Term=paste0(typeclus,"_Neutral"),Gene=filter(allnegs,Disease==disease&Typenet==typenet&NegType==paste0(typeclus,"_Neutral"))$Genes,Type=paste0(typeclus,"_Neutral"),Net=typenet))
  }
  return(outtib)
}

diffpos=read_tsv(paste0("../data/splits/",typediff,"_splits_string.tsv"))
colnames(diffpos)=c("ID","Name","Postrain","Postest")

bins=read_tsv(paste0("../data/study_bins/string_",typediff,"_bins.tsv"),col_types="cccc") %>% select(-"...1")
bins=bins %>% gather(Bin,Gene,1:ncol(bins))
highmedgenes = filter(bins,Bin!="low_genes")$Gene

typenets=c("strict_normal")
typecluss=c("domino","noprop")

outtib=tibble()
negfiles=c()
for(typeclus in typecluss){
  allnegs=tibble()
  for(typenet in typenets){
    clusters=read_tsv(paste0("../data/clusters/",typenet,"_",typediff,"_",typeclus,".tsv"),col_types="ccc")
    clusters = clusters %>% arrange(Disease)
    diseases=clusters$Disease%>% unique
    out=mclapply(diseases,get_negs,diffpos,highmedgenes,typenet,typeclus,clusters,mc.cores=detectCores()) %>% bind_rows
    out$Typeclus=typeclus
    allnegs = allnegs %>% bind_rows(out)
  }
  negfiles=c(negfiles,paste0("./../data/pos_neg_combinations/",typediff,"_",typeclus,"negs_tempout.tsv"))
  write_tsv(allnegs,paste0("./../data/pos_neg_combinations/",typediff,"_",typeclus,"negs_tempout.tsv"))
  allnegs=read_tsv(paste0("./../data/pos_neg_combinations/",typediff,"_",typeclus,"negs_tempout.tsv"))
}
alltypenegs=mclapply(negfiles,read_tsv,mc.cores=detectCores(length(negfiles))) %>% bind_rows
diseases=alltypenegs$Disease %>% unique
assignments=mclapply(diseases,get_assignment_multi,diffpos,alltypenegs,typecluss,mc.cores=detectCores()) %>% bind_rows

for(typeclus in typecluss){
  allnegs=read_tsv(paste0("./../data/pos_neg_combinations/",typediff,"_",typeclus,"negs_tempout.tsv"))
  newdiseases=allnegs$Disease %>% unique
  out=mclapply(newdiseases,get_split,diffpos,allnegs,typenets,typeclus,assignments,typediff,mc.cores=detectCores()) %>% bind_rows
  outtib=bind_rows(outtib,out)
}

for(typenet in typenets){
  curout=outtib %>% filter(Net == typenet | Net=="All") %>% select(-Net)
  print(nrow(curout))
  write_tsv(curout,paste0("./../data/pos_neg_combinations/",typenet,"_",typediff,"_neutral_diff_string_allclustype.tsv"))
}
