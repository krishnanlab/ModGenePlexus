library(tidyverse)
source("../moduleFunctions.R")
library(parallel)
library(reticulate)
use_condaenv("base")
args <- commandArgs(TRUE)
pyfile="/mnt/research/compbio/krishnanlab/projects/module_expansion/src/calc_negatives.py"
source_python(pyfile)

call_negative_calc<-function(posgenes,gsctouse,net,traingenes,trainsub){
  file_loc="/mnt/home/mckimale/play/plexpy/data_2022_05_17/"
  if(net=="string"){
    net_type="STRING"
  }else if(net=="string_exp"){
    net_type="STRING-EXP"
  }else if(net=="consensus"){
    net_type="ConsensusPathDB"
  }else if(net=="giant"){
    net_type="GIANT-TN"
  }

  if(gsctouse=="mondo" ){
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
  # the indexing at the end is for an edge case where floor+ceiling > length(neguni)
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
  print(length(aaneguni))
  acneguni=call_negative_calc(acpostrain,"mondo","string",highmedgenes,trainsub=T)
  print(length(acneguni))

  negout = tibble(Disease=disease,Typenet=typenet,NegType="AllAssign",Genes=aaneguni%>%paste0(collapse=", "))
  neutralgenes=c()
  for(cluster in unique(curclusters$Cluster)){
    curcluspos=filter(curclusters,Cluster==cluster)$Gene
    if(length(curcluspos)<CLUSTER_SIZE){
      next
    }
    clusneguni=call_negative_calc(curcluspos,"mondo","string",highmedgenes,trainsub=T)
    # genes in acneguni but not clusneguni are neutral for that cluster
    neutralgenes=c(neutralgenes,setdiff(acneguni,clusneguni)) %>% unique
  }
  subacneguni=acneguni[!(acneguni %in% neutralgenes)]

  neutralgenes=paste0(neutralgenes,collapse=", ")

  negout = negout %>% bind_rows(tibble(Disease=disease,Typenet=typenet,NegType=paste0(typeclus,"_AllClus"),Genes=subacneguni %>% paste0(collapse=", ")))
  negout = negout %>% bind_rows(tibble(Disease=disease,Typenet=typenet,NegType=paste0(typeclus,"_Neutral"),Genes=neutralgenes))
}


get_split<-function(disease,diffpos,allnegs,typenets,typeclus,assignments,typediff){
  if(disease=="CD4__Treg_39_73___NA_2015"){
    return(tibble())
  }

  curpos=diffpos%>%filter(Name==disease)
  aapostrain=curpos$Postrain %>% strsplit(", ") %>% unlist

  postest = curpos$Postest

  disneg=filter(allnegs,Disease==disease,NegType!=paste0(typeclus,"_Neutral"))$Genes %>% strsplit(", ") %>% unlist %>% unique
  negassigns = assignments%>%filter(Disease==disease) %>% filter(TypeNeg=="All" | TypeNeg==typeclus) %>% select(-TypeNeg)
  print(any(duplicated(negassigns$Neg_gene)))

  aaneguni=filter(allnegs,Disease==disease & NegType=="AllAssign")$Genes[1] %>% strsplit(", ") %>% unlist %>% unique
  outtib=get_output(disease,term=paste0(typediff,"_AllAssign"),aapostrain,postest,negassigns,aaneguni) %>% mutate(Net="All")
  for(typenet in typenets){
    clusters=read_tsv(paste0("../../results/clusters/",typenet,"_",typediff,"_mstrict_domino.tsv"),col_types="ccc")
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
      outtib = outtib %>% bind_rows(tibble(Trait=disease,Term=paste0(typediff,cl),
                                           Gene=c(curpostrain,postest,clusnegtrain,clusnegtest),
                                           Type=c("Postrain","Postest","Negtrain","Negtest"),
                                           Net=typenet))
    }
    outtib = outtib %>% bind_rows(tibble(Trait=disease,Term=paste0(typeclus,"_Neutral"),Gene=filter(allnegs,Disease==disease&Typenet==typenet&NegType==paste0(typeclus,"_Neutral"))$Genes,Type=paste0(typeclus,"_Neutral"),Net=typenet))
  }

    return(outtib)
}

make_diffpos<-function(magmatype){
  if(magmatype=="magma_1e2"){
    diffpos=read_tsv("../../data/splits/magma_1e2_mstrict_splits_string.tsv")
  }else if(magmatype=="magma_1e5"){
    diffpos=read_tsv("../../data/splits/magma_1e5_mstrict_splits_string.tsv")
  }else if(magmatype=="magma_5e2"){
    diffpos=read_tsv("../../data/splits/magma_5e2_mstrict_splits_string.tsv")
  }else if(magmatype=="magma_1e1"){
    diffpos=read_tsv("../../data/splits/magma_1e1_mstrict_splits_string.tsv")
  }
  colnames(diffpos)=c("ID","Name","Postrain","Postest")
  return(diffpos)
}

typenets=c("strict_normal")

typecluss=c("magma_1e2_mstrict_domino","magma_1e5_mstrict_domino","magma_5e2_mstrict_domino","magma_1e1_mstrict_domino")

outtib=tibble()
negfiles=c()
for(typeclus in typecluss){
  if(typeclus=="magma_1e2_mstrict_domino"){
    typediff="magma_1e2"
    diffpos=make_diffpos("magma_1e2")
  }else if(typeclus=="magma_1e5_mstrict_domino"){
    typediff="magma_1e5"
    diffpos=make_diffpos("magma_1e5")
  }else if(typeclus=="magma_5e2_mstrict_domino"){
    typediff="magma_5e2"
    diffpos=make_diffpos("magma_5e2")
  }else if(typeclus=="magma_1e1_mstrict_domino"){
    typediff="magma_1e1"
    diffpos=make_diffpos("magma_1e1")
  }

  bins=read_tsv(paste0("../../data/study_bins/string_",typediff,"_bins.tsv"),col_types="cccc") %>% select(-"...1")
  bins=bins %>% gather(Bin,Gene,1:ncol(bins))
  highmedgenes = filter(bins,Bin!="low_genes")$Gene

  allnegs=tibble()

  for(typenet in typenets){
    clusters=read_tsv(paste0("../../results/clusters/",typenet,"_",typeclus,".tsv"),col_types="ccc") %>% filter(Disease %in% diffpos$Name)
    clusters = clusters %>% arrange(Disease)
    diseases=clusters$Disease%>% unique

    out =mclapply(diseases,get_negs,diffpos,highmedgenes,typenet,typeclus,clusters,mc.cores=detectCores()) %>% bind_rows
    print(out)
    out  = out %>% bind_rows
    out$Typeclus=typeclus
    allnegs = allnegs %>% bind_rows(out)
  }
  negfiles=c(negfiles,paste0("./../../data/pos_neg_combinations/",typeclus,"negs_tempout_mstrict.tsv"))
  write_tsv(allnegs,paste0("./../../data/pos_neg_combinations/",typeclus,"negs_tempout_mstrict.tsv"))
  allnegs=read_tsv(paste0("./../../data/pos_neg_combinations/",typeclus,"negs_tempout_mstrict.tsv"))
}
alltypenegs=mclapply(negfiles,read_tsv,mc.cores=detectCores(length(negfiles))) %>% bind_rows
diseases=alltypenegs$Disease %>% unique

diffpos=make_diffpos("magma_1e5")
assignments=mclapply(diseases,get_assignment_multi,diffpos,alltypenegs,typecluss,mc.cores=detectCores()) %>% bind_rows


for(typeclus in typecluss){
  if(typeclus=="magma_1e2_mstrict_domino"){
    typediff="magma_1e2"
    diffpos=make_diffpos("magma_1e2")
  }else if(typeclus=="magma_1e5_mstrict_domino"){
    typediff="magma_1e5"
    diffpos=make_diffpos("magma_1e5")
  }else if(typeclus=="magma_5e2_mstrict_domino"){
    typediff="magma_5e2"
    diffpos=make_diffpos("magma_5e2")
  }else if(typeclus=="magma_1e1_mstrict_domino"){
    typediff="magma_1e1"
    diffpos=make_diffpos("magma_1e1")
  }

  allnegs=read_tsv(paste0("./../../data/pos_neg_combinations/",typeclus,"negs_tempout_mstrict.tsv"))
  newdiseases=allnegs$Disease %>% unique
  out=mclapply(newdiseases,get_split,diffpos,allnegs,typenets,typeclus,assignments,typediff,mc.cores=detectCores()) %>% bind_rows
  outtib=bind_rows(outtib,out)
}

typediff="mstrict"
for(typenet in typenets){
  curout=outtib %>% filter(Net == typenet | Net=="All") %>% select(-Net)
  print(nrow(curout))
  write_tsv(curout,paste0("./../../data/pos_neg_combinations/",typenet,"_",typediff,"_neutral_diff_string_allclustype.tsv"))

}

diseases=outtib$Trait %>% unique
bad=tibble()
for(trait in diseases){
  cur=outtib %>% filter(Trait==trait)
  terms=cur$Term %>% unique
  for(term in terms){
    curterm = cur %>% filter(Term==term)
    curnegtrain=filter(curterm,Type=="Negtrain")$Gene %>% strsplit(", ") %>% unlist
    # only care about AllAssigns for negtest
    aaterms=c("magma_1e2_AllAssign","magma_1e5_AllAssign","magma_5e2_AllAssign","magma_1e1_AllAssign")
    for(otherterm in aaterms){
      oterm=cur %>% filter(Term==otherterm)
      onegtest=filter(oterm,Type=="Negtest")$Gene %>% strsplit(", ") %>% unlist
      m=intersect(onegtest,curnegtrain)
      if(length(m)>0 & term!=otherterm){
        bad=bind_rows(bad,tibble(Term=term,Otherterm=otherterm))
      }
    }


  }
}
result=c()
for(disease in unique(assignments$Disease)){
  cur=filter(assignments,Disease==disease)
  for(tneg in unique(assignments$TypeNeg)){
    newcur=filter(cur,TypeNeg==tneg)
  has_duplicates <- any(duplicated(newcur$Neg_gene))

  result=c(result,has_duplicates)
  }
  print(has_duplicates)

}
