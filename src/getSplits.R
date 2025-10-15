#Genes in clusters are used to calculate posgenes. We then calculate negatives
#for each individual cluster and subset for neutrals

#output: trait | term (AllAssign,AllClus, clus num, neutral)| Gene | type

library(tidyverse)
source("moduleFunctions.R")
library(parallel)
library(reticulate)
use_condaenv("base")
args <- commandArgs(TRUE)
pyfile="calc_negatives.py"
source_python(pyfile)

#Getting arguments: Typediff
typediff=as.character(args[1])


call_negative_calc<-function(posgenes,gsctouse,net,traingenes,trainsub){
  
  #print("loaded")
  #Call function and put in arguments with it.
  
  #original function needed file_loc, net_type, hsc, and pos genes in network
  #1st is directory
  #2nd pretty sure is if STring, giant, etc
  #3rd is disgenet or not
  #4th is genes, so pass in vector posgenes
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
  
  #check if its mondo or not when making the decision
  if(gsctouse=="mondo" ){
    gsc="DisGeNet"
  }else{
    gsc=gsctouse
  }
  #I should have already subsetted for the network that I want to use
  pos_genes_in_net=posgenes
  out=get_negatives(file_loc,net_type,gsc,pos_genes_in_net,traingenes,trainsub) %>% as.vector
  
  
  return(out)
}

#Base splits on AllAssign for now. This gets called once
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

#To deal with having multiple cluster types
#Take genes in all - do the assign
#Then take genes in each individual - do the assignment
#This will get equal proportion and equality of AllAssign, and change negtrain depending on 
#contect

#I basically get the universe of negative genes for various cluster types. And then check
#"Okay are there any genes that are in one clustype (allassign/clus_allClus) but not in the other
get_assignment_multi<-function(disease,diffpos,alltypenegs,typecluss){
  #Loading negatives and getting genes common across all
  #and tibble recording uniques
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
  #Getting genes that are not commong
  untib = neggenetib %>% filter(!(Genes %in% commongenes))
  
  #Getting postest and postrain
  curpos=diffpos%>%filter(Name==disease)
  #Getting the AllAssigned postrain from the declared split
  postrain=curpos$Postrain %>% strsplit(", ") %>% unlist
  
  #Getting the postest. This currently gets postest from original declaration
  #and does not use clusters at all
  postest = curpos$Postest %>% strsplit(", ") %>% unlist

  pos_train_percent= length(postrain) / length(c(postrain,postest))
  pos_test_percent=1-pos_train_percent
  
  #Getting bins for common
  totalnegsplit=bin_gene_function(pos_train_percent,pos_test_percent,commongenes,length(commongenes),"All")
  
  #if(nrow(untib)>0){
  #  print(disease)
  #  return(totalnegsplit)
  #}
  
  #We now get splits for each clus type if there are genes that are in one but not the other
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

#negassigns come from AllAssign
#This function takes clus neg universe, and gets its traintest split based on
#how the genes were split in allassign
get_output<-function(disease,term,postrain,postest,negassigns,termneguni,aaneguni){
  

  #Join terms neggenes with negs tibble to get train/test assignments for that term
  
  termnegs=negassigns %>% filter(Neg_gene %in% termneguni)
    
  negtrain=filter(termnegs,Bin=="Negative_Train")$Neg_gene
  
  ###
  #Make negtest equal to allassign negtest
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
  #Filtering clusters for cur disease
  curclusters=clusters %>% filter(Disease==disease)
  
  #Filtering split for our disease
  curpos=diffpos%>%filter(Name==disease)
  #Getting the AllAssigned postrain from the declared split
  aapostrain=curpos$Postrain %>% strsplit(", ") %>% unlist %>% as.double
  #Getting the AllClus postrain from all cluster genes (genes made into clusters)
  acpostrain = curclusters$Gene
  #Getting the postest. This currently gets postest from original declaration
  #and does not use clusters at all
  postest = curpos$Postest %>% strsplit(", ") %>% unlist
  #Calculate negatives for the AllAssign
  aaneguni=call_negative_calc(aapostrain,"mondo","string",highmedgenes,trainsub=T)
  print(length(aaneguni))
  #Calculate negative for the AllClus
  acneguni=call_negative_calc(acpostrain,"mondo","string",highmedgenes,trainsub=T)
  print(length(acneguni))
  
  negout = tibble(Disease=disease,Typenet=typenet,NegType="AllAssign",Genes=aaneguni%>%paste0(collapse=", "))
  #negout = negout %>% bind_rows(tibble(Disease=disease,Typenet=typenet,NegType=paste0(typeclus,"_AllClus"),Genes=acneguni %>% paste0(collapse=", ")))
  neutralgenes=c()
  for(cluster in unique(curclusters$Cluster)){
    curcluspos=filter(curclusters,Cluster==cluster)$Gene
    if(length(curcluspos)<CLUSTER_SIZE){
      next
    }
    clusneguni=call_negative_calc(curcluspos,"mondo","string",highmedgenes,trainsub=T)
    #print(length(clusneguni))
    #Get genes in acneguni but not in clusneguni, as they are genes neutral in that cluster
    neutralgenes=c(neutralgenes,setdiff(acneguni,clusneguni)) %>% unique
  }
  #Have neutrals, so now can subset AllClusnegs
  subacneguni=acneguni[!(acneguni %in% neutralgenes)]
  
  #Add neutrals to output
  neutralgenes=paste0(neutralgenes,collapse=", ")
  
  negout = negout %>% bind_rows(tibble(Disease=disease,Typenet=typenet,NegType=paste0(typeclus,"_AllClus"),Genes=subacneguni %>% paste0(collapse=", ")))
  negout = negout %>% bind_rows(tibble(Disease=disease,Typenet=typenet,NegType=paste0(typeclus,"_Neutral"),Genes=neutralgenes))
}


get_split<-function(disease,diffpos,allnegs,typenets,typeclus,assignments,typediff){
  #print(disease)

  #Filtering split for our disease
  curpos=diffpos%>%filter(Name==disease)
  #print("yay")
  #Getting the AllAssigned postrain from the declared split
  aapostrain=curpos$Postrain %>% strsplit(", ") %>% unlist

  #Getting the postest. This currently gets postest from original declaration
  #and does not use clusters at all
  postest = curpos$Postest #%>% strsplit(", ") %>% unlist
  
  #Getting all seen negative genes for a disease
  disneg=filter(allnegs,Disease==disease,NegType!=paste0(typeclus,"_Neutral"))$Genes %>% strsplit(", ") %>% unlist %>% unique
  #print("yay2")
  #get negatives that are part of either all types, and just typeclus
  negassigns = assignments%>%filter(Disease==disease) %>% filter(TypeNeg=="All" | TypeNeg==typeclus) %>% select(-TypeNeg)
  print(any(duplicated(negassigns$Neg_gene)))
  #negassigns=get_assignment(aapostrain,postest,disneg)
  #print("yay3")

  #Just need to get these negatives one time as its the same no matter network
  aaneguni=filter(allnegs,Disease==disease & NegType=="AllAssign")$Genes[1] %>% strsplit(", ") %>% unlist %>% unique
  #Have assignments, can now set up split. Or load AllAssign for that trait calculated in "domino"
  if(typeclus=="domino"){
   outtib=get_output(disease,term="AllAssign",aapostrain,postest,negassigns,aaneguni) %>% mutate(Net="All")
  }else{
    outtib=tibble()
  }
  #Now each typenet gets unique things based on its AllClus. Double for loop with cluster assignments too
  for(typenet in typenets){

    #Read correct cluster file
    #clusters=read_tsv(paste0("../data/clusters/",typenet,"_",typeclus,".tsv"),col_types="ccc")
    clusters=read_tsv(paste0("../data/clusters/",typenet,"_",typediff,"_",typeclus,".tsv"),col_types="ccc")
    #Filtering clusters for cur disease
    curclusters=clusters %>% filter(Disease==disease)
    #If this typenet does not get this disease, skip
    if(nrow(curclusters)==0){
      print("Trait does not have big enough cluster in this typenet")
      next
    }
    
    #Getting the AllClus postrain from all cluster genes (genes made into clusters)
    acpostrain = curclusters$Gene
    
    #getting the allclus genes for this typenet
    acneguni=filter(allnegs,Disease==disease & Typenet==typenet & NegType==paste0(typeclus,"_AllClus"))$Genes[1] %>% strsplit(", ") %>% unlist %>% unique
    outtib = bind_rows(outtib,get_output(disease,term=paste0(typeclus,"_AllClus"),acpostrain,postest,negassigns,acneguni,aaneguni) %>% mutate(Net=typenet))
    
    clusnegtrain=filter(outtib,Term==paste0(typeclus,"_AllClus") & Net==typenet & Type=="Negtrain")$Gene
    clusnegtest=filter(outtib,Term==paste0(typeclus,"_AllClus") & Net==typenet &Type=="Negtest")$Gene
    
    #Now do cluster assignments
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
    #Add the neutrals
    outtib = outtib %>% bind_rows(tibble(Trait=disease,Term=paste0(typeclus,"_Neutral"),Gene=filter(allnegs,Disease==disease&Typenet==typenet&NegType==paste0(typeclus,"_Neutral"))$Genes,Type=paste0(typeclus,"_Neutral"),Net=typenet))
  }
  
    return(outtib)
}

#positive genes of traits. Must load separately to get test genes since not in cluster output
#This came from Chris stuff
#diffpos=read_tsv(paste0("../data/splits/diff_splits_string.tsv"))
diffpos=read_tsv(paste0("../data/splits/",typediff,"_splits_string.tsv"))
colnames(diffpos)=c("ID","Name","Postrain","Postest")

#Bins needed for negative calculations
#bins=read_tsv(paste0("../data/study_bins/string_diff_bins.tsv"),col_types="cccc") %>% select(-"...1")
bins=read_tsv(paste0("../data/study_bins/string_",typediff,"_bins.tsv"),col_types="cccc") %>% select(-"...1")
#Edit bins so I can join it by making "Bin" "Gene" columns
bins=bins %>% gather(Bin,Gene,1:ncol(bins))
highmedgenes = filter(bins,Bin!="low_genes")$Gene

#typenets=c("strict_normal","strict_wto","med_normal","med_wto")
typenets=c("strict_normal")

#Loop here for all types of runs. So we will have "domino","louvain","kmeans"
#ASSUMPTION: Only calculate AllAssign negs in domino
typecluss=c("domino","noprop")
#typecluss=c("domino","noprop","notest")
#typecluss=c("domino")

outtib=tibble()
negfiles=c()
for(typeclus in typecluss){
  
  allnegs=tibble()
  
  for(typenet in typenets){
    #clusters=read_tsv(paste0("../data/clusters/",typenet,"_",typeclus,".tsv"),col_types="ccc")
    clusters=read_tsv(paste0("../data/clusters/",typenet,"_",typediff,"_",typeclus,".tsv"),col_types="ccc")
    clusters = clusters %>% arrange(Disease)
    diseases=clusters$Disease%>% unique
    #diseases=diseases[1:20]
    #diseases="Diabetes_Mellitus__Experimental"
    #diseases=diseases[1:20]
    out =mclapply(diseases,get_negs,diffpos,highmedgenes,typenet,typeclus,clusters,mc.cores=detectCores()) %>% bind_rows
    print(out)
    out  = out %>% bind_rows
    out$Typeclus=typeclus
    allnegs = allnegs %>% bind_rows(out)
  }
  negfiles=c(negfiles,paste0("./../data/pos_neg_combinations/",typediff,"_",typeclus,"negs_tempout.tsv"))
  write_tsv(allnegs,paste0("./../data/pos_neg_combinations/",typediff,"_",typeclus,"negs_tempout.tsv"))
  #No loop for splits as need all typenets to do it
  allnegs=read_tsv(paste0("./../data/pos_neg_combinations/",typediff,"_",typeclus,"negs_tempout.tsv"))
}
#Getting negative assignment for all diseases
alltypenegs=mclapply(negfiles,read_tsv,mc.cores=detectCores(length(negfiles))) %>% bind_rows
diseases=alltypenegs$Disease %>% unique
#diseases=diseases[1:20]
assignments=mclapply(diseases,get_assignment_multi,diffpos,alltypenegs,typecluss,mc.cores=detectCores()) %>% bind_rows


#Getting splits for all diseases
for(typeclus in typecluss){ 
  allnegs=read_tsv(paste0("./../data/pos_neg_combinations/",typediff,"_",typeclus,"negs_tempout.tsv"))
  newdiseases=allnegs$Disease %>% unique
  #Pass outtib for access to AllAsisgn negatives assignment
  out=mclapply(newdiseases,get_split,diffpos,allnegs,typenets,typeclus,assignments,typediff,mc.cores=detectCores()) %>% bind_rows
  #lapply(newdiseases,get_split,diffpos,allnegs,typenets,typeclus,assignments,typediff)
  outtib=bind_rows(outtib,out)
}

#Split "out" up using the Net column and write multple files
for(typenet in typenets){
  #Filter for typenet and for "All" 
  curout=outtib %>% filter(Net == typenet | Net=="All") %>% select(-Net)
  print(nrow(curout))
  #print(curout$Term %>% table)
  #write
  write_tsv(curout,paste0("./../data/pos_neg_combinations/",typenet,"_",typediff,"_neutral_diff_string_allclustype.tsv"))
  
}
