
library(tidyverse)
library(parallel)
library(EGAD)
library(gdata)

maxeq<-function(vec){
  newvec=1-vec
  p=prod(newvec)
  p=1-p
  return(p)
}
#p=1- pi_i (1-p_i)
getmaxres<-function(res,type="Prob"){
  if(type=="Prob"){
    maxres=res %>% select(-Rank) %>% group_by(Entrez) %>% summarize_at(vars(Probability),list(MaxProb=maxeq))%>% arrange(MaxProb %>% desc)
  }else if(type=="normrank"){
    maxres=res %>% select(-Probability) %>% group_by(Entrez) %>% summarize_at(vars(PRank),list(MaxProb=maxeq))%>% arrange(MaxProb %>% desc)
  }else if(type=="maxmin"){
    #Do this twice, then choose the highest of MaxRank and MinRank. If MinRank is higher, the prob for that gene is 1-minrank. If MaxRank higher keep it
    maxres=res %>% select(-Probability) %>% group_by(Entrez) %>% summarize_at(vars(PRank),list(MaxRank=maxeq))%>% arrange(MaxRank %>% desc)
    minres=res %>% select(-Probability) %>% group_by(Entrez) %>% summarize_at(vars(InP),list(MinRank=maxeq))%>% arrange(MinRank %>% desc)
    maxres=maxres %>% left_join(minres,by="Entrez")
    maxres =maxres %>% mutate(MaxProb = case_when(MinRank>=MaxRank ~ 1-MinRank,
                                                  MaxRank>(MinRank) ~ MaxRank , TRUE ~ (1-MinRank))) %>% select(-MaxRank,-MinRank)
  }
  colnames(maxres)[2]="Probability"
  return(maxres)
}

#What gets passed in here is the maxres and then all results
getavgres<-function(maxres,clusres){
  badprobs = maxres %>% filter(Probability < .20)
  #Take just the bad prob genes, and get a score for each gene where its the average rather than max
  avggenes = clusres %>% select(-Rank) %>% group_by(Entrez)  %>% summarize_at(vars(Probability),list(MaxProb=min))%>% arrange(MaxProb %>% desc)
  #Only keep genes in avggenes that are in bad probs
  avgofbad = avggenes %>% filter(Entrez %in% badprobs$Entrez)
  colnames(avgofbad)[2]="Probability"
  #Binding and arranging in correct order
  avgres = maxres %>% filter(Probability >=.20) %>% bind_rows(avgofbad) %>% arrange(Probability %>% desc)
  
  
  return(avgres)
}


getlog2AUPRC<-function(plexres,postest,negtest,neutrals){
  plexres=plexres %>% mutate(Label=case_when(Entrez %in% postest~1,
                                             Entrez %in% negtest~ 0))
  #plexres=plexres %>% mutate(Label=case_when(Entrez %in% testgenes~1))
  plexres$Label= plexres$Label %>% replace_na(-1)
  #Get prior, the test pos size / (postest + negtest)
  prior=length(postest) / (length(postest)+length(negtest))
  #prior=1
  print(prior)
  #auprc is binary, so get rid of neutrals
  plexres=plexres %>% filter(Label!=-1)
  #Divide the probability over the prior
  score=auprc(plexres$Probability,plexres$Label)
  score=log2(score/prior)
  #Subset testgenes on neutrals and then calculate another score, then put scores
  #in collapsed string
  plexres = plexres %>% filter(!(Entrez %in% neutrals))
  newscore=auprc(plexres$Probability,plexres$Label)
  newscore=log2(newscore/prior)
  score=paste0(c(score,newscore),collapse=", ")
  return(score)
}

#Calculate tau using
#tau = sum(1-xi)/(n-1) where xi = xi/max(of all xis)
tau<-function(clusres){
  clusres <- clusres %>%
    group_by(Entrez) %>%
    mutate(normalized_prob = Probability / max(Probability)) %>%
    ungroup()
  
  # Step 2: Calculate τ for each Entrez
  clusres_tau <- clusres %>%
    group_by(Entrez) %>%
    summarise(
      tau = sum(1 - normalized_prob) / (n() - 1)
    )
  
  colnames(clusres_tau)=c("Entrez","Probability")
  return(clusres_tau)
}

#tau = sum(1-xi)/(n-1) where xi = xi/max(of all xis)
taurank<-function(clusres){
  clusres <- clusres %>%
    group_by(Entrez) %>%
    mutate(normalized_prob = Rank / max(Rank)) %>%
    ungroup()
  
  # Step 2: Calculate τ for each Entrez
  clusres_tau <- clusres %>%
    group_by(Entrez) %>%
    summarise(
      tau = sum(1 - normalized_prob) / (n() - 1)
    )
  #clusres_tau$tau = clusres_tau$tau * (-1)
  
  colnames(clusres_tau)=c("Entrez","Probability")
  return(clusres_tau)
}

tauminprob<-function(clusres,cursplit){
  

  
  clusres <- clusres %>% mutate(Normprob=Probability/Rank) %>%#%>% mutate(Normprob=log2(Probability/Prior)) %>%
    group_by(Entrez) %>%
    mutate(normalized_prob = Normprob / min(Normprob)) %>%
    ungroup()
  
  # Step 2: Calculate τ for each Entrez
  clusres_tau <- clusres %>%
    group_by(Entrez) %>%
    summarise(
      tau = sum(1 - normalized_prob) / (n() - 1)
    )
  
  colnames(clusres_tau)=c("Entrez","Probability")
  return(clusres_tau)
}

#Dealing with situation swhere if you have all of them be high, you get lower taus
taurank_normalized<-function(clusres,valtotry=1,ranktotry=10){
  clusres <- clusres %>%
    group_by(Entrez) %>%
    mutate(normalized_prob = Rank / max(Rank)) %>%
    ungroup()
  
  # Step 2: Calculate τ for each Entrez
  clusres_tau <- clusres %>%
    group_by(Entrez) %>%
    summarise(
      tau = (sum(1 - normalized_prob) / (n() - 1)) +
        ifelse(min(Rank) <= ranktotry, (valtotry/min(Rank))/(n()-1), 0) # Conditionally apply boosting
    )
  #clusres_tau$tau = clusres_tau$tau * (-1)
  
  
  
  colnames(clusres_tau)=c("Entrez","Probability")
  return(clusres_tau)
}

quantilerank<-function(clusres){
  clusres_quan <- clusres %>%
    group_by(Entrez) %>%
    summarize(Quantil = quantile(Rank,.05)) %>% 
    ungroup()
  colnames(clusres_quan)=c("Entrez","Probability")
  return(clusres_quan)
  
}

quantilemet<-function(clusres){
  clusres_quan <- clusres %>%
    group_by(Entrez) %>%
    summarize(Quantil = quantile(Probability,.8)) %>% 
    ungroup()
  colnames(clusres_quan)=c("Entrez","Probability")
  return(clusres_quan)
}

taurank_min<-function(clusres){
  clusres <- clusres %>%
    group_by(Entrez) %>%
    mutate(normalized_prob = (1/Rank) / min(Rank)) %>%
    ungroup()
  
  # Step 2: Calculate τ for each Entrez
  clusres_tau <- clusres %>%
    group_by(Entrez) %>%
    summarise(
      tau = sum(1 - normalized_prob) / (n() - 1)
    )
  #clusres_tau$tau = clusres_tau$tau * (-1)
  
  colnames(clusres_tau)=c("Entrez","Probability")
  return(clusres_tau)
  
  
  
}
taurank_trans<-function(clusres){
  N=length(unique(clusres$Entrez))
  
  clusres <- clusres %>%
    group_by(Entrez) %>%
    mutate(normalized_prob =((N-Rank+1)/N)) %>%
    ungroup()
  
  # Step 2: Calculate τ for each Entrez
  clusres_tau <- clusres %>%
    group_by(Entrez) %>%
    summarise(
      tau = sum(1 - normalized_prob) / (n() - 1)
    )
  #clusres_tau$tau = clusres_tau$tau * (-1)
  
  colnames(clusres_tau)=c("Entrez","Probability")
  return(clusres_tau)
  
  
  
}

meanrank<-function(clusres){
  
  clusres <- clusres %>%
    group_by(Entrez) %>%
    mutate(MeanRank=mean(Rank)) %>%
    ungroup() %>% select(Entrez,Rank,Cluster,MeanRank)
  
  clusres_rank=clusres %>% select(Entrez,MeanRank) %>% distinct
  colnames(clusres_rank)=c("Entrez","Probability")
  return(clusres_rank)
  
  
  
}


#Taus
taucalc<-function(file,typeclus){
  plexout=read_tsv(file,col_types="cccdccdccd") %>% filter(str_starts(Cluster,typeclus) | Cluster=="AllAssign")
  disease=plexout$Disease[1]
  cursplit=monsplits %>% filter(Trait==disease)

  outtib=tibble()
  
  #Auprc of AllAssign and AllClus
  wholes=c("AllAssign",paste0(typeclus,"_AllClus"))
  for(whole in wholes){
    #This is equal to AllAssign because I want just the postest/negtest. Not of allclus, but of AllAssign for evaluation
    curwhole=cursplit %>% filter(Term=="AllAssign"| Term==paste0(typeclus,"_Neutral"))
    cpostest=filter(curwhole,Type=="Postest")$Gene %>% strsplit(", ") %>% unlist
    cnegtest=filter(curwhole,Type=="Negtest")$Gene %>% strsplit(", ") %>% unlist
    cneutrals=filter(curwhole,Term==paste0(typeclus,"_Neutral"))$Gene %>% strsplit(", ") %>% unlist
    
    curres=plexout %>% filter(Cluster==whole)
    cauprcs=getlog2AUPRC(curres,cpostest,cnegtest,cneutrals)
    outtib=outtib %>% bind_rows(tibble(Trait=disease,Term=whole,AUPRCS=cauprcs))
    
  }
  #AUPRC of various max methods: Eval on AllClus like normal
  
  clusres=plexout %>% filter(!(Cluster %in% c("AllAssign",paste0(typeclus,"_AllClus"))))
  
  #Filter for where it has "typeclus" in Cluster
  clusres = clusres %>% filter(str_detect(Cluster,typeclus))
  
  maxres=getmaxres(clusres,"Prob")
  
  
  
  compostest=filter(cursplit,Term=="AllAssign" & Type=="Postest")$Gene %>% strsplit(", ") %>% unlist
  comnegtest=filter(cursplit,Term=="AllAssign" & Type=="Negtest")$Gene %>% strsplit(", ") %>% unlist
  maxauprcs=getlog2AUPRC(maxres,compostest,comnegtest,cneutrals)
  outtib=outtib %>% bind_rows(tibble(Trait=disease,Term="MaxProb",AUPRCS=maxauprcs))
  
  #Max res is default based on all maxes. Add version where replace those where max is too low
  #and take the average instead
  
  
  maxavgres=getavgres(maxres,clusres)
  maxavgauprcs=getlog2AUPRC(maxavgres,compostest,comnegtest,cneutrals)
  outtib=outtib %>% bind_rows(tibble(Trait=disease,Term="MaxAvgProb",AUPRCS=maxavgauprcs))
  
  
  #clusres<- clusres %>%
  #  group_by(Cluster) %>%
  #  mutate(Probability = (Probability - mean(Probability)) / sd(Probability)) %>%
  #  ungroup()
  
  
  #Getting tau
  taures=tau(clusres)
  tauauprcs=getlog2AUPRC(taures,compostest,comnegtest,cneutrals)
  outtib=outtib %>% bind_rows(tibble(Trait=disease,Term="Tau",AUPRCS=tauauprcs))
  
  taurankres=taurank(clusres)
  taurankauprcs=getlog2AUPRC(taurankres,compostest,comnegtest,cneutrals)
  outtib=outtib %>% bind_rows(tibble(Trait=disease,Term="Taurank",AUPRCS=taurankauprcs))
  
   #get number of clusters
  outtib$Numclus=(clusres$Cluster%>%unique %>% length)
  print(outtib)
  return(outtib)
}





typenets=c("strict_normal")
typecluss=c("domino","noprop")
nettoget="string"
typediffs=c("magma_1e2","magma_1e5","magma_1e8")

for(typediff in typediffs){
for(typeclus in typecluss){
  
  
  
  for(typenet in typenets){
    #split info
    #Filter this for clusters AllAssign and prefix is typeclus
    #  filter(str_starts(Term,typeclus) | Term=="AllAssign")
    monsplits=read_tsv(paste0("./../data/pos_neg_combinations/",typenet,"_",typediff,"_neutral_diff_string_allclustype.tsv"))%>%
      filter(str_starts(Term,typeclus) | Term=="AllAssign")
    
    monfiles=list.files(paste0("../plexout/",nettoget,"_",typenet),full.names=T)
    

    #Getting auprcs
    monauprcs=mclapply(monfiles,taucalc,typeclus,mc.cores=detectCores()) %>% bind_rows
    #testfiles=monfiles[seq(12,length(monfiles),by=12)]
    #test=lapply(testfiles,taucalc,typeclus)
    monauprcs=monauprcs  %>% separate(col="AUPRCS",sep= ", ", into = c("AUPRC", "Neu_AUPRC"))

    monall=monauprcs
    
    write_tsv(monall,paste0("../results/diff_auprcs/Tau_",typenet,"_diff_auprcs_allclustype_",typeclus,"_",typediff,".tsv"))
  }
}}


