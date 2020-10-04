get_extant <- function(tm,tree){
  extinct = tree$extinct[tree$extinct$brts<tm,]
  extant = tree$extant[tree$extant$brts<=tm,]
  extant$clade = NULL
  extant$t_ext = 0
  extinct$t_ext[tree$t_ext>tm] = 0 
  
  extended_tree = rbind(extant,extinct)
  extended_tree = extended_tree[order(extended_tree$brts),]
  
  extant_species = NULL
  for(i in 1:nrow(extended_tree)){
    if(extended_tree$t_ext[i]==0){
      extant_species = rbind(extant_species,data.frame(brts=extended_tree[i,"brts"],
                                                       parent=extended_tree[i,"parent"],
                                                       child=extended_tree[i,"child"]))
    }else{
      if(extended_tree$t_ext[i]!=999){
        kids = which((extended_tree$child[i]==extended_tree$parent))
        K=kids
        if(length(kids)!=0){
          while(  length(which((extended_tree$t_ext==0)&kids))==0 ){
            K2 = NULL
            for(i in 1:length(kids)){
              k = which((kids[i]==extended_tree$parent))
              K2 = c(K2,k)
            }
            kids = K2
          }
          mk = min(which((extended_tree$t_ext==0)&kids))
          parents = extended_tree$parent[mk]
          while(parents[1]!=extended_tree$child[i]){
            grandparent = extended_tree$par
            parents = c(1)
          }
          mk = min(kids)
          bt = extended_tree[mk,"brts"]
          extant_species =  rbind(extant_species,data.frame(brts=extended_tree[i,"brts"],
                                                            parent=extended_tree[i,"parent"],
                                                            child=extended_tree[mk,"child"]))
          extended_tree$t_ext[mk] = 999
        }
      }
    }
  }
  map_child = extended_tree$child
  
  return(extant_species)
}


transf <- function(name_spe,vec){
  which(vec==name_spe)
}

newick<- function(tree,CT){
  n<-nrow(tree)
  child.nms<-as.character(tree$child)
  parent.nms<-as.character(tree$parent)
  species.nms<-unique(child.nms,parent.nms)
  n.species<-length(species.nms)
  CT<-rep(CT,n.species)
  for (i in seq(n,1)){
    nw<-paste("(",parent.nms[i],":",as.character(CT[which(species.nms==parent.nms[i])]-tree$brts[i]),",",child.nms[i],":",as.character(CT[which(species.nms==child.nms[i])]-tree$brts[i]),")", sep = "")
    j<-which(parent.nms[i]==child.nms)
    rp<-which(parent.nms==child.nms[j])
    if (length(rp)>0){
      parent.nms[rp]<-nw
    }
    species.nms[which(species.nms==child.nms[j])]<-nw
    child.nms[j]<-nw
    CT[j]<-CT[j]-(CT[which(species.nms==parent.nms[i])]-tree$brts[i])
  }
  return(paste(child.nms[1],";",sep=""))
  #return(child.nms)
}

PDTT_plot <- function(tree){
  ct = max(tree$brts)
  times = seq(0,ct,length.out = 1000)
  times = sort(times,decreasing = T)
  PD=NULL
  for(i in 1:length(times)){
    G=GPD(times[i],tree[-1,])
    
    PD=rbind(PD,data.frame(time=rep(times[i],nrow(G)),
                         P=colSums(G)/(nrow(G)-1),
                         lineage=as.character(1:nrow(G))))
  }
  g1 = ggplot(PD)+
    geom_line(aes(x=time,y=P,colour=lineage,alpha=0.5))+
    theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
return(g1)
}

phylo2emph <- function(phylo){
  #transformation of ultrametric trees into data frame
  tree = DDD::phylo2L(phylo)
  brts_dd = tree[,1]
  brts = cumsum(-diff(c(brts_dd,0)))
  
  tree = list(extant = data.frame(brts = c(0,brts[-length(brts)]),
                                  parent=c(1,abs(tree[,2][-1])),
                                  child=abs(tree[,3])),
              extinct = data.frame(brts = numeric(),
                                   parent = numeric(),
                                   child = numeric(),
                                   t_ext = numeric()),
              
              ct=brts_dd[1])
  
  return(tree)
}


phylo2tree <- function(tree){
  # to map newick trees into ther xxxx format
  ltt = ape:::ltt.plot.coords(tree)
  t = diff(ltt[,1])
  ltt = ltt[-1,]
  n = ltt[,2]
  E = diff(n)
  E[E==-1] = 0
  return(list(brts=cumsum(t),to=E,t_e))
}
# emphasisLD
GPD2<-function(tm,tree){
  # input: an ultramedric tree defined by a data.frame
  # with columns brts, parent, 
  i1<-tree$brts<=tm 
  newtree<-tree[i1,]
  d<-nrow(newtree)
  gpd<-matrix(0,ncol=d+1,nrow=d+1)
  sets<-as.list(1:(d+1))
  #sets <-as.list(0:(d))
  map_child = 1:length(newtree$child)
  newtree$map_child = map_child
  map_parent = 1:length(newtree$parent)
  newtree$map_parent = map_parent
  for (i in d:1){
    s1<-lapply(sets,function(s,e){if(e%in%s) return(s) else return(-1)},e=newtree$map_parent[i])
    s2<-lapply(sets,function(s,e){if(e%in%s) return(s) else return(-1)},e=newtree$map_child[i])
    oldset1<-sapply(s1, function(x){x[1]==-1})
    oldset2<-sapply(s2, function(x){x[1]==-1})
    set1<-sets[!oldset1][[1]]
    set2<-sets[!oldset2][[1]]
    gpd[set1,set2]<-tm-newtree$brts[i]
    sets<-c(sets[oldset1&oldset2],list(c(set1,set2)))
  }
  return(gpd+t(gpd))
}


GPD<-function(tm,tree){
  # input: an ultramedric tree defined by a data.frame
  # with columns brts, parent, 
  i1<-(tree$brts<=tm)
  newtree<-tree[i1,]
  d<-nrow(newtree)
  gpd<-matrix(0,ncol=d+1,nrow=d+1)
  sets<-as.list(1:(d+1))
  for (i in d:1){
    s1<-lapply(sets,function(s,e){if(e%in%s) return(s) else return(-1)},e=newtree$parent[i])
    s2<-lapply(sets,function(s,e){if(e%in%s) return(s) else return(-1)},e=newtree$child[i])
    oldset1<-sapply(s1, function(x){x[1]==-1})
    oldset2<-sapply(s2, function(x){x[1]==-1})
    set1<-sets[!oldset1][[1]]
    set2<-sets[!oldset2][[1]]
    gpd[set1,set2]<-tm-newtree$brts[i]
    sets<-c(sets[oldset1&oldset2],list(c(set1,set2)))
  }
  return(gpd+t(gpd))
}

# more utilities  (emphasis)

n_from_time <- function(tm,tree){
  # return N at tm.
  extended_tree = extend_tree(tree)

  n = cumsum(extended_tree$event)+cumsum(extended_tree$event-1)+1
  brts = extended_tree$brts
  if(tm==0) tm = 0.000000000000001
  N = n[max(which(brts < tm))]
  return(N)
} 

n_for_all_bt <- function(tree){
  brts = c(tree$extant$brts,
           tree$extinct$brts,
           tree$extinct$t_ext)
  brts = c(0,sort(brts[brts!=0]))
  n = sapply(brts,n_from_time,tree)
  return(n)
}

extend_tree <- function(tree){
  if(is.null(tree$extinct)){
    tree = list(extant=tree,extinct=data.frame(brts=NULL,t_ext=NULL))
  } 
  extended_tree = data.frame(brts = c(tree$extant$brts,
                                      tree$extinct$brts,
                                      tree$extinct$t_ext),
                             event = c(rep(1,nrow(tree$extant)),
                                       rep(1,nrow(tree$extinct)),
                                       rep(0,nrow(tree$extinct))))
  extended_tree = extended_tree[order(extended_tree$brts),]
  extended_tree = rbind(extended_tree,data.frame(brts=tree$ct,event=2))
  if(extended_tree$brts[2]==0){
    extended_tree = extended_tree[-1,]
  }
  return(extended_tree)
}

#phylodiversity <- function(tm,tree,soc){
#  i1<-tree$brts<=tm 
#  i2<-tree$to==0&i1
#  i3<-tree$t_ext%in%tree$brts[i2]
#  dt<-diff(c(0,tree$brts[i1&!i2&!i3],tm))
#  return(sum(dt*(soc:(length(dt)+soc-1))))
#}

get_current_species <- function(tm,tree){
  species = c(tree$extant$child[tree$extant$brts<tm],
              tree$extinct$child[tree$extinct$t_ext>tm])
  return(species)
}

data_to_table <- function(df,replicant,left,right){
  df = df[df$rep==replicant,]
  df = df[df$iteration %in% left:right,]
  summ = data.frame(lfhat = mean(df$fhat),sd_fhat=sd(df$fhat),mad_fhat=mad(df$fhat),replicant=replicant,par1=median(df$par1),par2=median(df$par2),par3=median(df$par3),par4=median(df$par4),E_time = sum(df$E_time)/60, M_time = sum(df$M_time)/60, sample_size=mean(df$sample_size))
  return(summ)
}

AIC_llik <- function(LogLik,k){
  aic <- (2*k)-(2*LogLik)
  return(aic)
}

AICw <- function(l1,l2,k1,k2){
  IC <- AIC_llik(c(l1,l2),c(k1,k2))
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights[1])
}

vectors2phylo <- function(list){
  t=list$wt
  n=list$n
  E=list$E
  S=list$S
  ct=sum(t)
  newick = paste(sl[1],";",sep="")
  N=1
  identf = data.frame(Spec="aa",Time=0) # Labels of species
  for (i in 1:(length(t)-1)){
    # speciation
    sumt = sum(t[1:i])
    if( is.null(S)){
      BD = sample(1:N,1)
      species = as.character(identf[BD,1])
    }else{
      species = S[i]
    }
    if (E[i] == 1){
      ind = regexpr(species,newick)[1]-1
      atm = sumt-identf[which(identf[,1]==species),2]
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
      identf = rbind(identf,data.frame(Spec=substr(sl[i+1],1,2),Time=sumt))
      identf[identf$Spec == species,2] = sumt
      N = N+1
    }
    # extinction
    if (E[i]==0){
      ind = regexpr(species,newick)[1] + 2
      atm = sumt-identf[which(identf[,1]==species),2]
      identf = identf[!identf$Spec==species,]
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
      N=N-1
    }
  }
  newick = compphyl(newi=newick,identf=identf,ct=ct)
  newick = read.tree(text=newick)
  return(newick)
}



tree2phylo <- function(tree,initspec=1){
  wt = -diff(c(0,tree$brts))
  to = tree$to
  to[to==2] = 1
  ct = sum(wt)
  newick = paste(sl[1],";",sep="")
  N = 1
  identf = data.frame(Spec="a",Time=0) # Labels of species
  for (i in 1:(length(wt)-1)){
    # speciation
    bt = sum(wt[1:i])
    BD = sample(1:N,1)
    species = as.character(identf[BD,1])  
    if (to[i] == 1){
      ind = regexpr(species,newick)[1]-1
      atm = bt-identf[which(identf[,1]==species),2]
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
      identf = rbind(identf,data.frame(Spec=substr(sl[i+1],1,2),Time=bt))
      identf[identf$Spec == species,2] = bt
      N = N+1
    }
    # extinction
    if (to[i]==0){
      ind = regexpr(species,newick)[1] + 2
      atm = bt-identf[which(identf[,1]==species),2]
      identf = identf[!identf$Spec==species,]
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
      N=N-1
    }
  }
  newick = compphyl(newi=newick,identf=identf,ct=ct)
  newick = read.tree(text=newick)
  return(newick)
}

compphyl <- function(newi,identf,ct){
  #set to extant species to the present time
  identf[,1] = as.character(identf[,1])
  identf[,2] = ct-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}

# time calculation
get.time <- function(time,mode='sec'){
  dif = proc.time()-time
  ti = as.numeric(dif[3])
  if(mode == 'min')  ti = ti/60
  if(mode == 'hou') ti = ti/3600
  return(ti)
}


