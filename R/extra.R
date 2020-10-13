get_extant <- function(tm,tree){ #test1Ok, test2Ok
  # obtain extant tree from full tree
  extinct = tree$extinct[tree$extinct$brts<tm,]
  extant = tree$extant[tree$extant$brts<tm,]
  extant$clade = NULL
  if(nrow(extant)>0) extant$t_ext = 0
  extinct$t_ext[extinct$t_ext>tm] = 0 
  
  ext_tree = rbind(extant,extinct)
  ext_tree = ext_tree[order(ext_tree$brts),]
  for(i in 1:nrow(ext_tree)){
    # check if the species survived to the present, if not
    if((ext_tree$t_ext[i]!=0)&(ext_tree$t_ext[i]!=999)){
      # search for an extant lineage
      kids = (ext_tree$child[i]==ext_tree$parent)
      key = FALSE
      while(sum(kids)!=0 & !key){ # until there are no more descendants or...
        # if there is no child that survives to present, check next generation 
        if(sum((ext_tree$t_ext[kids]==0))==0){
          K2 = rep(F,nrow(ext_tree))
          for(j in which(kids)){
            k = (ext_tree$child[j]==ext_tree$parent)
            K2 = K2 | k 
          }
          kids = K2
        }else{
          key = TRUE
          # Search for species that survive to present
          mk = min(which(kids&(ext_tree$t_ext==0)))
          child = ext_tree[mk,"child"]
          parent = ext_tree[mk,"parent"]
          ext_tree$t_ext[mk] = 999
          ext_tree$parent[ext_tree$parent==child] = ext_tree$child[i]
          # check if parent survived to present, if not
          while(parent!=ext_tree$child[i]){
            mk = which(ext_tree$child==parent)
            child = parent
            parent = ext_tree[mk,"parent"]
            ext_tree$t_ext[mk] = 999
            ext_tree$parent[ext_tree$parent==child] = ext_tree$child[i]
            
          }
          ext_tree$t_ext[i] = 0
        }
      }
    }
  }
  ext_tree = ext_tree[ext_tree$t_ext==0,]
  ext_tree$t_ext = NULL
  if(nrow(ext_tree)>2){
    map_child = ext_tree$child[3:nrow(ext_tree)]
    for(i in 1:length(map_child)){
      ext_tree[ext_tree$parent==map_child[i],"parent"]=3+i
      ext_tree[i+2,"child"]=3+i
    }
  }
  return(ext_tree)
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



phylo2emph <- function(phylo){
  # test1Ok, test2Ok
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
  tree$extant[3:nrow(tree$extant),"parent"]=tree$extant[3:nrow(tree$extant),"parent"]+1
  tree$extant$child = tree$extant$child + 1 
  tree$extant$parent[1:2]=1
  tree$extinct$parent=tree$extinct$parent+1
  tree$extinct$child=tree$extinct$child+1
  return(tree)
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


GPD<-function(tree,ct){
  # input: an ultramedric tree defined by a data.frame
  # with columns brts, parent, 
  i1<-(tree$brts<=ct)
  tree<-tree[i1,]
  if(tree$brts[1]==tree$brts[2]){
    tree = tree[-1,]
    tree$parent = tree$parent - 1
    tree$child = tree$child -1
    tree$parent[1] = 1
  }
  d<-nrow(tree)
  gpd<-matrix(0,ncol=d+1,nrow=d+1)
  sets<-as.list(1:(d+1))
  for (i in d:1){
    s1<-lapply(sets,function(s,e){if(e%in%s) return(s) else return(-1)},e=tree$parent[i])
    s2<-lapply(sets,function(s,e){if(e%in%s) return(s) else return(-1)},e=tree$child[i])
    oldset1<-sapply(s1, function(x){x[1]==-1})
    oldset2<-sapply(s2, function(x){x[1]==-1})
    set1<-sets[!oldset1][[1]]
    set2<-sets[!oldset2][[1]]
    gpd[set1,set2]<-ct-tree$brts[i]
    sets<-c(sets[oldset1&oldset2],list(c(set1,set2)))
  }
  return(gpd+t(gpd))
}

# more utilities  (emphasis)

n_from_time <- function(tm,tree){
  # return N at tm.
  if(tm==0){
    N=2
  }else{
    extended_tree = extend_tree(tree)

    n = cumsum(extended_tree$event)+cumsum(extended_tree$event-1)+1
    brts = extended_tree$brts
    if(tm==0) tm = 0.000000000000001
    N = n[max(which(brts < tm))]
  }
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


foo <- function(phylo, metric = "colless") {
  if(!is.null(phylo$tree)) phylo = phylo$tree
  if (metric == "colless") {
    xx <- apTreeshape:::as.treeshape(phylo)  # convert to apTreeshape format
    apTreeshape:::colless(xx, "yule")  # calculate colless' metric
  } else if (metric == "gamma") {
    ape:::gammaStat(phylo)
  } else stop("metric should be one of colless or gamma")
}

#phylodiversity <- function(tm,tree,soc){
#  i1<-tree$brts<=tm 
#  i2<-tree$to==0&i1
#  i3<-tree$t_ext%in%tree$brts[i2]
#  dt<-diff(c(0,tree$brts[i1&!i2&!i3],tm))
#  return(sum(dt*(soc:(length(dt)+soc-1))))
#}

get_required_sampling_size <- function(M,tol=.05){
  n <- M$sample_size
  f<-  M$fhat
  hlp<-MASS:::rlm(f~I(1/n),weights = n)
  ab<-coef(hlp)
  
  f.r<-ab[1]-tol
  n.r<-ceiling(ab[2]/(f.r-ab[1]))
  return(n.r)
}

get_current_species <- function(tm,tree){
  species = c(tree$extant$child[tree$extant$brts<tm],
              tree$extinct$child[tree$extinct$t_ext>tm])
  return(species)
}

# time calculation
get_time <- function(time,mode='sec'){
  dif = proc.time()-time
  ti = as.numeric(dif[3])
  if(mode == 'min')  ti = ti/60
  if(mode == 'hou') ti = ti/3600
  return(ti)
}


