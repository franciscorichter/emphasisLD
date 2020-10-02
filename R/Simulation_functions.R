sample_tree_full <- function(diversification_model,ct){  
  
  cbt  = 0 
  tree = list(extant=data.frame(brts=c(0,0),
                                parent=c(0,1),
                                child=c(1,2),
                                clade=c(0,1)),
              extinct=data.frame(brts=NULL,
                                 parent=NULL,
                                 child=NULL,
                                 t_ext=NULL),
              ct=ct)
  next_bt = 0 
  while(cbt < ct & sum(tree$extant$clade==1)>0 & sum(tree$extant$clade==0)>0){
    
    tree_extant = get_extant(cbt,tree)
    brts = c(tree_extant$brts,ct)
    
    ### Draw speciation 
    event_time = draw_event_time(cbt,
                                 ct,
                                 diversification_model,
                                 tree=list(extant=tree_extant,
                                           extinct=data.frame(brts=NULL,
                                           parent=NULL,
                                           child=NULL,
                                           t_ext=NULL),
                                           ct = ct))
    if(event_time<ct){
      
      ## resolve allocation and update tree
      allocation = draw_allocation_full(event_time,
                                  diversification_model,
                                  tree)
      row_extant = which(tree$extant$child==allocation$species)
      if(allocation$event=="e"){
        to_add = tree$extant[tree$extant$child==allocation$species,]
        tree$extinct = rbind(tree$extinct,data.frame(brts=to_add$brts,
                                                     parent=to_add$parent,
                                                     child=to_add$child,
                                                     t_ext=event_time))
        tree$extant = tree$extant[-row_extant,]
      }
      
      if(allocation$event=="s"){
        next_child = max(tree$extant$child,tree$extinct$child)+1
        to_add = data.frame(brts=event_time,parent=allocation$species,child=next_child,clade=tree$extant$clade[row_extant])
        tree$extant = rbind(tree$extant,to_add)
      }

    }
    cbt = min(event_time, ct)
  }
  
  return(tree)
}

draw_allocation_full <- function(event_time,
                                 diversification_model,
                                 tree){
  
  # calculate rates
  l = speciation_rate(tm = event_time,
                      tree = tree,
                      diversification_model = diversification_model)
  m = extinction_rate(tm = event_time,
                      tree = tree,
                      diversification_model = diversification_model)
  mu = sum(m)
  lambda = sum(l)
  ## choose speciation or extinction 
  to = sample(c("e","s"),1,prob = c(mu/(mu+lambda),lambda/(mu+lambda))) 
  
  ## choose species
  current_species <- get_current_species(tm = event_time,
                                         tree = tree)
  if(to=="e") probs = m
  if(to=="s") probs = l
  species = sample(current_species,prob = probs,size=1)
  

  return(list(event=to,species=species))
}

