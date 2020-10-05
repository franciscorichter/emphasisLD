## testing augmentation function 
load("~/Dropbox/github/emphasisLD/data/FamilyBirdTrees.Rdata")
phylo = FamilyBirdTrees$Bucconidae$tree
tree = phylo2emph(phylo)
diversification_model = list(model="rpd1",pars=c(0.2,0.6,-0.01))
diversification_model = list(model="ldpd",pars=c(0.2,0.6,-0.02,0))

tree_extant = get_extant(tm,tree)

## check if they work 
sample_tree_full(diversification_model,ct)

sample_tree(diversification_model = diversification_model,phylo = phylo)


## check if emphasis and emphasisLD are equivalent
emphasis:::augment_tree(brts = as.numeric(ape:::branching.times(phylo)),
                        pars = diversification_model$pars,
                        model = "rpd1",
                        soc = 2,
                        sampler_spe = "rpd1")

## now check with GPD --- not working yet... GPD function failing

diversification_model = list(model="ldpd",
                             pars=c(0.1,0.5,-0.175,0))

lambda.rpd1(6,tree,diversification_model$pars)

lambda.ldpd(6,tree,c(diversification_model$pars,0.1))

## check if both models are equivalent for 