########################
#Main Tree Function
########################
treeft <- function(formula, f.var, s.var, s.var.type, target, data, method="ra", method.point="es", impurity.ft="mean_sse", MINDAT=5, error.method="cv", nfold=5, test.data=NULL) {
###formula: formula for fitting
###data: long or wide format

   print("SPLITTING########################################################")
   ###Splitting
   obj.split = split.treeft(formula=formula, f.var=f.var, s.var=s.var, s.var.type=s.var.type, target=target, data=data, 
                            method=method, method.point=method.point, impurity.ft=impurity.ft, MINDAT=MINDAT)

   if(error.method=="no") { #not pruned
      return(list(tree.model=obj.split$tree.model, coefs=obj.split$coefs, 
                  formula=obj.split$formula, f.var=obj.split$f.var, s.var=obj.split$s.var, s.var.type=obj.split$s.var.type, 
                  target=obj.split$target, method=obj.split$method, method.point=obj.split$method.point,
                  impurity.ft=obj.split$impurity.ft, n.samples=obj.split$n.samples, MINDAT=obj.split$MINDAT))
   }


   print("PRUNING##########################################################")
   ###Pruning
   obj.pruned = prune.treeft(obj=obj.split, target=target, data=data)


   print("SELECTING########################################################")
   ###Selecting
   if(error.method=="test") { 
      obj = select.treeft(obj=obj.pruned, target=target, data=test.data, error.method="test") #test data
   } else if(error.method=="cv") {
      obj = select.treeft(obj=obj.pruned, target=target, data=data, error.method="cv", nfold=nfold)
   } else {
      stop("The error method is not available!")
   }


   print("DONE!#############################################################")
   return(list(tree.model=obj$tree.model, coefs=obj$coefs, 
               formula=obj$formula, f.var=obj$f.var, s.var=obj$s.var, s.var.type=obj$s.var.type, 
               target=obj$target, method=obj$method, method.point=obj$method.point, impurity.ft=obj$impurity.ft, 
               n.samples=obj$n.samples, MINDAT=obj$MINDAT,
               ID.subtrees=obj$ID.subtrees, alpha=obj$alpha, error.subtrees=obj$error.subtrees,
               selected.tree.model=obj$selected.tree.model, error.selected.tree.model=obj$error.selected.tree.model, error.method=obj$error.method, nfold=obj$nfold,
               subtrees=obj$subtrees, coefs.subtrees=obj$coefs.subtrees))
}

############################################################################################################################################################# 