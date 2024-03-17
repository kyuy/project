###################################################
#Select an optimal tree by test error or CV error
###################################################

select.treeft <- function(obj, target, data, error.method, nfold=5){
    #obj: object from prune.treeft()
    #data: long or wide format
    
    print("select.treeft...")
    
    if(error.method=="test") { #Using independent test data
        
        nfold = NA
        n.trees = max(obj$subtrees$T)
        
        ii = which(obj$subtrees$T==0)
        jj = ncol(obj$subtrees)
        obj$tree.model = obj$subtrees[ii,-jj]
        terrors = test.error(obj=obj, target=target, data=data)

        
        for(t in 1:n.trees) {
            
            jj = ncol(obj$subtrees)
            ii = which(obj$subtrees$T==t)
            obj$tree.model = obj$subtrees[ii,-jj]
            
            jj = ncol(obj$coefs.subtrees)
            ii = which(obj$coefs.subtrees[,jj]==t)
            obj$coefs = obj$coefs.subtrees[ii,-c(1, jj)]
            
            rownames(obj$coefs) = obj$coefs.subtrees[ii,-jj]$rn
            tmp = test.error(obj=obj, target=target, data=data)
            terrors = c(terrors, tmp)
            
        }
        
    } else { #Using cross-validation (CV)
        
        terrors = cv.error(obj=obj, alt.var=alt.var, choice=choice, id=id, varying=varying, data=data, shape=shape, nfold=nfold)
        
    }
    
    
    #Select the final tree with the minimum error
    tt = which.min(terrors)
    test.error = terrors[tt]
    tt=tt-1
    
    jj = ncol(obj$subtrees)
    ii = which(obj$subtrees$T==tt)
    obj$tree.model = obj$subtrees[ii,-jj]
    
    jj = ncol(obj$coefs.subtrees)
    ii = which(obj$coefs.subtrees[,jj]==tt)
    obj$coefs = obj$coefs.subtrees[ii,-jj]
    rownames(obj$coefs) = obj$coefs$rn
    obj$coefs = obj$coefs[,-1]
    
    T = 1:length(obj$alpha); T= T-1 #Tree number
    ID.subtrees = paste("T",T, sep="")
    selected.tree.model = paste("T",tt, sep="") #Tree number selected
    
    print("select.treeft...ok")
    
    return(list(tree.model=obj$tree.model, coefs=obj$coefs,  
                formula=obj$formula, f.var=obj$f.var, s.var=obj$s.var, s.var.index=obj$s.var.index, s.var.type=obj$s.var.type, 
                target=obj$target,
                method=obj$method, method.point=obj$method.point,impurity.ft=obj$impurity.ft, n.samples=obj$n.samples, MINDAT=obj$MINDAT,
                ID.subtrees=ID.subtrees, alpha=obj$alpha, error.subtrees = terrors,
                selected.tree.model=selected.tree.model, error.selected.tree.model=test.error, error.method=error.method, nfold=nfold,
                subtrees=obj$subtrees, coefs.subtrees=obj$coefs.subtrees))
    
    
}
########################################################################