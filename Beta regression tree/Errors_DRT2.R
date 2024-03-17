###################################################
#Compute test or CV errors
###################################################

#####################################################################
# Compute test errors
test.error <- function(obj, target, data){
    ###obj: object from treeft() or prune.treeft()
    ###data: long/wide format
    
    print("test.error...")
    #print(dim(data)); print(head(data))
    
    x = data[,obj$s.var]                                              
    z = data[,obj$f.var] # n.classes = 1
    colnames(z) = obj$f.var
    
    
    
    #print(x); print(z)
    ###compute the predicted probability
    tmp = predict.treeft.prob(obj=obj, x, z)
    p = tmp$pred #check order!!
    
    
    ###compute the error
    n = nrow(data)
    sse = sum((data[,target]-p)^2)
    
    print("test.error...ok")
    return(sse)
    
}




######################################################################
#Compute Cross-Validation (CV) errors
cv.error <- function(obj, target, data, nfold) {
    ###obj: from prune.treeft()
    
    print("cv.error...")
    #print(dim(data)); print(head(data))
    
    ###Compute CV errors
    n.trees = 1+max(obj$subtrees$T) # number of T0, T1, T2, ...=length(unique(obj.subtrees$subtrees$T))
    CVerrors = matrix(NA, nfold, n.trees)
    
    V = sample(1:nfold, obj$n.samples, replace=TRUE) # V-fold selected index
    
    for(iv in 1:nfold) { #V-fold Cross-Validation(CV)
        
        print(paste("FOLD V =", iv, " ========================================="))
        
        ii = which(V==iv)
        data.train = data[-ii,] #train data
        data.test  = data[ii,] #test data
        
        print("SPLITTING##################CV")
        ###Splitting with train data
        obj.tree = split.treeft(formula=obj$formula, f.var=obj$f.var, s.var=obj$s.var, s.var.type=obj$s.var.type, 
                                target=target, data=data.train, method=obj$method, method.point=obj$method.point, 
                                impurity.ft=obj$impurity.ft, MINDAT=obj$MINDAT*(1-1/nfold)) #MINDAT=obj$MINDAT*(1-1/nfold)): considering sample size
        
        #print("OK: obj.tree:"); print(obj.tree)
        
        
        print("PRUNING#####################CV")
        ###Pruning with alpha' = g(t)
        ###Compute new alpha from 0=alpha_0 < alpha_1 < alpha_2 < ... < alpha_n.alpha-1: more precise alpha for cv
        ###                  for        T_0 >     T_1 >     T_2 > ... >     T_n.alpha-1
        n.alpha = length(obj$alpha) # 0, 1, ..., n.alpha-1
        alpha.prime = c() 
        for(k in 1:(n.alpha-1)) {
            alpha.prime = c(alpha.prime, sqrt(obj$alpha[k]*obj$alpha[k+1]))   
        }
        obj.subtrees = prune.treeft.cv(obj=obj.tree, alpha.prime=alpha.prime, target=target, data=data.train) #CHECK!!!
        
        #print("OK: obj.subtrees:"); print(obj.subtrees)
        
        
        print("SELECTING####################CV")
        ###Computing CV errors
        n.trees = max(obj.subtrees$subtrees$T) # number of T1, T2, ... but not T0
        ii = which(obj.subtrees$subtrees$T==0)
        jj = ncol(obj.subtrees$subtrees)
        obj.subtrees$tree.model = obj.subtrees$subtrees[ii,-jj]
        CVerrors[iv,1] = test.error(obj=obj.subtrees, target=target, data=data.test)
        
        for(t in 1:n.trees) {
            
            jj = ncol(obj.subtrees$subtrees)
            ii = which(obj.subtrees$subtrees$T==t)
            obj.subtrees$tree.model = obj.subtrees$subtrees[ii,-jj]
            
            jj = ncol(obj.subtrees$coefs.subtrees)
            ii = which(obj.subtrees$coefs.subtrees[,jj]==t)
            obj.subtrees$coefs = obj.subtrees$coefs.subtrees[ii,-c(1, jj)]
            
            rownames(obj.subtrees$coefs) = obj.subtrees$coefs.subtrees[ii,-jj]$rn
            tmp = test.error(obj=obj.subtrees, target=target, data=data.test)  
            CVerrors[iv,(t+1)] = tmp
            
        }
        
    }
    
    print("FOLD ends =========================================")
    #print("OK: CVerrors:"); print(CVerrors)
    
    ###CV error for T0, T1, T2, ....
    CVerrors = apply(CVerrors, 2, mean) #2=column mean
    
    #print(CVerrors)
    
    print("cv.error...ok")
    
    return(CVerrors)
}



#################################################################################
# Cut trees to get a series of subtrees for CV
prune.treeft.cv <- function(obj, alpha.prime, target, data) {
    ###obj: from treeft()
    ###alpha: from prune.treeft()
    
    print("prune.treeft.cv...")
    #print(dim(data)); print(head(data))
    
    T = rep(0, nrow(obj$tree.model))
    subtrees = cbind(obj$tree.model, T)
    
    T = rep(0, nrow(obj$coefs))
    coefs.subtrees = cbind(obj$coefs, T)
    coefs.subtrees = cbind(data.frame(rn = rownames(coefs.subtrees)), coefs.subtrees)
    
    obj.tmp = obj
    nodes = obj$tree.model[,1] #all the nodes in the maximal T
    tnodes = obj$tree.model[is.na(obj$tree.model[,2]),1] #all the terminal nodes the maximal T
    
    Max.T = rep(1, nrow(obj$tree.model)) #1= nodes in the maximal T
    obj$tree.model = cbind(obj$tree.model, Max.T)
    cut.nodes = c()
    t=0
    
    n.alpha= length(alpha.prime)+1 # number of 0, 1, 2, ...
    for(k in 2:n.alpha){ #if last alpha???
        
        #print("PRUNING>>>");print(t)
        
        if(k < n.alpha) {
            gt = compute.alpha(obj.tmp, data) #Compute alpha with train data
            ii = which(gt <= alpha.prime[k]) #if nothing exists???
            
            obj.tmp$tree.model[ii,2:3] = NA #cut and declare it as terminals
            nodeid.weakest = obj.tmp$tree.model[ii,1]
            cut.nodes = c(cut.nodes, nodeid.weakest)
            
            #print("nodeid.weakest="); print(nodeid.weakest) 
            #print("cut.nodes="); print(cut.nodes)
            
            
            if(length(nodeid.weakest) >= 1) { #if weak nodes exist, find their child nodes, and then cut them
                ii = rep(0, length(nodes))
                for(i in 1:length(nodeid.weakest)) {
                    jj = find.children(nodes, nodeid.weakest[i])#find its children
                    ii = ii + jj
                    #print(nodeid.weakest[i]); print(jj); print(ii)
                }
                
                nodes.Tt = nodes[ii >= 1]             
                ii2 = which2(obj.tmp$tree.model[,1], nodes.Tt)
                obj.tmp$tree.model = obj.tmp$tree.model[-ii2,] #NA error???
                
                #print("Of nodes:"); print(nodes); print("prune the following nodes:"); print(nodes.Tt)
            }
        } else { # gt < the largest alpha
            obj.tmp$tree.model = obj.tmp$tree.model[1,] #only root node
            obj.tmp$tree.model[,2:3] = NA #cut and declare it as terminals
        }
        
        T = rep(0, nrow(obj$tree.model))                          #0= no nodes in T
        T[which2(obj$tree.model[,1], obj.tmp$tree.model[,1])] = 1 #1= nodes in T
        obj$tree.model = cbind(obj$tree.model, T)
        
        #print(obj$tree.model)
        
        t=t+1
        T = rep(t, nrow(obj.tmp$tree.model))
        subtrees = rbind(subtrees, cbind(obj.tmp$tree.model, T))
        
        obj.coefs = list(tree.model=obj.tmp$tree.model,  
                         formula=obj$formula, f.var=obj$f.var, s.var=obj$s.var,  s.var.index=obj$s.var.index, s.var.type=obj$s.var.type, 
                         target=obj$target, method=obj$method, impurity.ft=obj$impurity.ft, n.samples=obj$n.samples, MINDAT=obj$MINDAT)
        coefs = get.coefs(obj.coefs, data)
        
        #print(coefs)
        
        T = rep(t, nrow(coefs))
        ttmp = cbind(coefs, T)
        ttmp = cbind(data.frame(rn = rownames(ttmp)), ttmp)
        coefs.subtrees = rbind(coefs.subtrees, ttmp)
        
        #print(coefs.subtrees)
    }
    
    
    n = ncol(obj$tree.model)
    colnames(obj$tree.model)[4:n] = paste("T", 0:(n-4), sep="")
    #colnames(obj$tree.model)[n] = "T.Root"
    
    
    print("prune.treeft.cv...ok")
    
    return(list(subtrees=subtrees, coefs.subtrees=coefs.subtrees, alpha=alpha.prime, nodes.pruned=cut.nodes,
                tree.model=obj$tree.model, coefs=obj$coefs, 
                formula=obj$formula, f.var=obj$f.var, s.var=obj$s.var, s.var.type=obj$s.var.type, 
                target=obj$target, method=obj$method, impurity.ft=obj$impurity.ft, n.samples=obj$n.samples, MINDAT=obj$MINDAT
    ))
    
}

#################################################################################################################################


