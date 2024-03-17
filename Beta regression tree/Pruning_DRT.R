######################################################
# Cut trees to get a series of subtrees
######################################################
prune.treeft <- function(obj, target, data) {
    #obj: from treeft()
    
    print("prune.treeft...")
    #print(dim(data)); print(head(data))
    
    # ###Data sizes for DC
    # if(shape=="long") {
    #     
    #     data = dataformat.long.to.wide(data=data, id=id, choice=choice, alt.var=alt.var, f.var=obj$f.var, s.var=obj$s.var) #from long to wide format if long 
    #     obj$varying = (1+length(obj$s.var)+1):ncol(data)
    #     
    # } else if(shape=="wide") {
    #     
    #     tmp = data[,varying]
    #     i = sort.list(names(tmp))
    #     
    #     choice.index = which(names(data)==choice)
    #     s.var.index = which2(names(data), obj$s.var)
    #     
    #     data = cbind(data[,c(choice.index,s.var.index)],tmp[,i])
    #     obj$varying = (1+length(obj$s.var)+1):ncol(data)
    #     
    # } else {
    #     stop("The shape is not available.")
    # }
    #print(dim(data)); print(head(data))
    
    T = rep(0, nrow(obj$tree.model))
    subtrees = cbind(obj$tree.model, T)
    
    T = rep(0, nrow(obj$coefs))
    coefs.subtrees = cbind(obj$coefs, T)
    coefs.subtrees = cbind(data.frame(rn = rownames(coefs.subtrees)), coefs.subtrees)
    #print("OK"); print(subtrees); print(coefs.subtrees)
    
    obj.tmp = obj
    nodes = obj$tree.model[,1] #all the nodes in the maximal T
    tnodes = obj$tree.model[is.na(obj$tree.model[,2]),1] #all the terminal nodes the maximal T
    
    Max.T = rep(1, nrow(obj$tree.model)) #1= nodes in the maximal T
    obj$tree.model = cbind(obj$tree.model, Max.T)
    alpha = 0 #c() CHECK!!!
    cut.nodes = c()
    t=0
    
    #print("OK1"); print(obj$tree.model)
    
    while(nrow(obj.tmp$tree.model) > 1) {
        
        #Compute 0=alpha_0 < alpha_1 < alpha_2 < ... < alpha_n.alpha 
        #for           T_0 >     T_1 >     T_2 > ... >     T_n.alpha
        gt = compute.alpha(obj.tmp, data)
        
        #print("g(t)="); print(gt)
        
        ii = which.min(gt) #the weakest node with min g(t)
        alpha = c(alpha, gt[ii])
        nodeid.weakest = obj.tmp$tree.model[ii,1]
        cut.nodes = c(cut.nodes, nodeid.weakest)
        
        #print("OK2"); print("cut:"); print(nodeid.weakest)
        
        obj.tmp$tree.model[ii,2:3] = NA
        
        ii = find.children(nodes, nodeid.weakest) #find its children
        nodes.Tt = nodes[ii==1]
        
        #print("OK3"); print(nodes); print(nodes.Tt)
        
        ii2 = which2(obj.tmp$tree.model[,1], nodes.Tt)
        obj.tmp$tree.model = obj.tmp$tree.model[-ii2,]
        
        T = rep(0, nrow(obj$tree.model))                          #0= no nodes in T
        T[which2(obj$tree.model[,1], obj.tmp$tree.model[,1])] = 1 #1= nodes in T
        obj$tree.model = cbind(obj$tree.model, T)
        
        t=t+1
        T = rep(t, nrow(obj.tmp$tree.model))
        subtrees = rbind(subtrees, cbind(obj.tmp$tree.model, T))
        
        #print("OK4")
        
        obj.coefs = list(tree.model=obj.tmp$tree.model,  
                         formula=obj$formula, f.var=obj$f.var, s.var=obj$s.var, s.var.index=obj$s.var.index, s.var.type=obj$s.var.type, 
                         target=obj$target, method=obj$method, method.point=obj$method.point,impurity.ft=obj$impurity.ft,
                         n.samples=obj$n.samples, MINDAT=obj$MINDAT)
        coefs = get.coefs(obj.coefs, data)
        
        T = rep(t, nrow(coefs))
        
        #print("OK5");print(coefs.subtrees)
        #print("OK6");print(cbind(coefs, T))
        #print("OK7");
        ttmp = cbind(coefs, T)
        ttmp = cbind(data.frame(rn = rownames(ttmp)), ttmp)
        coefs.subtrees = rbind(coefs.subtrees, ttmp) 
        
        #print(obj$tree.model)
        
    }
    
    
    n = ncol(obj$tree.model)
    colnames(obj$tree.model)[4:n] = paste("T", 0:(n-4), sep="")
    #colnames(obj$tree.model)[n] = "T.Root"
    
    
    print("prune.treeft...ok")
    
    return(list(subtrees=subtrees, coefs.subtrees=coefs.subtrees, alpha=alpha, nodes.pruned=cut.nodes,
                tree.model=obj$tree.model, coefs=obj$coefs, 
                formula=obj$formula, f.var=obj$f.var, s.var=obj$s.var, s.var.index=obj$s.var.index, s.var.type=obj$s.var.type, 
                target=obj$target, method=obj$method, method.point=obj$method.point,impurity.ft=obj$impurity.ft, n.samples=obj$n.samples, MINDAT=obj$MINDAT
    ))
    
}


#######################################################
#Compute alpha
compute.alpha <- function(obj, data) {
    
    #print("compute.alpha...")
    #print(dim(data));print(head(data))
    #print(obj$tree.model)
    
    
    x = data[,obj$s.var.index]     #data with split variables
    out.p = predict.treeft(obj, x) #prediction
    
    nr = nrow(obj$tree.model)
    alpha = rep(NA, nr)
    
    nodes = obj$tree.model[,1] #all the nodes
    tnodes = obj$tree.model[is.na(obj$tree.model[,2]),1] #all the terminal nodes
    #print(nodes); print(tnodes); print(out.p)
    for(j in 1:nr) {
        
        # print(c(j, obj$tree.model[j,1]))
        if(!is.na(obj$tree.model[j,2])) {#if intermediate nodes
            
            nodeid = obj$tree.model[j,1]
            ii = find.children(tnodes, nodeid) #find its children
            
            tnodes.Tt = tnodes[ii==1] #terminal nodes of nodeid
            n.Tt = length(tnodes.Tt)
            #for node t
            ii = which2(out.p, tnodes.Tt)
            w = (length(ii)/length(out.p))#weight P_t
            #print(obj)
            #if(obj$method=="nloglik") w=1 #???
            if(obj$impurity.ft == "nloglik"){R.t = impurity.DR(obj$formula, obj$target, data[ii,],  obj$impurity.ft)}
            else{
                impur = impurity.DR(obj$formula, obj$target, data[ii,],  obj$impurity.ft)/length(ii)
                R.t = impur*w
            }

            #for subnodes of node t
            R.Tt = 0 
            for(i in 1:n.Tt) {
                ii = which(out.p == tnodes.Tt[i])
                w = (length(ii)/length(out.p))#weight P_t
                #if(obj$method=="nloglik") w=1 #???
                if(obj$impurity.ft == "nloglik"){R.Tt = R.Tt + impurity.DR(obj$formula, obj$target, data[ii,], obj$impurity.ft)}
                else{
                    impur = impurity.DR(obj$formula, obj$target, data[ii,], obj$impurity.ft)/length(ii)
                    R.Tt = R.Tt + impur*w 
                }
            }
            alpha[j] = (R.t-R.Tt)/(n.Tt-1) 
            alpha[j] = abs(alpha[j]) #Should be non-negative!!!
            #print(R.t); print(R.Tt); print(alpha[j])
            
        }
        
        #print(c(j, obj$tree.model[j,1], alpha[j]))
        
    }
    
    print("compute.alpha...ok")
    
    return(alpha)
}


####################
# Find its children
find.children <- function(nodes, nodeid) {
    
    nr = length(nodes)
    children = rep(0, nr) #1=child, 0=not
    
    for(i in 1:nr) {
        
        ii = nodes[i]
        while(ii > 0) {
            ii = floor(ii / 2)
            if(ii==nodeid)  children[i] = 1
        }
    }  
    
    return(children)
}

################################################################################################################END




