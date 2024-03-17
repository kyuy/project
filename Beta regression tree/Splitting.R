########################
# Splitting
########################
split.treeft <- function(formula, f.var, s.var, s.var.type, target, data, method="ra", method.point="es", impurity.ft="mean_sse", MINDAT=5) {
    ###formula: formula for fitting
    ###data: long or wide format
    ###f.var: fit the regression model but not split the nodes
    ###s.var: split the nodes but not fit the regression model
    print("split.treeft...")
    #print(dim(data)); print(head(data))
    
    
    n.samples = nrow(data) #sample size (no missing is assumed!)
    s.var.index = which2(names(data), s.var)
    
    
    #print(s.var); print(s.var.index); print(s.var.type)
    #print(n.samples)
    
    ###Tree initiation
    tnodeid=rep(1, n.samples)
    nodeid=1
    nodes.all.cumulative = c()
    nodeid.cumulative = c()
    split.var.cumulative = c()
    split.set.cumulative = c()
    split.set.elements.cumulative = list()
    if(MINDAT < 1) MINDAT = round(n.samples*MINDAT) 
    
    ###Splitting recursively
    out = rsplit(formula=formula, s.var.index=s.var.index, s.var.type=s.var.type, target=target,
                 data=data, method=method, method.point=method.point, impurity.ft=impurity.ft, MINDAT=MINDAT, 
                 split.var.cumulative=split.var.cumulative, split.set.cumulative=split.set.cumulative, split.set.elements.cumulative=split.set.elements.cumulative, 
                 nodes.all.cumulative=nodes.all.cumulative, nodeid.cumulative=nodeid.cumulative, tnodeid=tnodeid, nodeid=nodeid)
    
    
    ###Tree model
    tree.model = data.frame(matrix(NA, length(out$nodes.all.cumulative), 4))
    colnames(tree.model) = c("Node.ID","Split.Var","Point/Set", " ")
    
    ii = which2(out$nodes.all.cumulative, out$nodeid.cumulative)
    tree.model[,1] = out$nodes.all.cumulative
    tree.model[ii,2] = s.var[out$split.var.cumulative] 
    tree.model[ii,3] = out$split.set.cumulative 
    tree.model[ii,4] = unlist(out$split.set.elements.cumulative)

    ii = which(is.na(tree.model[,4])==FALSE) #for categorical variables
    tree.model[ii,3] =  tree.model[ii,4]
    tree.model = tree.model[,-4]
    
    obj = list(tree.model=tree.model,  
               formula=formula, f.var=f.var, s.var=s.var, s.var.index=s.var.index, s.var.type=s.var.type, 
               target=target, method=method, impurity.ft=impurity.ft, n.samples=n.samples, MINDAT=MINDAT)

    coefs = get.coefs(obj, data)
    
    
    print("split.treeft...ok")
    
    return(list(tree.model=tree.model, coefs=coefs, formula=formula, f.var=f.var, s.var=s.var, s.var.index=s.var.index, s.var.type=s.var.type, 
                target=target, method=method, method.point=method.point, impurity.ft=impurity.ft, n.samples=n.samples, MINDAT=MINDAT))
    
}

############################################################################################################
#Split Recursively
rsplit <- function(formula, s.var.index, s.var.type, target,
                   data, method, method.point, impurity.ft, MINDAT, 
                   split.var.cumulative, split.set.cumulative, split.set.elements.cumulative, 
                   nodes.all.cumulative, nodeid.cumulative, tnodeid, nodeid) {
    
    
    print("rsplit...")
    print(paste("NODE ID #####", nodeid, "#####")); 
    #print(c(nodeid, split.var.cumulative, split.set.cumulative))
    #print(tnodeid)
    
    ###Splitting
    tmp = splitting(formula=formula, s.var.index=s.var.index, s.var.type=s.var.type, target=target,
                    data=data, method=method, method.point=method.point, impurity.ft=impurity.ft, MINDAT=MINDAT, tnodeid=tnodeid, nodeid=nodeid)
    #print(tmp)
    
    nodeid=tmp$nodeid
    nodes.all.cumulative = c(nodes.all.cumulative, tmp$nodeid) 
    #print(tmp$split.set)
    if(!is.na(tmp$split.set)) {  #if the split set is selected successfully
        
        nodeid.cumulative = c(nodeid.cumulative, tmp$nodeid)
        split.var.cumulative = c(split.var.cumulative, tmp$split.var)
        split.set.cumulative = c(split.set.cumulative, tmp$split.set)
        split.set.elements.cumulative = c.factor(split.set.elements.cumulative, tmp$split.set.elements) ##??? for characters c.factor() rather than c()
        
        tnodeid = tmp$tnodeid
        
        if(length(tnodeid == 2*nodeid) >= MINDAT)
        {  
            nodeid=2*tmp$nodeid
            tmp2 = rsplit(formula=formula, s.var.index=s.var.index,  s.var.type=s.var.type, target=target,
                          data=data, method=method, method.point=method.point, impurity.ft=impurity.ft, MINDAT=MINDAT, 
                          split.var.cumulative=split.var.cumulative, split.set.cumulative=split.set.cumulative, split.set.elements.cumulative=split.set.elements.cumulative,
                          nodes.all.cumulative=nodes.all.cumulative, nodeid.cumulative=nodeid.cumulative, tnodeid=tnodeid, nodeid=nodeid)
            
            tnodeid = tmp2$tnodeid
            nodeid.cumulative = tmp2$nodeid.cumulative
            nodes.all.cumulative = tmp2$nodes.all.cumulative
            split.var.cumulative = tmp2$split.var.cumulative 
            split.set.cumulative = tmp2$split.set.cumulative 
            split.set.elements.cumulative = tmp2$split.set.elements.cumulative 
            
        }
        
        
        if(length(tnodeid == 2*nodeid+1) >= MINDAT)
        {  
            nodeid=2*tmp$nodeid+1
            tmp3 = rsplit(formula=formula, s.var.index=s.var.index, s.var.type=s.var.type, target=target, 
                          data=data, method=method, method.point=method.point, impurity.ft=impurity.ft, MINDAT=MINDAT, 
                          split.var.cumulative=split.var.cumulative, split.set.cumulative=split.set.cumulative, split.set.elements.cumulative=split.set.elements.cumulative,
                          nodes.all.cumulative=nodes.all.cumulative, nodeid.cumulative=nodeid.cumulative, tnodeid=tnodeid, nodeid=nodeid)
            
            
            tnodeid = tmp3$tnodeid
            nodeid.cumulative = tmp3$nodeid.cumulative
            nodes.all.cumulative = tmp3$nodes.all.cumulative
            split.var.cumulative = tmp3$split.var.cumulative
            split.set.cumulative = tmp3$split.set.cumulative
            split.set.elements.cumulative = tmp3$split.set.elements.cumulative 
            
        }
        
    }
    
    
    print("rsplit...ok")
    
    return(list(split.var.cumulative=split.var.cumulative, split.set.cumulative=split.set.cumulative, split.set.elements.cumulative=split.set.elements.cumulative, 
                nodes.all.cumulative=nodes.all.cumulative, nodeid.cumulative=nodeid.cumulative, tnodeid=tnodeid, nodeid=nodeid))
    
}


################################################################
#Split by the selected split rule
splitting <- function(formula, s.var.index, s.var.type, target, data, 
                      method, method.point, impurity.ft, MINDAT, tnodeid, nodeid) {
    
    print("splitting...")
    
    split.var = NA
    split.set = NA
    split.set.elements = NA
    n.samples = nrow(data) #sample size (no missing is assumed!)
    
    ii = which(tnodeid==nodeid)
    #print("N.samples:"); print(length(ii))
    #print(nodeid); print(tnodeid)
    
    ###Split rule selection
    if(method=="es") { #Split rule selection by ES
        tmp = split.rule.es(formula=formula, s.var.index=s.var.index, s.var.type=s.var.type, 
                            target=target, data=data[ii,], impurity.ft=impurity.ft, MINDAT=MINDAT)
    }
    if(method=="ra") { #Split rule selection by RA
        tmp = split.rule.ra(formula=formula, s.var.index=s.var.index, s.var.type=s.var.type, 
                            target=target, data=data[ii,], method.point=method.point, impurity.ft=impurity.ft, MINDAT=MINDAT)
    }
    
    #print("The split variable and point/set are:"); print(tmp$split.var);print(tmp$split.set)
    
    X = data[,s.var.index]             #data with split variables
    if(!is.data.frame(X)) X=data.frame(X)
    
    if(!is.na(tmp$split.set)) {        #if the split set is selected successfully
        
        split.var = tmp$split.var
        n = nrow(X)
        
        if(s.var.type[split.var]=="n") { #for ordered X
            
            split.set = tmp$split.set
            
            ii2 = which(X[,split.var] <= tmp$split.set)
            ii3 = (1:n)[-ii2]
            
            ii2 = intersect(ii, ii2)
            ii3 = intersect(ii, ii3)
            
        } else {                        #for unordered X
            
            
            split.set = tmp$split.set
            
            i = tmp$split.set 
            L = length(tmp$split.set.elements)
            ii.tmp = as.character(binary(i)+10^L)
            aa = substring(ii.tmp, 2:(L+1), 2:(L+1))
            sts = which(aa=="0")
            
            split.set.elements = paste(tmp$split.set.elements[sts], collapse=",")
            
            ii2 = which2(X[, split.var], tmp$split.set.elements[sts])
            ii3 = (1:n)[-ii2]
            
            ii2 = intersect(ii, ii2)
            ii3 = intersect(ii, ii3)
            
        }
        
        tnodeid[ii2] = 2*tnodeid[ii2]    #left node
        tnodeid[ii3] = 2*tnodeid[ii3]+1  #right node 
        
    }
    
    print("splitting...ok")
    
    return(list(nodeid=nodeid, tnodeid=tnodeid, split.var=split.var, split.set=split.set, split.set.elements = split.set.elements))
    
}

################################################################
