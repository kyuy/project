########################################################
# Split Rule Selection by Exhaustive Search (ES) for DRT
########################################################
split.rule.es <- function(formula, s.var.index, s.var.type, target, data, impurity.ft, MINDAT) {
    
    print("split.rule.es...")
    
    X = data[,s.var.index]; if(!is.data.frame(X)) X=data.frame(X) #data with split variables
    nc = ncol(X)
    
    split.var = NA 
    split.set = NA 
    delta = NA 
    split.set.elements = NA
    
    tmp.set = c()
    tmp.delta = c()
    
    for(j in 1:nc)
    {
        if(s.var.type[j]=="n")  {   #X:ordered
            
            tmp = split.es.orderedX(formula=formula, target=target,
                                    x=X[,j], data=data, MINDAT=MINDAT, impurity.ft=impurity.ft)
            tmp.set = c(tmp.set, tmp$pt)
            tmp.delta = c(tmp.delta, tmp$max.delta)
            
        } else {                   #X:unordered
            
            tmp = split.es.unorderedX(formula=formula, target=target,
                                      x=X[,j], data=data, MINDAT=MINDAT, impurity.ft=impurity.ft)
            tmp.set = c(tmp.set, tmp$pt)
            tmp.delta = c(tmp.delta, tmp$max.delta)
            
        }
        
    }
    
    #print(tmp.delta)
    #print(tmp.set)
    
    ##########
    #Choose the split var with a split point/set
    jj <- which.max(tmp.delta)
    if(length(jj) >0) {
        delta = tmp.delta[jj] 
        split.var = jj 
        split.set = tmp.set[jj] 
        if(s.var.type[jj]=="c") split.set.elements = sort(unique(X[,jj])) #since some categories may be missing at this node
    }
    
    #print(split.var);print(split.set)
    
    print("split.rule.es...ok")
    
    return(list(split.var=split.var, split.set=split.set, split.set.elements=split.set.elements, delta=delta, impurity.ft=impurity.ft))
}


####################################################
#Selecting a split point by ES for ordered variable X
split.es.orderedX <- function(formula, target, x, data, MINDAT, impurity.ft) {
    
    print("split.es.orderedX...")
    
    
    pt = NA
    pts = sort(unique(x))
    n = length(pts)
    delta = c(rep(NA, n))
    max.delta=NA
    
    if(nrow(data) >= MINDAT) { 
        
        
        impur = impurity.DR(formula, target, data, impurity.ft)
        
        for(i in 1:n) {
            
            ii = which(x <= pts[i]) 
            
            if((nrow(data[ii,]) >= MINDAT) & (nrow(data[-ii,]) >= MINDAT)) { 
                impur.tL = impurity.DR(formula, target, data[ii,],  impurity.ft)
                impur.tR = impurity.DR(formula, target, data[-ii,], impurity.ft)
                delta[i] = impur - (impur.tL + impur.tR)
                
                #print(c(delta[i],impur, impur.tL,impur.tR, pts[i]))
                
            }
        }
        
    }
    
    
    ii = which.max(delta)
    if(length(ii) >0) {
        pt = pts[ii]
        max.delta = delta[ii]
    }
    
    #print("max.delta, pt"); print(c(max.delta, pt))
    
    print("split.es.orderedX...ok")
    
    
    return(list(pt=pt, max.delta=max.delta))
}


#####################################################
#Selecting a split set by ES for unordered variable X
split.es.unorderedX <- function(formula, target, x, data, MINDAT, impurity.ft) {
    
    print("split.es.unorderedX...")
    
    pt.best = NA
    set.best = NA
    set.elements = sort(unique(x))
    #print(x)
    #print(set.elements)
    L = length(set.elements)
    LL = 2^(L-1)-1
    max.delta = NA
    #print(set.elements)
    #print(L)
    #print(LL)
    if(nrow(data) >= MINDAT) { 
        impur = impurity.DR(formula, target, data, impurity.ft)
        
        for(i in 1:LL) {
            
            ii = as.character(binary(i)+10^L)
            aa = substring(ii, 2:(L+1), 2:(L+1))
            sts = which(aa=="0")
            ii = which2(x, set.elements[sts])
            
            if((nrow(data[ii,]) >= MINDAT) & (nrow(data[-ii,]) >= MINDAT)) { 
                impur.tL = impurity.DR(formula, target, data[ii,], impurity.ft)
                impur.tR = impurity.DR(formula, target, data[-ii,], impurity.ft)
                delta = impur - (impur.tL + impur.tR)
                
                if(is.na(max.delta)) {
                    set.best = set.elements[sts]
                    pt.best = i
                    max.delta = delta
                } else {
                    if(delta >= max.delta) {
                        set.best = set.elements[sts]
                        pt.best = i
                        max.delta = delta
                    }
                }
                
                
            }
        }
        
    }
    
    #print("max.delta, pt.best, set.best");print(c(max.delta, pt.best)); print(set.best)
    
    print("split.es.unorderedX...ok")
    
    return(list(pt=pt.best, set=set.best, set.elements=set.elements, max.delta=max.delta))
}

###########################################################

