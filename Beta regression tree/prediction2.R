##########################################################################################
#Prediction
##########################################################################################


##################################
#Get coefficients for a tree model
get.coefs <- function(obj, data){
    
    #print("get.coefs...")
    
    x = data[,obj$s.var.index]             #data with split variables: wide format
    pred.tnodes = predict.treeft(obj, x)   #predict terminal nodes
    uniq.tnodes = sort(unique(pred.tnodes))
    n.uniq.tnodes = length(uniq.tnodes)
    #print(obj)
    #print(pred.tnodes); print(table(pred.tnodes)); 
    
    require("betareg")
    n.f.var = length(obj$f.var) 
    n.s.var = length(obj$s.var)
    
    result_names = c('(Intercept)')
    for(i in obj$f.var){
        temp = data[,i]
        if(is.factor(temp)){
            result_names = c(result_names, paste0(i, levels(temp)[-1]))
        }
        else{
            result_names = c(result_names, i)
        }
    }
    n.parameters = length(result_names)
    
    ii = is.na(obj$tree.model[,2])
    id.tnodes =  obj$tree.model[ii,1]
    n.tnodes = length(id.tnodes) 
    
    tmp = data.frame(matrix(0, n.tnodes, n.parameters))
    colnames(tmp) = result_names
    rownames(tmp) = id.tnodes
    fit = NA
    for(i in 1:n.tnodes) {
        ii = which(pred.tnodes==id.tnodes[i])
        # if(length(ii) > n.parameters) {
        if(TRUE) {
            nvar = length(obj$f.var)
            check = FALSE
            # for(k in 1:nvar){ # when c-variable has only one level >> no betareg
            #     if(is.factor(data[,obj$f.var[k]]) && length(unique(data[ii,obj$f.var[k]])) == 1){
            #         check = TRUE
            #     }
            # }
            # 원래의 formula
            original_formula <- obj$formula
            
            # Formula에서 응답 변수와 예측 변수 분리
            response_var <- all.vars(original_formula)[1]
            predictor_vars <- all.vars(original_formula)[-1]
            
            # 단일 레벨 factor 변수 제거
            valid_vars <- sapply(data[predictor_vars], function(x) !(is.factor(x) && length(unique(x)) == 1))
            new_predictors <- names(valid_vars)[valid_vars]
            
            # 새로운 formula 생성
            new_formula <- reformulate(new_predictors, response = response_var)
            new_formula <- as.formula(paste(deparse(new_formula), "| 1"))
            
            
            # if(min(data[ii, 1])>0 && max(data[ii, 1])<1 && obj$impurity.ft != 'mean_sse' && !check){
            result <- tryCatch({
                fit = betareg(new_formula, data[ii,])
                fit$coefficients$mean
            }, warning = function(w) {
                mu = mean(data[ii, 1])
                cf = c(mu, rep(0, n.parameters-1))
                names(cf) = result_names
                cf
            })
            # }
            # else{
            #     mu = mean(data[ii, 1])
            #     result = c(mu, rep(0, n.parameters-1))
            #     names(result) = result_names
            # }
            for(j in names(result)){
                tmp[i, j] = result[j]
            }
        }
        
    }
    #print("get.coefs...ok")
    
    return(tmp)
}



###########################################
#Compute probabilities and predicted values
predict.treeft.prob <- function(obj, x, z, logit = TRUE){
    ###data z: long format, x: wide format   
    ###z=for fitting, x=for splitting
    
    #print("predict.treeft.prob...")
    
    if(!is.data.frame(z)) z = data.frame(z)
    if(!is.data.frame(x)) x = data.frame(x)
    
    n = nrow(x) #n = obj$n.samples
    n.f.var = ncol(z) #n.f.var= length(obj$f.var)
    n.tnodes = nrow(obj$coefs) 
    
    pred.tnodes = predict.treeft(obj, x) #predict terminal nodes
    uniq.tnodes = sort(unique(pred.tnodes))
    n.uniq.tnodes = length(uniq.tnodes) #can be != n.tnodes = nrow(obj$coefs)
    
    #print(pred.tnodes); print(uniq.tnodes)
    #print(obj$coefs)
    
    ###Intercepts alpha
    alpha0 = cbind(as.numeric(rownames(obj$coefs)), obj$coefs[,1])
    if(!is.data.frame(alpha0)) alpha0 = data.frame(alpha0) 
    #print(pred.tnodes)
    #print("alpha0="); print(alpha0);
    
    ###Coefficients beta
    beta0 = obj$coefs[,-1] #f.var not exist???
    if(!is.data.frame(beta0)) beta0 = data.frame(beta0) 
    
    #print("beta0="); print(beta0)
    #print(pred.tnodes)
    V = matrix(NA, n, 1)
    resid = matrix(NA, n, 1)
    for(i in 1:n) {
        kk = pred.tnodes[i]
        V[i,1] = alpha0[alpha0[ ,1] == kk, 2]
        idx = 1
        for(k in 1:n.f.var) {
            if(is.factor(z[ ,k])){
                pos = match(z[i,k], levels(z[, k]))
                V[i,1] = ifelse(pos == 1, V[i,1], V[i,1] + (pos != 1) * beta0[alpha0[ ,1] == kk,idx+pos-2])
                idx = idx + length(levels(z[, k])) - 1
            }
            else{
                V[i,1] = V[i,1] + beta0[alpha0[ ,1] == kk,idx] * z[i,k]
                idx = idx + 1
            }
        }
    }
    
    if(logit){
        V = exp(V)/(1+exp(V))
    }
    #print("predict.treeft.prob...ok")
    
    return(list(pred=V))
    
}



####################################
#Predict terminal nodes
predict.treeft <- function(obj, x) {
    ###data x: wide format
    ###x=for splitting
    
    print("predict.treeft...")
    
    #print(dim(x)); print(head(x))
    #print(obj$tree.model)
    
    if(!is.data.frame(x)) x = data.frame(x)
    n.samples = nrow(x) #sample size (no missing is assumed!)
    
    tnodeid = rep(1, n.samples)
    nodeid = 1
    
    if(nrow(obj$tree.model) > 1) { #if not a tree with only the root
        out = rsplit.predict(obj=obj, x=x, tnodeid=tnodeid, nodeid=nodeid)
        tnodeid=out$tnodeid
    }
    
    #print("predict.treeft...ok")
    
    return(tnodeid)
}


#####################################################
#Split Recursively
rsplit.predict <- function(obj, x, tnodeid, nodeid) {
    
    #print("rsplit.predict...")
    
    tmp = splitting.predict(obj=obj, x=x, tnodeid=tnodeid, nodeid=nodeid)
    #print(tmp)
    
    tnodeid = tmp$tnodeid
    nodeid = tmp$nodeid
    
    j = which(obj$tree.model[,1]==nodeid)
    split.var = obj$tree.model[j,2]
    split.set = obj$tree.model[j,3]
    #print(paste("nodeID", nodeid))
    if(!is.na(split.set))  #if intermediate node
    {   #print(paste("nodeID", nodeid))
        nodeid = 2*tmp$nodeid
        tmp2 = rsplit.predict(obj=obj, x=x, tnodeid=tnodeid, nodeid=nodeid)
        tnodeid = tmp2$tnodeid
        
        nodeid = 2*tmp$nodeid+1
        tmp3 = rsplit.predict(obj=obj, x=x, tnodeid=tnodeid, nodeid=nodeid)
        tnodeid = tmp3$tnodeid
        
    }
    
    #print("rsplit.predict...ok")
    
    return(list(tnodeid=tnodeid, nodeid=nodeid))
    
}



########################################################
#Split by the selected split rule
splitting.predict <- function(obj, x, tnodeid, nodeid) {
    
    #print("splitting.predict...")
    #print(dim(x)); print(head(x))
    
    ii = which(tnodeid==nodeid)
    jj = which(obj$tree.model[,1]==nodeid)
    split.var = obj$tree.model[jj,2]
    split.set = obj$tree.model[jj,3]
    if(!is.na(split.set))  #if intermediate node
    {
        kk = which(obj$s.var==split.var)
        #kk = which(names(x)==split.var)
        
        #print(obj$s.var); print(split.var);print(kk)
        if(obj$s.var.type[kk]=="n") { #for ordered variable
            
            split.set = as.numeric(obj$tree.model[jj,3])
            ii2 = which(x[,kk] <= split.set)
            ##ii2 = which(x[,split.var] <= split.set)
            if(length(ii2) == 0){ii3 = 1:nrow(x)}
            else{ii3 = (1:nrow(x))[-ii2]}
            
        } else { #for unordered variable
            split.set = unlist(strsplit(obj$tree.model[jj,3], ","))
            
            #print(x[,kk]); print(split.set)
            
            ii2 = which2(x[,kk], split.set)
            #ii2 = which2(x[, split.var], split.set)
            ii3 = (1:nrow(x))[-ii2]
            
        }
        #print(ii);print(ii2);print(ii3);
        ii2 = intersect(ii, ii2)
        ii3 = intersect(ii, ii3)
        
        tnodeid[ii2] = 2*tnodeid[ii2]    #left node
        tnodeid[ii3] = 2*tnodeid[ii3]+1  #right node 
        
    }
    
    #print("splitting.predict...ok")
    
    return(list(nodeid=nodeid, tnodeid=tnodeid))
    
}

####################################################################END
