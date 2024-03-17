#########################################################
# Split Rule Selection by Residual Analysis (RA) for DC
#########################################################
split.rule.ra <- function(formula, s.var.index, s.var.type, target, data, method.point, impurity.ft, MINDAT) {
    
    print("split.rule.ra...")
    #print(dim(data))
    
    
    ###Split variable selection
    split.var = 1
    if(length(s.var.index) > 1) {
        tmp = split.var.selection.ra(formula, s.var.index, s.var.type, target, data, impurity.ft, MINDAT)
        split.var = tmp$split.var #split variable selected 
    }
    #print("Split variable is... "); print(split.var)
    
    
    ###Split point/set selection for the selected split variable
    X = data[,s.var.index] #data with split variables
    if(!is.data.frame(X)) X=data.frame(X)
    
    if(s.var.type[split.var]=="n")  {   #X:ordered
        
        if(method.point=="median") {   #shortcut
            tmp = split.shortcut.orderedX(formula=formula, target=target, 
                                          x=X[,split.var], data=data, MINDAT=MINDAT, impurity.ft=impurity.ft)
        } else {                       #ES
            tmp = split.es.orderedX(formula=formula, target=target, 
                                    x=X[,split.var], data=data, MINDAT=MINDAT, impurity.ft=impurity.ft)
        }
        
        
    } else {                           #X:unordered
        
        #if(method.set=="shortcut") {  #shortcut
        #   tmp = split.shortcut.unorderedX(formula=formula, target=target, 
        #                                   x=X[,split.var], data=data, MINDAT=MINDAT, impurity.ft=impurity.ft)
        #} else {                      #ES
        tmp = split.es.unorderedX(formula=formula, target=target, 
                                  x=X[,split.var], data=data, MINDAT=MINDAT, impurity.ft=impurity.ft)
        #}   
        
    }
    
    delta = tmp$max.delta
    split.set = tmp$pt #split point/set selected 
    
    split.set.elements = NA    
    if(s.var.type[split.var]=="c") split.set.elements = sort(unique(X[,split.var])) 
    
    
    #print("Split point/set is... "); print(split.set); print(split.set.elements)
    
    print("split.rule.ra...ok")
    
    return(list(split.var=split.var, split.set=split.set, split.set.elements=split.set.elements, delta=delta, impurity.ft=impurity.ft))
}


###################################################
#Split variable selection by residual analysis (RA)
split.var.selection.ra <- function(formula, s.var.index, s.var.type, target, data, impurity.ft, MINDAT) {
    
    print("split.var.selection.ra...")
    
    #print(dim(data)); print(head(data))
    
    
    X = data[,s.var.index] #data with split variables
    if(!is.data.frame(X)) X=data.frame(X)
    
    ###Fit a model to the data (fitted or mean)
    require("betareg")
    nvar = ncol(data)-1
    if(impurity.ft == 'mean_sse'){
        y_hat = mean(data[,target])
    }
    else{
        # 원래의 formula
        original_formula <- formula
        
        # Formula에서 응답 변수와 예측 변수 분리
        response_var <- all.vars(original_formula)[1]
        predictor_vars <- all.vars(original_formula)[-1]
        
        # 단일 레벨 factor 변수 제거
        valid_vars <- sapply(data[predictor_vars], function(x) !(is.factor(x) && length(unique(x)) == 1))
        new_predictors <- names(valid_vars)[valid_vars]
        
        # 새로운 formula 생성
        new_formula <- reformulate(new_predictors, response = response_var)
        new_formula <- as.formula(paste(deparse(new_formula), "| 1"))
        check = FALSE
        
        # beta regression error 뜨거나 범주형 변수의 레벨이 1개면 mean으로 대체
        if(check){ # check = TRUE : 범주형 변수의 레벨이 1개
            mu = mean(data[,target])
            impur = mean((data[,target] - mu)^2)
        }
        else{ # optimization error
            result <- tryCatch({ # betareg 오류시 mean 사용
                fit = betareg(new_formula, data)
                y_hat = fitted(fit)
                y_hat
            }, warning = function(w) {
                y_hat = mean(data[,target])
                y_hat
            })
            impur = result
        }
        y_hat = result
    }
    
    ###Obtain the residuals
    r = data[,target] - y_hat
    r.sum = ifelse(r > 0, 1, -1)
    
    
    # stop("OKKKKKK")
    
    ###Compute uncertainty coefficient (entropy)
    split.var = 1
    assoc = 1
    p_val = c()
    for(j in 1:ncol(X)) {
        
        x = X[,j]
        if(s.var.type[j]=="n") { #Cut x at Q1,Q2,Q3 for ordered X
            x = ifelse(x <= median(x),
                       ifelse(x <= quantile(x)[2], 1, 2),
                       ifelse(x <= quantile(x)[4], 3, 4)) 
            tmp=chisq.test(table(x, r.sum))$p.value #Chisq test statistic
        }
        else{ # for c-var
            tab = table(x, r.sum)
            if(ncol(tab) == 1 || nrow(tab) == 1){
                tmp = 1
            }
            else{
                new_tab <- as.table(tab[rowSums(tab) != 0, , drop = FALSE])
                if(ncol(new_tab) == 1 || nrow(new_tab) == 1){tmp = 1}
                else{tmp=chisq.test(new_tab)$p.value}
            }
        }
        p_val = c(p_val, tmp)
        #print(table(x, r.sum))
        #print("statistic="); print(tmp)
        if(tmp <= assoc) {
            assoc = tmp #the largest chi-statistic = smallest p-value
            split.var = j #select the split variable with the largest statistic
        }
        
    }
    # Interaction term
    comb = combn(ncol(X), 2)
    for(i in 1:ncol(comb)){
        idx1 = comb[1, i]; idx2 = comb[2, i]
        #print(paste("idx:", idx1, idx2))
        d1 = X[ ,idx1]
        d2 = X[ ,idx2]
        # n-var & n-var
        if(is.numeric(d1) & is.numeric(d2)){
            tt = ifelse(d1 <= median(d1),
                        ifelse(d2 <= median(d2), 1, 2),
                        ifelse(d2 <= median(d2), 3, 4))
            tmp = chisq.test(table(tt, r.sum))$p.value
        }
        # c-var & c-var
        else if(!is.numeric(d1) & !is.numeric(d2)){
            itr = interaction(d1, d2)
            tab = table(itr, r.sum)
            if(ncol(tab) == 1 || nrow(tab) == 1){
                tmp = 1
            }
            else{
                new_tab <- as.table(tab[rowSums(tab) != 0, , drop = FALSE])
                #print(chisq.test(new_tab))
                if(ncol(new_tab) == 1 || nrow(new_tab) == 1){tmp = 1}
                else{tmp=chisq.test(new_tab)$p.value}
            }
        }
        # n-var vs c-var
        else{
            if(is.numeric(d1)){d1 = ifelse(d1 <= median(d1), 1, 2)}
            else{d2 = ifelse(d2 <= median(d1), 1, 2)}
            itr = interaction(d1, d2)
            tab = table(itr, r.sum)
            if(ncol(tab) == 1 || nrow(tab) == 1){
                tmp = 1
            }
            else{
                new_tab <- as.table(tab[rowSums(tab) != 0, , drop = FALSE])
                #print(chisq.test(new_tab))
                if(ncol(new_tab) == 1 || nrow(new_tab) == 1){tmp = 1}
                else{tmp=chisq.test(new_tab)$p.value}
            }
        }
        #print(tmp); print(assoc)
        if(tmp <= assoc) {
            assoc = tmp #the largest chi-statistic
            # n-var & n-var
            if(is.numeric(d1) & is.numeric(d2)){
                m1 = mean(d1); m2 = mean(d2)
                ii1 = which(d1 <= m1); ii2 = which(d2 <= m2)
                impur1 = impurity.DR(formula, target, data[ii1,], impurity.ft)+ impurity.DR(formula, target, data[-ii1,], impurity.ft)
                impur2 = impurity.DR(formula, target, data[ii2,], impurity.ft)+ impurity.DR(formula, target, data[-ii2,], impurity.ft)
                #print(impur1); print(impur2)
                split.var = ifelse(impur1 < impur2, idx1, idx2)
            }
            # c-var & c-var >> choose smallest p-value
            else if(is.numeric(d1) & is.numeric(d2)){
                split.var = ifelse(p_val[idx1] < p_val[idx2], idx1, idx2)
            }
            # n-var vs c-var >> choose c-var
            else{
                if(is.numeric(d1)){split.var = idx2}
                else{split.var = idx1}
            }
        }
    }
    
    print("split.var.selection.ra...ok")
    return(list(split.var=split.var))
}



####################################################
# split the node along the sample mean of each variable to choose the variable whose split yields the smaller SSE
split.es.orderedX2 <- function(formula, target, x, data, MINDAT, impurity.ft) {
    
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




##################################################################
# #Summary measure of association: Uncertainty coefficient: Entropy
# entropy = function(x,y){
#     
#     pxy = table(x,y)/sum(table(x,y))
#     px = table(x)/sum(table(x))
#     py= table(y)/sum(table(y))
#     
#     tmpy = 0 
#     for(j in 1:length(py)) {
#         tmpy = tmpy + py[j] * log(py[j])
#     }
#     
#     tmp = 0 
#     for(i in 1:length(px)){
#         for(j in 1:length(py)){
#             
#             if(pxy[i,j] > 0) {
#                 tmp = tmp + pxy[i,j] * log(pxy[i,j]/(px[i]*py[j]))
#             }
#             
#         }}
#     
#     ucoef = -1*tmp/tmpy
#     
#     return(ucoef)
# }


####################################################
#Selecting a split point by median for ordered variable X
split.shortcut.orderedX <- function(formula, choice, n.classes, varying, x, data, MINDAT, impurity.ft) {
    
    #print("split.shortcut.orderedX...")
    
    
    pt = median(x) #median of x as the split point
    max.delta=NA
    
    if(nrow(data) >= MINDAT) { 
        
        impur = impurity.DC(formula, choice, n.classes, varying, data, impurity.ft)
        
        ii = which(x <= pt)
        
        pL=1; pR=1
        if(impurity.ft == "rsquare") {pL = length(ii)/length(x); pR=1-pL}
        if(impurity.ft == "misclass") {pL = length(ii)/length(x); pR=1-pL}
        
        if((nrow(data[ii,]) >= MINDAT) & (nrow(data[-ii,]) >= MINDAT)) { 
            impur.tL = impurity.DC(formula, choice, n.classes, varying, data[ii,],  impurity.ft)
            impur.tR = impurity.DC(formula, choice, n.classes, varying, data[-ii,], impurity.ft)
            max.delta = impur - (pL*impur.tL + pR*impur.tR)
        }  
    }
    
    if(is.na(max.delta)) pt = NA
    
    #print("split.shortcut.orderedX...ok")
    
    return(list(pt=pt, max.delta=max.delta))
    
}


###################################################################################

