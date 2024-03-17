###############################
# Impurity for Regression Trees
###############################

###########################################################################################
# Compute impurity for DRT
impurity.DR <- function(formula, target, data, impurity.ft) {
    ###formula: formula for fitting
    ###data: wide format
    
    #print("imurity.DC...")
    
    #print(dim(data))
    
    
    #print(head(data))
    
    ###compute impurity
    require("betareg")
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
    if(impurity.ft=="nloglik") { # negative log-likelihood
        
        # beta regression error 뜨면 pass
        # beta regression error 뜨면 pass
        result <- tryCatch({
            fit = betareg(new_formula, data)
            impur = -fit$loglik
            return(impur)
        }, warning = function(w) {
            impur = 1e+10
            return(impur)
        })
        impur = result
        
    } else if(impurity.ft=="fitted_sse") { # (y - y_hat)^2
        # beta regression error 뜨거나 범주형 변수의 레벨이 1개면 mean으로 대체
        check = FALSE
        if(check){ # check = TRUE : 범주형 변수의 레벨이 1개
            mu = mean(data[,target])
            impur = sum((data[,target] - mu)^2)
        }
        else{ # optimization error
            result <- tryCatch({
                fit = betareg(new_formula, data)
                impur = sum((data[,target] - fitted(fit))^2)
                return(impur)
            }, warning = function(w) {
                mu = mean(data[,target])
                impur = sum((data[,target] - mu)^2)
                return(impur)
            })
            impur = result
        }
        
    }
    else if(impurity.ft=="linear_sse") { # (y - y_hat)^2
        nvar = ncol(data)-1
        check = FALSE
        for(i in 1:nvar){ # when c-variable has only one level >> no betareg
            if(is.factor(data[,(i+1)]) && length(unique(data[,(i+1)])) == 1){
                check = TRUE
            }
        }
        # beta regression error 뜨거나 범주형 변수의 레벨이 1개면 mean으로 대체
        if(check){ # check = TRUE : 범주형 변수의 레벨이 1개
            mu = mean(data[,target])
            impur = mean((data[,target] - mu)^2)
        }
        else{ # optimization error
            result <- tryCatch({
                fit = lm(formula, data)
                impur = mean((data[,target] - fitted(fit))^2)
                return(impur)
            }, warning = function(w) {
                mu = mean(data[,target])
                impur = mean((data[,target] - mu)^2)
                return(impur)
            })
            impur = result
        }
    }
    else if(impurity.ft=="mean_sse") { # (y - y_bar)^2
        
        mu = mean(data[,target])
        impur = mean((data[,target] - mu)^2)
        
    } else {
        
        stop("The impurity is not available.")
        
    }
    
    #print("impurity.DRT...ok")
    
    return(impur)
}

###########################################################################################

