# x1 has linear relationship with y when x2 is selected properly (X-shaped)
make_data = function(n, y_type, test_size = 100){
    if(y_type == "uniform"){
        y = rbeta(n, 1, 1)
    }
    else if(y_type == "normal"){
        y = rbeta(n, 2, 2)
    }
    else if(y_type == "tube"){
        y = rbeta(n, 0.5, 0.5)
    }
    else if(y_type == 'left skewed'){
        y = rbeta(n, 3, 1)
    }
    else{y = rbeta(n, 1, 3)} # right skewed
    
    ######### set independent variables ###########
    x2 = rnorm(n, 0, 1)
    x3 = rnorm(n, 0, 1)
    x4 = as.factor(sample(c('1','2','3','4','5'), size = n, replace = TRUE))
    x5 = as.factor(sample(c('a','b','c','d','e','f','g','h','i','j'), size = n, replace = TRUE))
    x1 = ifelse(x2 > 0,  2*(y-0.5), -2*(y-0.5)) + rnorm(n, 0, 0.4)
    ######### set independent variables ###########
    
    dat = data.frame(y=y, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5)
    n = nrow(dat)
    test = dat[1:test_size, ]
    train = dat[(test_size+1):n, ]
    return(list(train = train, test = test))
}

# simulation function
betatree_simul = function(formula, alpha = 0.05, MINDAT = 10, nrep, n, y_type, test_size, seed){
    set.seed(seed)
    result = data.frame()
    for(i in 1:nrep){
        print(paste("++++++++", "#reps", i, "++++++++"))
        dat = make_data(n, y_type, test_size)
        train = dat$train; test = dat$test
        bm = betatree(formula, data = train, alpha = alpha, minsize = MINDAT)
        ll = length(bm)-1
        model_output = capture.output(bm)
        split_line <- model_output[8]
        split_parts <- strsplit(split_line, " ")[[1]]
        split_variable <- split_parts[5] 
        split_point <- split_parts[7]   
        f1 = fitted(bm)
        f2 = predict(bm, test)
        result = rbind(result, c(ll, split_variable, split_point, mean((train$y-f1)^2), mean((test$y-f2)^2)))
    }
    colnames(result) = c("#nodes", "split_var", "split_pt", "train_mse", "test_mse")
    result$train_mse = as.numeric(result$train_mse); result$test_mse = as.numeric(result$test_mse)
    result[sapply(result, is.numeric)] <- round(result[sapply(result, is.numeric)], 4)
    return(result)
}



tree_simul = function(formula, f.var, s.var, s.var.type, target, method, method.point,
                      impurity.ft, MINDAT, error.method, nfold = 5, nrep, n, y_type, test_size, seed){
    set.seed(seed)
    result = data.frame()
    for(i in 1:nrep){
        print(paste("++++++++", "#reps", i, "++++++++"))
        dat = make_data(n, y_type, test_size)
        train = dat$train; test = dat$test
        tree = treeft(formula = formula, f.var = f.var, s.var = s.var, s.var.type = s.var.type, target = target,
                      data = train, method=method, method.point=method.point, impurity.ft=impurity.ft,
                      MINDAT=MINDAT, error.method=error.method, nfold=5, test.data=test)
        mse = tree$error.selected.tree.model/((n-test_size)/nfold); n.nodes = nrow(tree$coefs)
        split_var = tree$tree.model[1, 2]; split_pt = tree$tree.model[1, 3]
        p_test = predict.treeft.prob(tree, test[,s.var], test[,f.var], logit = TRUE)
        result = rbind(result, c(n.nodes, split_var, split_pt, mse, mean((test$y-p_test$pred)^2)))
    }
    colnames(result) = c("#nodes", "split_var", "split_pt", "train_mse", "test_mse")
    result$train_mse = as.numeric(result$train_mse); result$test_mse = as.numeric(result$test_mse)
    result[sapply(result, is.numeric)] <- round(result[sapply(result, is.numeric)], 6)
    return(result)
}

bm = betatree_simul(formula = y~x1+x2+x3|1|x1+x2+x3+x4+x5, alpha = 0.05, MINDAT = 20, nrep = 5, n = 300,
                      y_type = "uniform", test_size = 100, seed = 1)

ra = tree_simul(formula = y~x1+x2+x3|1, f.var = c('x1','x2','x3'),
                   s.var = c('x1','x2','x3','x4','x5'), 
                   s.var.type = c('n','n','n','c','c'),
                   target = 'y', method="ra", method.point="es",
                   impurity.ft="nloglik", MINDAT=20, error.method = 'cv', nfold = 5,
                   nrep = 5, n = 300,
                   y_type = "uniform", test_size = 100, seed = 1)

bm;ra

################## RESULT ########################

### comment: True first split variable is "x2" and point is 0. 
###          My model selects variables stably and consistently shows good performance.
### Result of beta regression tree (using betareg package)
#nodes split_var  split_pt train_mse test_mse
#1      4        x1  0.30328:    0.0419   0.0485
#2      4        x1 -0.14077:    0.0333   0.0487
#3      4        x1  -0.00839    0.0432   0.0451
#4      4        x1   0.09382    0.0419   0.0511
#5      2        x2 -0.01303:    0.0287   0.0322
### Result of beta regression tree (using my package)
#nodes split_var              split_pt train_mse test_mse
#1      2        x2    0.0048844495381889  0.032687 0.031596
#2      3        x2  -0.00385754983129943  0.030460 0.031401
#3      5        x2 -0.000745568823387066  0.032558 0.037369
#4      3        x2    0.0912471945496722  0.025026 0.024880
#5      4        x2   -0.0106012589033012  0.031623 0.027791

