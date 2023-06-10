library(glmnet)
library(tidyverse)
library(dplyr)
library("ggrepel")
library(lme4)

# linear rescale to interval [0,1] function, used to rescale real-valued constraints
# formula to map x \in [a,b] to [c,d] : f(x) = c + (d-c/b-a)(x-a)
scaling01 <- function(df, col){
  df[col] = round((df[col]-min(df[col]))/(max(df[col])-min(df[col])), digits=2)
  return (df)
}

# leave one out cross validation function for the glmnet regression algorithm (sum model)
leave_one_out_cv_penality <- function(X, Y, selected_regression) {
  # initializes a matrix to contain cross-validation loss for (L,A) pairs
  cv_matrix <- matrix(ncol=3, nrow=0)
  colnames(cv_matrix) <- c("minimum_cv_error", "lambda", "alpha") 
  # iterates over values of lambda (L)
  for(L in seq(0, 1, by=0.05)){
    #iterates over values of alpha (A)
    for(A in seq(0, 1, by=0.1)){
      # initialize cross-validation error for (L,A)
      cv_error <- 0
      # iterates over K items
      for(k in 1:nrow(X)){
        W = coef(glmnet(X[-k,], Y[-k], lambda = L, alpha = A, upper.limits = 0)) 
        # computes the cross-validation loss
        l = Y[k] - predict(t(X[k,]), W, selected_regression)
        cv_error <- cv_error + l**2
      }
      # updates cross-validation loss for (L,A)
      cv_error <- cv_error / nrow(X)
      # appends cross-validation loss for (L,A) to the cv_matrix 
      cv_matrix <- rbind(cv_matrix, c(cv_error, L, A))
    }
  }
  return(round(cv_matrix,2))
  
}

# leave one out cross validation function for the glmnet regression algorithm (product model)
leave_one_out_cv_penality_product <- function(X, Y, group_constraints,number_group) {
  cv_matrix <- matrix(ncol=3, nrow=0)
  colnames(cv_matrix) <- c("minimum_cv_error", "lambda", "alpha") 
  for(L in seq(0, 1, by=0.05)){
    for(A in seq(0, 1, by=0.1)){
      cv_error <- 0
      for(k in 1:nrow(X)){
        W = coef(glmnet(X[-k,], Y[-k], lambda = L, alpha = A, upper.limits = 0)) 
        l = Y[k] - predict_product(X[k,], W)
        cv_error <- cv_error + l**2
      }
      cv_error <- cv_error / nrow(X)
      cv_matrix <- rbind(cv_matrix, c(cv_error, L, A))
    }
  }
  return(round(cv_matrix,2))
  
}

# leave one out cross validation function for the lmer regression algorithm (sum model)
leave_one_out_cv_linear <- function(X, Y, selected_regression){
  cv_matrix <- matrix(ncol=1, nrow=0)
  colnames(cv_matrix) <- c("minimum_cv_error")
  for (k in 1:nrow(X)){
    cv_error <- 0
    W = coef(lm(Y[-k] ~ X[-k,]))
    W <- replace_na(W,0)
    # computes the cross-validation loss
    l = Y[k] - predict(X[k,], W, selected_regression)
    cv_error <- l**2
    cv_matrix <- rbind(cv_matrix, c(cv_error))
  }
  return(round(cv_matrix,2))
}

# leave one out cross validation function for the lmer regression algorithm (product model)
leave_one_out_cv_linear_product <- function(X, Y, W, group_constraints,number_group){
  cv_matrix <- matrix(ncol=1, nrow=0)
  colnames(cv_matrix) <- c("minimum_cv_error")
  for (k in 1:nrow(X)){
    cv_error <- 0
    W = coef(lm(Y[-k] ~ X[-k,]))
    W <- replace_na(W,0)
    # computes the cross-validation loss
    l = Y[k] - predict_product(X[k,], W)
    cv_error <- l**2
    cv_matrix <- rbind(cv_matrix, c(cv_error))
  }
  return(round(cv_matrix,2))
}

# acceptability prediction function  (sum model)
predict <- function(X, W, selected_regression){
  if (selected_regression == "penality"){
    return(round(W[1] + X %*% W[-c(1)], digits=2))
  }else{
    W <- replace_na(W,0)
    return(round(W[1] + X %*% W[-c(1)], digits=2)) 
  }
}

# acceptability prediction function (product model)
predict_product <- function(df, W){ 
  result <-  W[1]
  W <- W[-1]

  result <- result+df%*%W
  result <- sum(result)
  return(result)
}

# For mode of setting test: acceptability prediction function  (product model)
predict_set_product <- function(df, W, group_constraints,number_group){ 
  df1 <- data.frame(t(df$Constraints.violation))
  colnames(df1) <- c(df$Condition)
  
  # replace na value to 0
  W <- replace_na(W,0)
  X <- c()
  
  if(number_group == 2){
    g1 <- apply(df1[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df1[,group_constraints[[2]]], 2, as.numeric)
    
    for (c1 in g1){
      for (c2 in g2){
          X <- c(X, c1*c2)
      }
    }
    
  }else if(number_group == 3){
    g1 <- apply(df1[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df1[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df1[,group_constraints[[3]]], 2, as.numeric)
    
    for (c1 in g1){
      for (c2 in g2){
        for(c3 in g3){
          X <- c(X, c1*c2*c3)
        }
      }
    }
    
  }else if(number_group == 4){
    
    g1 <- apply(df1[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df1[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df1[,group_constraints[[3]]], 2, as.numeric)
    g4 <- apply(df1[,group_constraints[[4]]], 2, as.numeric)
    
    for (c1 in g1){
      for (c2 in g2){
        for(c3 in g3){
         for(c4 in g4){
           X <- c(X, c1*c2*c3*c4)
         }
        }
      }
    }
    
  }else if(number_group == 5){
    
    g1 <- apply(df1[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df1[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df1[,group_constraints[[3]]], 2, as.numeric)
    g4 <- apply(df1[,group_constraints[[4]]], 2, as.numeric)
    g5 <- apply(df1[,group_constraints[[5]]], 2, as.numeric)
    
    for (c1 in g1){
      for (c2 in g2){
        for(c3 in g3){
          for(c4 in g4){
            for(c5 in g5){
              X <- c(X, c1*c2*c3*c4*c5)
            }
          }
        }
      }
    }
    
  }
  
  result <-  W[1]
  W <- W[-1]
  result <- result+X%*%W
  result <- sum(result)

  return(result)
}


# For mode of uploading test file: acceptability prediction function  (product model)
predict_test_product <- function(df, W, group_constraints,number_group){
  # Remove na value to 0
  W <- replace_na(W,0)
  if(number_group == 2){
    g1 <- apply(df[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df[,group_constraints[[2]]], 2, as.numeric)
    X <- model.matrix(~g1:g2, data=df)
    colnames(X) <- c(colnames(X)[1], gsub("(g1|g2)", "", colnames(X)[-1])) 
  }else if(number_group == 3){
    g1 <- apply(df[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df[,group_constraints[[3]]], 2, as.numeric)
    X <- model.matrix(~g1:g2:g3, data=df)
    colnames(X) <- c(colnames(X)[1], gsub("(g1|g2|g3)", "", colnames(X)[-1])) 
  }else if(number_group == 4){
    g1 <- apply(df[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df[,group_constraints[[3]]], 2, as.numeric)
    g4 <- apply(df[,group_constraints[[4]]], 2, as.numeric)
    X <- model.matrix(~g1:g2:g3:g4, data=df)
    colnames(X) <- c(colnames(X)[1], gsub("(g1|g2|g3|g4)", "", colnames(X)[-1])) 
  }else if(number_group == 5){
    g1 <- apply(df[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df[,group_constraints[[3]]], 2, as.numeric)
    g4 <- apply(df[,group_constraints[[4]]], 2, as.numeric)
    g5 <- apply(df[,group_constraints[[5]]], 2, as.numeric)
    X <- model.matrix(~g1:g2:g3:g4:g5, data=df)
    colnames(X) <- c(colnames(X)[1], gsub("(g1|g2|g3|g4|g5)", "", colnames(X)[-1])) 
  }
  # Remove column (Intercept)
  X <- X[,-1]
  result <-  W[1,]
  W <- W[-1,]
  result <- result+X%*%W
  X <- rowSums(X)
  return(X)
} 

# For mode of setting test: acceptability prediction function  (sum model)
predict_set_sum <- function(X,W){
  list_values =  X$Constraints.violation
  result = 0.00
  W <- replace_na(W,0)
  i = 1 
  while (i <= length(list_values)){
    result = result + list_values[i]*W[i+1]
    i=i+1
  }
  return(result)
}


# Compute weighted product model
main_product <- function(df, group_constraints, list_scales, dep_val, selected_regression, number_group){
  # rescale real-valued constraints
  for(c in list_scales){
    df <- scaling01(df, c)
  }
  
  Y <- apply(df[dep_val], 2, as.numeric)
  
  if(number_group == 2){
    g1 <- apply(df[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df[,group_constraints[[2]]], 2, as.numeric)
    X <- model.matrix(~g1:g2)
    colnames(X) <- c(colnames(X)[1], gsub("(g1|g2)", "", colnames(X)[-1])) 
  }else if(number_group == 3){
    g1 <- apply(df[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df[,group_constraints[[3]]], 2, as.numeric)
    X <- model.matrix(~g1:g2:g3)
    colnames(X) <- c(colnames(X)[1], gsub("(g1|g2|g3)", "", colnames(X)[-1])) 
  }else if(number_group == 4){
    g1 <- apply(df[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df[,group_constraints[[3]]], 2, as.numeric)
    g4 <- apply(df[,group_constraints[[4]]], 2, as.numeric)
    X <- model.matrix(~g1:g2:g3:g4)
    colnames(X) <- c(colnames(X)[1], gsub("(g1|g2|g3|g4)", "", colnames(X)[-1])) 
  }else if(number_group == 5){
    g1 <- apply(df[,group_constraints[[1]]], 2, as.numeric)
    g2 <- apply(df[,group_constraints[[2]]], 2, as.numeric)
    g3 <- apply(df[,group_constraints[[3]]], 2, as.numeric)
    g4 <- apply(df[,group_constraints[[4]]], 2, as.numeric)
    g5 <- apply(df[,group_constraints[[5]]], 2, as.numeric)
    X <- model.matrix(~g1:g2:g3:g4:g5)
    colnames(X) <- c(colnames(X)[1], gsub("(g1|g2|g3|g4|g5)", "", colnames(X)[-1])) 
  }
  # Remove column (Intercept)
  X <- X[,-1]
  print(X)
  # Selon the tye of regression
  if (selected_regression == "penality"){
    cv_matrix <- matrix(ncol=3, nrow=0)
    colnames(cv_matrix) <- c("minimum_cv_error", "lambda", "alpha") 
    cv_matrix <- leave_one_out_cv_penality_product(X, Y, group_constraints,number_group)
    index_min <- which(cv_matrix[,1] == min(cv_matrix[,1]), arr.ind = TRUE)
    cv_error <- cv_matrix[index_min,'minimum_cv_error'][1]
    L <- cv_matrix[index_min, 'lambda'][1]
    A <- cv_matrix[index_min, 'alpha'][1]
    fit = glmnet(X,Y, lambda=L, alpha=A, upper.limits = 0)  
    W = as.matrix(round(coef(fit), digits=2))
    weights = W
    weights[,1] = W[,1]
    weights[1,1] = weights[1,1]
    output <- list(X=X, Y=Y, cv_matrix =  cv_matrix, L = L, A = A, W = W, cv_error = cv_error, weights=weights)
    
    print(W)
  } else {
    cv_matrix <- leave_one_out_cv_linear_product(X, Y, W, group_constraints,number_group)
    
    index_min <- which(cv_matrix[,1] == min(cv_matrix[,1]), arr.ind = TRUE)
    cv_error <- cv_matrix[index_min,'minimum_cv_error'][1]
    fit = lm(Y~X)
    W = as.matrix(round(coef(fit), digits=2))
    rownames(W) <- c(rownames(W)[1], gsub("X", "", rownames(W)[-1])) 
    weights = W
    weights[,1] = W[,1]
    weights[1,1] = weights[1,1]
    output <- list(X=X, Y=Y, cv_matrix =  cv_matrix, W = W, weights=weights)
  }
  return(output)
}

# Compute weighted sum model
main <- function(df, list_constraints, list_scales, dep_val, selected_regression){
  # rescale real-valued constraints
  for(c in list_scales){
    df <- scaling01(df, c)
  }
  
  X <- apply(df[,list_constraints], 2, as.numeric)
  Y <- apply(df[dep_val], 2, as.numeric)

  # Selon the tye of regression
  if (selected_regression == "penality"){
    cv_matrix <- leave_one_out_cv_penality(X,Y, selected_regression)
    index_min <- which(cv_matrix[,1] == min(cv_matrix[,1]), arr.ind = TRUE)
    cv_error <- cv_matrix[index_min,'minimum_cv_error'][1]
    L <- cv_matrix[index_min, 'lambda'][1]
    A <- cv_matrix[index_min, 'alpha'][1]
    fit = glmnet(X,Y, lambda=L, alpha=A, upper.limits = 0)  
    W = as.matrix(round(coef(fit), digits=2))
    weights = W
    weights[,1] = W[,1]
    weights[1,1] = weights[1,1]
    output <- list(X=X, Y=Y, cv_matrix =  cv_matrix, L = L, A = A, W = W, cv_error = cv_error, weights=weights)
  } else {
    cv_matrix <- leave_one_out_cv_linear(X,Y, selected_regression)
    index_min <- which(cv_matrix[,1] == min(cv_matrix[,1]), arr.ind = TRUE)
    cv_error <- cv_matrix[index_min,'minimum_cv_error'][1]
    fit = lm(Y~X)
    W = as.matrix(round(coef(fit), digits=2))
    rownames(W) <- c(rownames(W)[1], gsub("X", "", rownames(W)[-1])) 
    weights = W
    weights[,1] = W[,1]
    weights[1,1] = weights[1,1]
    output <- list(X=X, Y=Y, cv_matrix =  cv_matrix, W = W, weights=weights)
  }
  return (output)
}

# for plotting cv error 
plot_cv_error <- function(cv_matrix, selected_regression){
  # compute cv error selon the type of regression
  if (selected_regression == "penality"){
    list_mse_m <- c()
    list_stdev_mse <- c()
    list_se_mse <- c()
    list_upper_mse <- c()
    list_lower_mse <- c()
    
    result <- as.data.frame(cv_matrix)
    cv_error <- as.numeric(unlist(result["minimum_cv_error"]))
    lambda <- as.numeric(unlist(unique(result["lambda"])))
    for (i in lambda){
      error <- result[which(result$lambda == i), "minimum_cv_error"]
      mse_m = mean(error)
      stdev_mse= sd(error)
      se_mse = stdev_mse/sqrt(length(error))
      upper_mse= mse_m+se_mse*1.96
      lower_mse=mse_m-se_mse*1.96
      
      
      list_mse_m <- c(list_mse_m, mse_m)
      list_stdev_mse <- c(list_stdev_mse, stdev_mse)
      list_se_mse <- c(list_se_mse, se_mse)
      list_upper_mse <- c(list_upper_mse, upper_mse)
      list_lower_mse <- c(list_lower_mse, lower_mse)
    }
    
    result_error = data.frame('mse'=list_mse_m,'mse_upper'=list_upper_mse,'mse_lower'=list_lower_mse, 'lambda'=lambda)
    
    index_min <- which(cv_matrix[,1] == min(cv_matrix[,1]), arr.ind = TRUE)
    L <- cv_matrix[index_min, 'lambda'][1]
    
    return(
      ggplot(data = result_error, aes(x=log10(lambda), y=list_mse_m))+ 
        geom_point()+
        geom_errorbar(aes(ymin=list_lower_mse, ymax=list_upper_mse),  width=.2)+
        geom_line(color = 'red')+
        geom_vline(xintercept = log10(L), linetype="dotted", color = "gold4", size=1.5)+
        xlab (expression(log10(lambda)))+
        ylab ('Mean squared error')+
        theme_classic() +
        geom_tile()
    )
    
  } else {
    result <- as.data.frame(cv_matrix)
    cv_error <- as.numeric(unlist(result["minimum_cv_error"]))
    mse_m = mean(cv_error)
    stdev_mse= sd(cv_error)
    se_mse = stdev_mse/sqrt(length(cv_error))
    upper_mse= mse_m+se_mse*1.96
    lower_mse=mse_m-se_mse*1.96
    result_error = data.frame('mse'=mse_m,'mse_upper'=upper_mse,'mse_lower'=lower_mse)
    
    return(
      ggplot(data=result_error,aes(x="1", y=mse_m))+ 
        geom_bar(stat="identity", fill='steelblue')+
        geom_text(aes(label=mse_m), vjust=1.6, color="white", size=3.5)+
        geom_errorbar(aes(ymin=lower_mse, ymax=upper_mse), width=.2)+
        xlab ("Linear model")+
        ylab ('Mean squared error')+
        theme_classic()
    )
  }
  
}

# for plotting weights_histogram
plot_weights_histogram <- function(W){
  list_constraints = row.names(W)
  d <- data.frame(constraints= list_constraints[-1], weights=as.matrix(W[-c(1)]))
  print(W)
  return(ggplot(data=d, aes(x=constraints, y=weights, group=constraints,fill=constraints)) + 
           geom_bar(stat="identity") + guides(fill = FALSE) + 
           geom_text(aes(label=weights), vjust=1.6, color="black", size=3.5)+
           theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = rel(1.3)), axis.text.y = element_text(size = rel(1.3)), axis.title = element_text(size=13)))
}

# for plotting the result of prediction without dependent values (product model)
plot_predictions_test_no_gold_product <- function(df_test, group_constraints, W, list_scales, number_group){
  for(c in list_scales){
    df_test <- scaling01(df_test, c)
  }
  
  df_test['pred'] = predict_test_product(df_test, W, group_constraints, number_group)

  return(ggplot(df_test, aes(x=exp,y=pred, label=exp))+geom_bar(stat="identity") + guides(fill = FALSE)
         + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = rel(1.3)), axis.text.y = element_text(size = rel(1.3)),axis.title = element_text(size=13)))
  
}

# for plotting the result of prediction with dependent values (product model)
plot_predictions_test_product <- function(df_test, group_constraints, W, dep_val, list_scales, number_group){
  # rescale real-valued constraints
  for(c in list_scales){
   df_test <- scaling01(df_test, c)
   }

  Y_test <- apply(df_test[dep_val], 2, as.numeric)
  df_test['pred'] = predict_test_product(df_test, W, group_constraints, number_group)

  result <- data.frame(Y_test, "pred" = df_test['pred'])
  return(ggplot(result, aes(x=Y_test, y=pred, label=pred)) + geom_line()
           + theme(axis.title = element_text(size = 13), axis.text = element_text(size=13)) 
           + geom_text(size=4, hjust=-0.2, vjust=0) + guides(fill = FALSE))
}

# for plotting the result of prediction without dependent values (sum model)
plot_predictions_test_no_gold <- function(df_test, list_constraints, W, selected_regression, list_scales){
  # rescale real-valued constraints
  for(c in list_scales){
    df_test <- scaling01(df_test, c)
  }
  
  X_test <- apply(df_test[,list_constraints], 2, as.numeric)
  df_test['pred'] = predict(X_test, W, selected_regression)
  result <- data.frame("pred" = df_test['pred'], exp=df_test['exp'])
  return(ggplot(result, aes(x=exp,y=pred, label=pred))+geom_bar(stat="identity") + guides(fill = FALSE)
           + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = rel(1.3)), axis.text.y = element_text(size = rel(1.3)),axis.title = element_text(size=13)))
}

# computes predictions of the model during cv with (L,A) minimizing cv-error
# input: X, Y, W of type matrix; L, A of type FLOAT, upper_limit_zero of type BOOL
# output: Y_pred of type matrix
compute_predictions_loocv <- function(X, Y, W, L, A, upper_limit=TRUE) {
  # initializes prediction vector
  Y_pred = vector()
  # iterates over K items
  for(k in 1:nrow(X)){
    # fits without leaving out item k
    fit = glmnet(X[-k,], Y[-k], lambda=L, alpha=A, upper.limits = 0)
    # regression weights
    W = round(coef(fit), digits=2)
    # computes predictions 
    p = round(W[1] + X[k,] %*% W[-c(1)], digits=2)
    # apprends predition to prediction vector
    Y_pred <- c(Y_pred, p)  
  }
  return(Y_pred)
}

plot_predictions_test <- function(df_test, list_constraints, W, dep_val, selected_regression, list_scales, L, A){
  print(L)
  print(typeof(print(L)))
  # rescale real-valued constraints
  for(c in list_scales){
    df_test <- scaling01(df_test, c)
  }

  X_test <- apply(df_test[,list_constraints], 2, as.numeric)
  Y_test <- apply(df_test[dep_val], 2, as.numeric)
  predictions_loocv = compute_predictions_loocv(X_test, Y_test, W, L, A)
  # appends predictions vector to df to plot and export
  df_test['pred'] <- predictions_loocv

  write.csv(x = df_test,file = "pred.csv")
  
  ggplot(df_test, aes(x=mean, y=pred, label=exp, colour=extraction, shape=extracted_function)) + geom_point() + geom_abline(a=0, b=1, col = "black", lty="dashed") + xlim(0, 7) + ylim(0, 7) + geom_text(size=2.5, hjust=-0.2, vjust=0)
  ggsave(filename = "test2.png", 
         width = 4.5,
         height = 6, 
         units = "in",
         dpi = 700)
  
  return(ggplot(df_test, aes(x=mean, y=pred, label=exp, colour=extraction, shape=extracted_function)) + geom_point() + geom_abline(a=0, b=1, col = "black", lty="dashed") + xlim(0, 7) + ylim(0, 7) + geom_text(size=2.5, hjust=-0.2, vjust=0))
}

# # for plotting the result of prediction with dependent values (sum model)
# plot_predictions_test <- function(df_test, list_constraints, W, dep_val, selected_regression, list_scales){
#   # rescale real-valued constraints
#   for(c in list_scales){
#     df_test <- scaling01(df_test, c)
#   }
# 
#   X_test <- apply(df_test[,list_constraints], 2, as.numeric)
#   Y_test <- apply(df_test[dep_val], 2, as.numeric)
#   df_test['pred'] = predict(X_test, W, selected_regression)
#   result <- data.frame(Y_test, "pred" = df_test['pred'])
#   # compute error
#   test_error = mean(abs(result[,1] - result[,2]))
#   print(result)
#   print(Y_test)
#   print(test_error)
# 
#   return(ggplot(result, aes(x=Y_test, y=pred, label=pred)) + geom_line()
#          + theme(axis.title = element_text(size = 13), axis.text = element_text(size=13))
#          + geom_text(size=4, hjust=-0.2, vjust=0) + guides(fill = FALSE))
# }


