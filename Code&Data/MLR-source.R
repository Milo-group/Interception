# The following packages are installed and loaded 
# suppressMessages(install.packages(c('data.table','tidyr','ggplot2','scales',
#                    'reshape2','tibble','caret',
#                    'plyr','dplyr')))

library(data.table)
library(tidyr)
library(ggplot2)
library(scales)
library(reshape2)
library(tibble)
library(caret)
library(plyr)
library(dplyr)

data <- data.table::fread('pub_data.csv', header = T) # Import data

data <- data.frame(data)
row.names(data) <- data[,2]
data$position <- as.factor(data$position) # Changing variable class to fit
data$class <- as.factor(data$class)
data[,c(5:14)] <- as.numeric(as.matrix(data[,c(5:14)]))

data[,c(5:14)] <- data.frame(scale(data[,c(5:14)], T, T)) # Scaling of continuous variables 

kfold.data <- data[,-c(1:2,7,9)] # Remove name, O and Cring charge columns

# k.fold.multinom description:
# 
# Takes a formula, dataset, number of folds and the output vector. 
# Uses nnet::multinom (see nnet v7.3-14 documentation - multinom).
# Divides data into folds, each time using a single fold as the test set
# and the rest as training. 
# Creates a predicted probalities and classes dataframe for each run 
# which is then combined to form a complete table.
# Returns a classifictaion table for predicted outcome vs. experimental observations,
# an accuracy measure based on that table and the probabilities table. 

k.fold.multinom <- function(formula, data, folds, outcome.column) {
  models <- list()
  probalities <- list()
  acc <- list()
  class.pred <- list()
  split.assign <- sample(1:folds, nrow(data), replace = TRUE)
  new_dat <- cbind(data, split.assign)
  for (i in 1:folds) {
    train <- new_dat[split.assign != i,]
    test <- new_dat[split.assign == i,]
    models[[match(i,1:folds)]] <- nnet::multinom(formula = formula,
                                                 data = train,
                                                 maxit = 2000,
                                                 trace = FALSE)
    probalities[[match(i,1:folds)]] <- data.frame(predict(models[[match(i,1:folds)]], newdata = test , "probs")*100)
    class.pred[[match(i,1:folds)]] <- data.frame(predict(models[[match(i,1:folds)]], newdata = test , "class"))
    probalities[[match(i,1:folds)]] <- cbind(probalities[[match(i,1:folds)]],class.pred[[match(i,1:folds)]],test$flag)
  }
  probs <- data.frame(do.call(rbind, probalities))
  probs <- probs[order(probs$test.flag),]
  colnames(probs) <- c('cond1','cond2','cond3','prediction')
  probs[,1:3] <- round(probs[,1:3],digits = 0)
  pred <- probs[,4]
  pred <- plyr::mapvalues(pred, 
                    from=c("1","2","3"), 
                    to=c("SIPr&BA","SIPr",'TAC'))
  actual <- kfold.data[[outcome.column]]
  actual <- plyr::mapvalues(actual, 
                      from=c("1","2","3"), 
                      to=c("SIPr&BA","SIPr",'TAC'))
  ct <- table(actual, pred)
  acc <- round((sum(diag(ct))/sum(ct))*100,2)
  return(list(acc, ct, probs))
}

# kf.iter description:
# 
# Takes a formula, dataset, number of folds and the output vector along with 
# number of iterations. 
# Runs k.fold.multinom by iterations numbers and returns the average 
# accuracy over all iterations, the best and worst accuracy and 
# best and worst classifictaion tables.

kf.iter <- function(formula, data, folds, out.col, iterations) {
  iter.list <- list()
  ct.list <- list()
  for (i in 1:iterations) {
    mod <- k.fold.multinom(formula, data, folds, out.col)
    iter.list[[match(i,1:iterations)]] <- mod[[1]]
    ct.list[[match(i,1:iterations)]] <- mod[[2]]
  }
  over.all.accuracy <- round(Reduce(`+`,iter.list)/iterations,digits = 2)
  best <- iter.list[which.max(iter.list)]
  worst <- iter.list[which.min(iter.list)]
  Accuracies <- knitr::kable(cbind(over.all.accuracy, best,worst))
  tab <- cbind(dcast(data.frame(ct.list[which.max(iter.list)]),actual~pred),
               rep("***",3),
               dcast(data.frame(ct.list[which.min(iter.list)]),actual~pred))
  names(tab)[5] <- ''
  cts <- knitr::kable(tab,
                      caption = "\n\nBest (left) and Worst (right) Classification Tables")
  print(Accuracies)
  print(cts)
  
}

# mod.info description:
# 
# Takes a nnet::multinom model object.
# using the full data set as both training and testing sets, model information is extracted.
# Creates a class.table object to be used in confusion matrix plotting,
# computes and returnes varibale importance as the sum of each coefficient's absolute value
# over all predictied classes (See caret v6.0-86 documentation - varImp),
# accuracy based on classification table, 
# coeffcient's exponent as another measure of importance 
# and McFadden's pseudo R squared for multinomial logistic regression (T. Domencich and D. L. McFadden, 1975).

mod.info <- function(model) {
  pred <- predict(test,newdata = data, 'class')
  pred <- plyr::mapvalues(pred, 
                    from=c("1","2","3"), 
                    to=c("SIPr&BA","SIPr",'TAC'))
  actual <- data$class
  actual <- plyr::mapvalues(actual, 
                      from=c("1","2","3"), 
                      to=c("SIPr&BA","SIPr",'TAC'))
  class.table <<- table(actual, pred)
  
  Accuracy <- paste(round((sum(diag(class.table))/sum(class.table))*100,2),"%",sep = '')
  test.0 <- nnet::multinom(class ~ 1,data = data, maxit = 2000, trace=F)
  test.1 <- test
  
  loglik.0 <- as.numeric(nnet:::logLik.multinom(test.0))
  loglik.1 <- as.numeric(nnet:::logLik.multinom(test.1))
  
  McFadden_R2 <- round((1 - (loglik.1/loglik.0)),digits = 3)
  st <- cbind(Accuracy,McFadden_R2)
  
  stats <- knitr::kable(st)
  var_imp <- knitr::kable(varImp(test),caption = "\n\nVariable Importance")
  ce <- knitr::kable(coef(test),caption ="\n\nCoefficients")
  CE <- knitr::kable(exp(coef(test)),caption = "\n\nexp(coefficients)")
  print(stats)
  
  print(var_imp)
  
  print(ce)
  
  print(CE)
  
}


# ct.plot description:
# 
# Takes a class table and computes;
# * Recall - TP/(TP+FN) & complimentary percentages of wrong prdicitions in each cell.
# * Precision - TP/(TP+FP) for each prediction class.
# * Accuracy - all.TP/total.
# Classifies into Right and Wrong predictions and presents class size.

ct.plot <- function(class.table) {
  ct <- as.matrix(class.table)
  total <- rep(NA, nrow(ct))
  ct <- cbind(ct,total)
  for (i in 1:dim(ct)[1]) {
    ct[i,4] <- sum(ct[i,1:3])
  }
  total <- rep(NA, ncol(ct))
  ct <- rbind(ct,total)
  for (i in 1:dim(ct)[2]) {
    ct[4,i] <- sum(ct[1:3,i])
  }
  ct <- melt(ct)
  names(ct) <- c('Exp','Pred','Freq')
  ct$Exp <- as.factor(ct$Exp)
  ct$Pred <- as.factor(ct$Pred)
  
  ct <- dplyr::mutate(ct, Right.Wrong=rep(NA, nrow(ct)))
  ct <- dplyr::mutate(ct,prop = rep(NA, nrow(ct)))
  
  
  
  for (i in as.numeric(row.names(ct[ct$Exp!='total' & ct$Pred!='total', ]))) {
    if (ct$Exp[[i]]=="SIPr&BA") {
      ct$prop[[i]] <- ct$Freq[[i]]/sum(ct$Freq[ct$Exp=="SIPr&BA" & ct$Pred!='total'])
    }
    if (ct$Exp[[i]]=="SIPr") {
      ct$prop[[i]] <- ct$Freq[[i]]/sum(ct$Freq[ct$Exp=="SIPr" & ct$Pred!='total'])
    }
    if (ct$Exp[[i]]=='TAC') {
      ct$prop[[i]] <- ct$Freq[[i]]/sum(ct$Freq[ct$Exp=='TAC' & ct$Pred!='total'])
    }
    if (ct$Exp[i] == ct$Pred[i]) {
      ct$Right.Wrong[i] <- 'Right'
    } else {
      ct$Right.Wrong[i] <- 'Wrong'
    }
  } 
  ct$prop <- round(ct$prop,digits = 3)*100
  
  for (i in as.numeric(row.names(ct[ct$Exp =='total' | ct$Pred =='total', ]))) {
    if (ct$Pred[i]=='total' & ct$Exp[i]!='total') {
      ct$prop[i] <- (ct$Freq[i]/ct$Freq[ct$Exp == 'total' & ct$Pred == 'total'])*100
      ct$Right.Wrong[i] <- 'Size'
    }
    if (ct$Pred[i]!='total' & ct$Exp[i]=='total') {
      ct$prop[[i]] <- (ct$Freq[ct$Exp == ct$Pred[i] & ct$Pred == ct$Pred[i]]/ct$Freq[i])*100
      ct$Right.Wrong[i] <- 'Precision'
    }
    if (ct$Pred[i]=='total' & ct$Exp[i]=='total') {
      ct$prop[[i]] <- (sum(diag(class.table))/ct$Freq[i])*100
      ct$Right.Wrong[i] <- 'Accuracy'
    }
  } 
  ct$prop <- round(ct$prop,digits = 1)
  
  ct <- mutate(ct,value=ct$prop)
  
  for (i in 1:nrow(ct)) {
    if (ct$Right.Wrong[i]!='Right' & ct$Right.Wrong[i]!='Wrong') {
      ct$prop[i] <- 45
    }
  }
  
  base <- ggplot(data = ct, 
                 mapping = aes(x =ordered(Pred, levels =sort(unique(Pred))),
                               y=ordered(Exp, levels =rev(sort(unique(Exp)))),
                               fill = Right.Wrong,
                               alpha = prop))+
    geom_tile(color = 'black',size = 1.2)+
    coord_fixed()+
    geom_text(aes(label=paste(Freq,"\n",
                              '(',value,'%',')',sep = '')),
              size=5,vjust = .5, fontface  = "bold", alpha = 1)+
    scale_fill_manual(values = c(Right = '#66a182', Wrong = '#d1495b',
                                 Size = 'lightgrey', Precision = '#FFCC33',
                                 Accuracy = 'lightblue'))+
    theme(axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size = 11,face = 'bold'),
          axis.text.y = element_text(size = 11,face = 'bold'),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank())
  base
}

# prob.heatmap description:
# 
# Takes a nnet::multinom model object.
# Represents the final probabilities table as a heat map with a color gradient based
# on probability percentage. The correctly predicted case names are colored green,
# wrong but second in percentage are colored orange and absolute wrong is colored red.
# The experimental outcome is presented in its own column, color coded for convenience.

prob.heatmap <- function(model) {
  pred <- predict(model,newdata = data, 'class')
  r.w <- pred == data$class
  probs <- fitted(model)*100
  verif <- data.frame(cbind(data$class, pred, r.w, probs, rep(NA, nrow(probs))))
  row.names(verif) <- row.names(probs)
  
  for (i in 1:dim(verif)[1]) {
    if (verif$r.w[i] == 1) { 
      verif$V7[i] <- "#66a182" # green
    } else {
      if (verif[i,verif[i,1]+3] == min(verif[i,4:6])) {
        verif$V7[i] <- '#d1495b'
      } else {
        if (is.na(verif$V7[i])) {
          verif$V7[i] <-  'tan1'
        }
      }
    }
  }
  
  pro.df <- data.frame(probs)
  pro.df <- rownames_to_column(pro.df)
  pro.df[,1] <- factor(pro.df[,1],levels = pro.df[,1])
  pro.df[,5] <- as.numeric(data$class)
  colnames(pro.df) <- c('Aldehyde',"SIPr&BA",'SIPr','TAC','Exp')
  row.names(pro.df) <- row.names(probs)
  
  long <- melt(pro.df,id.vars = 'Aldehyde')
  long[,3] <- round(long[,3],digits = 2)
  long <- mutate(long, exp_shape=rep(NA,nrow(long)))
  
  
  for (i in 1:nrow(long)) {
    if (long$variable[i] == 'Exp' & long$value[i] == 1) {
      long$exp_shape[i] <- "SIPr&BA"
    }
    if (long$variable[i] == 'Exp' & long$value[i] == 2) {
      long$exp_shape[i] <- 'SIPr'
    }
    if (long$variable[i] == 'Exp' & long$value[i] == 3) {
      long$exp_shape[i] <- 'TAC'
    }
  }
  shape_vec <- long[long$variable == 'Exp',4]
  col_vec <- shape_vec
  
  for (i in 1:length(col_vec)) {
    if (col_vec[i] == "SIPr&BA") {
      col_vec[i] <- 'darkorchid'
    }
    if (col_vec[i] == 'SIPr') {
      col_vec[i] <- "darkslategray"
    }
    if (col_vec[i] == 'TAC') {
      col_vec[i] <- "darkgoldenrod4"
    }
  }
  
  
  
  prob.heatmap <- ggplot(mapping = aes(x = variable,
                                       y =ordered(Aldehyde, 
                                                  levels =rev(factor(pro.df$Aldehyde, 
                                                                     levels = pro.df$Aldehyde)))))+
    geom_tile(data = long[long$variable != 'Exp',], 
              color = 'black',aes(fill = value))+
    coord_fixed(ratio = 0.2)+
    geom_text(data = long[long$variable != 'Exp',], 
              aes(label=value))+
    scale_fill_gradient(name = "% probability",
                        low = "#FFFFFF",
                        high = "dodgerblue3",
                        guide = guide_colorbar(frame.colour = "black", 
                                               ticks.colour = "black"))+
    theme(axis.text.x = element_text(size = 10, face = 'bold',),
          axis.text.y = element_text(size = 10, face = 'bold', 
                                     colour = rev(verif$V7)),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank())+
    scale_x_discrete(position = "top",limits = levels(long$variable))+
    geom_tile(data = long[long$variable == 'Exp',],
              alpha = 0,inherit.aes = F,
              aes(x=rev(variable),
                  y=ordered(Aldehyde, levels =rev(factor(pro.df$Aldehyde, 
                                                         levels = pro.df$Aldehyde)))))+
    geom_text(data = long[long$variable == 'Exp',], label = shape_vec,
              size = 2.4, color = col_vec, fontface = 'bold')
  
  prob.heatmap
}

## Example - Published Model

kf.iter(formula = class ~ C.O + L + B5 + diff_Ca_O,data = kfold.data ,
          folds = 3,out.col = 11,iterations = 1000)

test <- nnet::multinom(class ~ C.O + L + B5 + diff_Ca_O,
                       data = data,
                       maxit = 2000)

mod.info(test)

ct.plot(class.table)

prob.heatmap(test)


