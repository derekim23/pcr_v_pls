#This script involves five main steps:
#Stratified k-fold data split to ensure label balance across folds.
#Nested x-validation for optimizing PLSDA.
#Emulating LASSO stability selection for variable selection.
#Testing the optimized model.
#And finally visualizing some outputs of the model.

#The original dataset has been removed because of its sensitive nature.
#However, any "typical" dataset with input and output matrices should work.


library(mixOmics)
library(tidyverse)
library(caret)
library(CrossValidate)
library(kernlab)
library(boot)
library(glmnet)
library(splitstackshape)
library(Boruta)

#Relevant functions.
avg <- function(x, indices){
  mean(x[indices], trim=0, na.rm=T)
}

stratified_k_fold <- function(data, n_folds = 5){
  #Stratified k-folds generation.
  folds=list()
  data1=data
  
  for(k in 1:n_folds){
    j=n_folds-(k-1)
    if(j>1){
      a=stratified(data1, group = "y", size = 1/j)
      folds[[k]]=a$id
      data1=data1%>%filter(id %in% setdiff(data1$id, a$id))
    } else{
      folds[[k]]=data1$id
    }
  }
  return(folds)
}

get_mode <- function(x) {
  #Mode function
  unique_x <- unique(x)
  tabulate_x <- tabulate(match(x, unique_x))
  unique_x[tabulate_x == max(tabulate_x)][1] #Break ties.
}

#Load data.
load(file=file.path("WPtestDataSet.RData"),verbose=TRUE)
y <- dat0[,OUTCOME_METRIC]
x <- dat0[,phenotypes]

#Ignore y instances that have NA.
idx_not_na <- !is.na(y)
y_clean <- as.factor(y[idx_not_na])
x_clean <- scale(x[idx_not_na,])
subjects <- dat0$UID[idx_not_na]
names(y_clean) <- subjects
rownames(x_clean) <- subjects
dat <- data.frame(1:length(y_clean),y_clean,x_clean)
colnames(dat)[1:2] <- c("id","y")

#Check whether the dataset is balanced.
sum(y_clean==1)
sum(y_clean==0)

#Quick overview of the input x_clean.
glimpse(x_clean)
summary(x_clean)

#Key dims.
p = ncol(x_clean)
n = nrow(x_clean)

set.seed(933933)

#Nested x-validation.
#Please note that the loop below will take a while!
#A greater number of trials leads to memory running out on a 64GB laptop.
gridX <- 1:75
n_folds <- 5
f1 <- c()
acc <- c()
ntrials <- 5
opt_ncomp <- c()

#Max number of components for the PLSR algo without causing singularity.
n_comp = 7
opt_keepX <- matrix(nrow=n_comp, ncol=n_folds*ntrials)

for (j in 1:ntrials){
  #Create balanced folds
  folds = stratified_k_fold(dat, n_folds)
  r = 1
  for (i in folds){
    test_idx <- i
    train_idx <- setdiff(1:n, i)
    
    y_test <- y_clean[test_idx]
    x_test <- x_clean[test_idx,]
    
    y_train <- y_clean[train_idx]
    x_train <- x_clean[train_idx,]
    
    #Inner validation for finding the optimal #s of comps and variables to keep.
    #Any higher ncomp causes problems with matrix inversions.
    mod <- mixOmics::tune.splsda(x_train, y_train, ncomp = n_comp, 
                                 test.keepX = gridX,
                                 folds = n_folds, nrepeat = 10, cpu = 8)
    keepX <- mod$choice.keepX
    ncomp <- mod$choice.ncomp$ncomp
    if (is.null(ncomp)) ncomp = 1
    
    opt_keepX[,(r) + n_folds*(j-1)] <- keepX 
    opt_ncomp <- c(opt_ncomp, ncomp)
    
    mod_train <- mixOmics::splsda(x_train, y_train, ncomp = ncomp, 
                            keepX = keepX[1:ncomp], max.iter = 300)
    pred_test <- predict(mod_train, data.frame(x_test))
    cm <- confusionMatrix(y_test, factor(pred_test$class$max.dist[,ncomp],
                                         levels = c(0,1)))
    f1 <- c(f1, cm$byClass["F1"])
    acc <- c(acc, cm$overall["Accuracy"])
    r = r + 1
  }
}

f1 #Some instances of 1s.
acc[f1 == 1] #Matches.
opt_ncomp[f1 == 1]
opt_keepX[1,f1 == 1]
#Max # of variables chosen = 9
#Choose that to be on the safe side.

#Now, check which variables get picked up consistently in the first comp.
#This is a process emulating the LASSO stability selection.
nvar = 9
ncomp_ = 1
f12 <- c()
acc2 <- c()
var_sel <- matrix(nrow=nvar, ncol=ntrials)
ntrials = 100

for (i in 1:ntrials){
  #50% subsampling as in the LASSO stability criterion.
  sample_mask <- balancedSplit(y_clean, 0.5)
  x_train <- x_clean[sample_mask,]
  y_train <- y_clean[sample_mask]
  x_test <- x_clean[!sample_mask,]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- mixOmics::splsda(x_train, y_train, ncomp = ncomp_, 
                              keepX = nvar, max.iter = 300)
  pred_test <- predict(mod_train, data.frame(x_test))
  cm <- confusionMatrix(y_test, factor(pred_test$class$max.dist[,ncomp_], 
                                       levels = c(0,1)))
  f12 <- c(f12, cm$byClass["F1"])
  acc2 <- c(acc2, cm$overall["Accuracy"])
  var_sel[,i] = names(mod_train$loadings[[1]][,1])[mod_train$loadings[[1]][,1]!=0]
}

#There are four variables selected at least 80% of the time.
tab <- table(var_sel)
sum(tab>=80)
freq_selected <- names(tab)[tab>=80]

#Try fitting a full PLSDA with just those variables.
f13 <- c()
acc3 <- c()
ntrials = 1000

#Test its performance.
for (i in 1:ntrials){
  #80:20 split.
  sample_mask <- balancedSplit(y_clean, 0.8)
  x_train <- x_clean[sample_mask,freq_selected]
  y_train <- y_clean[sample_mask]
  x_test <- x_clean[!sample_mask,freq_selected]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- mixOmics::plsda(x_train, y_train, ncomp = ncomp_, max.iter = 300)
  pred_test <- predict(mod_train, data.frame(x_test))
  cm <- confusionMatrix(y_test, factor(pred_test$class$max.dist[,ncomp_], 
                                       levels=c(0,1)))
  f13 <- c(f13, cm$byClass["F1"])
  acc3 <- c(acc3, cm$overall["Accuracy"])
}

summary(f13)
summary(acc3)

sd(f13)
sd(acc3)

boot.ci(boot(f13, R = 2000, statistic = avg))
boot.ci(boot(acc3, R = 2000, statistic = avg))

#Try fitting a full PLSDA with all variables.
f14 <- c()
acc4 <- c()

#Test its performance.
for (i in 1:ntrials){
  #80:20 split.
  sample_mask <- balancedSplit(y_clean, 0.8)
  x_train <- x_clean[sample_mask,]
  y_train <- y_clean[sample_mask]
  x_test <- x_clean[!sample_mask,]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- mixOmics::plsda(x_train, y_train, ncomp = ncomp_, max.iter = 300)
  pred_test <- predict(mod_train, data.frame(x_test))
  cm <- confusionMatrix(y_test, factor(pred_test$class$max.dist[,ncomp_], 
                                       levels=c(0,1)))
  f14 <- c(f14, cm$byClass["F1"])
  acc4 <- c(acc4, cm$overall["Accuracy"])
}

summary(f14)
summary(acc4)

sd(f14)
sd(acc4)

boot.ci(boot(f14, R = 2000, statistic = avg))
boot.ci(boot(acc4, R = 2000, statistic = avg))

#Model with the entire dataset with three components and selected variables.
mod_full <- mixOmics::plsda(x_clean[,freq_selected], 
                            y_clean, ncomp = ncomp_, max.iter=300)
bg_full <- background.predict(mod_full, comp.predicted=1, dist = "max.dist")

mod_full$prop_expl_var

#Plot decision boundary.
plotIndiv(mod_full, comp = 1:2, group = y_clean, ind.names = T, 
          title = "Decision Boundary [All Data]",
          legend = T, ellipse = T, background = bg_full)

#Save some outputs.
write.csv(mod_full$variates$X,"/Users/Derek/Downloads/x2.csv")
write.csv(mod_full$variates$Y,"/Users/Derek/Downloads/y2.csv")
write.csv(mod_full$Y,"/Users/Derek/Downloads/outcome2.csv")
write.csv(subjects,"/Users/Derek/Downloads/subjects2.csv")
write.csv(mod_full$loadings.star[[1]],"/Users/Derek/Downloads/var_loadings.csv")

#Some correlation calculations.
cor(mod_full$variates$X[,1], mod_full$variates$Y[,1], method = "spearman")
cor(mod_full$variates$X[,2], mod_full$variates$Y[,2], method = "spearman")
cor(mod_full$variates$X[,3], mod_full$variates$Y[,3], method = "spearman")
cor(mod_full$variates$X[,4], mod_full$variates$Y[,4], method = "spearman")

cor(mod_full$variates$X[,1], as.numeric(mod_full$Y), method = "spearman")
cor(mod_full$variates$X[,2], as.numeric(mod_full$Y), method = "spearman")
cor(mod_full$variates$X[,3], as.numeric(mod_full$Y), method = "spearman")
cor(mod_full$variates$X[,4], as.numeric(mod_full$Y), method = "spearman")

#Variable loading plots.
plotLoadings(mod_full, comp = 1, title = 'Variable Loadings of \n Comp. 1', 
             size.title = rel(1.25),
             contrib = 'max', method = 'mean')
plotLoadings(mod_full, comp = 2, title = 'Variable Loadings of \n Comp. 2', 
             size.title = rel(1.25),
             contrib = 'max', method = 'mean')
plotLoadings(mod_full, comp = 3, title = 'Variable Loadings of \n Comp. 3', 
             size.title = rel(1.25),
             contrib = 'max', method = 'mean')
plotLoadings(mod_full, comp = 4, title = 'Variable Loadings of \n Comp. 4', 
             size.title = rel(1.25),
             contrib = 'max', method = 'mean')

#################################
#What about kernelized versions?#
#################################
rbf <- kernelMatrix(x_clean, kernel = rbfdot(sigma = 1))
pol2 <- kernelMatrix(x_clean, kernel = polydot(degree = 2, scale = 1, 
                                               offset = 1))
pol3 <- kernelMatrix(x_clean, kernel = polydot(degree = 3, scale = 1,
                                               offset = 1))
pol5 <- kernelMatrix(x_clean, kernel = polydot(degree = 5, scale = 1, 
                                               offset = 1))

f1_rbf <- c()
acc_rbf <- c()

#RBF kernel first.
for (i in 1:ntrials){
  #80:20 split.
  tryCatch({
    sample_mask <- balancedSplit(y_clean, 0.8)
    x_train <- rbf[sample_mask,]
    y_train <- y_clean[sample_mask]
    x_test <- rbf[!sample_mask,]
    y_test <- y_clean[!sample_mask]
    
    mod_train <- mixOmics::plsda(x_train, y_train, ncomp = ncomp_, max.iter = 300)
    pred_test <- predict(mod_train, x_test)
    cm <- confusionMatrix(y_test, factor(pred_test$class$max.dist[,ncomp_], 
                                         levels=c(0,1)))
    f1_rbf <- c(f1_rbf, cm$byClass["F1"])
    acc_rbf <- c(acc_rbf, cm$overall["Accuracy"])
  }, error = function(e){})
}

summary(f1_rbf)
summary(acc_rbf)
sd(f1_rbf)
sd(acc_rbf)

boot.ci(boot(f1_rbf, R = 2000, statistic = avg))
boot.ci(boot(acc_rbf, R = 2000, statistic = avg))

f1_pol2 <- c()
acc_pol2 <- c()

#Quadratic.
for (i in 1:ntrials){
  #80:20 split.
  sample_mask <- balancedSplit(y_clean, 0.8)
  x_train <- pol2[sample_mask,]
  y_train <- y_clean[sample_mask]
  x_test <- pol2[!sample_mask,]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- mixOmics::plsda(x_train, y_train, ncomp = ncomp_, max.iter = 300)
  pred_test <- predict(mod_train, x_test)
  cm <- confusionMatrix(y_test, factor(pred_test$class$max.dist[,ncomp_], 
                                       levels=c(0,1)))
  f1_pol2 <- c(f1_pol2, cm$byClass["F1"])
  acc_pol2 <- c(acc_pol2, cm$overall["Accuracy"])
}

summary(f1_pol2)
summary(acc_pol2)
sd(f1_pol2,na.rm=T)
sd(acc_pol2)

boot.ci(boot(f1_pol2, R = 2000, statistic = avg))
boot.ci(boot(acc_pol2, R = 2000, statistic = avg))

f1_pol3 <- c()
acc_pol3 <- c()

#Cubic.
for (i in 1:ntrials){
  #80:20 split.
  sample_mask <- balancedSplit(y_clean, 0.8)
  x_train <- pol3[sample_mask,]
  y_train <- y_clean[sample_mask]
  x_test <- pol3[!sample_mask,]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- mixOmics::plsda(x_train, y_train, ncomp = ncomp_, max.iter = 300)
  pred_test <- predict(mod_train, x_test)
  cm <- confusionMatrix(y_test, factor(pred_test$class$max.dist[,ncomp_], 
                                       levels=c(0,1)))
  f1_pol3 <- c(f1_pol3, cm$byClass["F1"])
  acc_pol3 <- c(acc_pol3, cm$overall["Accuracy"])
}

summary(f1_pol3)
summary(acc_pol3)
sd(f1_pol3,na.rm=T)
sd(acc_pol3)

boot.ci(boot(f1_pol3, R = 2000, statistic = avg))
boot.ci(boot(acc_pol3, R = 2000, statistic = avg))

f1_pol5 <- c()
acc_pol5 <- c()

#Quintic.
for (i in 1:ntrials){
  #80:20 split.
  sample_mask <- balancedSplit(y_clean, 0.8)
  x_train <- pol5[sample_mask,]
  y_train <- y_clean[sample_mask]
  x_test <- pol5[!sample_mask,]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- mixOmics::plsda(x_train, y_train, ncomp = ncomp_, max.iter = 300)
  pred_test <- predict(mod_train, x_test)
  cm <- confusionMatrix(y_test, factor(pred_test$class$max.dist[,ncomp_], 
                                       levels=c(0,1)))
  f1_pol5 <- c(f1_pol5, cm$byClass["F1"])
  acc_pol5 <- c(acc_pol5, cm$overall["Accuracy"])
}

summary(f1_pol5)
summary(acc_pol5)
sd(f1_pol5,na.rm=T)
sd(acc_pol5)

boot.ci(boot(f1_pol5, R = 2000, statistic = avg))
boot.ci(boot(acc_pol5, R = 2000, statistic = avg))

#Try regressing on PCs.
pcs <- read.csv("pa_measures.csv")
pcs = pcs[!is.na(y),-1]

lambda = c(0, 0.05, 0.1, 0.125, 0.175, 0.18, 0.2, 0.205)
f1_glm <- c()
acc_glm <- c()
ntrials_glm <- 1000
opt_lambda <- c()

#Run this by sampling instead of nested x-val, because glmnet complains
#of data imbalance otherwise.
for (i in 1:ntrials_glm){
  test_mask <- balancedSplit(y_clean, 0.2)
  train_mask <- !test_mask
  
  test_idx <- (1:n)[test_mask]
  train_idx <- (1:n)[train_mask]
  
  y_test <- y_clean[test_idx]
  x_test <- pcs[test_idx,]
  
  y_train <- y_clean[train_idx]
  x_train <- pcs[train_idx,]
  
  #Inner validation for finding the optimal hyperparameter.
  mod <- cv.glmnet(as.matrix(x_train), y_train, family = "binomial", 
                       alpha = 1, lambda = lambda, n_foldss = 5, 
                       type.measure = "class")
  opt_lambda <- c(opt_lambda, mod$lambda.min)
  
  mod_train <- glmnet(as.matrix(x_test), y_test, family = "binomial",
                    alpha = 1, lambda = mod$lambda.min)
  pred_test <- predict(mod_train, newx = as.matrix(x_test), s = mod$lambda.min,
                       type = "class")
  cm <- confusionMatrix(y_test, factor(pred_test,
                                       levels = c(0,1)))
  f1_glm <- c(f1_glm, cm$byClass["F1"])
  acc_glm <- c(acc_glm, cm$overall["Accuracy"])
}
f1_glm
acc_glm
ld <- get_mode(opt_lambda[acc_glm == 1 & f1_glm == 1])

#Stability
f1_glm2 <- c()
acc_glm2 <- c()
nvar = 7
n_trials_glm_stability = 100
var_sel_glm <- matrix(nrow=nvar, ncol=n_trials_glm_stability)

for (i in 1:n_trials_glm_stability){
  #50% subsampling as in the LASSO stability criterion.
  sample_mask <- balancedSplit(y_clean, 0.5)
  x_train <- pcs[sample_mask,]
  y_train <- y_clean[sample_mask]
  x_test <- pcs[!sample_mask,]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- glmnet(x_train, y_train, family = "binomial", alpha = 1, lambda = ld)
  pred_test <- predict(mod_train, newx = as.matrix(x_test), s = ld, type = "class")
  cm <- confusionMatrix(y_test, factor(pred_test, levels = c(0,1)))
  f1_glm2 <- c(f1_glm2, cm$byClass["F1"])
  acc_glm2 <- c(acc_glm2, cm$overall["Accuracy"])
  selected <- names(mod_train$beta[,1])[mod_train$beta[,1]!=0]
  var_sel_glm[,i] = c(selected, rep(NA, 7 - length(selected)))
}
tab_glm <- table(var_sel_glm)
sum(tab_glm>=90)
freq_selected_glm <- names(tab_glm)[tab_glm>=90]

#Try fitting a full PLSDA with just those variables.
f1_glm3 <- c()
acc_glm3 <- c()

#Test its performance.
for (i in 1:ntrials_glm){
  #80:20 split.
  sample_mask <- balancedSplit(y_clean, 0.8)
  x_train <- pcs[sample_mask,freq_selected_glm]
  y_train <- y_clean[sample_mask]
  x_test <- pcs[!sample_mask,freq_selected_glm]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- glmnet(x_train, y_train, family = "binomial", alpha = 1, lambda = 0)
  pred_test <- predict(mod_train, newx = as.matrix(x_test), s = ld, type = "class")
  cm <- confusionMatrix(y_test, factor(pred_test, levels = c(0,1)))
  f1_glm3 <- c(f1_glm3, cm$byClass["F1"])
  acc_glm3 <- c(acc_glm3, cm$overall["Accuracy"])
}

summary(f1_glm3)
summary(acc_glm3)

sd(f1_glm3, na.rm = T)
sd(acc_glm3, na.rm = T)

boot.ci(boot(f1_glm3, R = 2000, statistic = avg))
boot.ci(boot(acc_glm3, R = 2000, statistic = avg))

#What about a model with all variables and no regularization?
f1_glm4 <- c()
acc_glm4 <- c()

for (i in 1:ntrials_glm){
  #80:20 split.
  sample_mask <- balancedSplit(y_clean, 0.8)
  x_train <- pcs[sample_mask,]
  y_train <- y_clean[sample_mask]
  x_test <- pcs[!sample_mask,]
  y_test <- y_clean[!sample_mask]
  
  mod_train <- glmnet(x_train, y_train, family = "binomial", alpha = 1, lambda = 0)
  pred_test <- predict(mod_train, newx = as.matrix(x_test), s = ld, type = "class")
  cm <- confusionMatrix(y_test, factor(pred_test, levels = c(0,1)))
  f1_glm4 <- c(f1_glm4, cm$byClass["F1"])
  acc_glm4 <- c(acc_glm4, cm$overall["Accuracy"])
}

summary(f1_glm4)
summary(acc_glm4)

sd(f1_glm4, na.rm = T)
sd(acc_glm4, na.rm = T)

boot.ci(boot(f1_glm4, R = 2000, statistic = avg))
boot.ci(boot(acc_glm4, R = 2000, statistic = avg))
