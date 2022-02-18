#Candidate algorithms to test:
#PCR, sPCR, PLS-DA, sPLS-DA

#=============================================================================#

#PLSR
library(mixOmics)
library(sparsepca)
library(nnet)

#=============================================================================#

#Functions.
train_dev_test_split <- function(label, train_ratio=0.6, dev_ratio=0.2){
  #Function for creating data splits.
  #Ad-hoc for the lymphoma dataset.
  
  #Count the number of target class instances.
  ews <- sum(label=="EWS") #Again, ews
  bl <- sum(label=="BL") #bl
  nb <- sum(label=="NB") #nb
  rms <- sum(label=="RMS") #rms
  
  #Use train/dev/test split of 60/20/20.
  train_lab_ews <- sample(which(label == "EWS"), floor(train_ratio*ews))
  train_lab_bl <- sample(which(label == "BL"), floor(train_ratio*bl))
  train_lab_nb <- sample(which(label == "NB"), floor(train_ratio*nb))
  train_lab_rms <- sample(which(label == "RMS"), floor(train_ratio*rms))
  
  #trainaining indices.
  train_idx <- c(train_lab_ews, train_lab_bl, train_lab_nb, train_lab_rms)
  all_idx = 1:length(label)
  
  #Test candidate indices.
  test_cand <- sort(all_idx[!(all_idx%in%train_idx)])
  test_ratio = 1-dev_ratio-train_ratio
  
  test_lab_ews <- sample(test_cand[which(label[test_cand] == "EWS")], ceiling(test_ratio*ews))
  test_lab_bl <- sample(test_cand[which(label[test_cand] == "BL")], ceiling(test_ratio*bl))
  test_lab_nb <- sample(test_cand[which(label[test_cand] == "NB")], ceiling(test_ratio*nb))
  test_lab_rms <- sample(test_cand[which(label[test_cand] == "RMS")], ceiling(test_ratio*rms))
  
  #Test indices.
  test_idx <- c(test_lab_ews, test_lab_bl, test_lab_nb, test_lab_rms)
  
  #Validation indices.
  dev_idx <- test_cand[!test_cand%in%test_idx]
  
  output <- list(train_idx, dev_idx, test_idx)
  names(output) <- c("train_idx", "dev_idx", "test_idx")
  
  return(output)
}

f1 <- function(cm){
  #Macro-/micro-F1 calculator.
  
  #True positives.
  tp <- diag(cm)
  
  #Row and column sums of the confusion matrix.
  row_sum <- apply(cm,1,sum)
  col_sum <- apply(cm,2,sum)
  
  #Precision and recall.
  prec <- tp/col_sum
  rec <-  tp/row_sum
  
  #Macro- and micro-f1.
  macro_f1 <- mean(2 * prec * rec/(prec + rec))
  micro_f1 <- sum(tp)/sum(row_sum)
  out <- c(macro_f1, micro_f1)
  names(out) <- c("Macro-F1", "Micro-F1")
  return(out)
}

get_mode <- function(x) {
  #Mode function
  unique_x <- unique(x)
  tabulate_x <- tabulate(match(x, unique_x))
  unique_x[tabulate_x == max(tabulate_x)][1] #Break ties.
}

#=============================================================================#

#Load data from spls.
data("srbct")
?srbct

x <- scale(srbct$gene)
y <- as.factor(srbct$class)

#Number of gene expression features.
ncol(x)

#Sample size.
n <- nrow(x)

#Explore the dataset.
#Quite a bit of class imbalance.
bp <- barplot(table(y), names = levels(y),
              xlab = "Tumor Type", ylab = "Freqency",
              main = "Bar-plot of Tumor Types (Target)",
              yaxt="n", cex.lab = 1.4, cex.main = 1.4)

#Count the number of target class instances.
(ews <- sum(y=="EWS"))
(bl <- sum(y=="BL")) 
(nb <- sum(y=="NB"))
(rms <- sum(y=="RMS"))

text(bp, 0, c(ews, bl, nb, rms), cex=1.4, pos=3, col="white", font = 2)

#=============================================================================#

#Can't arbitrarily cross validate, because of the severe class imbalance and
#small n, so the train_dev_test_split function has been implemented to
#conduct stratified (by target class) sampling.
macro_pcr_best_n_comps <- NULL
macro_pcr_test_scores <- NULL
micro_pcr_best_n_comps <- NULL
micro_pcr_test_scores <- NULL

macro_spcr_best_n_comps <- NULL
macro_spcr_test_scores <- NULL
micro_spcr_best_n_comps <- NULL
micro_spcr_test_scores <- NULL

macro_plsda_best_n_comps <- NULL
macro_plsda_test_scores <- NULL
micro_plsda_best_n_comps <- NULL
micro_plsda_test_scores <- NULL

macro_splsda_best_n_comps <- NULL
macro_splsda_test_scores <- NULL
micro_splsda_best_n_comps <- NULL
micro_splsda_test_scores <- NULL

#Will need to loop starting from here.

for (j in 1:40){ 
  #Run tests 40 times for each method and append the optimal # of components
  #and test scores to the variables declared above.
  #First get the random, balanced indices for the data split.
  idx <- train_dev_test_split(y)
  train_idx <- idx$train_idx
  dev_idx <- idx$dev_idx
  test_idx <- idx$test_idx
  
  #Create the data split.
  x_train <- x[train_idx,]
  x_dev <- x[dev_idx,]
  x_test <- x[test_idx,]
  
  y_train <- y[train_idx]
  y_dev <- y[dev_idx]
  y_test <- y[test_idx]
  
  #Create the PCs.
  x_train_pca <- pca(x_train, ncomp = 12)
  x_dev_pca <- pca(x_dev, ncomp = 12)
  x_test_pca <- pca(x_test, ncomp = 12)
  
  #Run the validation grid search.
  best_macro_pcr <- NULL
  best_micro_pcr <- NULL
  best_macro_pcr_dev <- 0
  best_micro_pcr_dev <- 0
  best_macro_pcr_n_comp <- 0
  best_micro_pcr_n_comp <- 0
  
  for (i in 2:10){
    pcr <- multinom(y_train ~., data = data.frame(y_train,x_train_pca$variates$X[,1:i]), 
                    maxit = 250)
    pred_dev_pcr <- predict(pcr, newdata=data.frame(x_dev_pca$variates$X[,1:i]))
    
    #Create confusion matrix and calculate F1 scores.
    cm_dev_pcr <- as.matrix(confusionMatrix(y_dev, pred_dev_pcr))
    f1_dev_pcr <- f1(cm_dev_pcr)
    if (!is.na(f1_dev_pcr[1]) && best_macro_pcr_dev < f1_dev_pcr[1]){
      best_macro_pcr_dev <- f1_dev_pcr[1]
      best_macro_pcr <- pcr
      best_macro_pcr_n_comp <- i
    }
    if (best_micro_pcr_dev < f1_dev_pcr[2]){
      best_micro_pcr_dev <- f1_dev_pcr[2]
      best_micro_pcr <- pcr
      best_micro_pcr_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  #In some cases, macro F1 comes out NA because no predictions are made for a certain class.
  if (!is.null(best_macro_pcr)){
    pred_test_macro_pcr <- predict(best_macro_pcr, newdata = x_test_pca$variates$X[,1:best_macro_pcr_n_comp])
    cm_test_macro_pcr <- as.matrix(confusionMatrix(y_test, pred_test_macro_pcr))
    f1_test_macro_pcr <- f1(cm_test_macro_pcr)[1]
  } else {
    f1_test_macro_pcr <- NA
  }
  
  pred_test_micro_pcr <- predict(best_micro_pcr, newdata = x_test_pca$variates$X[,1:best_micro_pcr_n_comp])
  cm_test_micro_pcr <- as.matrix(confusionMatrix(y_test, pred_test_micro_pcr))
  f1_test_micro_pcr <- f1(cm_test_micro_pcr)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_pcr_best_n_comps <- c(macro_pcr_best_n_comps, best_macro_pcr_n_comp)
  macro_pcr_test_scores <- c(macro_pcr_test_scores, f1_test_macro_pcr)
  
  micro_pcr_best_n_comps <- c(micro_pcr_best_n_comps, best_micro_pcr_n_comp)
  micro_pcr_test_scores <- c(micro_pcr_test_scores, f1_test_micro_pcr)
  
  #############
  
  #Create the sPCs; keep default alpha of 0.0001 in the interest of time.
  #Run the validation grid search.
  best_macro_spcr <- NULL
  best_micro_spcr <- NULL
  best_macro_spcr_dev <- 0
  best_micro_spcr_dev <- 0
  best_macro_spcr_n_comp <- 0
  best_micro_spcr_n_comp <- 0
  
  for (i in 2:10){
    #alpha value set s.t. the median number of selected variables is roughly
    #200 across components.
    x_train_spca <- spca(x_train, k = i, alpha = 0.000693)
    x_dev_spca <- spca(x_dev, k = i, alpha = 0.000693)
    x_test_spca <- spca(x_test, k = i, alpha = 0.000693)
    spcr <- multinom(y_train ~., data = data.frame(y_train,x_train_spca$scores[,1:i]), 
                    maxit = 250)
    pred_dev_spcr <- predict(spcr, newdata=data.frame(x_dev_spca$scores[,1:i]))
    
    #Create confusion matrix and calculate F1 scores.
    cm_dev_spcr <- as.matrix(confusionMatrix(y_dev, pred_dev_spcr))
    f1_dev_spcr <- f1(cm_dev_spcr)
    if (!is.na(f1_dev_spcr[1]) && best_macro_spcr_dev < f1_dev_spcr[1]){
      best_macro_spcr_dev <- f1_dev_spcr[1]
      best_macro_spcr <- spcr
      best_macro_spcr_n_comp <- i
    }
    if (best_micro_spcr_dev < f1_dev_spcr[2]){
      best_micro_spcr_dev <- f1_dev_spcr[2]
      best_micro_spcr <- spcr
      best_micro_spcr_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  if (!is.null(best_macro_spcr)){
    pred_test_macro_spcr <- predict(best_macro_spcr, 
                                    newdata = data.frame(x_test_spca$scores[,1:best_macro_spcr_n_comp]))
    cm_test_macro_spcr <- as.matrix(confusionMatrix(y_test, pred_test_macro_spcr))
    f1_test_macro_spcr <- f1(cm_test_macro_spcr)[1]
  } else {
    f1_test_macro_spcr <- NA
  }
  
  pred_test_micro_spcr <- predict(best_micro_spcr, 
                                  newdata = data.frame(x_test_spca$scores[,1:best_micro_spcr_n_comp]))
  cm_test_micro_spcr <- as.matrix(confusionMatrix(y_test, pred_test_micro_spcr))
  f1_test_micro_spcr <- f1(cm_test_micro_spcr)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_spcr_best_n_comps <- c(macro_spcr_best_n_comps, best_macro_spcr_n_comp)
  macro_spcr_test_scores <- c(macro_spcr_test_scores, f1_test_macro_spcr)
  
  micro_spcr_best_n_comps <- c(micro_spcr_best_n_comps, best_micro_spcr_n_comp)
  micro_spcr_test_scores <- c(micro_spcr_test_scores, f1_test_micro_spcr)
  
  best_macro_plsda <- NULL
  best_micro_plsda <- NULL
  best_macro_plsda_dev <- 0
  best_micro_plsda_dev <- 0
  best_macro_plsda_n_comp <- 0
  best_micro_plsda_n_comp <- 0
  
  for (i in 2:10){
    plsda <- mixOmics::plsda(data.frame(x_train), y_train, ncomp = i, max.iter = 250)
    pred_dev_plsda <- predict(plsda, data.frame(x_dev))
    #Create confusion matrix and calculate F1 scores.
    cm_dev_plsda <- as.matrix(confusionMatrix(y_dev, 
                                              factor(pred_dev_plsda$class$max.dist[,i],
                                                        levels = levels(y_dev))))
    f1_dev_plsda <- f1(cm_dev_plsda)
    if (!is.na(f1_dev_plsda[1]) && best_macro_plsda_dev < f1_dev_plsda[1]){
      best_macro_plsda_dev <- f1_dev_plsda[1]
      best_macro_plsda <- plsda
      best_macro_plsda_n_comp <- i
    }
    if (best_micro_plsda_dev < f1_dev_plsda[2]){
      best_micro_plsda_dev <- f1_dev_plsda[2]
      best_micro_plsda <- plsda
      best_micro_plsda_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  pred_test_macro_plsda <- predict(best_macro_plsda, newdata = x_test)
  cm_test_macro_plsda <- as.matrix(confusionMatrix(y_test, 
                                                    factor(pred_test_macro_plsda$
                                                             class$max.dist[,best_macro_plsda$ncomp],
                                                           levels = levels(y_test))))
  f1_test_macro_plsda <- f1(cm_test_macro_plsda)[1]
  
  pred_test_micro_plsda <- predict(best_micro_plsda, newdata = x_test)
  cm_test_micro_plsda <- as.matrix(confusionMatrix(y_test, 
                                                   factor(pred_test_micro_plsda$
                                                            class$max.dist[,best_micro_plsda$ncomp],
                                                          levels = levels(y_test))))
  f1_test_micro_plsda <- f1(cm_test_micro_plsda)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_plsda_best_n_comps <- c(macro_plsda_best_n_comps, best_macro_plsda_n_comp)
  macro_plsda_test_scores <- c(macro_plsda_test_scores, f1_test_macro_plsda)
  
  micro_plsda_best_n_comps <- c(micro_plsda_best_n_comps, best_micro_plsda_n_comp)
  micro_plsda_test_scores <- c(micro_plsda_test_scores, f1_test_micro_plsda)
  
  best_macro_splsda <- NULL
  best_micro_splsda <- NULL
  best_macro_splsda_dev <- 0
  best_micro_splsda_dev <- 0
  best_macro_splsda_n_comp <- 0
  best_micro_splsda_n_comp <- 0
  
  for (i in 2:10){
    #Keeping 200 variables for component 1, roughly in line with the approach for sPCA.
    splsda <- mixOmics::splsda(data.frame(x_train), y_train, keepX = 200,
                               ncomp = i, max.iter = 250)
    pred_dev_splsda <- predict(splsda, data.frame(x_dev))
    #Create confusion matrix and calculate F1 scores.
    cm_dev_splsda <- as.matrix(confusionMatrix(y_dev, 
                                              factor(pred_dev_splsda$class$max.dist[,i],
                                                     levels = levels(y_dev))))
    f1_dev_splsda <- f1(cm_dev_splsda)
    if (!is.na(f1_dev_splsda[1]) && best_macro_splsda_dev < f1_dev_splsda[1]){
      best_macro_splsda_dev <- f1_dev_splsda[1]
      best_macro_splsda <- splsda
      best_macro_splsda_n_comp <- i
    }
    if (best_micro_splsda_dev < f1_dev_splsda[2]){
      best_micro_splsda_dev <- f1_dev_splsda[2]
      best_micro_splsda <- splsda
      best_micro_splsda_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  pred_test_macro_splsda <- predict(best_macro_splsda, newdata = x_test)
  cm_test_macro_splsda <- as.matrix(confusionMatrix(y_test, 
                                                    factor(pred_test_macro_splsda$
                                                             class$max.dist[,best_macro_splsda$ncomp],
                                                           levels = levels(y_test))))
  f1_test_macro_splsda <- f1(cm_test_macro_splsda)[1]
  
  pred_test_micro_splsda <- predict(best_micro_splsda, newdata = x_test)
  cm_test_micro_splsda <- as.matrix(confusionMatrix(y_test, 
                                                    factor(pred_test_micro_splsda$
                                                             class$max.dist[,best_micro_splsda$ncomp],
                                                           levels = levels(y_test))))
  f1_test_micro_splsda <- f1(cm_test_micro_splsda)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_splsda_best_n_comps <- c(macro_splsda_best_n_comps, best_macro_splsda_n_comp)
  macro_splsda_test_scores <- c(macro_splsda_test_scores, f1_test_macro_splsda)
  
  micro_splsda_best_n_comps <- c(micro_splsda_best_n_comps, best_micro_splsda_n_comp)
  micro_splsda_test_scores <- c(micro_splsda_test_scores, f1_test_micro_splsda)
}

#Check best stats from tests
get_mode(macro_pcr_best_n_comps[macro_pcr_best_n_comps!=0])
get_mode(micro_pcr_best_n_comps)
mean(macro_pcr_best_n_comps[macro_pcr_best_n_comps!=0])
mean(micro_pcr_best_n_comps)
sd(macro_pcr_best_n_comps[macro_pcr_best_n_comps!=0])
sd(micro_pcr_best_n_comps)
mean(macro_pcr_test_scores, na.rm = TRUE)
mean(micro_pcr_test_scores, na.rm = TRUE)
sd(macro_pcr_test_scores, na.rm = TRUE)
sd(micro_pcr_test_scores, na.rm = TRUE)

get_mode(macro_spcr_best_n_comps[macro_spcr_best_n_comps!=0])
get_mode(micro_spcr_best_n_comps)
mean(macro_spcr_best_n_comps[macro_spcr_best_n_comps!=0])
mean(micro_spcr_best_n_comps)
sd(macro_spcr_best_n_comps[macro_spcr_best_n_comps!=0])
sd(micro_spcr_best_n_comps)
mean(macro_spcr_test_scores, na.rm = TRUE)
mean(micro_spcr_test_scores, na.rm = TRUE)
sd(macro_spcr_test_scores, na.rm = TRUE)
sd(micro_spcr_test_scores, na.rm = TRUE)

get_mode(macro_plsda_best_n_comps)
get_mode(micro_plsda_best_n_comps)
mean(macro_plsda_best_n_comps)
mean(micro_plsda_best_n_comps)
sd(macro_plsda_best_n_comps)
sd(micro_plsda_best_n_comps)
mean(macro_plsda_test_scores, na.rm = TRUE)
mean(micro_plsda_test_scores, na.rm = TRUE)
sd(macro_plsda_test_scores, na.rm = TRUE)
sd(micro_plsda_test_scores, na.rm = TRUE)

get_mode(macro_splsda_best_n_comps)
get_mode(micro_splsda_best_n_comps)
mean(macro_splsda_best_n_comps)
mean(micro_splsda_best_n_comps)
sd(macro_splsda_best_n_comps)
sd(micro_splsda_best_n_comps)
mean(macro_splsda_test_scores, na.rm = TRUE)
mean(micro_splsda_test_scores, na.rm = TRUE)
sd(macro_splsda_test_scores, na.rm = TRUE)
sd(micro_splsda_test_scores, na.rm = TRUE)

#Visualize the PCs.
par(mar=c(5.1, 4.1, 4.1, 7), xpd=TRUE)
pca_full <- pca(x)
plot(pca_full$variates$X[,1:2], col = y, 
     main = "First Two PCs", 
     cex.lab = 1.4, cex.main = 1.4,
     xlab = paste("PC1 (", round(pca_full$prop_expl_var[[1]][1]*100,1),"%)", sep =""),
     ylab = paste("PC2 (", round(pca_full$prop_expl_var[[1]][2]*100,1),"%)", sep =""))
legend("right", inset=c(-0.3,0), legend = levels(y), col = 1:4, pch = 1)

#Visualize the sPCs.
spca_full <- spca(x, k = 2, alpha = 0.000693)
plot(spca_full$scores[,1:2], col = y, 
     main = "First Two sPCs",
     xlab = "sPC1", ylab = "sPC2",
     cex.lab = 1.4, cex.main = 1.4)
legend("right", inset=c(-0.3,0), legend = levels(y), col = 1:4, pch = 1)

#Visualize PLSDA solution.
plsda_full <- mixOmics::plsda(x, y, ncomp = 3, max.iter = 250)
background <- background.predict(plsda_full, comp.predicted=2, dist = "max.dist")
plotIndiv(plsda_full, comp = 1:2, group = y, ind.names = F, 
          title = "PLSDA Decision Boundaries",
          legend = T, ellipse = T, background = background)


#Visualize PLSDA solution.
splsda_full <- mixOmics::splsda(x, y, ncomp = 3, keepX = 200, max.iter = 250)
background <- background.predict(splsda_full, comp.predicted=2, dist = "max.dist")
plotIndiv(splsda_full, comp = c(1,2), group = y, ind.names = F, 
          title = "sPLSDA Decision Boundaries",
          legend = T, ellipse = T, background = background)

#Visualize spca loading factors.
library(ggplot2)
spca_full_loadings <- spca_full$loadings[,1][abs(spca_full$loadings[,1]) > 0.049]
genes <- colnames(x)[abs(spca_full$loadings[,1]) > 0.049]
group <- as.character(spca_full_loadings)
group[spca_full_loadings<0] = "Negative Loadings"
group[spca_full_loadings>0] = "Positive Loadings"

spca_full_df <- data.frame(Group = group, Loadings = spca_full_loadings, Genes = genes)

#Order the dataset by the absolute value of loadings.
spca_full_df <- spca_full_df[abs(spca_full_df$Loadings)>0.049,]
spca_full_df <- spca_full_df[order(abs(spca_full_df$Loadings), decreasing = T),]

ggplot(spca_full_df, aes(x = rownames(spca_full_df), y = Loadings, fill = Group)) +
  coord_flip() +
  geom_bar(stat = "identity", position = "identity", width = 0.525) +
  ggtitle("sPCA Loadings of Comp. 1 for \nFirst 200 Genes by Absolute Loadings") + xlab("Genes") +
  theme(legend.position="top", 
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
        plot.title = element_text(hjust = 0.5))
plot(spca_full$loadings[,1])

#Plot loading factors of sPLSDA.
plotLoadings(splsda_full, comp = 1, title = 'sPLSDA Loadings of \n Comp. 1 for Selected 200', 
             size.title = rel(1.25),
             contrib = 'max', method = 'mean')

#Plot loading factors of sPLSDA
plotLoadings(splsda_full, comp = 2, title = 'sPLSDA Loadings of \n Comp. 2', 
             size.title = rel(1.25),
             contrib = 'max', method = 'mean')

#PLSDA on the selected variables.
x_selected <- x[,selectVar(splsda_full)[1]$name]

macro_plsda_s_best_n_comps <- NULL
macro_plsda_s_test_scores <- NULL
micro_plsda_s_best_n_comps <- NULL
micro_plsda_s_test_scores <- NULL

#Try testing plsda_s with selected variables.
for (j in 1:40){ 
  #Run tests 40 times for each method and append the optimal # of components
  #and test scores to the variables declared above.
  #First get the random, balanced indices for the data split.
  idx <- train_dev_test_split(y)
  train_idx <- idx$train_idx
  dev_idx <- idx$dev_idx
  test_idx <- idx$test_idx
  
  #Create the data split.
  x_train <- x_selected[train_idx,]
  x_dev <- x_selected[dev_idx,]
  x_test <- x_selected[test_idx,]
  
  y_train <- y[train_idx]
  y_dev <- y[dev_idx]
  y_test <- y[test_idx]
  
  best_macro_plsda_s <- NULL
  best_micro_plsda_s <- NULL
  best_macro_plsda_s_dev <- 0
  best_micro_plsda_s_dev <- 0
  best_macro_plsda_s_n_comp <- 0
  best_micro_plsda_s_n_comp <- 0
  
  for (i in 2:10){
    plsda_s <- mixOmics::plsda(data.frame(x_train), y_train, ncomp = i, max.iter = 250)
    pred_dev_plsda_s <- predict(plsda_s, data.frame(x_dev))
    #Create confusion matrix and calculate F1 scores.
    cm_dev_plsda_s <- as.matrix(confusionMatrix(y_dev, 
                                              factor(pred_dev_plsda_s$class$max.dist[,i],
                                                     levels = levels(y_dev))))
    f1_dev_plsda_s <- f1(cm_dev_plsda_s)
    if (!is.na(f1_dev_plsda_s[1]) && best_macro_plsda_s_dev < f1_dev_plsda_s[1]){
      best_macro_plsda_s_dev <- f1_dev_plsda_s[1]
      best_macro_plsda_s <- plsda_s
      best_macro_plsda_s_n_comp <- i
    }
    if (best_micro_plsda_s_dev < f1_dev_plsda_s[2]){
      best_micro_plsda_s_dev <- f1_dev_plsda_s[2]
      best_micro_plsda_s <- plsda_s
      best_micro_plsda_s_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  pred_test_macro_plsda_s <- predict(best_macro_plsda_s, newdata = x_test)
  cm_test_macro_plsda_s <- as.matrix(confusionMatrix(y_test, 
                                                   factor(pred_test_macro_plsda_s$
                                                            class$max.dist[,best_macro_plsda_s$ncomp],
                                                          levels = levels(y_test))))
  f1_test_macro_plsda_s <- f1(cm_test_macro_plsda_s)[1]
  
  pred_test_micro_plsda_s <- predict(best_micro_plsda_s, newdata = x_test)
  cm_test_micro_plsda_s <- as.matrix(confusionMatrix(y_test, 
                                                   factor(pred_test_micro_plsda_s$
                                                            class$max.dist[,best_micro_plsda_s$ncomp],
                                                          levels = levels(y_test))))
  f1_test_micro_plsda_s <- f1(cm_test_micro_plsda_s)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_plsda_s_best_n_comps <- c(macro_plsda_s_best_n_comps, best_macro_plsda_s_n_comp)
  macro_plsda_s_test_scores <- c(macro_plsda_s_test_scores, f1_test_macro_plsda_s)
  
  micro_plsda_s_best_n_comps <- c(micro_plsda_s_best_n_comps, best_micro_plsda_s_n_comp)
  micro_plsda_s_test_scores <- c(micro_plsda_s_test_scores, f1_test_micro_plsda_s)
}

#Check best stats from tests.
get_mode(macro_plsda_s_best_n_comps)
get_mode(micro_plsda_s_best_n_comps)
mean(macro_plsda_s_best_n_comps)
mean(micro_plsda_s_best_n_comps)
sd(macro_plsda_s_best_n_comps)
sd(micro_plsda_s_best_n_comps)
mean(macro_plsda_s_test_scores, na.rm = TRUE)
mean(micro_plsda_s_test_scores, na.rm = TRUE)
sd(macro_plsda_s_test_scores, na.rm = TRUE)
sd(micro_plsda_s_test_scores, na.rm = TRUE)

#==================================================================#
#Try another dataset: nutrimouse.
#==================================================================#

#Function.
train_dev_test_split2 <- function(label, train_ratio=0.6, dev_ratio=0.2){
  #Function for creating data splits.
  #Ad-hoc for the lymphoma dataset.
  
  #Count the number of target class instances.
  coc <- sum(label=="coc") 
  fish <- sum(label=="fish")
  lin <- sum(label=="lin") 
  ref <- sum(label=="ref") 
  sun <- sum(label=="sun") 
  
  #Use train/dev/test split of 60/20/20.
  train_lab_coc <- sample(which(label == "coc"), floor(train_ratio*coc))
  train_lab_fish <- sample(which(label == "fish"), floor(train_ratio*fish))
  train_lab_lin <- sample(which(label == "lin"), floor(train_ratio*lin))
  train_lab_ref <- sample(which(label == "ref"), floor(train_ratio*ref))
  train_lab_sun <- sample(which(label == "sun"), floor(train_ratio*sun))
  
  #trainaining indices.
  train_idx <- c(train_lab_coc, train_lab_fish, train_lab_lin, train_lab_ref, 
                 train_lab_sun)
  all_idx = 1:length(label)
  
  #Test candidate indices.
  test_cand <- sort(all_idx[!(all_idx%in%train_idx)])
  test_ratio = 1-dev_ratio-train_ratio
  
  test_lab_coc <- sample(test_cand[which(label[test_cand] == "coc")], ceiling(test_ratio*coc))
  test_lab_fish <- sample(test_cand[which(label[test_cand] == "fish")], ceiling(test_ratio*fish))
  test_lab_lin <- sample(test_cand[which(label[test_cand] == "lin")], ceiling(test_ratio*lin))
  test_lab_ref <- sample(test_cand[which(label[test_cand] == "ref")], ceiling(test_ratio*ref))
  test_lab_sun <- sample(test_cand[which(label[test_cand] == "sun")], ceiling(test_ratio*sun))
  
  #Test indices.
  test_idx <- c(test_lab_coc, test_lab_fish, test_lab_lin, test_lab_ref, 
                test_lab_sun)
  
  #Validation indices.
  dev_idx <- test_cand[!test_cand%in%test_idx]
  
  output <- list(train_idx, dev_idx, test_idx)
  names(output) <- c("train_idx", "dev_idx", "test_idx")
  
  return(output)
}

#===========================================================================#
data("nutrimouse")
?nutrimouse

x2 <- list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
y2 <-  nutrimouse$diet

bp2 <- barplot(table(y2), names = levels(y2),
              xlab = "Tumor Type", ylab = "Freqency",
              main = "Bar-plot of Tumor Types (Target)",
              yaxt="n", cex.lab = 1.4, cex.main = 1.4)

#Count the number of target class instances.
(ews <- sum(y=="EWS"))
(bl <- sum(y=="BL")) 
(nb <- sum(y=="NB"))
(rms <- sum(y=="RMS"))

text(bp, 0, c(ews, bl, nb, rms), cex=1.4, pos=3, col="white", font = 2)

#===========================================================================#

#Can't arbitrarily cross validate, because of the severe class imbalance and
#small n, so the train_dev_test_split function has been implemented to
#conduct stratified (by target class) sampling.
macro_pcr2_best_n_comps <- NULL
macro_pcr2_test_scores <- NULL
micro_pcr2_best_n_comps <- NULL
micro_pcr2_test_scores <- NULL

macro_spcr2_best_n_comps <- NULL
macro_spcr2_test_scores <- NULL
micro_spcr2_best_n_comps <- NULL
micro_spcr2_test_scores <- NULL

macro_plsda2_best_n_comps <- NULL
macro_plsda2_test_scores <- NULL
micro_plsda2_best_n_comps <- NULL
micro_plsda2_test_scores <- NULL

macro_splsda2_best_n_comps <- NULL
macro_splsda2_test_scores <- NULL
micro_splsda2_best_n_comps <- NULL
micro_splsda2_test_scores <- NULL

#Will need to loop starting from here.

for (j in 1:40){ 
  #Run tests 40 times for each method and append the optimal # of components
  #and test scores to the variables declared above.
  #First get the random, balanced indices for the data split.
  idx <- train_dev_test_split2(y2)
  train_idx <- idx$train_idx
  dev_idx <- idx$dev_idx
  test_idx <- idx$test_idx
  
  #Create the data split.
  x_train <- list(gene = x2[[1]][train_idx,], lipid = x2[[2]][train_idx,])
  x_dev <- list(gene = x2[[1]][dev_idx,], lipid = x2[[2]][dev_idx,]) 
  x_test <- list(gene = x2[[1]][test_idx,], lipid = x2[[2]][test_idx,]) 
  
  y_train <- y2[train_idx]
  y_dev <- y2[dev_idx]
  y_test <- y2[test_idx]
  
  #Create the PCs.
  x_train_pca <- pca(data.frame(x_train), ncomp = 10)
  x_dev_pca <- pca(data.frame(x_dev), ncomp = 10)
  x_test_pca <- pca(data.frame(x_test), ncomp = 10)
  
  #Run the validation grid search.
  best_macro_pcr2 <- NULL
  best_micro_pcr2 <- NULL
  best_macro_pcr2_dev <- 0
  best_micro_pcr2_dev <- 0
  best_macro_pcr2_n_comp <- 0
  best_micro_pcr2_n_comp <- 0
  
  for (i in 2:10){
    pcr2 <- multinom(y_train ~., data = data.frame(y_train,x_train_pca$variates$X[,1:i]), 
                    maxit = 250)
    pred_dev_pcr2 <- predict(pcr2, newdata=data.frame(x_dev_pca$variates$X[,1:i]))
    
    #Create confusion matrix and calculate F1 scores.
    cm_dev_pcr2 <- as.matrix(confusionMatrix(y_dev, pred_dev_pcr2))
    f1_dev_pcr2 <- f1(cm_dev_pcr2)
    if (!is.na(f1_dev_pcr2[1]) && best_macro_pcr2_dev < f1_dev_pcr2[1]){
      best_macro_pcr2_dev <- f1_dev_pcr2[1]
      best_macro_pcr2 <- pcr2
      best_macro_pcr2_n_comp <- i
    }
    if (best_micro_pcr2_dev < f1_dev_pcr2[2]){
      best_micro_pcr2_dev <- f1_dev_pcr2[2]
      best_micro_pcr2 <- pcr2
      best_micro_pcr2_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  #In some cases, macro F1 comes out NA because no predictions are made for a certain class.
  if (!is.null(best_macro_pcr2)){
    pred_test_macro_pcr2 <- predict(best_macro_pcr2, newdata = x_test_pca$variates$X[,1:best_macro_pcr2_n_comp])
    cm_test_macro_pcr2 <- as.matrix(confusionMatrix(y_test, pred_test_macro_pcr2))
    f1_test_macro_pcr2 <- f1(cm_test_macro_pcr2)[1]
  } else {
    f1_test_macro_pcr2 <- NA
  }
  
  pred_test_micro_pcr2 <- predict(best_micro_pcr2, newdata = x_test_pca$variates$X[,1:best_micro_pcr2_n_comp])
  cm_test_micro_pcr2 <- as.matrix(confusionMatrix(y_test, pred_test_micro_pcr2))
  f1_test_micro_pcr2 <- f1(cm_test_micro_pcr2)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_pcr2_best_n_comps <- c(macro_pcr2_best_n_comps, best_macro_pcr2_n_comp)
  macro_pcr2_test_scores <- c(macro_pcr2_test_scores, f1_test_macro_pcr2)
  
  micro_pcr2_best_n_comps <- c(micro_pcr2_best_n_comps, best_micro_pcr2_n_comp)
  micro_pcr2_test_scores <- c(micro_pcr2_test_scores, f1_test_micro_pcr2)
  
  #############
  
  #Create the sPCs; keep default alpha of 0.0001 in the interest of time.
  #Run the validation grid search.
  best_macro_spcr2 <- NULL
  best_micro_spcr2 <- NULL
  best_macro_spcr2_dev <- 0
  best_micro_spcr2_dev <- 0
  best_macro_spcr2_n_comp <- 0
  best_micro_spcr2_n_comp <- 0
  
  for (i in 2:10){
    #alpha value set s.t. the median number of selected variables is roughly
    #200 across components.
    x_train_spca <- spca(data.frame(x_train), k = i, alpha = 0.000015)
    x_dev_spca <- spca(data.frame(x_dev), k = i, alpha = 0.000015)
    x_test_spca <- spca(data.frame(x_test), k = i, alpha = 0.000015)
    spcr2 <- multinom(y_train ~., data = data.frame(y_train,x_train_spca$scores[,1:i]), 
                     maxit = 250)
    pred_dev_spcr2 <- predict(spcr2, newdata=data.frame(x_dev_spca$scores[,1:i]))
    
    #Create confusion matrix and calculate F1 scores.
    cm_dev_spcr2 <- as.matrix(confusionMatrix(y_dev, pred_dev_spcr2))
    f1_dev_spcr2 <- f1(cm_dev_spcr2)
    if (!is.na(f1_dev_spcr2[1]) && best_macro_spcr2_dev < f1_dev_spcr2[1]){
      best_macro_spcr2_dev <- f1_dev_spcr2[1]
      best_macro_spcr2 <- spcr2
      best_macro_spcr2_n_comp <- i
    }
    if (best_micro_spcr2_dev < f1_dev_spcr2[2]){
      best_micro_spcr2_dev <- f1_dev_spcr2[2]
      best_micro_spcr2 <- spcr2
      best_micro_spcr2_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  if (!is.null(best_macro_spcr2)){
    pred_test_macro_spcr2 <- predict(best_macro_spcr2, 
                                    newdata = data.frame(x_test_spca$scores[,1:best_macro_spcr2_n_comp]))
    cm_test_macro_spcr2 <- as.matrix(confusionMatrix(y_test, pred_test_macro_spcr2))
    f1_test_macro_spcr2 <- f1(cm_test_macro_spcr2)[1]
  } else {
    f1_test_macro_spcr2 <- NA
  }
  
  pred_test_micro_spcr2 <- predict(best_micro_spcr2, 
                                  newdata = data.frame(x_test_spca$scores[,1:best_micro_spcr2_n_comp]))
  cm_test_micro_spcr2 <- as.matrix(confusionMatrix(y_test, pred_test_micro_spcr2))
  f1_test_micro_spcr2 <- f1(cm_test_micro_spcr2)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_spcr2_best_n_comps <- c(macro_spcr2_best_n_comps, best_macro_spcr2_n_comp)
  macro_spcr2_test_scores <- c(macro_spcr2_test_scores, f1_test_macro_spcr2)
  
  micro_spcr2_best_n_comps <- c(micro_spcr2_best_n_comps, best_micro_spcr2_n_comp)
  micro_spcr2_test_scores <- c(micro_spcr2_test_scores, f1_test_micro_spcr2)
  
  best_macro_plsda2 <- NULL
  best_micro_plsda2 <- NULL
  best_macro_plsda2_dev <- 0
  best_micro_plsda2_dev <- 0
  best_macro_plsda2_n_comp <- 0
  best_micro_plsda2_n_comp <- 0
  
  for (i in 2:10){
    plsda2 <- mixOmics::block.plsda(x_train, y_train, ncomp = i, max.iter = 250)
    pred_dev_plsda2 <- predict(plsda2, x_dev)
    #Create confusion matrix and calculate F1 scores.
    cm_dev_plsda2 <- as.matrix(confusionMatrix(y_dev, 
                                              factor(pred_dev_plsda2$WeightedPredict.class$max.dist[,i],
                                                     levels = levels(y_dev))))
    f1_dev_plsda2 <- f1(cm_dev_plsda2)
    if (!is.na(f1_dev_plsda2[1]) && best_macro_plsda2_dev < f1_dev_plsda2[1]){
      best_macro_plsda2_dev <- f1_dev_plsda2[1]
      best_macro_plsda2 <- plsda2
      best_macro_plsda2_n_comp <- i
    }
    if (best_micro_plsda2_dev < f1_dev_plsda2[2]){
      best_micro_plsda2_dev <- f1_dev_plsda2[2]
      best_micro_plsda2 <- plsda2
      best_micro_plsda2_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  if (!is.null(best_macro_plsda2)){
    pred_test_macro_plsda2 <- predict(best_macro_plsda2, newdata = x_test)
    temp <- pred_test_macro_plsda2$WeightedPredict.class$max.dist 
    cm_test_macro_plsda2 <- as.matrix(confusionMatrix(y_test, 
                                                     factor(temp[,ncol(temp)],
                                                            levels = levels(y_test))))
    f1_test_macro_plsda2 <- f1(cm_test_macro_plsda2)[1]
  } else {
    f1_test_macro_plsda2 <- NA
  }
  
  pred_test_micro_plsda2 <- predict(best_micro_plsda2, newdata = x_test)
  temp <- pred_test_micro_plsda2$WeightedPredict.class$max.dist
  cm_test_micro_plsda2 <- as.matrix(confusionMatrix(y_test, 
                                                   factor(temp[,ncol(temp)],
                                                          levels = levels(y_test))))
  f1_test_micro_plsda2 <- f1(cm_test_micro_plsda2)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_plsda2_best_n_comps <- c(macro_plsda2_best_n_comps, best_macro_plsda2_n_comp)
  macro_plsda2_test_scores <- c(macro_plsda2_test_scores, f1_test_macro_plsda2)
  
  micro_plsda2_best_n_comps <- c(micro_plsda2_best_n_comps, best_micro_plsda2_n_comp)
  micro_plsda2_test_scores <- c(micro_plsda2_test_scores, f1_test_micro_plsda2)
  
  best_macro_splsda2 <- NULL
  best_micro_splsda2 <- NULL
  best_macro_splsda2_dev <- 0
  best_micro_splsda2_dev <- 0
  best_macro_splsda2_n_comp <- 0
  best_micro_splsda2_n_comp <- 0
  
  for (i in 2:10){
    #Keeping 200 variables for component 1, roughly in line with the approach for sPCA.
    list.keepX = list(gene = rep(10, 2), lipid = rep(5,2))
    splsda2 <- mixOmics::block.splsda(x_train, y_train, keepX = list.keepX,
                               ncomp = i, max.iter = 250)
    pred_dev_splsda2 <- predict(splsda2, x_dev)
    #Create confusion matrix and calculate F1 scores.
    cm_dev_splsda2 <- as.matrix(confusionMatrix(y_dev, 
                                               factor(pred_dev_splsda2$WeightedPredict.class$max.dist[,i],
                                                      levels = levels(y_dev))))
    f1_dev_splsda2 <- f1(cm_dev_splsda2)
    if (!is.na(f1_dev_splsda2[1]) && best_macro_splsda2_dev < f1_dev_splsda2[1]){
      best_macro_splsda2_dev <- f1_dev_splsda2[1]
      best_macro_splsda2 <- splsda2
      best_macro_splsda2_n_comp <- i
    }
    if (best_micro_splsda2_dev < f1_dev_splsda2[2]){
      best_micro_splsda2_dev <- f1_dev_splsda2[2]
      best_micro_splsda2 <- splsda2
      best_micro_splsda2_n_comp <- i
    }  
  }
  
  #Run the best models on the test set.
  if (!is.null(best_macro_splsda2)){
    pred_test_macro_splsda2 <- predict(best_macro_splsda2, newdata = x_test)
    temp <- pred_test_macro_splsda2$WeightedPredict.class$max.dist 
    cm_test_macro_splsda2 <- as.matrix(confusionMatrix(y_test, 
                                                      factor(temp[,ncol(temp)],
                                                             levels = levels(y_test))))
    f1_test_macro_splsda2 <- f1(cm_test_macro_splsda2)[1]
  } else {
    f1_test_macro_splsda2 <- NA
  }
  
  pred_test_micro_splsda2 <- predict(best_micro_splsda2, newdata = x_test)
  temp <- pred_test_micro_splsda2$WeightedPredict.class$max.dist 
  cm_test_micro_splsda2 <- as.matrix(confusionMatrix(y_test, 
                                                     factor(temp[,ncol(temp)],
                                                            levels = levels(y_test))))
  f1_test_micro_splsda2 <- f1(cm_test_micro_splsda2)[2]
  
  #Pick out the best hyper-parameter and test score at each turn.
  macro_splsda2_best_n_comps <- c(macro_splsda2_best_n_comps, best_macro_splsda2_n_comp)
  macro_splsda2_test_scores <- c(macro_splsda2_test_scores, f1_test_macro_splsda2)
  
  micro_splsda2_best_n_comps <- c(micro_splsda2_best_n_comps, best_micro_splsda2_n_comp)
  micro_splsda2_test_scores <- c(micro_splsda2_test_scores, f1_test_micro_splsda2)
}

#Check best stats from tests
get_mode(macro_pcr2_best_n_comps[macro_pcr2_best_n_comps!=0])
get_mode(micro_pcr2_best_n_comps)
mean(macro_pcr2_best_n_comps[macro_pcr2_best_n_comps!=0])
mean(micro_pcr2_best_n_comps)
sd(macro_pcr2_best_n_comps[macro_pcr2_best_n_comps!=0])
sd(micro_pcr2_best_n_comps)
mean(macro_pcr2_test_scores, na.rm = TRUE)
mean(micro_pcr2_test_scores, na.rm = TRUE)
sd(macro_pcr2_test_scores, na.rm = TRUE)
sd(micro_pcr2_test_scores, na.rm = TRUE)

get_mode(macro_spcr2_best_n_comps[macro_spcr2_best_n_comps!=0])
get_mode(micro_spcr2_best_n_comps)
mean(macro_spcr2_best_n_comps[macro_spcr2_best_n_comps!=0])
mean(micro_spcr2_best_n_comps)
sd(macro_spcr2_best_n_comps[macro_spcr2_best_n_comps!=0])
sd(micro_spcr2_best_n_comps)
mean(macro_spcr2_test_scores, na.rm = TRUE)
mean(micro_spcr2_test_scores, na.rm = TRUE)
sd(macro_spcr2_test_scores, na.rm = TRUE)
sd(micro_spcr2_test_scores, na.rm = TRUE)

get_mode(macro_plsda2_best_n_comps)
get_mode(micro_plsda2_best_n_comps)
mean(macro_plsda2_best_n_comps)
mean(micro_plsda2_best_n_comps)
sd(macro_plsda2_best_n_comps)
sd(micro_plsda2_best_n_comps)
mean(macro_plsda2_test_scores, na.rm = TRUE)
mean(micro_plsda2_test_scores, na.rm = TRUE)
sd(macro_plsda2_test_scores, na.rm = TRUE)
sd(micro_plsda2_test_scores, na.rm = TRUE)

get_mode(macro_splsda2_best_n_comps)
get_mode(micro_splsda2_best_n_comps)
mean(macro_splsda2_best_n_comps)
mean(micro_splsda2_best_n_comps)
sd(macro_splsda2_best_n_comps)
sd(micro_splsda2_best_n_comps)
mean(macro_splsda2_test_scores, na.rm = TRUE)
mean(micro_splsda2_test_scores, na.rm = TRUE)
sd(macro_splsda2_test_scores, na.rm = TRUE)
sd(micro_splsda2_test_scores, na.rm = TRUE)