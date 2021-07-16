rm(list = ls())
library(rsample)
library(iRF)
library(randomForest)

load("~/Desktop/RData/Omics_genus.RData")
load("~/Desktop/Rdata/omics_genus_Class.Rdata")
omics <- as.data.frame(omics)
#sampling
set.seed(70)
X <- list()
for (time in 1:20){
  X[[time]] <- initial_split(omics, prop = 0.8, strata = "Class")
}

#sun iRF
set.seed(70)
irf <- list()
for (time in 1:20){
  irf[[time]] <- iRF(x=X_omics[X[[time]]$in_id,], y=Y_omics[X[[time]]$in_id], 
                     xtest=X_omics[-X[[time]]$in_id,], ytest=Y_omics[-X[[time]]$in_id], 
                     n.iter=5,
                     ntree=500,
                     n.core=1,
                     n.bootstrap=30,
                     interactions.return=5,
                     bootstrap.forest=T,
                     rit.param=list(depth=5, ntree=500, nchild=2, 
                                    class.id=1, class.cut=NULL))
  print(head(irf[[time]]$interaction[5]))
}

#check interactions
for (time in 1:20){ 
  
  toplot_omics <- irf[[time]]$interaction[[5]]
  dotchart(rev(toplot_omics[1:min(20, length(toplot_omics))]), 
           xlab='Prevalence enrichment', xlim=c(0, 1),
           main='Prevalent Features/Interactions \n on Decision paths')
  
}

#ROC curve
library(AUC)
plot(0:1, 0:1, type='l', lty = 2, xlab = 'FPR', ylab = 'TPR', main='ROC Curve')
for (iter in 1:5){
  # performance on test set
  cat(paste('iter = ', iter, ':: '))
  roc.info <- roc(irf[[11]]$rf.list[[iter]]$test$votes[,2], Y_omics[-X[[11]]$in_id])
  lines(roc.info$fpr, roc.info$tpr, type='l', col=iter, lwd=2)
  cat(paste('AUROC: ', round(100*auc(roc.info), 2), '%\n', sep=''))
} 
legend('bottomright', legend=paste('iter:', 1:iter), col=1:iter, lwd=2, bty='n')


#PR curve
library(ROCR)
plot(1:0, xlim=c(0,1), ylim=c(0.5,1), col = "white",
     xlab = 'Recall', ylab = 'Precision', main='PR Curve')
for (iter in 1:5){
  # performance on test set
  cat(paste('iter = ', iter, ':: '))
  pred <- prediction(as.numeric(irf[[18]]$rf.list[[iter]]$test$votes[,2]), 
                     Y_omics[-X[[18]]$in_id])
  perf <- performance(pred,"prec","rec")
  lines(as.numeric(unlist(perf@x.values)), as.numeric(unlist(perf@y.values)), 
        type='l',col=iter,lwd=2)
}
legend('bottomleft', legend=paste('iter:', 1:iter), col=1:iter, lwd=2, bty='n')

#variable importance
varImpPlot(irf[[11]]$rf.list[[1]],n.var=15,
           main='Variable Importance (iter:1)')
varImpPlot(irf[[11]]$rf.list[[2]],n.var=15,
           main='Variable Importance (iter:2)')
varImpPlot(irf[[11]]$rf.list[[3]],n.var=15,
           main='Variable Importance (iter:3)')
varImpPlot(irf[[11]]$rf.list[[4]],n.var=15,
           main='Variable Importance (iter:4)')
varImpPlot(irf[[11]]$rf.list[[5]],n.var=15,
           main='Variable Importance (iter:5)')

save(irf,X,file = "irf20times.RData")


