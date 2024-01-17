##机器学习，如下：
###基于Lasso、Ridge、Enet、StepCox、survivalSVM、CoxBoost、SuperPC、plsRcox、RSF、GBM这10种机器学习模型，两两相互交联使用，寻找关键基因并构建预后模型，评估C-index
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(dplyr)
library(tibble)
library(BART)
library(caret)
library(glmnet)
library(pec)
library(UpSetR)
library(limma)
library(Matrix)
library(data.table)

my<-read.table("surSigExp.txt",header=T,row.names=1)

rm(result)
result <- data.frame()
set.seed(1234)
seed=123
nfold=10 #启用10折交叉验证

#全10折交叉验证

#### 1.RSF ####

rf_nodesize <- 5
fold_indices <- sample(1:nfold, size = nrow(my), replace = TRUE)
c_index_values <- numeric(nfold)
gene_importance <- rep(0, ncol(my) - 2)
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
  # 拟合 RSF 模型
fit <- rfsrc(Surv(OS.time,OS)~.,data = train_data,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整,nodesize决定树深  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
  gene_importance_fold <- fit$importance
  gene_importance <- gene_importance + gene_importance_fold# 累加基因重要性
val_data_list <- list(test_data)
est_data <- test_data
pre_var <- colnames(test_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',pre_var)]})
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%rownames_to_column('ID')
c_index_values[fold] <- cc$Cindex
}
average_gene_importance <- gene_importance / nfold
fit$importance<-average_gene_importance
attr(fit$importance,"names") <-colnames(my[-c(1:2)])
rftop<-data.frame(Feature=var.select(fit)$topvars,vimp=var.select(fit)$varselect[var.select(fit)$topvars,2])

mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- 'RSF'
result <- rbind(result,mean_c_index)


keygene<-rftop$Feature
my2 <- cbind(my[,1:2],my[,-c(1,2)][,keygene])
val_dd_list<-list(my2)

est_data <- my2
pre_var <- colnames(my2)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]


#### 1.1.RSF+Lasso ####

fold_indices <- sample(1:nfold, size = nrow(my2), replace = TRUE)
c_index <- numeric(nfold)
la<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
  # 拟合 RSF 模型
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
cvfit = cv.glmnet(x1, x2,nfold=10,family = "cox",alpha=1,  type.measure = 'C')#10折
best_model <- glmnet(x1,x2, alpha = 1, lambda = cvfit$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
c_index[fold] <- cc$Cindex
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
geneids <- rownames(coef.min)[active.min]
index.min = coef.min[active.min]
combine <- cbind(geneids, index.min)
lasso.result<-as.data.frame(combine)
lasso.result$index.min<-as.numeric(lasso.result$index.min)
la[fold]<-lasso.result
}
mean_c_index <- data.frame(mean(c_index))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('RSF+','Lasso')
result <- rbind(result,mean_c_index)

#### 1.2.RSF+Ridge ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
ridge_model <- cv.glmnet(x1, x2, alpha = 0, family = "cox", nfolds = 10,type.measure = 'C')
best_model <- glmnet(x1,x2, alpha = 0, lambda = ridge_model$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('RSF+','Ridge')
result <- rbind(result,mean_c_index)



#### 1.3.RSF+Enet ####

result2<-data.frame()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
for (alpha in seq(0.1,0.9,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10,type.measure = 'C')
best_model <- glmnet(x1,x2, alpha = alpha, lambda = ridge_model$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
colnames(cc)<-'C'
  cc$Model <- paste0('RSF+','Enet','[α=',alpha,']')
  result2 <- rbind(result2,cc)
}}

split_data<-split(result2,result2$Model)
avg_values<-c()
for(i in 1:length(split_data)){avg_values[i]<-mean(split_data[[i]]$C)}
index<-as.data.frame(avg_values)
for (i in 1:9){index[i,2]<-names(split_data)[i]}
colnames(index)<-c("C-index","Model")
result <- rbind(result,index)

#### 1.4.RSF+StepCox ####

ccc<-list()
cccc<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,train_data),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values <- cc$Cindex
cc$Model<-paste0('rLasso+','StepCox','[',direction,']')
ccc[direction]<-cc
}
cccc[[fold]]<-ccc
}
ccx<-data.frame()
ccx<-as.data.frame(rbindlist( cccc, fill = FALSE, idcol = NULL))
ccccx<-data.frame()
for(i in 1:3){ccccx[i,1]<-sum(ccx[i])/nfold}
ccccx[1,2]<-paste0('RSF+','StepCox','[','both',']')
ccccx[2,2]<-paste0('RSF+','StepCox','[','backward',']')
ccccx[3,2]<-paste0('RSF+','StepCox','[','forward',']')
colnames(ccccx)<-c('C-index','Model')
result<-rbind(result,ccccx)

#### 1.5.RSF+GBM ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)
fit <- gbm(formula = Surv(OS.time,OS)~.,data =train_data,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('RSF+','GBM')
result <- rbind(result,mean_c_index)


#### 1.6.RSF+survival-SVM ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)

fit = survivalsvm(Surv(OS.time,OS)~., data= train_data, gamma.mu = 1)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('RSF+','survival-SVM')
result <- rbind(result,mean_c_index)


#### 1.7.RSF+superPC ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)

data<-list(x=t(train_data[,-c(1,2)]),y=train_data$OS.time,censoring.status=train_data$OS,featurenames=colnames(train_data)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs<-lapply(val_dd_list,function(w){test<-list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
ff<-superpc.predict(fit,data,test,threshold=cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('RSF+','superPC')
result <- rbind(result,mean_c_index)


#### 1.8.RSF+plsRcox ####
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_data_list <- list(test_data)
est_data <- train_data
pre_var <- colnames(train_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(list(train_data),function(x){x[,c('OS.time','OS',pre_var)]})
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,pre_var],time=train_data$OS.time,status=train_data$OS),nt=10,verbose = FALSE)
fit<-plsRcox(train_data[,pre_var],time=train_data$OS.time,event=train_data$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs<-lapply(list(test_data),function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('RSF+','plsRcox')
result <- rbind(result,mean_c_index)



#### 1.9.RSF+coxBoost ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]

est_data <- train_data
pre_var <- colnames(train_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(list(train_data),function(x){x[,c('OS.time','OS',pre_var)]})

set.seed(seed)
pen<-optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(list(test_data),function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('RSF+','coxBoost ')
result <- rbind(result,mean_c_index)

#### 2. Lasso ####

nfold=10
fold_indices <- sample(1:nfold, size = nrow(my), replace = TRUE)
c_index <- numeric(nfold)
la<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
  # 拟合 RSF 模型
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
alpha=1
cvfit = cv.glmnet(x1, x2,nfold=10,family = "cox",alpha=alpha,  type.measure = 'C')#10折
best_model <- glmnet(x1,x2, alpha = alpha, lambda = cvfit$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
c_index[fold] <- cc$Cindex
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
geneids <- rownames(coef.min)[active.min]
index.min = coef.min[active.min]
combine <- cbind(geneids, index.min)
lasso.result<-as.data.frame(combine)
lasso.result$index.min<-as.numeric(lasso.result$index.min)
la[fold]<-lasso.result
}
mean_c_index <- data.frame(mean(c_index))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- 'Lasso'
union_set <- Reduce(union, la)#intersect则是并集
intersect_set <- Reduce(intersect, la)
result <- rbind(result,mean_c_index)

keygene<-union_set
my2 <- cbind(my[,1:2],my[,-c(1,2)][,keygene])
val_dd_list<-list(my2)

est_data <- my2
pre_var <- colnames(my2)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]

result
#### 2. Lasso for union_set (rLasso)####

c_index_lasso <- data.frame()
x1 <- as.matrix(my2[, 3:ncol(my2)])
x2 <- as.matrix(Surv(my2$OS.time, my2$OS))
alpha=1
cvfit = cv.glmnet(x1, x2,nfold=10,family = "cox",alpha=alpha,  type.measure = 'C')#10折
best_model <- glmnet(x1,x2, alpha = alpha, lambda = cvfit$lambda.min, family = "cox")
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
geneids <- rownames(coef.min)[active.min]
index.min = coef.min[active.min]
combine <- cbind(geneids, index.min)
lasso.result<-as.data.frame(combine)
lasso.result$index.min<-as.numeric(lasso.result$index.min)
lasso.result
rs<-cbind(my2[,1:2],RS=predict(cvfit, newx = as.matrix(my2[,-c(1:2)]), s = cvfit$lambda.min))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
c_index_lasso <- as.data.frame(cc$Cindex)
colnames(c_index_lasso)<-'C-index'
c_index_lasso$Model <- 'rLasso'

keygene<-lasso.result$geneids
result <- rbind(result,c_index_lasso)

my2 <- cbind(my[,1:2],my[,-c(1,2)][,keygene])
val_dd_list<-list(my2)

est_data <- my2
pre_var <- colnames(my2)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]

ggplot(lasso.result, aes(x = geneids, y = index.min)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Gene ID", y = "Coefficient", title = "Weights of Key Genes in Lasso Model") +
  coord_flip()  # Flips the axes for a horizontal plot


#### 2.1 rLasso+RSF ####

fold_indices <- sample(1:nfold, size = nrow(my2), replace = TRUE)
c_index_values <- numeric(nfold)
gene_importance <- rep(0, ncol(my2) - 2)# -2 是因为去除了 OS.time 和 OS
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
  # 拟合 RSF 模型
fit <- rfsrc(Surv(OS.time,OS)~.,data = train_data,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整,nodesize决定树深  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
  gene_importance_fold <- fit$importance# 提取基因重要性
  gene_importance <- gene_importance + gene_importance_fold# 累加基因重要性
val_data_list <- list(test_data)
est_data <- test_data
pre_var <- colnames(test_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(val_data_list,function(x){x[,c('OS.time','OS',pre_var)]})
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%rownames_to_column('ID')
c_index_values[fold] <- cc$Cindex
}
average_gene_importance <- gene_importance / nfold
fit$importance<-average_gene_importance
attr(fit$importance,"names") <-colnames(my2[-c(1:2)])
rftop<-data.frame(Feature=var.select(fit)$topvars,vimp=var.select(fit)$varselect[var.select(fit)$topvars,2])
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('rLasso+','RSF')
result <- rbind(result,mean_c_index)

#### 2.2 rLasso+Ridge ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
ridge_model <- cv.glmnet(x1, x2, alpha = 0, family = "cox", nfolds = 10,type.measure = 'C')

best_model <- glmnet(x1,x2, alpha = 0, lambda = ridge_model$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('rLasso+','Ridge')
result <- rbind(result,mean_c_index)

#### 2.3 rLasso+Enet ####

result2<-data.frame()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
for (alpha in seq(0.1,0.9,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10,type.measure = 'C')
best_model <- glmnet(x1,x2, alpha = alpha, lambda = ridge_model$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
colnames(cc)<-'C'
  cc$Model <- paste0('rLasso+','Enet','[α=',alpha,']')
  result2 <- rbind(result2,cc)
}}

split_data<-split(result2,result2$Model)
avg_values<-c()
for(i in 1:length(split_data)){avg_values[i]<-mean(split_data[[i]]$C)}
index<-as.data.frame(avg_values)
for (i in 1:9){index[i,2]<-names(split_data)[i]}
colnames(index)<-c("C-index","Model")
result <- rbind(result,index)


#### 2.4 rLasso+StepCox ####

ccc<-list()
cccc<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,train_data),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values <- cc$Cindex
cc$Model<-paste0('rLasso+','StepCox','[',direction,']')
ccc[direction]<-cc}
cccc[[fold]]<-ccc}
ccx<-data.frame()
ccx<-as.data.frame(rbindlist( cccc, fill = FALSE, idcol = NULL))
ccccx<-data.frame()
for(i in 1:3){ccccx[i,1]<-sum(ccx[i])/nfold}
ccccx[1,2]<-paste0('rLasso+','StepCox','[','both',']')
ccccx[2,2]<-paste0('rLasso+','StepCox','[','backward',']')
ccccx[3,2]<-paste0('rLasso+','StepCox','[','forward',']')
colnames(ccccx)<-c('C-index','Model')
result<-rbind(result,ccccx)

#### 2.5 rLasso+GBM ####


for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)
fit <- gbm(formula = Surv(OS.time,OS)~.,data =train_data,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('rLasso+','GBM')
result <- rbind(result,mean_c_index)


#### 2.6 rLasso+survival-SVM ####


for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)

fit = survivalsvm(Surv(OS.time,OS)~., data= train_data, gamma.mu = 1)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('rLasso+','survival-SVM')
result <- rbind(result,mean_c_index)

#### 2.7.rLasso+superPC ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)

data<-list(x=t(train_data[,-c(1,2)]),y=train_data$OS.time,censoring.status=train_data$OS,featurenames=colnames(train_data)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs<-lapply(val_dd_list,function(w){test<-list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
ff<-superpc.predict(fit,data,test,threshold=cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('rLasso+','superPC')
result <- rbind(result,mean_c_index)

#### 2.8.rLasso+plsRcox ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_data_list <- list(test_data)
est_data <- train_data
pre_var <- colnames(train_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(list(train_data),function(x){x[,c('OS.time','OS',pre_var)]})
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,pre_var],time=train_data$OS.time,status=train_data$OS),nt=10,verbose = FALSE)
fit<-plsRcox(train_data[,pre_var],time=train_data$OS.time,event=train_data$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs<-lapply(list(test_data),function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('rLasso+','plsRcox')
result <- rbind(result,mean_c_index)

#### 2.9.rLasso+coxBoost####


for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]

est_data <- train_data
pre_var <- colnames(train_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(list(train_data),function(x){x[,c('OS.time','OS',pre_var)]})

set.seed(seed)
pen<-optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(list(test_data),function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('rLasso+','coxBoost ')
result <- rbind(result,mean_c_index)



#### 3.Enet ####

rm(ccc)
rm(cccc)
nfold=10
ccc<-list()
cccc<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
result2<-data.frame()
for (alpha in seq(0.1,0.9,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10,type.measure = 'C')
rs<-lapply(list(test_data),function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type='link',newx=as.matrix(x[,-c(1,2)]),s=fit$lambda.min)))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
cc$Model<-paste0('Enet','[α=',alpha,']')
ccc[alpha*10]<-cc
}
cccc[[fold]]<-ccc
}
ccx<-as.data.frame(rbindlist( cccc, fill = FALSE, idcol = NULL))
ccccx<-data.frame()
for(i in 1:9){ccccx[i,1]<-sum(ccx[i])/nfold}
for(i in 1:9){ccccx[i,2]<-paste0('Enet','[α=',i/10,']')}
colnames(ccccx)<-c('C-index','Model')
result<-rbind(result,ccccx)


#### 4.Ridge ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
ridge_model <- cv.glmnet(x1, x2, alpha = 0, family = "cox", nfolds = 10,type.measure = 'C')
best_model <- glmnet(x1,x2, alpha = 0, lambda = ridge_model$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('Ridge')
result <- rbind(result,mean_c_index)

#### 5.GBM ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
val_dd_list<-list(test_data)
fit <- gbm(formula = Surv(OS.time,OS)~.,data =train_data,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('GBM')
result <- rbind(result,mean_c_index)

#### 6.Stepwise-Cox ####这个真的慢

ccc<-list()
cccc<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
val_dd_list<-list(test_data)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,train_data),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values <- cc$Cindex
cc$Model<-paste0('rLasso+','StepCox','[',direction,']')
ccc[direction]<-cc
}
cccc[[fold]]<-ccc
}
ccx<-data.frame()
ccx<-as.data.frame(rbindlist( cccc, fill = FALSE, idcol = NULL))
ccccx<-data.frame()
for(i in 1:3){ccccx[i,1]<-sum(ccx[i])/nfold}
ccccx[1,2]<-paste0('StepCox','[','both',']')
ccccx[2,2]<-paste0('StepCox','[','backward',']')
ccccx[3,2]<-paste0('StepCox','[','forward',']')
colnames(ccccx)<-c('C-index','Model')
result<-rbind(result,ccccx)

#### 7.survival-SVM ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
val_dd_list<-list(test_data)

fit = survivalsvm(Surv(OS.time,OS)~., data= train_data, gamma.mu = 1)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('survival-SVM')
result <- rbind(result,mean_c_index)

#### 8.superPC ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
val_dd_list<-list(test_data)

data<-list(x=t(train_data[,-c(1,2)]),y=train_data$OS.time,censoring.status=train_data$OS,featurenames=colnames(train_data)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs<-lapply(val_dd_list,function(w){test<-list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
ff<-superpc.predict(fit,data,test,threshold=cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('superPC')
result <- rbind(result,mean_c_index)



#### 9.plsRcox ####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]
val_data_list <- list(test_data)
est_data <- train_data
pre_var <- colnames(train_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(list(train_data),function(x){x[,c('OS.time','OS',pre_var)]})
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,pre_var],time=train_data$OS.time,status=train_data$OS),nt=10,verbose = FALSE)
fit<-plsRcox(train_data[,pre_var],time=train_data$OS.time,event=train_data$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs<-lapply(list(test_data),function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('plsRcox')
result <- rbind(result,mean_c_index)

#### 10.coxBoost ####

la<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my[fold_indices != fold, ]
  test_data <- my[fold_indices == fold, ]

est_data <- train_data
pre_var <- colnames(train_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(list(train_data),function(x){x[,c('OS.time','OS',pre_var)]})

set.seed(seed)
pen<-optimCoxBoostPenalty(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),trace=TRUE,start.penalty=500,parallel = T)
cv.res <- cv.CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                      maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
fit <- CoxBoost(est_dd[,'OS.time'],est_dd[,'OS'],as.matrix(est_dd[,-c(1,2)]),
                stepno=cv.res$optimal.step,penalty=pen$penalty)
rs <- lapply(list(test_data),function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,newdata=x[,-c(1,2)], newtime=x[,1], newstatus=x[,2], type="lp")))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
coefficients_matrix <- fit$coefficients
non_zero_counts_per_column <- diff(coefficients_matrix@p)
non_zero_columns <- which(non_zero_counts_per_column > 0)
important_genes <- fit$xnames[non_zero_columns]
la[fold]<-as.data.frame(important_genes)
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('coxBoost ')
result <- rbind(result,mean_c_index)

union_set <- Reduce(union, la)#intersect则是并集
keygene<-union_set
my2 <- cbind(my[,1:2],my[,-c(1,2)][,keygene])
val_dd_list<-list(my2)

est_data <- my2
pre_var <- colnames(my2)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]

#### 10.1.coxBoost+RSF ####

fold_indices <- sample(1:nfold, size = nrow(my2), replace = TRUE)
c_index_values <- numeric(nfold)
gene_importance <- rep(0, ncol(my2) - 2)# -2 是因为去除了 OS.time 和 OS
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
  # 拟合 RSF 模型
fit <- rfsrc(Surv(OS.time,OS)~.,data = train_data,
             ntree = 1000,nodesize = rf_nodesize,##该值建议多调整,nodesize决定树深  
             splitrule = 'logrank',
             importance = T,
             proximity = T,
             forest = T,
             seed = seed)
  gene_importance_fold <- fit$importance# 提取基因重要性
  gene_importance <- gene_importance + gene_importance_fold# 累加基因重要性
rs <- lapply(list(test_data),function(x){cbind(x[,1:2],RS=predict(fit,newdata = x)$predicted)})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))%>%rownames_to_column('ID')
c_index_values[fold] <- cc$Cindex
}
average_gene_importance <- gene_importance / nfold
fit$importance<-average_gene_importance
attr(fit$importance,"names") <-colnames(my2[-c(1:2)])
rftop<-data.frame(Feature=var.select(fit)$topvars,vimp=var.select(fit)$varselect[var.select(fit)$topvars,2])
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('coxBoost+','RSF')
result <- rbind(result,mean_c_index)


#### 10.2.coxBoost+Rideg####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
ridge_model <- cv.glmnet(x1, x2, alpha = 0, family = "cox", nfolds = 10,type.measure = 'C')
best_model <- glmnet(x1,x2, alpha = 0, lambda = ridge_model$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('coxBoost+','Ridge')
result <- rbind(result,mean_c_index)

#### 10.3.coxBoost+Enet####

result2<-data.frame()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
for (alpha in seq(0.1,0.9,0.1)) {
  set.seed(seed)
  fit = cv.glmnet(x1, x2,family = "cox",alpha=alpha,nfolds = 10,type.measure = 'C')
best_model <- glmnet(x1,x2, alpha = alpha, lambda = ridge_model$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
colnames(cc)<-'C'
  cc$Model <- paste0('coxBoost+','Enet','[α=',alpha,']')
  result2 <- rbind(result2,cc)
}}

split_data<-split(result2,result2$Model)
avg_values<-c()
for(i in 1:length(split_data)){avg_values[i]<-mean(split_data[[i]]$C)}
index<-as.data.frame(avg_values)
for (i in 1:9){index[i,2]<-names(split_data)[i]}
colnames(index)<-c("C-index","Model")
result <- rbind(result,index)


#### 10.4.coxBoost+StepCox####

ccc<-list()
cccc<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)
for (direction in c("both", "backward", "forward")) {
  fit <- step(coxph(Surv(OS.time,OS)~.,train_data),direction = direction)
  rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=predict(fit,type = 'risk',newdata = x))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values <- cc$Cindex
cc$Model<-paste0('rLasso+','StepCox','[',direction,']')
ccc[direction]<-cc
}
cccc[[fold]]<-ccc
}
ccx<-data.frame()
ccx<-as.data.frame(rbindlist( cccc, fill = FALSE, idcol = NULL))
ccccx<-data.frame()
for(i in 1:3){ccccx[i,1]<-sum(ccx[i])/nfold}
ccccx[1,2]<-paste0('coxBoost+','StepCox','[','both',']')
ccccx[2,2]<-paste0('coxBoost+','StepCox','[','backward',']')
ccccx[3,2]<-paste0('coxBoost+','StepCox','[','forward',']')
colnames(ccccx)<-c('C-index','Model')
result<-rbind(result,ccccx)

#### 10.5.coxBoost+GBM####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)
fit <- gbm(formula = Surv(OS.time,OS)~.,data =train_data,distribution = 'coxph',
           n.trees = 10000,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 6)
# find index for number trees with minimum CV error
best <- which.min(fit$cv.error)
set.seed(seed)
fit <- gbm(formula = Surv(OS.time,OS)~.,data = train_data,distribution = 'coxph',
           n.trees = best,
           interaction.depth = 3,
           n.minobsinnode = 10,
           shrinkage = 0.001,
           cv.folds = 10,n.cores = 8)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,x,n.trees = best,type = 'link')))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('coxBoost+','GBM')
result <- rbind(result,mean_c_index)

#### 10.6.coxBoost+survival-SVM####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)

fit = survivalsvm(Surv(OS.time,OS)~., data= train_data, gamma.mu = 1)
rs <- lapply(val_dd_list,function(x){cbind(x[,1:2],RS=as.numeric(predict(fit, x)$predicted))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('coxBoost+','survival-SVM')
result <- rbind(result,mean_c_index)

#### 10.7.coxBoost+superPC####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_dd_list<-list(test_data)

data<-list(x=t(train_data[,-c(1,2)]),y=train_data$OS.time,censoring.status=train_data$OS,featurenames=colnames(train_data)[-c(1,2)])
set.seed(seed)
fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default 
                     n.fold = 10,
                     n.components=3,
                     min.features=5,
                     max.features=nrow(data$x),
                     compute.fullcv= TRUE,
                     compute.preval=TRUE)
rs<-lapply(val_dd_list,function(w){test<-list(x=t(w[,-c(1,2)]),y=w$OS.time,censoring.status=w$OS,featurenames=colnames(w)[-c(1,2)])
ff<-superpc.predict(fit,data,test,threshold=cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)
  rr <- as.numeric(ff$v.pred)
  rr2 <- cbind(w[,1:2],RS=rr)
  return(rr2)
})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('coxBoost+','superPC')
result <- rbind(result,mean_c_index)

#### 10.8.coxBoost+plsRcox####

for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
val_data_list <- list(test_data)
est_data <- train_data
pre_var <- colnames(train_data)[-c(1:2)]
est_dd <- est_data[,c('OS.time','OS',pre_var)]
val_dd_list <- lapply(list(train_data),function(x){x[,c('OS.time','OS',pre_var)]})
set.seed(seed)
cv.plsRcox.res=cv.plsRcox(list(x=train_data[,pre_var],time=train_data$OS.time,status=train_data$OS),nt=10,verbose = FALSE)
fit<-plsRcox(train_data[,pre_var],time=train_data$OS.time,event=train_data$OS,nt=as.numeric(cv.plsRcox.res[5]))
rs<-lapply(list(test_data),function(x){cbind(x[,1:2],RS=as.numeric(predict(fit,type="lp",newdata=x[,-c(1,2)])))})
cc<-data.frame(Cindex=sapply(rs,function(x){as.numeric(summary(coxph(Surv(OS.time,OS)~RS,x))$concordance[1])}))
c_index_values[fold] <- cc$Cindex
}
mean_c_index <- data.frame(mean(c_index_values))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('coxBoost+','plsRcox')
result <- rbind(result,mean_c_index)


#### 10.9.coxBoost+lasso####

nfold=10
fold_indices <- sample(1:nfold, size = nrow(my2), replace = TRUE)
c_index <- numeric(nfold)
la<-list()
for (fold in 1:nfold) {
  # 划分训练集和测试集
  train_data <- my2[fold_indices != fold, ]
  test_data <- my2[fold_indices == fold, ]
  # 拟合 RSF 模型
x1 <- as.matrix(train_data[, 3:ncol(train_data)])
x2 <- as.matrix(Surv(train_data$OS.time, train_data$OS))
alpha=1
cvfit = cv.glmnet(x1, x2,nfold=10,family = "cox",alpha=alpha,  type.measure = 'C')#10折
best_model <- glmnet(x1,x2, alpha = alpha, lambda = cvfit$lambda.min, family = "cox")
rs<-cbind(test_data[,1:2],RS=predict(best_model, newx = as.matrix(test_data[,-c(1:2)]), s = best_model$lambda))
colnames(rs)<-c("OS.time","OS","RS")
cc<-data.frame(Cindex=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
c_index[fold] <- cc$Cindex
coef.min = coef(cvfit, s = "lambda.min")
active.min = which(coef.min != 0)
geneids <- rownames(coef.min)[active.min]
index.min = coef.min[active.min]
combine <- cbind(geneids, index.min)
lasso.result<-as.data.frame(combine)
lasso.result$index.min<-as.numeric(lasso.result$index.min)
la[fold]<-lasso.result
}
mean_c_index <- data.frame(mean(c_index))
colnames(mean_c_index)<-'C-index'
mean_c_index$Model <- paste0('coxBoost+','Lasso')
result <- rbind(result,mean_c_index)

write.table(result,"result.txt",sep="\t")






#rs
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(survival)
library(survminer)
library(timeROC)
library(pROC)
library(export)
library(ggsci)

my<-cbind(rs[,1:3],my2[,3:length(my2)])

###桑基图展示la
# 假设la是你的10个基因列表组成的列表
library(ggplot2)
library(ggalluvial)
# 将列表转换为数据框，每个基因对应其列表编号
la_data <- do.call(rbind,
                   lapply(seq_along(la), function(i) {
                     data.frame(list_id = paste("List", i, sep=""), gene = la[[i]])
                   }))

# 假设la_data是你的数据框
library(dplyr)

# 找到所有列表共有的基因
common_genes <- la_data %>%
  group_by(gene) %>%
  summarise(count = n_distinct(list_id)) %>%
  filter(count == max(count)) %>%
  pull(gene)

# 过滤出只在所有列表中出现的基因
la_data_common <- la_data %>% 
  filter(gene %in% common_genes)

# 因为每个基因都在所有列表中出现，所以我们可以给它们相同的频率
la_data_common$Freq <- 1

# 绘制桑基图
library(ggplot2)
library(ggalluvial)

ggplot(la_data_common,
       aes(axis1 = list_id, axis2 = gene, y = Freq)) +
  geom_alluvium(aes(fill = list_id)) +   # 桑基流动部分
  geom_stratum() +                       # 分层显示每个列表
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +  # 添加标签
  theme_minimal() +
  ggtitle("Sankey Diagram of Common Genes Across Lists")




library(ggplot2)
library(dplyr)
data<-read.table("all-result.txt",header = T)
df<-data

# Calculate the mean C-index across the three studies
df$Mean_C_Index <- rowMeans(df[, c('GSE84437', 'STAD', 'GSE26257')])

# Reshape the data for plotting with ggplot2
df_long <- tidyr::pivot_longer(df, cols = c('GSE84437', 'STAD', 'GSE26257', 'Mean_C_Index'), 
                               names_to = 'Dataset', values_to = 'C_Index')

# Plot
ggplot(df_long, aes(x = Model, y = C_Index, group = Dataset, color = Dataset)) + 
  geom_line() +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "C-index Values for Each Model", x = "Model", y = "C-Index")







#RSF

library(ggplot2)

keygene_importance <- average_gene_importance[keygene]
keygene_data <- data.frame(Gene = keygene, Importance = keygene_importance)
ggplot(keygene_data, aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() +  # 翻转坐标轴，使基因名称在y轴上显示
  labs(title = "Average Importance of Key Genes for RSF Model", x = "Importance", y = "Gene")

#GBM
#在代码中，我使用了10折交叉。我想得到所有基因交叉后的平均权重，并以此作为GBM的模型，并绘图展示重要基因。
library(gbm)
library(dplyr)
library(ggplot2)

df<-my2
# 初始化一个向量来存储所有基因的累积重要性
gene_importance_cumulative <- rep(0, ncol(df) - 2)  # 减去2是因为前两列是生存时间和事件

# 进行交叉验证
for (fold in 1:nfold) {
  # 分割数据为训练集和测试集
  train_data <- df[fold_indices != fold, ]
  test_data <- df[fold_indices == fold, ]
  
  # 拟合GBM模型
  fit <- gbm(Surv(OS.time, OS) ~ ., data = train_data, distribution = 'coxph',
             n.trees = 10000, interaction.depth = 3, n.minobsinnode = 10,
             shrinkage = 0.001, cv.folds = 10, n.cores = NULL)

  gbm_imp <- summary(fit, n.trees = fit$n.trees, plot.it = FALSE)
  
gene_importance_cumulative <- gene_importance_cumulative + gbm_imp$rel.inf
}

average_gene_importance <- gene_importance_cumulative / nfold


gene_importance_df <- data.frame(
Gene = colnames(df)[-c(1:2)],
Importance = average_gene_importance
)

ordered_genes <- gene_importance_df %>%
arrange(desc(Importance))

top_genes <- head(ordered_genes, 20)

ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance)) +
geom_bar(stat = "identity", fill = "steelblue") +
coord_flip() + # 翻转坐标，使基因名在y轴
labs(title = "Average Importance of Genes across All Folds", x = "Importance", y = "Gene") +
theme_light()

top_features<-top_genes[1:7,]$Gene

train_data_reduced <- my2[, c("OS.time", "OS", top_features)]

fit_reduced <- gbm(Surv(OS.time, OS) ~ ., data = train_data_reduced, distribution = 'coxph',
                   n.trees = 10000, interaction.depth = 3, n.minobsinnode = 10,
                   shrinkage = 0.001, cv.folds = 10, n.cores = NULL)
risk_scores <- predict(fit_reduced, newdata = train_data_reduced , n.trees = 10000, type = "response")

train_data_reduced$RiskScore <- risk_scores




