# Information --------------------------------------------------------------------------------------------------------
# Version : 0.1.0
# Date : 15-November-2020
# Discription : Treatment response to bumetanide is better characterized by an immuno-behavioural covariation in young children with autism spectrum disorder

# Working Directory --------------------------------------------------------------------------------------------------
path_wd <- "D:/lqy/cytokine/final_final"
setwd(path_wd)

# Load Packages
library(readxl)
library(sva)
library(psych)
library(reshape2)
library(ggplot2)
library(factoextra)
library(sRDA)
library(corrplot)
library(plyr)
library(ggm)
library(cluster)
library(NbClust)
library(effectsize)
library(fmsb) 
library(ggpubr)
library(gridExtra)
library(caret)
library(pROC)
source(file.path("model_function.R"))

# Data I/O -----------------------------------------------------------------------------------------------------------
# get data
mydata <- as.data.frame(read_xlsx("data.xlsx"))
bldata <- mydata[mydata$bt==1,]    # the baseline levels of data
deltadata <- mydata[mydata$bt==3,]   #the change levels of data

# difference analyses of change levels of cytokines
diff_2batch_ccyto <- list()
w_cyto <- data.frame()
diff_2batch_ccyto[[1]] <- aggregate(deltadata[,c(25:59)],by=list(deltadata$batch),FUN = mean)
diff_2batch_ccyto[[2]] <- aggregate(deltadata[,c(25:59)],by=list(deltadata$batch),FUN = sd)
for (j in c(25:59)) {
  w_cyto[1,j-24] <- wilcox.test(deltadata[deltadata$batch==1,j],deltadata[deltadata$batch==2,j])$statistic
  w_cyto[2,j-24] <- wilcox.test(deltadata[deltadata$batch==1,j],deltadata[deltadata$batch==2,j])$p.value
  w_cyto[3,j-24] <- t.test(deltadata[deltadata$batch==1,j])$statistic
  w_cyto[4,j-24] <- t.test(deltadata[deltadata$batch==1,j])$p.value
  w_cyto[5,j-24] <- t.test(deltadata[deltadata$batch==2,j])$statistic
  w_cyto[6,j-24] <- t.test(deltadata[deltadata$batch==2,j])$p.value
}
for(i in 1:3){
  w_cyto[2*i,] <- p.adjust(w_cyto[2*i,],method = "fdr",35)
}
names(w_cyto) <- colnames(deltadata)[25:59]
diff_2batch_ccyto[[3]] <- w_cyto

# Data Preprocessing ------------------------------------------------------------------------------------------------
# PCA plot before adjustment for batch effect
pca <- prcomp(bldata[,25:59],scale=T)
fviz_pca_ind(pca, col.ind=as.factor(bldata$batch), geom = c("point"),mean.point=F, addEllipses = T, legend.title="batch") 
pca <- prcomp(deltadata[,25:59],scale=T)
fviz_pca_ind(pca, col.ind=as.factor(deltadata$batch), mean.point=F, addEllipses = T, legend.title="batch")

# min-max normalization and log transformation for the baseline and changes levels of cytokines each batch
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
bldata[,25:59] <- rbind(as.data.frame(lapply(bldata[bldata$batch==1,25:59],normalize)),
                  as.data.frame(lapply(bldata[bldata$batch==2,25:59],normalize)))
deltadata[,25:59] <- rbind(as.data.frame(lapply(deltadata[deltadata$batch==1,25:59],normalize)),
                     as.data.frame(lapply(deltadata[deltadata$batch==2,25:59],normalize)))
bldata[,25:59] <- log2(bldata[,25:59]+1)
deltadata[,25:59] <- log2(deltadata[,25:59]+1)

# adjustment for batch effect and PCA 
mydata1 <- bldata[,25:59]    #get baseline levels of cytokine data
cdata <- as.matrix(mydata1)
cdata <- t(cdata)
colnames(cdata) <- 1:79
combat_mdata = ComBat(dat=cdata, batch=bldata$batch,  par.prior=T,prior.plots=F,ref.batch=2)
write.csv(combat_mdata,"combat_mdata.csv",quote = F)
bldata[,25:59] <- t(combat_mdata)
pca <- prcomp(bldata[,25:59],scale=T)
fviz_pca_ind(pca, col.ind=as.factor(bldata$batch), geom = c("point"),mean.point=F, addEllipses = T, legend.title="batch") 
#PCA plot after adjustment for batch effect

bldata1 <- bldata

mydata1 <- deltadata[,25:59]   #get changes levels of cytokines data
cdata <- as.matrix(mydata1)
cdata <- t(cdata)
colnames(cdata) <- 1:79
combat_mdata = ComBat(dat=cdata, batch=deltadata$batch,  par.prior=F,prior.plots=T)
deltadata[,25:59] <- t(combat_mdata)
pca <- prcomp(deltadata[,25:59],scale=T)
fviz_pca_ind(pca, col.ind=as.factor(deltadata$batch), mean.point=F, addEllipses = T, legend.title="batch")

deltadata1<-deltadata

# get different type of datasets
deltacy <- deltadata
deltacy[,6:24] <- bldata[,6:24]
deltaca <- deltadata
deltaca[,25:59] <- bldata[,25:59]
dataset <- list(bldata,deltadata,deltacy,deltaca)

# Difference analyses on two batches -------------------------------------------------------------------------------
# difference analyses for baseline levels of cytokines after adjustment for batch effect
diff_2batch_bcyto <- list()
w_cyto <- data.frame()
diff_2batch_bcyto[[1]] <- aggregate(bldata[,c(25:59)],by=list(bldata$batch),FUN = mean)
diff_2batch_bcyto[[2]] <- aggregate(bldata[,c(25:59)],by=list(bldata$batch),FUN = sd)
for (j in c(25:59)) {
  w_cyto[1,j-24] <- wilcox.test(bldata[bldata$batch==1,j],bldata[bldata$batch==2,j])$statistic
  w_cyto[2,j-24] <- wilcox.test(bldata[bldata$batch==1,j],bldata[bldata$batch==2,j])$p.value
}
w_cyto[2,] <- p.adjust(w_cyto[2,],method = "fdr",35)
names(w_cyto) <- colnames(deltadata)[25:59]
diff_2batch_bcyto[[3]] <- w_cyto


# difference analyses of baseline and changes levels of CARS
diff_2batch_cars <- list()
w_cyto<-data.frame()
diff_2batch_cars[[1]] <- aggregate(bldata[,c(21:24)],by=list(bldata$batch),FUN = mean)
diff_2batch_cars[[2]] <- aggregate(bldata[,c(21:24)],by=list(bldata$batch),FUN = sd)
diff_2batch_cars[[3]] <- aggregate(deltadata[,c(21:24)],by=list(deltadata$batch),FUN = mean)
diff_2batch_cars[[4]] <- aggregate(deltadata[,c(21:24)],by=list(deltadata$batch),FUN = sd)
for (j in c(21:24)) {
  w_cyto[1,j-20] <- t.test(bldata[bldata$batch==1,j],bldata[bldata$batch==2,j])$statistic
  w_cyto[2,j-20] <- t.test(bldata[bldata$batch==1,j],bldata[bldata$batch==2,j])$p.value
  w_cyto[3,j-20] <- wilcox.test(deltadata[deltadata$batch==1,j],deltadata[deltadata$batch==2,j])$statistic
  w_cyto[4,j-20] <- wilcox.test(deltadata[deltadata$batch==1,j],deltadata[deltadata$batch==2,j])$p.value
}
w_cyto[4,] <- p.adjust(w_cyto[4,],method = "fdr",4)
names(w_cyto) <- colnames(deltadata)[21:24]
diff_2batch_cars[[5]] <- w_cyto

# difference analyses for baseline levels of ADOS and SRS
diff_2batch_as <- list()
w_cyto <- data.frame()
diff_2batch_as[[1]] <- aggregate(bldata[,c(61:71)],by=list(bldata$batch),FUN = mean)
diff_2batch_as[[2]] <- aggregate(bldata[,c(61:71)],by=list(bldata$batch),FUN = sd)
for (j in c(61:71)) {
  w_cyto[1,j-60] <- wilcox.test(bldata[bldata$batch==1,j],bldata[bldata$batch==2,j])$statistic
  w_cyto[2,j-60] <- wilcox.test(bldata[bldata$batch==1,j],bldata[bldata$batch==2,j])$p.value
}
names(w_cyto) <- colnames(deltadata)[61:71]
diff_2batch_as[[3]] <- w_cyto

# difference analyses for baseline levels of sex,age and bmi
diff_2batch_other <- list()
w_cyto<-data.frame()
diff_2batch_other[[1]] <- aggregate(bldata[,c(3:5)],by=list(bldata$batch),FUN = mean)
diff_2batch_other[[2]] <- aggregate(bldata[,c(3:5)],by=list(bldata$batch),FUN = sd)
for (j in c(3:5)) {
  w_cyto[1,j-2] <- wilcox.test(bldata[bldata$batch==1,j],bldata[bldata$batch==2,j])$statistic
  w_cyto[2,j-2] <- wilcox.test(bldata[bldata$batch==1,j],bldata[bldata$batch==2,j])$p.value
}
gender <- matrix(table(bldata$gender,bldata$batch),nrow = 2)
w_cyto[1,1] <- chisq.test(gender)$statistic
w_cyto[2,1] <- chisq.test(gender)$p.value
names(w_cyto) <- colnames(deltadata)[3:5]
diff_2batch_other[[3]] <- w_cyto

# Nonparametric correlations with FDR between the CARS_total score and 35 cytokines --------------------------------
col1 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                           "#FFFFFF", "#FDDBC7", "#F4A582","#D6604D", "#B2182B","#67001F"))
for(i in 1:4){
  case <- dataset[[i]]
  case <- data.frame(case$CARS_total,case[,25:59])
  colnames(case)[1] <- "CASR_total"
  corcase <- corr.test(case,use = "complete")[["r"]]
  corrplot(corcase,tl.srt = 45,method = "color",tl.cex = 0.8,tl.col = "black",col = col1(20),type = "upper")
}


# sCCA -------------------------------------------------------------------------------------------------------------
# find significant sparse canonical component
ccol <- colnames(mydata)[22:24]
fcol <- colnames(mydata)[25:59]

for(i in 1:4){
  case_tr <- dataset[[i]][dataset[[i]]$batch==2,]
  case_te <- dataset[[i]][dataset[[i]]$batch==1,]
  y <- case_tr[,ccol]
  x <- case_tr[,fcol]
  cca.res <- sCCA(x,y,nonzero = 3,ridge_penalty = 1,penalization = "ust",multiple_LV  = TRUE,nr_LVs = 3)
  for(j in 1:3){
    x1_coef <- cca.res$ALPHA[[j]]
    y1_coef <- cca.res$BETA[[j]]
    x_s <- as.matrix(case_te[,fcol])%*%x1_coef
    y_s <- as.matrix(case_te[,ccol])%*%y1_coef
    x_n <- as.matrix(case_tr[,fcol])%*%x1_coef
    y_n <- as.matrix(case_tr[,ccol])%*%y1_coef
    x1 <- data.frame(x_n,y_n,rep(1,37))
    names(x1) = c("x","y","group")
    x2 <- data.frame(x_s,y_s,rep(2,42))
    names(x2) = c("x","y","group")
    x_y <- rbind(x1,x2)
    x_y$group <- as.factor(x_y$group)
    p1 <- cor.test(x_n,y_n,method = "spearman")$p.value
    p2 <- cor.test(x_s,y_s,method = "spearman")$p.value
    if(p1<0.05 & p2<0.05){
      print(i)
      print(j)
    }
  }
}   # i = 2, j = 1
    #only deltadata on mode 1 was significant

# get significant sparse canonical component
case_tr <- deltadata[deltadata$batch==2,]
case_te <- deltadata[deltadata$batch==1,]
y <- case_tr[,ccol]
x <- case_tr[,fcol]
cca.res <- sCCA(x,y,nonzero = 3,ridge_penalty = 1,penalization = "ust",multiple_LV  = TRUE,nr_LVs = 3)
rho <- rep(NA,6)
p_c <- rep(NA,6)
for(i in 1:3){
  x1_coef <- cca.res$ALPHA[[i]]
  y1_coef <- cca.res$BETA[[i]]
  x_s <- as.matrix(case_te[,fcol])%*%x1_coef
  y_s <- as.matrix(case_te[,ccol])%*%y1_coef
  x_n <- as.matrix(case_tr[,fcol])%*%x1_coef
  y_n <- as.matrix(case_tr[,ccol])%*%y1_coef
  x1 <- data.frame(x_n,y_n,rep(1,37))
  names(x1) = c("x","y","group")
  x2 <- data.frame(x_s,y_s,rep(2,42))
  names(x2) = c("x","y","group")
  x_y <- rbind(x1,x2)
  x_y$group <- as.factor(x_y$group)
  cor.test(x_n,y_n,method = "spearman")
  cor.test(x_s,y_s,method = "spearman")
  rho[2*i-1] <- cor.test(x_n,y_n,method = "spearman")$estimate
  rho[2*i] <- cor.test(x_s,y_s,method = "spearman")$estimate
  p_c[2*i-1] <- cor.test(x_n,y_n,method = "spearman")$p.value
  p_c[2*i] <- cor.test(x_s,y_s,method = "spearman")$p.value
}

# run sCCA with permuted CARS domains
nPerms <- 5000
set.seed(1111)
cars_perm_tr <- rlply(nPerms, y[sample(nrow(y)), ])
rho_perm <- data.frame()
for(i in 1:nPerms){
  y <- cars_perm_tr[[i]]
  cca.res <- sCCA(x,y,nonzero = 3,ridge_penalty = 1,penalization = "ust",multiple_LV  = TRUE,nr_LVs = 3)
  for(j in 1:3){
    x1_coef <- cca.res$ALPHA[[j]]
    y1_coef <- cca.res$BETA[[j]]
    x_s <- as.matrix(case_te[,fcol])%*%x1_coef
    y_s <- as.matrix(case_te[,ccol])%*%y1_coef
    x_n <- as.matrix(case_tr[,fcol])%*%x1_coef
    y_n <- as.matrix(case_tr[,ccol])%*%y1_coef
    x1 <- data.frame(x_n,y_n,rep(1,37))
    names(x1) = c("x","y","group")
    x2 <- data.frame(x_s,y_s,rep(2,42))
    names(x2) = c("x","y","group")
    x_y <- rbind(x1,x2)
    x_y$group <- as.factor(x_y$group)
    cor.test(x_n,y_n,method = "spearman")
    cor.test(x_s,y_s,method = "spearman")
    rho_perm[2*j-1,i] <- cor.test(x_n,y_n,method = "spearman")$estimate
    rho_perm[2*j,i] <- cor.test(x_s,y_s,method = "spearman")$estimate
  }
}

p_perm <- rep(NA,6)
for(i in 1:6){
  p_perm[i] <- sum(rho_perm[i,]>=rho[i])/nPerms
}
                 #only deltadata on mode 1 was significant through permutation

# Sensitivity analyses for canonical components
ccol <- colnames(mydata)[22:24]
fcol <- colnames(mydata)[25:59]
case_tr <- deltadata[deltadata$batch==2,]
case_te <- deltadata[deltadata$batch==1,]
y <- case_tr[,ccol]
x <- case_tr[,fcol]
cca.res <- sCCA(x,y,nonzero = 3,ridge_penalty = 1,penalization = "ust")
x1_coef <- cca.res$ALPHA
y1_coef <- cca.res$BETA
x_s <- as.matrix(case_te[,fcol])%*%x1_coef
y_s <- as.matrix(case_te[,ccol])%*%y1_coef
x_n <- as.matrix(case_tr[,fcol])%*%x1_coef
y_n <- as.matrix(case_tr[,ccol])%*%y1_coef
x1 <- data.frame(x_n,y_n,rep(1,37))
names(x1) = c("x","y","group")
x2 <- data.frame(x_s,y_s,rep(2,42))
names(x2) = c("x","y","group")
x_y <- rbind(x1,x2)
x_y$group <- as.factor(x_y$group)
# sCCA plot
p1 <- ggplot(x_y[x_y$group==1,],aes(x,y))+geom_point()+geom_smooth(method = "lm")+
      theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,face="bold",color="black"))+ 
      guides(color=guide_legend(title="data set"))+xlab("Cytokines Canonical Component Scores")+ylab("CARS Canonical Components Scores")
p2 <- ggplot(x_y[x_y$group==2,],aes(x,y))+geom_point()+geom_smooth(method = "lm")+
      theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,face="bold",color="black")) +
      guides(color=guide_legend(title="data set"))+xlab("Cytokines Canonical Component Scores")+ylab("CARS Canonical Components Scores")
grid.arrange(p1,p2,ncol = 2)

covariates <- c(3:5,22:24,28,35,37)
bl_cov <- bldata[bldata$batch==1,covariates]
p_cov <- data.frame()
for(i in 1:9){
  x <- data.frame(x_s,y_s,bl_cov[,i])
  names(x) <- c("x","y","covariates")
  p <- pcor(c("x","y","covariates"),var(x))
  p_cov[1,i] <- pcor.test(p,1,42)$tval
  p_cov[2,i] <- pcor.test(p,1,42)$pvalue
}
colnames(p_cov) <- colnames(bldata)[covariates]

# Cluster -------------------------------------------------------------------------------------------------------
# clustering based on the immuno-behavioural plane
case_x <- x_y[x_y$group==1,1:2]
set.seed(1234)
devAskNewPage(ask = T)
nc <- NbClust(scale(case_x),min.nc = 2,max.nc = 15,method = "kmeans")  #k=3
set.seed(1234)
fit.km <- kmeans(scale(case_x),3,nstart = 25)
fviz_cluster(fit.km, data = scale(case_x))
fit.km$size
km.mean <- data.frame()
km.mean <- aggregate(case_x,by=list(cluster=fit.km$cluster),mean)

x_y1 <- x_y[x_y$group==2,1:2]
for(i in 1:nrow(x_y1)){
  x_y1[i,3] <- sqrt((x_y1[i,1]-km.mean[1,2])**2+(x_y1[i,2]-km.mean[1,3])**2)
  x_y1[i,4] <- sqrt((x_y1[i,1]-km.mean[2,2])**2+(x_y1[i,2]-km.mean[2,3])**2)
  x_y1[i,5] <- sqrt((x_y1[i,1]-km.mean[3,2])**2+(x_y1[i,2]-km.mean[3,3])**2)
  
}
x_y1$clustering <- 0
for(i in 1:nrow(x_y1)){
  x_y1[i,]$clustering <- which.min(x_y1[i,3:5])
}
x_y$clustering <- c(fit.km$cluster,x_y1$clustering)
x_y$clustering <- as.factor(x_y$clustering)

# plot three groups with both Data Set 1 & 2
par(mfrow=c(1,2))
b <- data.frame()
b <- as.data.frame(fit.km$cluster)
plot(case_x,col = (fit.km$cluster + 1),pch = (fit.km$cluster + 15), cex = 2,xlim=c(5,35),ylim=c(-30,60),
                  main="Data Set 1",xlab="Cytokines Canonical Component Scores",ylab="CARS Canonical Component Scores")
points(km.mean[1, 2], km.mean[1, 3], pch = 10, col = "red", cex = 2)
points(km.mean[2, 2], km.mean[2, 3], pch = 10, col = "green", cex = 2)
points(km.mean[3, 2], km.mean[3, 3], pch = 10, col = "blue", cex = 2)
plot(case_x,col = "grey",pch = (fit.km$cluster + 15), cex = 2,xlim=c(5,35),ylim=c(-30,60),
                  main="Data Set 2",xlab="Cytokines Canonical Component Scores",ylab="CARS Canonical Component Scores")
points(km.mean[1, 2], km.mean[1, 3], pch = 10, col = "red", cex = 2)
points(km.mean[2, 2], km.mean[2, 3], pch = 10, col = "green", cex = 2)
points(km.mean[3, 2], km.mean[3, 3], pch = 10, col = "blue", cex = 2)
for(i in 1:nrow(x_y1)){
  points(x_y1[i,1],x_y1[i,2], col = (x_y1[i,]$clustering+1),
         pch = (x_y1[i,]$clustering + 15), cex = 2)
}
a = (km.mean[3,3]-km.mean[1,3])/(km.mean[3,2]-km.mean[1,2])
b = km.mean[3,3]-a*km.mean[3,2]
a1 = -1/a
b1 = (km.mean[3,3]+km.mean[1,3])/2-(km.mean[3,2]+km.mean[1,2])/2*a1
a = (km.mean[3,3]-km.mean[2,3])/(km.mean[3,2]-km.mean[2,2])
b = km.mean[3,3]-a*km.mean[3,2]
a2 = -1/a
b2 = (km.mean[3,3]+km.mean[2,3])/2-(km.mean[3,2]+km.mean[2,2])/2*a2
x1 <- c(0,-b1/a1)
y1 <- c(b1,0)
lines(x1,y1,col="black")
x2 <- c(0,40,-b2/a2)
y2 <- c(b2,40*a2+b2,0)
lines(x2,y2,col="black")

# Difference analyses on three immuno-behavioural groups ------------------------------------------------------
deltadata <- mydata[mydata$bt==3,]
deltadata$clustering <- c(x_y1$clustering,fit.km$cluster)
# difference analyses for change levels of 3 cytokines
diff_3group_ccyto <- list()
w_cyto <- data.frame()
p <- NA
case <- deltadata[,c(28,35,37,95)]
diff_3group_ccyto[[1]] <- aggregate(case[,-4],by=list(case$clustering),FUN = mean)
diff_3group_ccyto[[2]] <- aggregate(case[,-4],by=list(case$clustering),FUN = sd)
for (i in 1:3) {
  for(j in 1:3){
    w_cyto[i,3*j-2] <- t.test(case[case$clustering==j,i])$statistic
    w_cyto[i,3*j-1] <- t.test(case[case$clustering==j,i])$p.value
    w_cyto[i,3*j] <- hedges_g(case[case$clustering==j,i])$Hedges_g
    p <- c(p,w_cyto[i,3*j-1])
  }
  w_cyto[i,10] <- kruskal.test(case[,i]~case$clustering)$statistic
  w_cyto[i,11] <- kruskal.test(case[,i]~case$clustering)$p.value
}
p <- p[-1]
p <- p.adjust(p,method = "fdr",9)
p <- matrix(p,nrow = 3,byrow = T)
for(j in 1:3){
  w_cyto[,3*j-1] <- p[,j]
}
rownames(w_cyto) <- c("IFNgamma","MIG","IFNalpha2")
colnames(w_cyto) <- c(paste0(c("t","p","g"),c(rep(1,3),rep(2,3),rep(3,3))),"chi-squared","pchi")
w_cyto[,11] <- p.adjust(w_cyto[,11],method = "fdr")
diff_3group_ccyto[[3]] <- w_cyto


# difference analyses for change levels of 3 CARS subscales 
diff_3group_ccars <- list()
w_cyto <- data.frame()
p <- NA
case <- deltadata[,c(22:24,95)]
diff_3group_ccars[[1]] <- aggregate(case[,-4],by=list(case$clustering),FUN = mean)
diff_3group_ccars[[2]] <- aggregate(case[,-4],by=list(case$clustering),FUN = sd)
for (i in 1:3) {
  for(j in 1:3){
    w_cyto[i,3*j-2] <- t.test(case[case$clustering==j,i])$statistic
    w_cyto[i,3*j-1] <- t.test(case[case$clustering==j,i])$p.value
    w_cyto[i,3*j] <- hedges_g(case[case$clustering==j,i])$Hedges_g
    p <- c(p,w_cyto[i,3*j-1])
  }
  w_cyto[i,10] <- kruskal.test(case[,i]~case$clustering)$statistic
  w_cyto[i,11] <- kruskal.test(case[,i]~case$clustering)$p.value
}
p <- p[-1]
p <- p.adjust(p,method = "fdr",9)
p <- matrix(p,nrow = 3,byrow = T)
for(j in 1:3){
  w_cyto[,3*j-1] <- p[,j]
}
rownames(w_cyto) <- c("CARS_S","CARS_N","CARS_D")
colnames(w_cyto) <- c(paste0(c("t","p","g"),c(rep(1,3),rep(2,3),rep(3,3))),"chi-squared","pchi")
w_cyto[,11] <- p.adjust(w_cyto[,11],method = "fdr")
diff_3group_ccars[[3]] <- w_cyto

# difference analyses for CARS_component_scores and Cytokine_component_scores
a <- x1_coef[4]+x1_coef[11]+x1_coef[13]
deltadata$Cytokines_component_score <- x1_coef[4]/a*deltadata$IFNgamma+
                                       x1_coef[11]/a*deltadata$MIG+x1_coef[13]/a*deltadata$IFNalpha2
b <- y1_coef[1]+y1_coef[2]+y1_coef[3]
deltadata$CARS_component_score <- y1_coef[1]/b*deltadata$CARS_socialimpairment+
                                y1_coef[2]/b*deltadata$CARS_negativeemotionality+y1_coef[3]/b*deltadata$CARS_distortedsensoryresponse

# Cytokine_component_scores
p <- rep(NA,3)
case <- deltadata[,c(95,96)]
for (i in c(2)) {
  p[i-1] <- t.test(case[case$clustering==1,i])$p.value
  p[i] <- t.test(case[case$clustering==2,i])$p.value
  p[i+1] <- t.test(case[case$clustering==3,i])$p.value
}
p <- p.adjust(p,method = "fdr",3)

# CARS_component_scores
p <- rep(NA,3)
case <- deltadata[,c(95,97)]
for (i in c(2)) {
  p[i-1] <- t.test(case[case$clustering==1,i])$p.value
  p[i] <- t.test(case[case$clustering==2,i])$p.value
  p[i+1] <- t.test(case[case$clustering==3,i])$p.value
}
p <- p.adjust(p,method = "fdr",3)

# Radar chart
bldata <- mydata[mydata$bt==1,]
mycaco_c <- deltadata[,6:59]/bldata[,6:59]
mycaco_c <- data.frame(mycaco_c,bldata$clustering)
mean1 <- data.frame(median(mycaco_c[mycaco_c$bldata.clustering=="1",]$CARS1))
mean2 <- data.frame(median(mycaco_c[mycaco_c$bldata.clustering=="2",]$CARS1))
mean3 <- data.frame(median(mycaco_c[mycaco_c$bldata.clustering=="3",]$CARS1))
for(i in 2:54){
  mean1 <- cbind(mean1,median(mycaco_c[mycaco_c$bldata.clustering=="1",i]))
  mean2 <- cbind(mean2,median(mycaco_c[mycaco_c$bldata.clustering=="2",i]))
  mean3 <- cbind(mean3,median(mycaco_c[mycaco_c$bldata.clustering=="3",i]))
}
mean1 <- as.numeric(as.character(as.matrix(mean1[1,])))
mean2 <- as.numeric(as.character(as.matrix(mean2[1,])))
mean3 <- as.numeric(as.character(as.matrix(mean3[1,])))
mean_t <- as.data.frame(rbind(mean1,mean2,mean3))
ccol <- colnames(bldata)[6:59]
colnames(mean_t) <- ccol
rownames(mean_t) <- c("SS-responding group","E-responding group","A-responding group")
max(mean_t)
min(mean_t)
mean_ma <- rep(0.5,54)
mean_mi <- rep(-0.5,54)
mean_w <- rbind(mean_ma,mean_mi,mean_t)
mean_w <- mean_w[,16:54]
cccol <- c("IL6","IL1alpha","IL7","IL17","IL9","IFNalpha2","IL1beta","IFNgamma","MIG","IP10","Eotaxin","IL16","IL18",
           "TNFalpha","TNFbeta","TRAIL","LIF","MIF","PDGFBB","MCP1MCAF","GROalpha","IL8","SDF1alpha","MIP1alpha","MIP1beta",
          "RANTES","CTACK","CARS_total","CARS_socialimpairment","CARS_negativeemotionality","CARS_distortedsensoryresponse",
          "GCSF","SCF","HGF","MCSF","SCGFbeta","IL4","IL13","IL2Ralpha")

par(mfrow=c(1,1))
set.seed(1234)
colors_border = c( "#CD5C5CE6" ,  "#2E8B57E6", "#4169E1E6")
colors_in = c( "#CD5C5C26" , "#2E8B5726" ,"#4169E126")
radarchart(mean_w[,cccol],axistype=1 , pty=32, centerzero = T,
           pcol=colors_border , pfcol=colors_in , plwd=3 , plty=1,
           cglcol="grey" ,cglty=1, axislabcol="black", cglwd=1.0,caxislabels=c('-0.5','-0.25','0','0.25','0.5'),
           vlcex=1.0)
legend(x=1.0, y=1.2, legend = rownames(mean_w[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.0, pt.cex=3)

# boxplot
names(deltadata)[95] <- "cluster"
deltadata$cluster <- as.factor(deltadata$cluster)
my_comparisons <- list(c(1,2),c(1,3),c(2,3))
p1 <- ggplot(deltadata,aes(cluster,CARS_socialimpairment,col=cluster))+geom_boxplot()+geom_jitter(width = 0.2,alpha=0.5) +
       theme(legend.position="none",panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA),axis.text=element_text(size=10,color="black")，
      axis.title=element_text(size=12,face="bold",color="black"),axis.title.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 8)
p2 <- ggplot(deltadata,aes(cluster,CARS_negativeemotionality,col=cluster))+geom_boxplot()+geom_jitter(width = 0.2,alpha=0.5) +
       theme(legend.position="none",panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA),axis.text=element_text(size=10,color="black") , 
      axis.title=element_text(size=12,face="bold",color="black"),axis.title.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 3.5)
p3 <- ggplot(deltadata,aes(cluster,CARS_distortedsensoryresponse,col=cluster))+geom_boxplot()+geom_jitter(width = 0.2,alpha=0.5)+
       theme(legend.position="none",panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA),axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,face="bold",color="black"),axis.title.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 3)
p4 <- ggplot(deltadata,aes(cluster,CARS_component_score,col=cluster))+geom_boxplot()+geom_jitter(width = 0.2,alpha=0.5)+
       theme(legend.position="none",panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA),axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,face="bold",color="black"),axis.title.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 7.2)
p5 <- ggplot(deltadata,aes(cluster,IFNgamma,col=cluster))+geom_boxplot()+geom_jitter(width = 0.2,alpha=0.5)+
       theme(legend.position="none",panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA),axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,face="bold",color="black"),axis.title.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 60)
p6 <- ggplot(deltadata,aes(cluster,MIG,col=cluster))+geom_boxplot()+geom_jitter(width = 0.2,alpha=0.5)+
       theme(legend.position="none",panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA),axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,face="bold",color="black"),axis.title.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 4500)
p7 <- ggplot(deltadata,aes(cluster,IFNalpha2,col=cluster))+geom_boxplot()+geom_jitter(width = 0.2,alpha=0.5)+
       theme(legend.position="none",panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA),axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,face="bold",color="black"),axis.title.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 75)
p8 <- ggplot(deltadata,aes(cluster,Cytokines_component_score,col=cluster))+geom_boxplot()+geom_jitter(width = 0.2,alpha=0.5)+
       theme(legend.position="none",panel.background=element_blank(),panel.border=element_rect(linetype="solid",fill=NA),axis.text=element_text(size=10,color="black"),
      axis.title=element_text(size=12,face="bold",color="black"),axis.title.x = element_blank())+
      stat_compare_means(comparisons = my_comparisons)+stat_compare_means(label.y = 2000)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol = 4)

# Prediction of the treatment response to bumetanide using the baseline information -----------------------------------------
# Prediction for SS-responding group
bldata <- bldata1
bldata$clustering <- c(x_y1$clustering,fit.km$cluster)
bldata[bldata$clustering==2,]$clustering <- 3
bldata$clustering <- make.names(bldata$clustering,unique = FALSE, allow_ = TRUE)
bldata$clustering <- as.factor(bldata$clustering)
newdata1 <- bldata[,c(1,3:4,21:71,95)]
Process = preProcess(newdata1,method = 'knnImpute')
newdata = predict(Process, newdata1)
newdata$batch <- bldata$batch
case_tr <- newdata[newdata$batch==2,]
case_te <- newdata[newdata$batch==1,]
# including baseline levels of cytokines 
# feature selection
subsets = c(3,4,5,6,7,8)
trctrl <- trainControl(method = "repeatedcv",number = 3,repeats = 5)
ctrl = rfeControl(functions = rfFuncs, method = "cv",verbose = FALSE, returnResamp = "final")
Profile = rfe(case_tr[,c(4:42,44:54)], case_tr[,55], sizes = subsets, rfeControl = ctrl)
print(Profile)   # variables = 3 
fcol <- Profile$optVariables
i <- length(fcol)
fcol[i+1] <- "clustering"
fcol[i+2] <- "gender"
fcol[i+3] <- "age_month"
fcol
datTrain <- case_tr[,fcol]
datTest <- case_te[,]
datTest1 <- datTest
datTrain <- upSample(datTrain[,-(i+1)],datTrain[,i+1],yname = "clustering")
table(datTrain$clustering)
# modeling and plot ROC curve
model_ss_c <- models(datTrain)
ROC_X1(model_ss_c,datTrain,datTest)

# without baseline levels of cytokines 
# feature selection
subsets = 3
Profile = rfe(case_tr[,c(4:7,44:54)], case_tr[,55], sizes = subsets, rfeControl = ctrl)
print(Profile)   # variables = 3 
fcol <- Profile$optVariables
i <- length(fcol)
fcol[i+1] <- "clustering"
fcol[i+2] <- "gender"
fcol[i+3] <- "age_month"
fcol
datTrain <- case_tr[,fcol]
datTest <- case_te[,]
datTrain <- upSample(datTrain[,-(i+1)],datTrain[,i+1],yname = "clustering")
table(datTrain$clustering)
# modeling and plot ROC curve
model_ss_nc <- models(datTrain)
ROC_X1(model_ss_nc,datTrain,datTest)

# delta AUC of models with and without using the baseline cytokines
deltaauc_ss <- auc(model_ss_c,model_ss_nc,datTest)

# Prediction for E-responding group
bldata <- bldata1
bldata$clustering <- c(x_y1$clustering,fit.km$cluster)
bldata[bldata$clustering==1,]$clustering <- 3
bldata$clustering <- make.names(bldata$clustering,unique = FALSE, allow_ = TRUE)
bldata$clustering <- as.factor(bldata$clustering)
newdata1 <- bldata[,c(1,3:4,21:71,95)]
Process = preProcess(newdata1,method = 'knnImpute')
newdata = predict(Process, newdata1)
newdata$batch <- bldata$batch
case_tr <- newdata[newdata$batch==2,]
case_te <- newdata[newdata$batch==1,]
# including baseline levels of cytokines 
# feature selection
subsets = c(3,4,5,6,7,8)
Profile = rfe(case_tr[,c(4:42,44:54)], case_tr[,55], sizes = subsets, rfeControl = ctrl)
print(Profile)   # variables = 3 
fcol <- Profile$optVariables
i <- length(fcol)
fcol[i+1] <- "clustering"
fcol[i+2] <- "gender"
fcol[i+3] <- "age_month"
fcol
datTrain <- case_tr[,fcol]
datTest <- case_te[,]
datTrain <- upSample(datTrain[,-(i+1)],datTrain[,i+1],yname = "clustering")
table(datTrain$clustering)
# modeling and plot ROC curve
model_e_c <- models(datTrain)
ROC_X2(model_e_c,datTrain,datTest)

# without baseline levels of cytokines 
# feature selection
subsets = 3
Profile = rfe(case_tr[,c(4:7,44:54)], case_tr[,55], sizes = subsets, rfeControl = ctrl)
print(Profile)   # variables = 3 
fcol <- Profile$optVariables
i <- length(fcol)
fcol[i+1] <- "clustering"
fcol[i+2] <- "gender"
fcol[i+3] <- "age_month"
fcol
datTrain <- case_tr[,fcol]
datTest <- case_te[,]
datTrain <- upSample(datTrain[,-(i+1)],datTrain[,i+1],yname = "clustering")
table(datTrain$clustering)
# modeling and plot ROC curve
model_e_nc <- models(datTrain)
ROC_X2(model_e_nc,datTrain,datTest)

# delta AUC of models with and without using the baseline cytokines
deltaauc_e <- auc(model_e_c,model_e_nc,datTest)

# Performance comparison between SS-responding group and behaviourally defined responder group
# Prediction for behaviourally defined responder group(bumetanide response rate is 30%)
bldata <- bldata1
newdata <- bldata[,c(1,3:4,21:65)]
newdata$clustering <- "X1"
newdata[deltadata$CARS_total<2.5,]$clustering <- "X3"
newdata$clustering <- as.factor(newdata$clustering)
table(newdata$clustering)

subsets = c(3,4,5,6,7,8)
Profile = rfe(newdata[,c(4:7,44:48)], newdata[,49], sizes = subsets, rfeControl = ctrl)
print(Profile)
fcol <- Profile$optVariables
i <- length(fcol)
fcol[i+1] <- "clustering"
fcol[i+2] <- "gender"
fcol[i+3] <- "age_month"
fcol
case_tr <- newdata[newdata$batch==2,]
case_te <- newdata[newdata$batch==1,]
datTrain <- case_tr[,fcol]
datTest <- case_te[,]
datTest2 <- datTest
datTrain <- upSample(datTrain[,-(i+1)],datTrain[,i+1],yname = "clustering")
table(datTrain$clustering)

# modeling and plot ROC curve
model_r30_nc <- models(datTrain)
ROC_X1(model_r30_nc,datTrain,datTest)

# comparison between SS-responding group and behaviourally defined responder group
deltaauc_r30 <- auc_r(model_ss_c,model_r30_nc,datTest1,datTest2)


# Prediction for behaviourally defined responder group(bumetanide response rate is 40%)
bldata <- bldata1
newdata <- bldata[,c(1,3:4,21:65)]
newdata$clustering <- "X1"
newdata[deltadata$CARS_total<2,]$clustering <- "X3"
newdata$clustering <- as.factor(newdata$clustering)
table(newdata$clustering)

subsets = c(3,4,5,6,7,8)
Profile = rfe(newdata[,c(4:7,44:48)], newdata[,49], sizes = subsets, rfeControl = ctrl)
print(Profile)
fcol <- Profile$optVariables
i <- length(fcol)
fcol[i+1] <- "clustering"
fcol[i+2] <- "gender"
fcol[i+3] <- "age_month"
fcol
case_tr <- newdata[newdata$batch==2,]
case_te <- newdata[newdata$batch==1,]
datTrain <- case_tr[,fcol]
datTest <- case_te[,]
datTest2 <- datTest
datTrain <- upSample(datTrain[,-(i+1)],datTrain[,i+1],yname = "clustering")
table(datTrain$clustering)

# modeling and plot ROC curve
model_r40_nc <- models(datTrain)
ROC_X1(model_r40_nc,datTrain,datTest)

# comparison between SS-responding group and behaviourally defined responder group
deltaauc_r40 <- auc_r(model_ss_c,model_r40_nc,datTest1,datTest2)

