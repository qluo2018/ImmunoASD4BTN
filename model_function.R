models <- function(Train_data){
  models <- list()
  #svm
  set.seed(123)
  svm_fit <- train(clustering ~ ., data = Train_data,
                   method = "svmLinear2",
                   preProcess = c("center", "scale"),
                   trControl = trctrl,         
                   tuneLength = 3
  )
  
  #pls
  set.seed(123)
  pls_fit <- train(clustering ~ ., data = Train_data,
                   method = "pls",
                   preProcess = c("center", "scale"),
                   trControl = trctrl,         
                   tuneLength = 3
  )
  
  #nn
  set.seed(123)
  nn_fit <- train(clustering ~ ., data = Train_data,
                  method = "nnet",
                  preProcess = c("center", "scale"),
                  trControl = trctrl,         
                  tuneLength = 3
  )
  
  #slda
  set.seed(123)
  slda_fit <- train(clustering ~ ., data = Train_data,
                    method = "sparseLDA",
                    preProcess = c("center", "scale"),
                    trControl = trctrl,         
                    tuneLength = 3
  )
  
  #orf
  set.seed(123)
  orf_fit <- train(clustering ~ ., data = Train_data,
                   method = "ORFpls",
                   preProcess = c("center", "scale"),
                   trControl = trctrl,  
                   tuneLength = 3
  )
  
  models[[1]] <- svm_fit
  models[[2]] <- pls_fit
  models[[3]] <- nn_fit
  models[[4]] <- slda_fit
  models[[5]] <- orf_fit
  models
}

ROC_X1 <- function(models,Train_data,Test_data){
  #ROC curve plot
  prob <- predict(models[[5]], Test_data,type = "prob")
  nn.ROC = roc(Test_data$clustering,predictor = prob$X1,levels = c("X1","X3"),smooth="T",plot = T,col= "#EE3B3B")
  prob <- predict(models[[4]], Test_data,type = "prob")
  nn.ROC = roc(Test_data$clustering,predictor = prob$X1,levels = c("X1","X3"),smooth="T")
  plot(nn.ROC,col= "#458B74",add=T)
  prob <- predict(models[[3]], Test_data,type = "prob")
  nn.ROC = roc(Test_data$clustering,predictor = prob$X1,levels = c("X1","X3"),smooth="T")
  plot(nn.ROC,col= "#FFB90F",add=T)
  prob <- predict(models[[2]], Test_data,type = "prob")
  nn.ROC = roc(Test_data$clustering,predictor = prob$X1,levels = c("X1","X3"),smooth="T")
  plot(nn.ROC,col= "#B23AEE",add=T)
  svm_fit <- svm(clustering ~ ., data = Train_data,probability=T)
  prob <- predict(svm_fit, Test_data,probability=T)
  prob <- as.data.frame(attr(prob, "probabilities"))
  nn.ROC = roc(Test_data$clustering,predictor = prob$X1,levels = c("X1","X3"),smooth="T")
  plot(nn.ROC,col= "#1E90FF",add=T)
  legend("bottomright", legend=c("ORF", "sLDA","NN","PLS","SVM"),
         col=c("#EE3B3B","#458B74","#FFB90F","#B23AEE","#1E90FF"), lwd=2,cex=0.7)
}

ROC_X2 <- function(models,Train_data,Test_data){
  #ROC curve plot
  prob <- predict(models[[5]], Test_data,type = "prob")
  nn.ROC = roc(Test_data$clustering,predictor = prob$X2,levels = c("X2","X3"),smooth="T",plot = T,col= "#EE3B3B")
  prob <- predict(models[[4]], Test_data,type = "prob")
  nn.ROC = roc(Test_data$clustering,predictor = prob$X2,levels = c("X2","X3"),smooth="T")
  plot(nn.ROC,col= "#458B74",add=T)
  prob <- predict(models[[3]], Test_data,type = "prob")
  nn.ROC = roc(Test_data$clustering,predictor = prob$X2,levels = c("X2","X3"),smooth="T")
  plot(nn.ROC,col= "#FFB90F",add=T)
  prob <- predict(models[[2]], Test_data,type = "prob")
  nn.ROC = roc(Test_data$clustering,predictor = prob$X2,levels = c("X2","X3"),smooth="T")
  plot(nn.ROC,col= "#B23AEE",add=T)
  svm_fit <- svm(clustering ~ ., data = Train_data,probability=T)
  prob <- predict(svm_fit, Test_data,probability=T)
  prob <- as.data.frame(attr(prob, "probabilities"))
  nn.ROC = roc(Test_data$clustering,predictor = prob$X2,levels = c("X2","X3"),smooth="T")
  plot(nn.ROC,col= "#1E90FF",add=T)
  legend("bottomright", legend=c("ORF", "sLDA","NN","PLS","SVM"),
         col=c("#EE3B3B","#458B74","#FFB90F","#B23AEE","#1E90FF"), lwd=2,cex=0.7)
}

auc <- function(model1,model2,Test_data){
  index <- createResample(Test_data$clustering,100,list = F)
  auc <- data.frame()
  for(i in 1:100){
    datTest1 <- Test_data[index[,i],]
    auc[1,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[1]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[2,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[1]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[3,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[2]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[4,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[2]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[5,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[3]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[6,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[3]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[7,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[4]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[8,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[4]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[9,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[5]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[10,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[5]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
  }
  for(i in 1:5){
    auc[10+i,] <- auc[2*i-1,]-auc[2*i,]
  }
  t_test <- list()
  for(j in 1:5){
    t_test[[j]] <- t.test(auc[j+10,])
  }
  t_test[[6]] <- auc
  t_test
}

auc_r <- function(model1,model2,Test_data1,Test_data2){
  index <- createResample(Test_data1$clustering,100,list = F)
  auc <- data.frame()
  for(i in 1:100){
    datTest1 <- Test_data1[index[,i],]
    datTest2 <- Test_data2[index[,i],]
    auc[1,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[1]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[2,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[1]], newdata = datTest2)),datTest2$clustering),measure = "auc")@y.values[[1]]
    auc[3,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[2]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[4,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[2]], newdata = datTest2)),datTest2$clustering),measure = "auc")@y.values[[1]]
    auc[5,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[3]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[6,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[3]], newdata = datTest2)),datTest2$clustering),measure = "auc")@y.values[[1]]
    auc[7,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[4]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[8,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[4]], newdata = datTest2)),datTest2$clustering),measure = "auc")@y.values[[1]]
    auc[9,i] <- ROCR::performance(prediction(as.numeric(predict(model1[[5]], newdata = datTest1)),datTest1$clustering),measure = "auc")@y.values[[1]]
    auc[10,i] <- ROCR::performance(prediction(as.numeric(predict(model2[[5]], newdata = datTest2)),datTest2$clustering),measure = "auc")@y.values[[1]]
  }
  for(i in 1:5){
    auc[10+i,] <- auc[2*i-1,]-auc[2*i,]
  }
  t_test <- list()
  for(j in 1:5){
    t_test[[j]] <- t.test(auc[j+10,])
  }
  t_test[[6]] <- auc
  t_test
}






