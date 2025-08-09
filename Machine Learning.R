library(caret)#数据拆分及交叉验证
library(randomForest)#随机森林
library(glmnet)#lasso
library(xgboost)#xgboost
library(easyTCGA)
library(IOBR)
library(caret)
library(viridis)
gene<- read.table("/Volumes/WD/ccRCC_Epi/Bulk/TCGA_ICGC/Gene2.csv",header = T)
eset_com <- rmb_cmb
eset_com<-as.matrix(eset_com[gene$Gene,])
eset_com<-as.data.frame(t(eset_com))

pheno_new <- combined_datapd %>% filter(combined_datapd$Sample %in% rownames(eset_com))
#pheno_new<-pheno_new %>% dplyr::select(Sample,group)
pheno_new<-pheno_new %>% dplyr::rename(Event=group)
rownames(pheno_new) <-pheno_new$Sample
#pheno_new<-pheno_new %>% dplyr::select(Event)
encoded_df <- cbind(pheno_new$Event,eset_com)
encoded_df<-encoded_df %>% dplyr::rename(Event=`pheno_new$Event`)
encoded_df$Event <- as.factor(encoded_df$Event)

# 将数据分为训练集和测试集
set.seed(123)
#trainIndex <- createDataPartition(encoded_df$Event, p = .7,
#                                               list = FALSE,
#                                               times = 1)
#trainData<-encoded_df[ trainIndex,]
#testData<-encoded_df[-trainIndex,]

trainsample<-subset(pheno_new,pheno_new$dataset=="TCGA"|pheno_new$dataset=="GTEX")
testsample<-subset(pheno_new,pheno_new$dataset!="TCGA"|pheno_new$dataset!="GTEX")
trainData<-encoded_df %>% dplyr::filter(row.names(encoded_df) %in% trainsample$Sample)
testData<-encoded_df %>% dplyr::filter(row.names(encoded_df) %in% testsample$Sample)

fit.rf <-randomForest(Event ~ ., data =trainData, importance = TRUE, ntree = 500)#10折交叉验证
set.seed(123)
ctrl <- trainControl(method = "cv", number = 10)
rf_model_cv <- train(Event ~ ., data =trainData, method = "rf", trControl = ctrl)
#模型评价
rf_pred <- predict(rf_model_cv, newdata =testData,type="prob")[,2]
rf_roc<-roc(testData$Event,rf_pred)
auc_value <- auc(rf_roc)
p1<-ggroc(rf_roc,legacy.axes = TRUE,size=1,color="#69b3a2")+
  ggtitle("ROC curve of Randomforest")+
  geom_abline(intercept = 0,slope = 1,linetype="dashed",color="grey")+
  theme_bw()+
  annotate(geom = "text",label=paste0("AUC in test set:",round(auc_value,3)),
           x=0.7,y=0.05)
p1
ggsave("RFroc.pdf",plot =p1,width = 4,height = 4,units = 'in',dpi = 600)
varImpPlot(fit.rf)
pre_lasso <- predict(cv_lasso, s = cv_lasso$lambda.min, newx = x_test, type = "class")
pre_lasso = ifelse(pre_lasso=="0","normal","tumor")
tab = table(pre_lasso,testData$Event,dnn=c("true","pre"))
result_matrix<-caret::confusionMatrix(tab)
draw_confusion_matrix(result_matrix)

library(yardstick)
library(ggplot2)
data <- data.frame(Actual = testData$Event,
                   Prediction = pre_lasso
)
data <- data %>% dplyr::rename(Prediction=s1)
data$Prediction<-as.factor(data$Prediction)
table(data$Prediction, data$Actual)
cm <- conf_mat(data, Actual, Prediction)
autoplot(cm, type = "heatmap") +
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1")
ggsave(filename = "confusematrix_lasso.pdf",width = 3,height = 3)
