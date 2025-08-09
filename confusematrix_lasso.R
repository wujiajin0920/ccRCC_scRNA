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
