
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
#names(getModelInfo())

set.seed(123)      #????????
inputFile="merge.normalize.txt"      #?????????ļ?
geneFile="interGenes.txt"             #?????б??ļ?
setwd("")      #???ù???Ŀ¼

#??ȡ?????????ļ?
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#??ȡ?????б??ļ?,??ȡ?????????????ı???��
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
row.names(data)=gsub("-", "_", row.names(data))

#??ȡ??Ʒ??????Ϣ(????????ʵ????)
data=t(data)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
data=as.data.frame(data)
data$Type=group

#?????ݽ??з???(ѵ��???Ͳ?????)
inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

#??ȡ????Ԥ????????
yTestClass=test$Type
yTest=ifelse(yTestClass=="Control", 0, 1)
test<-test[,-ncol(test)]

#??ȡ????ѧϰ???????ļ?
control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)

#?Ի???ѧϰ????????ѭ??,????????ѧϰ??ģ??
modelList=list()
AUC=c()
ROCcolor=rainbow(nrow(methodRT))
for(i in 1:nrow(methodRT)){
	name=methodRT[i,"Name"]
	method=methodRT[i,"Method"]
	#????????ѧϰ??ģ??
	model=train(Type ~ ., data = train, method=method, trControl = control)
	if(name=="SVM"){
		model=train(Type ~ ., data = train, method=method, prob.model=TRUE, trControl = control)
	}
	#?õ?ÿ????ƷԤ???Ľ???
	pred=predict(model, newdata=test, type="prob")
	#????ROC????
	roc=roc(yTest, as.numeric(pred[,2]))
	AUC=c(AUC, paste0(name, ': ', sprintf("%.03f",roc$auc)))
	modelList[[method]]=as.numeric(roc$auc)
	if(i==1){
		pdf(file="ROC.pdf", width=5.5, height=5)
		plot(roc, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3)
	}else{
		plot(roc, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3, add=T)
	}
}
legend('bottomright', AUC, col=ROCcolor, lwd=3, bty = 'n', cex=0.9)
dev.off()

#????ROC?????µ???????ȡ???ŵ?ģ??
aucValue=unlist(modelList)
bestMethod=names(which(aucValue==max(aucValue)))
bestMethod
train$Type=ifelse(train$Type=="Control", 0, 1)
model=train(Type ~., data = train, method = bestMethod, trControl=control)

#????SHAPֵ
fit=permshap(model, train[,-ncol(train)])
#fit=kernelshap(model, train[,-ncol(train)])      #?????Ƚ϶???ʱ????????????
shp <- shapviz(fit, X_pred = train[,-ncol(train)], X=train[,-ncol(train)], interactions=T)
#???ݹ??׳̶ȶԻ???????????
important=sort(colMeans(abs(shp$S)), decreasing=T)
showVars=names(important)

#??????״ͼ
pdf(file="barplot.pdf", width=6, height=6)
sv_importance(shp, kind="bar", show_numbers=TRUE)+theme_bw()
dev.off()

#???Ʒ?Ⱥͼ
pdf(file="bee.pdf", width=7, height=6)
sv_importance(shp, kind = "bee", show_numbers=TRUE)+theme_bw()
dev.off()

#ɢ??ͼ(??��ͼ)
pdf(file="dependence.pdf", width=9, height=6)
sv_dependence(shp, v = showVars)+theme_bw()
dev.off()

#?ٲ?ͼ
pdf(file="waterfall.pdf", width=7, height=5)
sv_waterfall(shp, row_id = 1)
dev.off()

#????Ʒ��ͼ
pdf(file="force.pdf", width=9, height=5)
sv_force(shp, row_id = 1)
dev.off()


######??????ѧ??