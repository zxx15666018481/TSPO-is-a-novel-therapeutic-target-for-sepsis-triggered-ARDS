
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)


#???ù???Ŀ¼
setwd("")
source("refer.ML.R")

#??ȡѵ��?????????ļ?
Train_data <- read.table("data.train.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr=Train_data[,1:(ncol(Train_data)-1),drop=F]
Train_class=Train_data[,ncol(Train_data),drop=F]

#??ȡ??֤?????????ļ?
Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors = F)
Test_expr=Test_data[,1:(ncol(Test_data)-1),drop=F]
Test_class=Test_data[,ncol(Test_data),drop=F]
Test_class$Cohort=gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
Test_class=Test_class[,c("Cohort", "Type")]

#??ȡѵ��??????֤???Ľ???????
comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[,comgene])
Test_expr <- as.matrix(Test_expr[,comgene])
Train_set = scaleData(data=Train_expr, centerFlags=T, scaleFlags=T) 
names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort))
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)

#??ȡ????ѧϰ???????ļ?
methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)
methods=methodRT$Model
methods <- gsub("-| ", "", methods)


#׼??????ѧϰģ?͵Ĳ???
classVar = "Type"         #???÷????ı?��??
min.selected.var = 5      #??????Ŀ????ֵ
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))


######################????ѵ��?????ݹ???????ѧϰģ??######################
#????ģ?????ϵ?һ?ֻ???ѧϰ????ɸѡ??��
preTrain.var <- list()       #???ڱ??????㷨ɸѡ?ı?��
set.seed(seed = 123)         #????????
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method,              #????ѧϰ????
                                 Train_set = Train_set,        #ѵ��???Ļ???????????
                                 Train_label = Train_class,    #ѵ��???ķ???????
                                 mode = "Variable",            #ѡ??????ģʽ(ɸѡ??��)
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)

#????ģ?????ϵڶ??ֻ???ѧϰ????????ģ??
model <- list()            #??ʼ??ģ?ͽ????б?
set.seed(seed = 123)       #????????
Train_set_bk = Train_set
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)
  Variable = preTrain.var[[method[1]]]
  Train_set = Train_set_bk[, Variable]
  Train_label = Train_class
  model[[method_name]] <- RunML(method = method[2],           #????ѧϰ????
                                Train_set = Train_set,        #ѵ��???Ļ???????????
                                Train_label = Train_label,    #ѵ��???ķ???????
                                mode = "Model",               #ѡ??????ģʽ(????ģ??)
                                classVar = classVar)
  
  #????ĳ?ֻ???ѧϰ????ɸѡ???ı?��С????ֵ?????÷???????Ϊ??
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk; rm(Train_set_bk)
#???????л???ѧϰģ?͵Ľ???
saveRDS(model, "model.MLmodel.rds")

#????????��?߼??ع?ģ??
FinalModel <- c("panML", "multiLogistic")[2]
if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){    #?????߼??ع?ģ?ͼ???ÿ???????ķ???????
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial", 
               data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit)
    return(tmp)
  })
}
#?????????Զ???��?߼??ع?ģ??
saveRDS(logisticmodel, "model.logisticmodel.rds")


#???ݻ???????��????ÿ???????ķ????÷?
model <- readRDS("model.MLmodel.rds")            #ʹ?ø???????ѧϰģ?͵????????Ϻ????????÷?
#model <- readRDS("model.logisticmodel.rds")     #ʹ?ö???��?߼??ع?ģ?ͼ????÷?
methodsValid <- names(model)                     #??????????????Ŀ??ȡ??Ч??ģ??
#???ݻ???????��Ԥ???????ķ??յ÷?
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], new_data = rbind.data.frame(Train_set,Test_set))
}
riskTab=as.data.frame(t(do.call(rbind, RS_list)))
riskTab=cbind(id=row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)

#???ݻ???????��Ԥ????Ʒ?ķ???
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = rbind.data.frame(Train_set,Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) # ??Ҫ?ϲ????Լ?????????????Ϣ?ļ??????д???
classTab=cbind(id=row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)

#??ȡÿ?ֻ???ѧϰ????ɸѡ???ı?��
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file="model.genes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#????ÿ??ģ?͵?AUCֵ
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],      #????ѧϰģ??
                                Test_set = Test_set,        #??֤???ı???????
                                Test_label = Test_class,    #??֤???ķ???????
                                Train_set = Train_set,      #ѵ��???ı???????
                                Train_label = Train_class,  #ѵ��???ķ???????
                                Train_name = "Train",       #ѵ��???ı?ǩ
                                cohortVar = "Cohort",       #GEO??id
                                classVar = classVar)        #??????��
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab=cbind(Method=row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)


##############################????AUC??ͼ##############################
#׼??ͼ?ε?????
AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)

#????AUC?ľ?ֵ?Ի???ѧϰģ?ͽ???????
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = T)
AUC_mat <- AUC_mat[names(avg_AUC),]
#??ȡ????ģ??(ѵ��??+????????AUC??ֵ????)
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

#??????ͼע?͵???ɫ
CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
names(CohortCol) <- colnames(AUC_mat)

#????ͼ??
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = AUC_mat,       #AUC????
                    avg_Cindex = avg_AUC,       #AUC??ֵ
                    CohortCol = CohortCol,      #???ݼ?????ɫ
                    barCol = "steelblue",       #?Ҳ???״ͼ????ɫ
                    cellwidth = cellwidth, cellheight = cellheight,    #??ͼÿ?????ӵĿ??Ⱥ͸߶?
                    cluster_columns = F, cluster_rows = F)      #?Ƿ??????ݽ??о???

#??????ͼ
pdf(file="model.AUCheatmap.pdf", width=cellwidth * ncol(AUC_mat) + 6, height=cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź##Rr: Windows 11 x64 (build 26200)

#Matrix products: default


#locale:
# [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
#[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
#[3] LC_MONETARY=Chinese (Simplified)_China.utf8
#[4] LC_NUMERIC=C                               
#[5] LC_TIME=Chinese (Simplified)_China.utf8    


#attached base packages:
# [1] stats     graphics  grDevices utils    
#[5] datasets  methods   base     

#loaded via a namespace (and not attached):
#  [1] KEGGREST_1.44.1         fastmatch_1.1-4        
#[3] gtable_0.3.6            ggplot2_4.0.0          
#[5] Biobase_2.64.0          lattice_0.22-6         
#[7] vctrs_0.6.5             tools_4.4.1            
#[9] generics_0.1.4          stats4_4.4.1           
#[11] parallel_4.4.1          tibble_3.3.0           
#[13] AnnotationDbi_1.66.0    RSQLite_2.3.7          
#[15] blob_1.2.4              pkgconfig_2.0.3        
#[17] Matrix_1.7-4            data.table_1.17.8      
#[19] forestploter_1.1.2      RColorBrewer_1.1-3     
#[21] S7_0.2.0                S4Vectors_0.42.0       
#[23] GenomeInfoDbData_1.2.12 lifecycle_1.0.4        
#[25] compiler_4.4.1          farver_2.1.2           
#[27] Biostrings_2.72.1       fgsea_1.30.0           
#[29] codetools_0.2-20        GenomeInfoDb_1.40.1    
#[31] preprocessCore_1.66.0   crayon_1.5.3           
#[33] pillar_1.11.0           tidyr_1.3.1            
#[35] BiocParallel_1.38.0     affy_1.82.0            
#[37] cachem_1.1.0            genefilter_1.86.0      
#[39] tidyselect_1.2.1        dplyr_1.1.4            
#[41] purrr_1.1.0             splines_4.4.1          
#[43] forcats_1.0.0           cowplot_1.2.0          
#[45] fastmap_1.2.0           grid_4.4.1             
#[47] cli_3.6.3               magrittr_2.0.3         
#[49] patchwork_1.3.1         survival_3.7-0         
#[51] XML_3.99-0.17           UCSC.utils_1.0.0       
#[53] scales_1.4.0            bit64_4.6.0-1          
#[55] XVector_0.44.0          httr_1.4.7             
#[57] affyio_1.74.0           matrixStats_1.5.0      
#[59] bit_4.6.0               gridExtra_2.3          
#[61] png_0.1-8               memoise_2.0.1          
#[63] IRanges_2.38.0          rlang_1.1.6            
#[65] Rcpp_1.1.0              xtable_1.8-4           
#[67] glue_1.8.0              DBI_1.2.3              
#[69] BiocManager_1.30.26     BiocGenerics_0.50.0    
#[71] jsonlite_2.0.0          rstudioapi_0.16.0      
#[73] annotate_1.82.0         R6_2.6.1               
#[75] MatrixGenerics_1.16.0   zlibbioc_1.50.0     

