### 读入Targets
library(limma)
#Targets <- readTargets("affymetrix/GSE29450/Targets.txt")  # 不指定行名所在的列
Targets <- readTargets("affymetrix/GSExxx/Targets.txt",
                       row.names = "sample_id")            # 指定行名所在的列

### 读入CEL
library(affy)
cel_files <- affy::list.celfiles("affymetrix/GSExxx/GSExxx_RAW/", 
                                 full.name = T)
cel_files

GSExxx_cel <- ReadAffy(filenames = Targets$FileName,
                         celfile.path = "affymetrix/GSExxx/GSExxx_RAW",
                         phenoData = Targets)

cel <- GSExxx_cel
# 表达矩阵 exprs
head(exprs(cel))[, 1:5]
# 样本信息 pData
GSExxx_targets <- pData(cel)
# 样本名称 sampeNames
sampleNames(cel)
# ScanDate
cel@protocolData@data$ScanDate
# expression matrix
GSExxx_expr_cel <- exprs(cel)

### 保存数据，备下一步分析用
save(GSExxx_cel, file = "affymetrix/GSExxx/GSExxx_cel.Rdata")
load_input <- load("affymetrix/GSExxx/GSExxx_cel.Rdata")
load_input

library(affy)
GSExxx_expr_cel <- exprs(GSExxx_cel)#提取表达矩阵
head(GSExxx_expr_cel)[, 1:5]

GSExxx_targets <- pData(GSExxx_cel) #提取targets（样本分组）数据，在phenoData中的data


###  RMA ----
GSExxx_rma <- affy::rma(GSExxx_cel)#用的是原始cel数据，不是背景矫正后的数据
GSExxx_expr_rma <- exprs(GSExxx_rma)
write.table(GSExxx_expr_rma,file="GSExxx_expr_rma.txt",sep="\t",quote=T,row.names = T)

