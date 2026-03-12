rsFile="model.classMatrix.txt"      #?????ľ????ļ?
method="glmBoost+RF"                #ѡ??????ѧϰ?ķ???(??Ҫ??????ͼ?????޸?)
setwd("=")     #???ù???Ŀ¼

#??ȡ?????ľ????ļ?
riskRT=read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)
#??ȡ???ݼ???ID
CohortID=gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID=gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort=CohortID

#?????????????ĺ???
draw_confusion_matrix <- function(cm=null, titleName=null) {
#	layout(matrix(c(1,1,2)))
	par(mar=c(2,2,3,2))
	plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
	title(paste0('CONFUSION MATRIX (', titleName,')'), cex.main=1.5)
	
	# create the matrix 
	rect(150, 430, 240, 370, col='#3F97D0')
	text(195, 435, 'Control', cex=1.2)
	rect(250, 430, 340, 370, col='#F7AD50')
	text(295, 435, 'Treat', cex=1.2)
	text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
	text(245, 450, 'Actual', cex=1.3, font=2)
	rect(150, 305, 240, 365, col='#F7AD50')
	rect(250, 305, 340, 365, col='#3F97D0')
	text(140, 400, 'Control', cex=1.2, srt=90)
	text(140, 335, 'Treat', cex=1.2, srt=90)
	  
	# add in the cm results 
	res <- as.numeric(cm)
	text(195, 400, res[1], cex=1.6, font=2, col='white')
	text(195, 335, res[2], cex=1.6, font=2, col='white')
	text(295, 400, res[3], cex=1.6, font=2, col='white')
	text(295, 335, res[4], cex=1.6, font=2, col='white')
}

#?????ݼ?????ѭ??
for(Cohort in unique(riskRT$Cohort)){
	#??ȡ??Ʒ?ķ?????Ϣ(????????ʵ????)
	rt=riskRT[riskRT$Cohort==Cohort,]
	y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
	y=ifelse(y=="Control", 0, 1)
	
	#????????
	result_matrix=table(rt[,method], y)
	pdf(file=paste0("confusion.", Cohort, ".pdf"), width=6, height=5)
	draw_confusion_matrix(cm=result_matrix, titleName=Cohort)
	dev.off()
}

#sessionInfo()
#R version 4.4.1 (2024-06-14 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 26200)

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


