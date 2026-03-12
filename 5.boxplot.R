

#???ð?
library(reshape2)
library(ggplot2)

setwd("")     #???ù???Ŀ¼

#????????ͼ?ĺ???
bioBoxplot=function(inputFile=null, outFile=null, titleName=null){
	#??ȡ?????ļ?,??ȡ????
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(rt)
	Project=gsub("(.*?)\\_.*", "\\1", rownames(data))            #??ȡ?????о???id
	Sample=gsub("(.+)\\_(.+)\\_(.+)", "\\2", rownames(data))     #??ȡ??Ʒ??????
	data=cbind(as.data.frame(data), Sample, Project)
	
	#??????ת????ggplot2?????ļ?
	rt1=melt(data, id.vars=c("Project", "Sample"))
	colnames(rt1)=c("Project","Sample","Gene","Expression")

	#????????ͼ
	pdf(file=outFile, width=10, height=5)
	p=ggplot(rt1, mapping=aes(x=Sample, y=Expression))+
  		geom_boxplot(aes(fill=Project), notch=T, outlier.shape=NA)+
  		ggtitle(titleName)+ theme_bw()+ theme(panel.grid=element_blank())+ 
  		theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5,size=2), plot.title=element_text(hjust = 0.5))
	print(p)
	dev.off()
}

#???ú???, ???????ν???ǰ??????ͼ
bioBoxplot(inputFile="merge.preNorm.txt", outFile="boxplot.preNorm.pdf", titleName="Before batch correction")
#???ú???, ???????ν???????????ͼ
bioBoxplot(inputFile="merge.normalize.txt", outFile="boxplot.normalzie.pdf", titleName="After batch correction")

#sessionInfo()
#R version 4.4.1 (2024-06-14 ucrt)
#Platform: x86_64-w64-mingw32/x64
#Running under: Windows 11 x64 (build 26200)
#locale:
#[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
#[2] LC_CTYPE=Chinese (Simplified)_China.utf8   
#[3] LC_MONETARY=Chinese (Simplified)_China.utf8
#[4] LC_NUMERIC=C                               
#[5] LC_TIME=Chinese (Simplified)_China.utf8    
#tzcode source: internal

#attached base packages:
# [1] stats     graphics  grDevices utils    
#[5] datasets  methods   base     

#loaded via a namespace (and not attached):
# [1] KEGGREST_1.44.1         fastmatch_1.1-4        
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



