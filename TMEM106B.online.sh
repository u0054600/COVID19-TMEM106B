
###############################################################################################################
library(Seurat)
library(cowplot)
library(grid)
library(gridExtra)
library(ggplot2)
library(lattice)


median.stat <- function(x){
   out <- quantile(x, probs = c(0.5))
   names(out) <- c("ymed")
   return(out) 
}

median.stat75 <- function(x){
   out <- quantile(x, probs = c(0.75))
   names(out) <- c("ymed")
   return(out) 
}

`%ni%`<- Negate(`%in%`)

###############################################################################################################

### Data are available at http://covid19.lambrechtslab.org/
data <- readRDS(file= "OnlineTMEM106B.rds")


###############################################################################################################

pdf("Figure6A.pdf")
DimPlot(data, group.by="Domain")
FeaturePlot(data, features="TMEM106B", max.cutoff=2.75,col=c(rgb (0.8,0.8,0.8,0.11),rgb (1,0,0,1)))
dev.off()


###############################################################################################################


	

pdf("Figure_6B_CovidInfections_spikeProtein.pdf")
data$Domain_FIN  <- factor(data$Domain_FIN , levels=c("Epithelial_control", "Epithelial_COVID19", "Epithelial_COVID19_Infected","Epithelial_COVID19_Infected_S",
"Lymphoid_control", "Lymphoid_COVID19", "Lymphoid_COVID19_Infected","Lymphoid_COVID19_Infected_S",
"Myeloid_control",  "Myeloid_COVID19", "Myeloid_COVID19_Infected" , "Myeloid_COVID19_Infected_S" ))
Idents(data ) <- 'Domain_FIN'
VlnPlot(data, features="TMEM106B", pt.size=0.1)+ stat_summary(fun.y = mean, geom='point', size = 10, colour = "green", shape = 95) + NoLegend() + stat_summary(fun.y = median.stat75, geom='point', size = 10, colour = "red", shape = 95)
dev.off()

######### Start Statistics related to Figure 6B


DF <- data@meta.data

TYPE <- "Lymphoid"
Epi <- DF [grep(TYPE, DF$ Domain_FIN),]
Epi$ Domain_FIN<- factor(Epi$ Domain_FIN , levels=c( paste0(TYPE,"_control"),  paste0(TYPE,"_COVID19"),  paste0(TYPE,"_COVID19_Infected"),  paste0(TYPE,"_COVID19_Infected_S")))
kruskal.test(TMEM106B1 ~ Domain_FIN, data = Epi)
pairwise.wilcox.test(DF$TMEM106B1, DF$Domain_FIN,
                  p.adjust.method = "none")

TYPE <- "Myeloid"
Epi <- DF [grep(TYPE, DF$ Domain_FIN),]
Epi$ Domain_FIN<- factor(Epi$ Domain_FIN , levels=c( paste0(TYPE,"_control"),  paste0(TYPE,"_COVID19"),  paste0(TYPE,"_COVID19_Infected"),  paste0(TYPE,"_COVID19_Infected_S")))
kruskal.test(TMEM106B1 ~ Domain_FIN, data = Epi)
pairwise.wilcox.test(DF$TMEM106B1, DF$Domain_FIN,
                  p.adjust.method = "none")

TYPE <- "Epithelial"
Epi <- DF [grep(TYPE, DF$ Domain_FIN),]
Epi$ Domain_FIN<- factor(Epi$ Domain_FIN , levels=c( paste0(TYPE,"_control"),  paste0(TYPE,"_COVID19"),  paste0(TYPE,"_COVID19_Infected"),  paste0(TYPE,"_COVID19_Infected_S")))
kruskal.test(TMEM106B1 ~ Domain_FIN, data = Epi)
pairwise.wilcox.test(DF$TMEM106B1, DF$Domain_FIN,
                  p.adjust.method = "none")

######### Stop Statistics related to Figure 6B ################################################################################################################################

 
 


 pdf("Figure6C_Epithelial.pdf")
  data$Infec<- factor (data$Infec, levels= c("AT2_control","AT2_COVID19","AT2_COVID19_Infected","Basal_control","Basal_COVID19", "Basal_COVID19_Infected",
 "Ciliated_control","Ciliated_COVID19","Ciliated_COVID19_Infected", "Hillock_control","Hillock_COVID19","Hillock_COVID19_Infected",
  "Inflammatory_control","Inflammatory_COVID19","Inflammatory_COVID19_Infected", 
 "Ionocyte_control","Ionocyte_COVID19","Ionocyte_COVID19_Infected",
  "Secretory_control","Secretory_COVID19","Secretory_COVID19_Infected"))

 Idents(data)<- 'Infec'
 VlnPlot(data, features="TMEM106B", pt.size=0.1)+ stat_summary(fun.y = mean, geom='point', size = 10, colour = "green", shape = 95) + NoLegend()+ stat_summary(fun.y = median.stat75, geom='point', size = 10, colour = "red", shape = 95)
 #### "NA" are the non-epithelial cells and are allready shown in Figure 6B
 dev.off()

 
 


 ########## Start Statistics Figure 6C
  Pln<-c()
 for (TYPE in c("AT2","Basal","Ciliated","Hillock", "Inflammatory","Ionocyte","Secretory","Lymphoid","Myeloid"))
{ 
Epi <- DF [grep(TYPE, DF$ Infec),]
Epi$ Infec<- factor(Epi$ Infec , levels=c( paste0(TYPE,"_control"),  paste0(TYPE,"_COVID19"),  paste0(TYPE,"_COVID19_Infected")))

# Compute the analysis of variance
res.kkn <- kruskal.test(TMEM106B1 ~ as.numeric(Infec), data = Epi)

Pln<- cbind(Pln,  res.kkn$p.value)
}
colnames(Pln) <- c("AT2","Basal","Ciliated","Hillock", "Inflammatory","Ionocyte","Secretory","Lymphoid","Myeloid")


#             AT2     Basal     Ciliated     Hillock Inflammatory  Ionocyte    Secretory   
# [1,] 0.06809637 0.1619966 1.501615e-11 0.009225837    0.2454293 0.9191453 4.268949e-05 


pairwise.wilcox.test(DF$TMEM106B1, DF$Infec,
                  p.adjust.method = "BH")

pairwise.wilcox.test(DF$TMEM106B1, DF$Infec,
                  p.adjust.method = "none")
########## Stop Statistics Figure 6C
 
 
 
 

 R version 3.6.0 (2019-04-26)
 Platform: x86_64-apple-darwin15.6.0 (64-bit)
 Running under: macOS  10.16

 Matrix products: default
 BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
 LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

 locale:
 [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

 attached base packages:
 [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

 other attached packages:
 [1] ggplot2_3.3.3 gridExtra_2.3 cowplot_1.1.1 Seurat_3.2.3 

 loaded via a namespace (and not attached):
   [1] Rtsne_0.15                  colorspace_2.0-0            deldir_0.2-3                ellipsis_0.3.1              ggridges_0.5.3              XVector_0.24.0             
   [7] GenomicRanges_1.36.1        spatstat.data_1.7-0         farver_2.0.3                leiden_0.3.6                listenv_0.8.0               ggrepel_0.9.0              
  [13] codetools_0.2-18            splines_3.6.0               polyclip_1.10-0             jsonlite_1.7.2              ica_1.0-2                   cluster_2.1.0              
  [19] png_0.1-7                   uwot_0.1.10                 shiny_1.5.0                 sctransform_0.3.2           compiler_3.6.0              httr_1.4.2                 
  [25] Matrix_1.3-2                fastmap_1.0.1               lazyeval_0.2.2              later_1.1.0.1               htmltools_0.5.0             tools_3.6.0                
  [31] rsvd_1.0.3                  igraph_1.2.6                gtable_0.3.0                glue_1.4.2                  GenomeInfoDbData_1.2.1      RANN_2.6.1                 
  [37] reshape2_1.4.4              dplyr_1.0.2                 Rcpp_1.0.5                  spatstat_1.64-1             scattermore_0.7             Biobase_2.44.0             
  [43] vctrs_0.3.6                 nlme_3.1-151                lmtest_0.9-38               stringr_1.4.0               globals_0.14.0              mime_0.9                   
  [49] miniUI_0.1.1.1              lifecycle_0.2.0             irlba_2.3.3                 goftest_1.2-2               future_1.21.0               zlibbioc_1.30.0            
  [55] MASS_7.3-53                 zoo_1.8-8                   scales_1.1.1                promises_1.1.1              spatstat.utils_1.20-2       parallel_3.6.0             
  [61] SummarizedExperiment_1.14.1 RColorBrewer_1.1-2          reticulate_1.18             pbapply_1.4-3               rpart_4.1-15                stringi_1.5.3              
  [67] S4Vectors_0.22.1            BiocGenerics_0.30.0         BiocParallel_1.18.1         GenomeInfoDb_1.20.0         rlang_0.4.10                pkgconfig_2.0.3            
  [73] matrixStats_0.57.0          bitops_1.0-6                lattice_0.20-41             ROCR_1.0-11                 purrr_0.3.4                 tensor_1.5                 
  [79] labeling_0.4.2              patchwork_1.1.1             htmlwidgets_1.5.3           tidyselect_1.1.0            parallelly_1.23.0           RcppAnnoy_0.0.18           
  [85] plyr_1.8.6                  magrittr_2.0.1              R6_2.5.0                    IRanges_2.18.3              generics_0.1.0              DelayedArray_0.10.0        
  [91] withr_2.3.0                 pillar_1.4.7                mgcv_1.8-33                 fitdistrplus_1.1-3          survival_3.2-7              abind_1.4-5                
  [97] RCurl_1.98-1.2              tibble_3.0.4                future.apply_1.7.0          crayon_1.3.4                KernSmooth_2.23-18          plotly_4.9.3               
 [103] data.table_1.13.6           digest_0.6.27               xtable_1.8-4                tidyr_1.1.2                 httpuv_1.5.4                stats4_3.6.0               
 [109] munsell_0.5.0               viridisLite_0.3.0          
 