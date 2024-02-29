devtools::install_github("saeyslab/synthspot")

library(Seurat)
library(synthspot)
library(dplyr)

# scrnaseq data
DimPlot(seurat_obj, group.by = "subclass", label = T)

set.seed(1)
synthetic_visium_data = generate_synthetic_visium(seurat_obj = seurat_obj, dataset_type = "artificial_diverse_distinct", 
                                                  clust_var = "subclass", n_regions = 5, n_spots_min = 5, n_spots_max = 20, 
                                                  visium_mean = 30000, visium_sd = 8000)

synthetic_visium_data$counts %>% as.matrix() %>% .[1:5,1:5]
##      priorregion1_spot_1 priorregion1_spot_2 priorregion1_spot_3 priorregion1_spot_4 priorregion1_spot_5
## Vip                 5287               11155                6689                4996                4294
## Sst                    0                   0                   0                   6                   0
## Npy                    2                   0                   0                   0                   1
## Tac2                   0                1822                 815                   0                 154
## Crh                   11                 920                1300                  83                 163

synthetic_visium_data$spot_composition %>% .[1:10,]
##    L2.3.IT L4 L5.IT L5.PT L6.CT L6.IT L6b Lamp5 NP Pvalb Sst Vip                name       region
## 1        0  0     0     0     0     0   0     0  7     0   0   1 priorregion1_spot_1 priorregion1
## 2        0  0     0     0     0     0   0     0  1     0   0   1 priorregion1_spot_2 priorregion1
## 3        0  0     0     0     0     0   0     0  7     0   0   1 priorregion1_spot_3 priorregion1
## 4        0  0     0     0     0     0   0     0  3     0   0   1 priorregion1_spot_4 priorregion1
## 5        0  0     0     0     0     0   0     0  6     0   0   3 priorregion1_spot_5 priorregion1
## 6        0  0     0     0     0     0   0     0  2     0   0   0 priorregion1_spot_6 priorregion1
## 7        0  0     0     0     0     9   0     0  0     0   0   0 priorregion2_spot_1 priorregion2
## 8        0  0     0     0     0     8   0     0  0     0   0   0 priorregion2_spot_2 priorregion2
## 9        0  0     0     0     0     8   0     0  0     0   0   0 priorregion2_spot_3 priorregion2
## 10       0  0     0     0     0     8   0     0  0     0   0   0 priorregion2_spot_4 priorregion2

synthetic_visium_data$relative_spot_composition %>% .[1:10,]
##    L2.3.IT L4 L5.IT L5.PT L6.CT L6.IT L6b Lamp5        NP Pvalb Sst       Vip                name       region
## 1        0  0     0     0     0     0   0     0 0.8750000     0   0 0.1250000 priorregion1_spot_1 priorregion1
## 2        0  0     0     0     0     0   0     0 0.5000000     0   0 0.5000000 priorregion1_spot_2 priorregion1
## 3        0  0     0     0     0     0   0     0 0.8750000     0   0 0.1250000 priorregion1_spot_3 priorregion1
## 4        0  0     0     0     0     0   0     0 0.7500000     0   0 0.2500000 priorregion1_spot_4 priorregion1
## 5        0  0     0     0     0     0   0     0 0.6666667     0   0 0.3333333 priorregion1_spot_5 priorregion1
## 6        0  0     0     0     0     0   0     0 1.0000000     0   0 0.0000000 priorregion1_spot_6 priorregion1
## 7        0  0     0     0     0     1   0     0 0.0000000     0   0 0.0000000 priorregion2_spot_1 priorregion2
## 8        0  0     0     0     0     1   0     0 0.0000000     0   0 0.0000000 priorregion2_spot_2 priorregion2
## 9        0  0     0     0     0     1   0     0 0.0000000     0   0 0.0000000 priorregion2_spot_3 priorregion2
## 10       0  0     0     0     0     1   0     0 0.0000000     0   0 0.0000000 priorregion2_spot_4 priorregion2

synthetic_visium_data$gold_standard_priorregion %>% head()
## # A tibble: 6 x 4
##   prior_region celltype  freq present
##   <chr>        <chr>    <dbl> <lgl>  
## 1 priorregion1 Vip      0.143 TRUE   
## 2 priorregion1 NP       0.857 TRUE   
## 3 priorregion1 Sst      0     FALSE  
## 4 priorregion1 L4       0     FALSE  
## 5 priorregion1 L5 IT    0     FALSE  
## 6 priorregion1 L6 IT    0     FALSE

synthetic_visium_data$dataset_properties 
## # A tibble: 1 x 10
##   dataset_id                   dataset_type                real_artificial uniform_diverse distinct_overlap dominant_celltype missing_celltype rare_celltype missing_celltype_sc real_region_var
##   <chr>                        <chr>                       <chr>           <chr>           <chr>            <lgl>             <lgl>            <lgl>         <lgl>               <lgl>          
## 1 artificial_diverse_distinct1 artificial_diverse_distinct artificial      diverse         distinct         FALSE             FALSE            FALSE         NA                  NA



