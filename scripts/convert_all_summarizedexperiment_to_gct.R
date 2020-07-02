library(SummarizedExperiment)
library(cmapR)
library(here)
library(tidyverse)

all_se_obj_pths = fs::dir_ls(here('final_output/'), glob = '*.rds', recurse = TRUE)

all_se_obj_pths %>%
    walk(function(pth){
        se = readRDS(pth)
        out_pth = fs::path_ext_set(pth, 'gct')
        g = cmapR::GCT(
            mat = assay(se),
            rdesc = rowData(se) %>% as.data.frame(),
            cdesc = colData(se) %>% as.data.frame()
        )
        cmapR::write_gct(g, out_pth, precision = 6, appenddim = FALSE)
    })



