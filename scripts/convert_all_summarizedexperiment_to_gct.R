library(future.apply)
plan(multicore(workers = 6L))  # parallel cores

library(SummarizedExperiment)
library(cmapR)
library(here)
library(fs)
library(dplyr)

all_se_obj_pths = fs::dir_ls(here('final_output/'), glob = '*.rds', recurse = TRUE)

# Read SummarizedExperiment object from pth and export as GCT
se_to_gct = function(pth) {
    se = readRDS(pth)
    out_pth = fs::path_ext_set(pth, 'gct')
    g = cmapR::GCT(
        mat = assay(se),
        rdesc = rowData(se) %>% as.data.frame(),
        cdesc = colData(se) %>% as.data.frame()
    )
    cmapR::write_gct(g, out_pth, precision = 6, appenddim = FALSE)
}

y = future_lapply(all_se_obj_pths, se_to_gct)
