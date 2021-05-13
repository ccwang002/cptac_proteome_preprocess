require(SummarizedExperiment)
require(here)
require(fs)
require(dplyr)

DATA_FREEZE_VER = 'v3.0'
DATA_FREEZE_ROOT = here(glue::glue(
    '~/Box/Ding_Lab/Projects_Current/PTMcosmos/Data_freeze/ptmcosmos_data_freeze_{DATA_FREEZE_VER}'
))

ORDERED_COHORTS = c(
    "TCGA_BRCA_retrospective",
    "TCGA_OV_retrospective",
    "CPTAC2_BRCA_prospective",
    "CPTAC2_CRC_prospective",
    "CPTAC2_OV_prospective",
    "CPTAC3_CCRCC_discovery",
    "CPTAC3_GBM_discovery",
    "CPTAC3_HNSCC_discovery",
    "CPTAC3_LSCC_discovery",
    "CPTAC3_LUAD_discovery",
    "CPTAC3_UCEC_discovery"
)

ORDERED_DATA_TYPES = c(
    'proteome',
    'phosphoproteome',
    'acetylome',
    'ubiquitylome'
)

locate_all_proteome_rds = function() {
    cptac_proteome_tbl = tibble(
        rds_pth = fs::dir_ls(
            fs::path(DATA_FREEZE_ROOT, 'cptac_proteome'),
            recurse = TRUE,
            glob = '*.rds'
        )
    ) %>%
        extract(
            rds_pth, c('pipeline', 'data_type'), '([^_]+)_([^_]+ome)\\.v[0-9.]+\\.rds$',
            remove = FALSE
        ) %>%
        mutate(
            cohort = factor(rds_pth %>% fs::path_dir() %>% fs::path_file(),
                            ORDERED_COHORTS),
            data_type = factor(data_type, ORDERED_DATA_TYPES)
        ) %>%
        select(cohort, data_type, pipeline, rds_pth) %>%
        arrange(cohort, data_type)

    cptac_proteome_tbl
}


# Read in all the SummarizedExperiment RDS objects of one data type
read_one_data_type_rds = function(cptac_proteome_tbl, data_type) {
    one_dt_tbl = cptac_proteome_tbl %>%
        filter(data_type == {{ data_type }})

    se_per_cohort = one_dt_tbl$rds_pth %>%
        set_names(one_dt_tbl$cohort) %>%
        map(readRDS)

    se_per_cohort
}


# Read all CPTAC proteome data
read_all_cptac_proteome = function(cptac_proteome_tbl) {
    all_cptac_proteome = list(
        global = read_one_data_type_rds(cptac_proteome_tbl, 'proteome'),
        phospho = read_one_data_type_rds(cptac_proteome_tbl, 'phosphoproteome'),
        acetyl = read_one_data_type_rds(cptac_proteome_tbl, 'acetylome'),
        ubiquityl = read_one_data_type_rds(cptac_proteome_tbl, 'ubiquitylome')
    )
    all_cptac_proteome
}
