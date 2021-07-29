rename_gbm_samples = function(sample_type) {
    case_when(
        sample_type == 'GTEx Normal' ~ 'normal',
        sample_type == 'Primary Tumor' ~ 'tumor'
    )
}

rename_ucec_and_ov_samples = function(sample_type) {
    case_when(
        sample_type == 'Tumor' ~ 'tumor',
        endsWith(sample_type, 'normal') ~ 'normal'
    )
}

rename_default_samples = function(sample_type) {
    case_when(
        sample_type == 'Primary Tumor' ~ 'tumor',
        endsWith(sample_type, 'Solid Tissue Normal') ~ 'normal'
    )
}

rename_chop_samples = function(sample_type) {
    case_when(
        endsWith(sample_type, 'Tumor') ~ 'tumor',
        endsWith(sample_type, 'Normal') ~ 'normal'
    )
}

rename_sample_type = function(se_list) {
    se_list %>% imap(function(se, cohort) {
        if (cohort == "CPTAC3_GBM_discovery") {
            colData(se)$simplified_sample_type = rename_gbm_samples(colData(se)$sample_type)
        } else if (cohort == "CPTAC3_UCEC_discovery" || cohort == "TCGA_OV_retrospective") {
            colData(se)$simplified_sample_type = rename_ucec_and_ov_samples(colData(se)$sample_type)
        } else if (cohort == "CPTAC3_PBTA_discovery" || cohort == "CPTAC3_HOPEAYA_discovery") {
            colData(se)$simplified_sample_type = rename_chop_samples(colData(se)$sample_type)
        } else {
            colData(se)$simplified_sample_type = rename_default_samples(colData(se)$sample_type)
        }
        se
    })
}

rename_all_proteome = function(all_proteome) {
    all_proteome = map(all_proteome, rename_sample_type)
    all_proteome
}
