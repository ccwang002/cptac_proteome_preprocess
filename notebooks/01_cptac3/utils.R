read_refseq_to_uniprot_coord_mapping = function(refseq_ver, uniprot_ver) {
    read_tsv(
        here(glue::glue(
            'annotation/mapping_uniprot_{uniprot_ver}/tracked_results/mappings/refseq_{refseq_ver}_to_uniprot_{uniprot_ver}_mappable.coord_mapping.tsv.gz'
        )),
        col_types = cols(
            .default = col_character(),
            source_segment_start = col_integer(),
            source_segment_end = col_integer(),
            target_segment_start = col_integer(),
            target_segment_end = col_integer(),
            target_to_source_distance = col_integer()
        )) %>%
        select(refseq_prot_id, uniprot_acc, target_uniparc_id, mapping_approach,
               source_segment_start:target_to_source_distance)
}

annotate_uniprot_mapping = function(site_row_tbl, site_col, coord_mapping_tbl) {
    mapped_to_uniprot_tbl = site_row_tbl %>%
        # Mapping to the target protein entries
        # The full peptide range (plus padding) must be included in one consecutive mapping segment
        select(unique_site_id, refseq_prot_id, orig_sites = {{ site_col }}, peptide_start, peptide_end) %>%
        left_join(coord_mapping_tbl, by = 'refseq_prot_id') %>%
        filter(case_when(
            mapping_approach == 'identical_sequence' ~ TRUE,
            mapping_approach == 'global_seq_align' ~ (
                peptide_start - 1L >= source_segment_start
                & peptide_end + 1L <= source_segment_end
            )
        )) %>%
        # Re-annotate the site
        mutate(
            mapped_sites = case_when(
                mapping_approach == 'identical_sequence' ~ {
                    glue::glue("{uniprot_acc}:{orig_sites}")
                },
                mapping_approach == 'global_seq_align' ~ {
                    # Shift all the sites position
                    # Ex. S20;S26 ----(+4)---> S24;S30
                    new_sites = map2_chr(
                        str_split(orig_sites, fixed(';')),
                        target_to_source_distance,
                        function(sites, d) {
                            new_site_locs = as.integer(str_sub(sites, start = 2L)) + d
                            str_c(str_sub(sites, end = 1L), new_site_locs, collapse = ';')
                        }
                    )
                    glue::glue("{uniprot_acc}:{new_sites}")
                }
            ),
            mapping_coord_change = case_when(
                mapping_approach == 'identical_sequence' ~ FALSE,
                mapping_approach == 'global_seq_align' ~ target_to_source_distance != 0
            )
        ) %>%
        # Only keep these columns in the site metadata
        select(
            unique_site_id,
            mapped_sites,
            mapped_uniparc_id = target_uniparc_id,
            mapping_approach,
            mapping_coord_change
        ) %>%
        # One site can be mapped to mulitple entries. Collapse their records
        group_by(unique_site_id) %>%
        summarize(
            mapped_to_multi_entries = n() > 1,
            across(
                c(mapped_sites, mapping_approach, mapped_uniparc_id, mapping_coord_change),
                str_c,
                collapse = '|'
            ),
            .groups = "drop"
        )

    site_row_tbl %>%
        left_join(mapped_to_uniprot_tbl, by = 'unique_site_id')
}