## CPTAC proteome data preprocess
Preprocessing the CPTAC proteome data (global/phospho/acetyl/ubiquityl). The final outputs become `cptac_proteome` of the PTMcosmos data freeze.

Currently it covers (`-` = not available):

|    CPTAC stage     | Cohort | Global | Phospho | Acetyl | Ubiquityl | Peptide search database |
| :----------------- | ------ | ------ | ------- | ------ | --------- | ----------------------- |
| CPTAC2/TCGA        | BRCA   | v      | v       | -      | -         | RefSeq 20130727         |
|                    | OV     | v      | v       | -      | -         | RefSeq 20111201         |
| CPTAC2 prospective | BRCA   | v      | v       | v      | -         | RefSeq 20160914         |
|                    | CRC    | v      | v       | -      | -         | RefSeq 20171003         |
|                    | OV     | v      | v       | -      | -         | RefSeq 20160914         |
| CPTAC3 discovery   | CCRCC  | v      | v       | -      | -         | RefSeq 20180629         |
|                    | HNSCC  | v      | v       | -      | -         | RefSeq 20180629         |
|                    | LUAD   | v      | v       | v      | -         | RefSeq 20180629*        |
|                    | LSCC   | v      | v       | v      | v         | RefSeq 20180629*        |
|                    | GBM    | v      | v       | v      | -         | RefSeq 20180629         |
|                    | UCEC   | v      | v       | v      | -         | RefSeq 20180629         |

Notes for the peptide databases (bold ones are recommended):
- RefSeq 20111201: CDAP (RefSeq release 37); Incomplete annotation; hg19
- RefSeq 20130727: Incomplete annotation; hg19
- RefSeq 20171003: Currently overwrite the annotation using RefSeq 20180629; hg38
- **RefSeq 20160914**: CDAP (Refseq 2016); hg19
- **RefSeq 20180629**: CDAP (RefSeq 2018); hg38
- RefSeq 20180629*: database includes smORFs


### Folder structure
```
annotation      Annotation and their processing scripts
notebooks       All the notebooks to preprocess the data
intermediates   Intermediate/temporary files
final_output    Final output files
```


### Upstream data sources
Upstream data sources were collected in `~/Box/MyCPTAC/CPTAC_proteome_v3.0`, which is compressed at `~/Box/Ding_Lab/Projects_Current/PTMcosmos/Data_freeze/ptmcosmos_data_freeze_v3.0/Bobo_CPTAC_proteome_v3.0.zip`.


### Setup
```
conda create -n cptac \
    python=3.8 \
    htslib=1.10 samtools=1.10 \
    snakemake=5.19 \
    notebook=6.0 pandas=1.1 \
    zstd zstandard=0.14 \
    xsv xmlschema orjson loguru \
    biopython=1.79 pysam=0.16 crc64iso=0.0.2 \
    "'anndata>=0.7,<0.8'" "'h5py=*=nompi*'" "'hdf5>=1.10=nompi*'" "'pyarrow>=4.0'"
```

Additionally it needs R 4.0 and Bioconductor 3.11.


### Run the workflow
For the first run, one should run the R and Python notebooks in order.

For subsequence runs, where intermediates files are available,
one can update all the final output files by running all the R notebooks in parallel:

```bash
parallel -j4 --bar \
    "R -e \"rmarkdown::render('{}')\"" \
    ::: ??_*.Rmd
```


### Processed data format
Column (sample) metadata:
- (column names): preferred sample names that matches the corresponding flagship paper. (eg C3L-00104)
- `case_id`, `aliquot_id`: CPTAC IDs. (eg C3L-00104, CPT0092440003)
- `tmt_experiment`, `tmt_channel`: TMT experiment run and channel (eg 11, 129C)
- `sample_type`: Sample type (tumor or normal). (eg Primary Tumor)

Row (feature; protein) metadata:
- (row names): gene symbol for protein (eg EGFR)
- `refseq_prot_id`: RefSeq Protein ID and version (eg NP_005219.2)
- `hgnc_id`: HGNC gene ID (eg HGNC:3236)
- `entrez_gene_id`: Entrez gene ID (eg 1956)
- `ensembl_gene_id`: Ensembl gene ID (eg ENSG00000146648)

Row (feature; PTM) metadata:
- (row names): a unique standardized naming of the PTM sites.
  It follows the format `<symbol>:<refseq_prot_id>:<sites>[.duplication]`
  (eg EGFR:NP_005219.2:S1064 and EGFR:NP_005219.2:S1064.1)
- `original_id`: the original site ID that matches the corresponding flagship paper.
  (eg NGLQSCPIKEDS*FLQR)
- `refseq_prot_id`: (same as protein metadata above) (eg NP_005219.2)
- `symbol`: gene symbol (eg EGFR)
- `phosphosites`, `acetylsites`, or `ubiquitylsites`: the PTM sites.
   (eg S1064 or S618;T620;S624)
- `peptide`: the underlying peptide sequence. PTM sites are in letter case.
   (eg EDsFLQR)
- `peptide_start`, `peptide_end`: the start and end amino acid position of the peptide on the protein.
   (eg 1062, 1068)
- `refseq_tx_id`: Corresponding RefSeq transcript ID of the protein.
   (eg NM_005228)
- `uniparc_id`: UniParc ID of the protein.
    This ID is the best way to check if the underlying protein sequence is identical to other IDs from a different annotation systems.
    (eg UPI000003E750)
- `hgnc_id`, `entrez_gene_id`, `ensembl_gene_id`: (same as protein metadata above)
   (eg HGNC:3236, 1956, ENSG00000146648)
- `multi_genomic_loci`: whether the protein can be mapped to multiple genomic locations.
   (eg FALSE)
- `num_exons`, `aa_len`, `tx_len`, `cds_len`: Details about the transcript and protein.

- `uniparc_crc64_checksum`: UniParc sequence CRC64 checksum. (eg D8A2A50B4EFB6ED2)
- `uniparc_xref_refseq_prot_ids`: All RefSeq protein IDs with the identical protein sequence cross-referenced using UniParc.
   (eg NP_005219.2)
- `uniparc_xref_ensembl_prot_ids`: All Ensembl protein IDs with the identical protein sequence cross-referenced using UniParc.
   (eg ENSP00000275493.1)
- `uniparc_xref_uniprot_ids`: All UniProt IDs with the identical protein sequence cross-referenced using UniParc.
   (eg P00533-1.1;P00533.2)

- `mapped_to_multi_entries`: Whether the peptide sequence can be mapped to multiple UniProt entries.
   (eg FALSE)
- `mapped_sites`: The mapped PTM sites to UniProt. It has the format `<uniprot_id>:<new_site_locs>`.
   If `mapped_to_multi_entries` is TRUE, corresponding results are concatenated by `|`.
   (eg P00533:S1064 or E9PAV3:T2037|Q13765:T174)
- `mapping_approach`: How the peptide sequence is mapped from RefSeq to UniProt.
   Possible options:

   - `identical_sequence`: Two sequences are identical
   - `global_seq_align`: Two sequences are globally aligned

- `mapped_uniparc_id`: The UniParc IDs of the mapped UniProt entry.
   (eg UPI000003E750)
- `mapping_coord_change`: Whether the PTM site coordinates are changed.
   (eg FALSE)


### Ding lab internal
The project is mirrored on katmai:
- This repo: `/diskmnt/Projects/PTMcosmos/cptac_proteome_preprocess`
- Upstream data sources: `/diskmnt/Projects/PTMcosmos/Bobo_CPTAC_data_collection/CPTAC_proteome_v3.0`

Mirror script (from local to katmai):
```bash
# Add -n for dry run
rsync -arv --delete --files-from=remote_rsync_files.list --exclude='.DS_Store' \
    ./ katmai:/diskmnt/Projects/PTMcosmos/cptac_proteome_preprocess/
```

Mirror script (from katmai to local):
```bash
# Add -n for dry run
rsync -arv --files-from=remote_rsync_files.list \
    katmai:/diskmnt/Projects/PTMcosmos/cptac_proteome_preprocess/ .
```
