## CPTAC proteome data preprocess
Preprocessing the CPTAC proteome data (global/phospho/acetyl/ubiquityl).

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
|                    | LSCC   | v      | v       | -      | v         | RefSeq 20180629*        |
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
Upstream data sources were collected in `~/Box/MyCPTAC/CPTAC_data_collection_v1`, which is compressed at `~/Box/Ding_Lab/Projects_Current/PTMcosmos/Data_freeze/ptmcosmos_data_freeze_v1.0/Bobo_CPTAC_data_collection_v1.zip`.


### Setup
```
conda create -n cptac \
    python=3.8 \
    htslib=1.10 samtools=1.10 \
    snakemake=5.19 \
    notebook=6.0 pandas=1.0 \
    zstd zstandard=0.14 \
    xsv xmlschema

    pysam=0.16 # not in bioconda yet
    crc64iso
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


### Ding lab internal
The project is mirrored on katmai:
- This repo: `/diskmnt/Projects/PTMcosmos/cptac_proteome_preprocess`
- Upstream data sources: `/diskmnt/Projects/PTMcosmos/Bobo_CPTAC_data_collection/CPTAC_data_collection_v1`

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