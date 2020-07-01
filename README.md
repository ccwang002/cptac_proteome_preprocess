## CPTAC proteome data preprocess

### Folder structure
```
notebooks       All the notebooks to preprocess the data
intermediates   Intermediate/temporary files
final_output    Final output files
annotation      Annotation and their processing scripts
```


### Upstream data sources
Upstream data sources were collected in `~/Box/MyCPTAC/CPTAC_data_collection_v1`.


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
- Upstream data sources: `/diskmnt/Projects/PTMcosmos/CPTAC_data/CPTAC_data_collection_v1`

Mirror script (from local to katmai):
```bash
rsync -a --info=progress2 ~/Box/MyCPTAC/ katmai:/diskmnt/Projects/PTMcosmos/CPTAC_data/

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