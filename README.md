## CPTAC proteome data preprocess

### Folder structure
```
intermediates
final_output
annotation
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
    zstd zstandard=0.14

    pysam=0.16 # not in bioconda yet
```

Additionally it needs R 4.0 and Bioconductor 3.11.


### Ding lab internal
The project is mirrored on katmai:
- This repo: `/diskmnt/Projects/PTMcosmos/cptac_proteome_preprocess`
- Upstream data sources: `/diskmnt/Projects/PTMcosmos/CPTAC_data/CPTAC_data_collection_v1`

Mirror script (from local to katmai):
```
rsync -a --info=progress2 ~/Box/MyCPTAC/ katmai:/diskmnt/Projects/PTMcosmos/CPTAC_data/
```