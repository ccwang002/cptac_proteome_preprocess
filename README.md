## CPTAC proteome data preprocess

### Folder structure
```
intermediates
final_output
annotation
```

Data and upstream sources were collected in `~/Box/Projects/2020_CPTAC_data_collection`


### Setup

R 4.0 and Bioconductor 3.11.

```
conda create -n cptac \
    python=3.8 \
    htslib=1.10 samtools=1.10
```