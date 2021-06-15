## UniProt release 2021_02 reviewed human proteome
The entries in use are defined by this query:

    reviewed:yes AND organism:"Homo sapiens (Human) [9606]" AND proteome:up000005640

See PTMcosmos preprocessing for the initial steps


## Create mapping from RefSeq annotations
Run `notebooks/{02-05}_map_refseq*.Rmd` to create the mappings.

Run sequence alignment per RefSeq annotation version:

```bash
# RefSeq 2018
python 3_run_alignment.py \
    tracked_results/mappings/refseq_20180629_to_uniprot_2021_02_mappable.tsv.gz \
    ~/Box/MyCPTAC/CPTAC_proteome_v2.1/DCC/RefSeq_20180629/RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2021_02.fasta.gz  \
    tracked_results/mappings/refseq_20180629_to_uniprot_2021_02_mappable.coord_mapping.tsv.gz


# RefSeq 2016
python 3_run_alignment.py \
    tracked_results/mappings/refseq_20160914_to_uniprot_2021_02_mappable.tsv.gz \
    ~/Box/MyCPTAC/CPTAC_proteome_v2.1/DCC/RefSeq_20160914/RefSeq.20160914_Human_ucsc_hg19_customProDBnr_mito_150contams.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2021_02.fasta.gz  \
    tracked_results/mappings/refseq_20160914_to_uniprot_2021_02_mappable.coord_mapping.tsv.gz


# RefSeq 2013
python 3_run_alignment.py \
    tracked_results/mappings/refseq_20130727_to_uniprot_2021_02_mappable.tsv.gz \
    ~/Box/MyCPTAC/CPTAC_proteome_v2.1/DCC/RefSeq_20130727/RefSeq.20130727-Human.contams.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2021_02.fasta.gz  \
    tracked_results/mappings/refseq_20130727_to_uniprot_2021_02_mappable.coord_mapping.tsv.gz \
    --remove-source-id-version


# RefSeq 2011
python 3_run_alignment.py \
    tracked_results/mappings/refseq_20111201_to_uniprot_2021_02_mappable.tsv.gz \
    ~/Box/MyCPTAC/CPTAC_proteome_v2.1/DCC/RefSeq_20111201/20111201_RefSeq_Human_37-Mouse_37_Trypsin.renamed.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2021_02.fasta.gz  \
    tracked_results/mappings/refseq_20111201_to_uniprot_2021_02_mappable.coord_mapping.tsv.gz
```


## Create mapping from GENCODE annotations
Run `notebooks/{06}_map_gencode*.Rmd` to create the mappings.

Run sequence alignment per GENCODE annotation version:

```bash
# GENCODE 34
python 4_run_alignment_for_gencode.py \
    tracked_results/mappings/gencode_34_to_uniprot_2021_02_mappable.tsv.gz \
    ../../annotation/gencode_34/tracked_results/gencode_34_unique_protein_seq.fasta.gz \
    tracked_results/uniprot_reviewed_human_proteome_hgnc_only.v2021_02.fasta.gz  \
    tracked_results/mappings/gencode_34_to_uniprot_2021_02_mappable.coord_mapping.tsv.gz
```
