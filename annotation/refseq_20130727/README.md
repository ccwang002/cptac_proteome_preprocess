## RefSeq 20160914
This database has relatively incomplete information. We use UniParc to retrieve the full version ID.

The genomic information and gene level annotation is missing.


### Workflow
Extract the list of protein IDs

```bash
gzcat ~/Box/MyCPTAC/CPTAC_data_collection_v1/DCC/RefSeq_20130727/RefSeq.20130727-Human.contams.categories.tsv.gz \
    | tail -n+2 \
    | cut -f1 \
    > ../intermediates/refseq_20130727/refseq_protein_ids.list
```