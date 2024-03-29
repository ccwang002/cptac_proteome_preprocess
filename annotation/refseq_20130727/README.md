## RefSeq 20160914
This database has relatively incomplete information. We use UniParc to retrieve the full version ID.

The genomic information and gene level annotation is missing.


### Workflow
Extract the list of protein IDs

```bash
gzcat ~/Box/MyCPTAC/CPTAC_proteome_v1/DCC/RefSeq_20130727/RefSeq.20130727-Human.contams.categories.tsv.gz \
    | tail -n+2 \
    | cut -f1 \
    > ../intermediates/refseq_20130727/refseq_protein_ids.no_ver.list
```

Parse the XMLs:

```bash
python scripts/parse_uniparc_xmls.py \
    intermediates/refseq_20130727/refseq_protein_ids.no_ver.list \
    intermediates/refseq_20130727/xmls \
    intermediates/refseq_20130727/refseq_uniparc_mapping.tsv.gz \
    2> intermediates/refseq_20130727/refseq_uniparc_mapping.log
```

Compress the XMLs:

```bash
cd xmls
zstd -T8 -9 *.xml
rm -rf *.xml
cd ..
zip -r -0 xmls.zip xmls/
```
