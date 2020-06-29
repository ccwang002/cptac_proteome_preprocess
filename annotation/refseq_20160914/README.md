## RefSeq 20160914

### Workflow
Use `notebooks/01_cptac_refseq_granges_annotation.Rmd` to obtain the full RefSeq protein ID list.

Run `notebooks/02_make_uniparc_xml_query_tsv.ipynb` to generate the CRC64 checksum for the missing RefSeq IDs.

Use the checksum to query UniProt to get the UniParc XMLs:

```bash
aria2c -c -j10 --download-result=hide --console-log-level=warn -i uniparc_aria2_downloads.links
```

Parse the XMLs:

```bash
python scripts/parse_uniparc_xmls.py \
    intermediates/refseq_20160914/refseq_protein_ids.list \
    intermediates/refseq_20160914/xmls \
    intermediates/refseq_20160914/refseq_uniparc_mapping.tsv.gz \
    2> intermediates/refseq_20160914/refseq_uniparc_mapping.log
```

Run `notebooks/01_cptac3_refseq_granges_annotation.Rmd` to merge the RefSeq annotations with UniParc.