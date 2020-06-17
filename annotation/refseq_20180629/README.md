## RefSeq 20180629 

### Workflow
Use EBI protein API to get UniParc IDs (but with many missing IDs):

```bash
python scripts/query_uniparc_by_seq.py \
    ../intermediates/refseq_20180629/refseq_protein_ids.list \
    <data_root>/DCC/RefSeq_20180629/RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams.fasta.gz \
    ../intermediates/refseq_20180629/refseq_uniparc_mapping.tsv.gz \
    --json-dir ../intermediates/refseq_20180629/uniparc_jsons
    2> ../intermediates/refseq_20180629/refseq_uniparc_mapping.log
```

Run `playground/query_uniparc_xml.ipynb` to generate the CRC64 checksum for the missing RefSeq IDs.

Use the checksum to query UniProt to get the UniParc XMLs:

```bash
tail -n+2 refseq_checksums.tsv \
    | parallel --bar -j5 --colsep \t \
        "curl -L 'https://www.uniprot.org/uniparc/?query=checksum%3A{2}&format=xml' > xmls/{1}.xml"
```

Parse the XMLs:

```bash
python scripts/parse_uniparc_xmls.py \
    custom_sources/refseq_20180629_missing_uniparc_ids.list \
    ../intermediates/refseq_20180629/ebi_api_missing/refseq_uniparc_mapping.ebi_missing.tsv.gz \
    ../intermediates/refseq_20180629/ebi_api_missing/xmls \
    2> intermediates/refseq_20180629/ebi_api_missing/refseq_uniparc_mapping.ebi_missing.log
```

Run `01_cptac3_refseq_granges_annotation.Rmd` to merge the RefSeq annotations with UniParc.