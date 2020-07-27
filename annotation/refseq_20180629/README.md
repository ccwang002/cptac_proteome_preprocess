## RefSeq 20180629

### Workflow
Get the list of RefSeq protein IDs by:

```bash
gzip -cd custom_sources/cptac3_refseq_20180629/refseq_20180629_human_only_unique_loci.tsv.gz \
    | xsv select -d "\t" refseq_prot_id \
    | tail -n+2 \
    > intermediates/refseq_20180629/refseq_protein_ids.list
```

Run `notebooks/02_make_uniparc_xml_query_tsv.ipynb` to generate the CRC64 checksum for the missing RefSeq IDs.

Download the UniParc XMLs by:

```bash
aria2c -c -j10 --download-result=hide --console-log-level=warn -i uniparc_aria2_downloads.links
```

Parse the XMLs:

```bash
python scripts/parse_uniparc_xmls.py \
    intermediates/refseq_20180629/refseq_protein_ids.list \
    intermediates/refseq_20180629/xmls \
    intermediates/refseq_20180629/refseq_uniparc_mapping.tsv.gz \
    2> intermediates/refseq_20180629/refseq_uniparc_mapping.log
```

Compress the XMLs:

```bash
cd xmls
zstd -T8 -9 *.xml
rm -rf *.xml
cd ..
zip -r -0 xmls.zip xmls/
```

Run `notesbooks/01_cptac3_refseq_granges_annotation.Rmd` to merge the RefSeq annotations with UniParc.